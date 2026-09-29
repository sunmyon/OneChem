# Implicit chemistry solver manual

## Highest-priority numerical invariants

1. **No negative abundance is admissible.** This applies to every species,
   including conservation-dependent and equilibrium-intermediate species, at
   initialization, every accepted Newton iterate, and export. Do not use a
   negative-abundance tolerance or silently clip a negative reservoir to zero;
   reject the trial and repair its independent corrections or change the
   conservation chart.
2. **Never shrink the entire correction by one common factor.** Independent
   species are bounded component by component. If the combined correction
   makes an algebraic reservoir inadmissible, repair only responsible
   independent components, leaving unrelated valid components (such as H2)
   unchanged. This applies to the explicit predictor and Newton corrections.
3. A repaired trial must be re-evaluated with the coupled algebraic equations
   and pass the all-species nonnegativity check before it can be accepted.

These invariants take precedence over convergence heuristics and recovery
policies. In particular, a tiny negative dependent abundance is a rejected
state, not a roundoff exception.

This document describes the full-implicit chemistry path implemented in
`src/chemistry/solver/chemistry_solver.c`.  It is intended as a maintenance
manual: it explains the numerical sequence, the responsibility of each helper,
and where compile-time options alter that sequence.

The implementation is divided by responsibility:

- `chemistry_solver.c`: physical preparation and Newton control flow;
- `chemistry_newton_control.c`: residual, bounded-step, and convergence rules;
- `chemistry_newton_reset.c`: deterministic reset compositions;
- `chemistry_newton_diagnostics.c`: printing, failure dumps, and diagnostic-only
  state capture;
- `chemistry_profile.c`: timing and aggregate profiling counters.

Printing and profiling implementations must not be added to
`chemistry_solver.c`.  The Newton driver exposes only one diagnostic event at
each observed stage; formatting, conditional debug output, diagnostic-only
state, and counter storage belong to the dedicated modules above.  RHS and
bounded-step timing calls remain directly adjacent to the operations they
measure.  Thin wrappers that merely forward arguments to one or two functions
are intentionally avoided.

## 1. Public entry points

### `chemistry_noneq_chem_implicit`

Production entry point for one implicit chemistry step.  It calls
`chemistry_noneq_chem_implicit_impl` with:

- no externally supplied Newton guess;
- the production iteration limit (`CHEMISTRY_MAX_ITERATIONS`);
- deterministic and random restart recovery enabled;
- terminal failure reporting enabled.

The input `ChemState` is restored if chemistry cannot produce an acceptable
state.  The timestep is never shortened inside this routine.

### `chemistry_noneq_chem_implicit_trial`

Internal entry point used by energy iteration.  Its numerical chemistry path
is the same as the production entry point, but terminal failure details are not
printed.  The caller can therefore try another energy state without generating
a full failure dump for every rejected trial.

### `chemistry_noneq_chem_implicit_with_guess`

Runs the same implicit equation from a caller-supplied Newton guess.  Restarts
are disabled and the caller supplies the iteration limit.  This is a diagnostic
and controlled-initialization interface; it does not add a previous-step state
to the production chemistry API.

### `chemistry_noneq_chem_implicit_impl`

Common implementation behind the three entry points.  Its control flow has
four phases:

1. Build the normalized chemistry problem and initial Newton guess.
2. Iterate Newton corrections, changing conservation charts or basins only
   through the explicit recovery policy.
3. Certify convergence and report a terminal failure when requested.
4. Apply output floors, reconstruct algebraic species, export the state, and
   evaluate chemical cooling.

The implementation function should remain a control-flow function.  Detailed
linear algebra, diagnostics, and reset recipes belong in the helpers described
below.

## 2. Equation being solved

For each independent non-equilibrium abundance, the backward-Euler equation is

\[
F_i(y) = y_i-y_i^0-\Delta t\,R_i(y)=0,
\]

where `y` is normalized by the hydrogen mass fraction, `y0` is the composition
at the start of the chemistry call, and `R` is the reaction-network right-hand
side.

The reduced Newton system is

\[
\left(I-\Delta t\,\frac{\partial R}{\partial y}\right)\delta y
=\Delta t\,R-(y-y^0).
\]

Conservation-dependent species are not independent columns of this system.
They are reconstructed from the selected conservation chart.  Equilibrium
species are likewise reconstructed algebraically when the selected scheme
requires it.

## 3. Problem construction and initialization

### `ChemImplicitNewtonProblem`

Owns the arrays that are fixed or initialized before Newton iteration:

- `abundance`: current normalized iterate;
- `initial_abundance`: normalized composition at the beginning of the call;
- `cooling_initial_abundance`: initial composition in the units expected by
  chemical cooling;
- `initialization_floor`: floor used only to construct a usable initial state
  and reapplied when exporting a converged state;
- `newton_lower_bound` and `upper_bound`: bounds enforced during Newton steps;
- `total_step`: current `abundance - initial_abundance`;
- `solver_species`: independent variables in the current conservation chart;
- `conservation`: current algebraic conservation system.

An initialization floor is not a lower bound on the mathematical Newton root.
Newton may resolve a nonnegative abundance below the initialization floor.

### `chemistry_prepare_implicit_newton_problem`

Performs the preprocessing that is independent of the initial-guess policy:

1. converts the caller's abundances to normalized abundances;
2. updates equilibrium intermediates required before initialization;
3. constructs initialization floors, Newton bounds, and scale estimates from
   one visible species-parameter table; C/O bounds are assigned from the
   species element metadata;
4. calls `repair_initial_abundances`;
5. constructs the conservation system;
6. initializes or reconstructs algebraic species;
7. records the independent solver-species list and initial composition.

Failure here means that a valid algebraic starting state could not be built.
No Newton iteration has yet been attempted.

### `chemistry_initialize_implicit_newton_guess`

Chooses exactly one initial-guess path:

- without `newton_guess`, build the established explicit predictor;
- with `newton_guess`, copy only independent species from the supplied array
  and reconstruct algebraic species.

It then applies the optional full-implicit intermediate correction and computes
the initial total step.  The supplied input state is never modified in order to
transport additional history into this function.

### `chemistry_correct_full_implicit_intermediate_guess`

Used only by the full-implicit equilibrium-initialization experiment.  It solves
the selected intermediate equilibrium system to improve the initial guess,
then reconstructs conservation-dependent species.  If this produces an invalid
dependent abundance, the complete pre-correction guess is restored.

### `repair_initial_abundances`

Preprocessing that repairs grossly inconsistent input guesses before the
conservation system is created. H, He, and D soft normalization uses one
table-driven loop while retaining the historical policies: H repairs either
sign of a large total error, whereas He and D are reduced only when their totals
are too large. Charge repair and final dominant-reservoir reconstruction remain
explicitly ordered numerical policies. This routine is not the Newton
conservation solver and should not be used to hide an error accumulated during
Newton iteration.

## 4. Conservation and algebraic species

### `chemistry_build_conservation_system`

Selects dependent reservoirs for the enabled conserved quantities and records
their totals. H, He, D, C, and O use one table-driven construction loop; the
compile-time species-family branches only add rows and reservoir candidates to
that table. The chart determines which species are removed from the reduced
Newton vector.

C and O retain one necessary special case. A carbon-bearing molecule contributes
to both conservation rows, so the selected carbon and oxygen reservoirs form a
coupled 2-by-2 block. When most oxygen is in CO or CO2, restricting the oxygen
reservoir to an atomic or oxygen-only species can select a numerically tiny
pivot. `chemistry_coupled_oxygen_reservoir` therefore searches all stored C/O
species for the largest oxygen reservoir whose C/O composition vector is
linearly independent of the selected carbon reservoir. The nonzero determinant
test prevents a singular conservation block.

### `chemistry_update_implicit_algebraic_species`

Reconstructs all algebraic species for the current independent composition.  In
the reduced equilibrium scheme this includes equilibrium intermediates; in all
schemes it reconstructs conservation-dependent reservoirs.

### `chemistry_newton_reselect_conservation_chart`

Changes the dependent reservoir when the current chart reaches an unsuitable
boundary.  It rebuilds the independent-species list without changing the
physical timestep.  Chart reselection is the first recovery attempted after an
invalid bounded Newton step.

### `chemistry_prepare_deterministic_reset_basin`

Applies one named deterministic reset, selects a compatible conservation chart,
optionally corrects the intermediate guess, and recomputes the total step.
Reset recipes themselves live in `chemistry_newton_reset.c`; they do not belong
in the main Newton control flow.

## 5. Residual and Jacobian construction

### `chemistry_network_rhs_gross`

Evaluates the compensated reaction RHS together with gross production and loss
rates.  Production and loss are retained separately because convergence of a
stiff species cannot be assessed reliably from the small net rate alone.

### `chemistry_newton_evaluate_residual`

Builds the convergence diagnostics for the current backward-Euler residual.
This is evaluated before applying the next Newton correction; the subsequent
step norm provides the second part of convergence certification.

The certified residual is normalized so that 1 means 1e-10 of the gross
reaction activity.  It is a sanity veto, not the convergence criterion:
`finish_newton_step` counts a confirmation only when the step norm is within
tolerance, no relevant step was clipped, the clipping hold has expired, and the
certified residual is at most `CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT` (1e3).  A limit
near 1 cannot be met by a stiff species whose Newton correction has fallen below
one ulp of its abundance, so it would reject a converged solution and trigger the
basin-restart cascade.

### `chemistry_build_newton_jacobian`

Owns the Jacobian policy:

1. use the analytic reaction Jacobian when enabled;
2. reduce its response through algebraic species;
3. if algebraic reduction fails, switch this basin to the numerical Jacobian;
4. skip rebuilding only when the tightly guarded LU-reuse condition allows it.

This keeps compile-time analytic/numerical selection out of the main solver.

### `chemistry_build_numerical_newton_jacobian`

Finite-differences each independent solver species.  It prefers a centered
difference, accepts a one-sided difference near a physical boundary, and halves
the perturbation when neither perturbed state can be reconstructed.  Derivatives
below the relative roundoff threshold are set to zero.

### `chemistry_assemble_reduced_newton_system`

Constructs `I - dt*dR/dy` in the column-major layout expected by LAPACK and
builds the right-hand side `dt*R - total_step`.  Algebraic species never enter
this reduced vector as independent unknowns.

### `chemistry_solve_reduced_newton_system`

Solves the dense reduced system:

- without LU reuse, calls `dgesv`;
- with LU reuse, calls `dgetrf` when rebuilding and `dgetrs` to solve.

The accepted dense intermediate solve and dense reduced Newton solve are
intentional.  This helper owns the LAPACK-specific storage and status handling.

### `chemistry_should_rebuild_newton_matrix`

Encapsulates the conservative LU-reuse rule.  A factorization is reused only
for the mandatory second convergence confirmation, after the previous Newton
correction is already below `1e-12`.  Otherwise the Jacobian and factorization
are rebuilt.

### `chemistry_unpack_reduced_newton_step`

Maps the LAPACK solution back into a full `N_SP` array.  Algebraic entries are
left at zero because they are reconstructed rather than independently stepped.

## 6. Bounded step and convergence

### `chemistry_newton_apply_step`

Applies a bounded Newton correction and reconstructs algebraic species.  It
returns both validity and a per-species clipping record. Independent
corrections are clipped or repaired individually; there is no common
line-search factor. Every accepted species abundance must be nonnegative.

### `chemistry_newton_material_step_norm`

Computes the dimensionless step norm over species relevant to convergence and
reports the limiting species.  Trace abundances and optionally initialized
equilibrium intermediates follow their configured convergence policy.

### `chemistry_step_has_relevant_clip`

Distinguishes a material boundary clip from clipping of a numerically negligible
species.  In the full-implicit equilibrium-initialization scheme, initialized
intermediates can be excluded from this test when
`CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE=0`.

### `chemistry_update_newton_convergence_state`

Updates four related pieces of state in one place:

- the post-clipping convergence hold;
- the consecutive boundary-step count;
- the request to reselect the conservation chart;
- the number of consecutive convergence confirmations.

A materially clipped step cannot itself confirm convergence.  After clipping,
the configured number of clean iterations must pass before confirmation resumes.

### `chemistry_recompute_total_steps`

Recomputes `abundance - initial_abundance` after every accepted correction and
zeros algebraic entries.  Recomputing from the two states avoids accumulating
roundoff through a long sum of Newton corrections.

## 7. Recovery policy

### `chemistry_prepare_newton_basin_entry`

Normalizes a proposed basin entry within the allowed bounds and reconstructs
algebraic species.  An invalid entry is rejected before Newton iteration.  A
valid entry is passed once to the diagnostics module and its algebraic total
steps are cleared.

### `chemistry_newton_apply_deterministic_reset`

Applies a reset selected by `ChemNewtonResetIndex`.  The enum declaration order
is the deterministic search order.  Named strategies replace the former raw
numbers such as 2, 3, and 18.

The C/O build and non-C/O build intentionally use different names for the
second and third strategy because those slots have different physical meanings.

### Random Newton restart

The Newton driver calls `chemistry_find_random_newton_guess` directly after
deterministic basins are exhausted and increments the reset sequence only when
a usable random guess is found.  The random search changes the guess at the
same physical timestep.

### Basin iteration budget

Each deterministic basin currently receives 70 Newton iterations.  Reaching
that budget without confirmation advances to the next named reset.  Random
search follows the deterministic sequence.  This policy is separate from the
global chemistry iteration limit.

## 8. Diagnostics, profiling, and output

### `ChemNewtonDiagnostics`

Holds diagnostic-only state for one Newton solve.  In a normal build it stores
no abundance array and occupies only one byte.  The `N_SP` basin snapshot is
compiled in only for `CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC`.

### `chemistry_newton_diagnostics_enter_basin`

Is the single basin-entry diagnostic event.  It prints the compact H/electron
composition when runtime diagnostics are enabled and, only in a random-restart
diagnostic build, records the validated basin entry.

### `chemistry_newton_diagnostics_record_iteration`

Is the single per-iteration diagnostic event.  It owns both the compact
iteration output and optional Jacobian output.  Additional output is compiled
only for `CHEMISTRY_DEBUG=1`; the Newton flow does not contain separate print calls.

### `chemistry_newton_diagnostics_finish`

Is the single normal-completion diagnostic event.  It emits a successful random
basin only when that diagnostic is compiled.

### `chemistry_newton_diagnostics_report_failure`

Emits a reproducible terminal-failure record: physical conditions, complete
input composition, final abundances, Newton steps, residuals, and algebraic
status.  This output is deliberately outside the main control-flow function.

### `chemistry_update_implicit_timescale`

Computes the minimum abundance-change timescale over independent species for
the caller's timestep control diagnostics.

### `chemistry_export_implicit_solution`

Applies output floors to independent species, reconstructs algebraic species,
rejects any negative dependent reservoir, converts abundances back to external
units, evaluates chemical cooling, and optionally updates the
continuum-opacity level populations.

### Profiling API

`chemistry_profile.h` exposes the same functions in every build.  With
`CHEMISTRY_PROFILE` disabled they are inline no-ops.  This removes profiling
`#ifdef` blocks from the numerical control flow.  Timer start/stop calls are
placed immediately around the RHS, Jacobian, linear solve, and bounded Newton
step; no profile-only forwarding wrappers are used.  Disabled profile calls are
optimized away.

## 9. Compile-time options

The Makefile options and their principal effects are:

| Make option | Values | Effect |
|---|---|---|
| `CARBON_OXYGEN_CHEMISTRY` | `0`, `1` | Includes C/O species and reactions. |
| `CHEMISTRY_CARBON_OXYGEN_CONSERVATION` | `0`, `1` | Includes C and O constraints in the Newton conservation system. |
| `CHEMISTRY_EQUILIBRIUM_INITIALIZATION` | `0`, `1` | Enables equilibrium-based treatment or initialization of intermediates. |
| `CHEMISTRY_NEWTON_SYSTEM` | `FULL`, `REDUCED` | Solves initialized intermediates in the full Newton system or eliminates equilibrium intermediates algebraically. |
| `CHEMISTRY_NEWTON_LU_REUSE` | `0`, `1` | Enables the guarded reuse of a reduced Newton factorization. |
| `CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE` | `0`, `1` | Includes or excludes initialized equilibrium species from convergence clipping checks. |
| `CHEMISTRY_EXPLICIT_UPDATE` | `SEQUENTIAL`, `SIMULTANEOUS` | Selects ExplicitSolver update order; default SIMULTANEOUS. Does not select the Newton predictor. |

The Makefile translates these options into internal preprocessor definitions.
Callers should use the Makefile options, not define the internal macros directly.

## 10. Which `#ifdef` branches should remain

The main implicit solver no longer contains profiling or random-diagnostic
preprocessor branches.  Remaining branches fall into three categories:

1. **Species-layout branches** such as `HELIUM`, `HD`, and `H2_PLUS`. These
   change the valid species indices. They are concentrated in the X-macro species
   catalog and in optional rows of visible data tables. `sSTORED_END` is counted
   from `CHEM_SPECIES_STORED_LIST`; it is no longer maintained by a parallel
   preprocessor ladder. Full-implicit intermediates are stored, while reduced
   equilibrium intermediates are internal.
2. **Numerical-policy branches** such as analytic versus numerical Jacobian,
   full versus reduced equilibrium treatment, and LU reuse.  These are confined
   to the helper that owns the affected numerical operation.
3. **Optional diagnostics or physics output** such as detailed debug printing
   and continuum-opacity populations.  These are confined to dedicated helper
   functions and do not interrupt the main Newton flow.

Charge neutrality in the explicit update is reconstructed by one loop over the
stored species updated by that scheme, using charge metadata. Initialization
parameters and element-conservation rows are likewise data-driven. The remaining
species-family directives decide which data exists; they do not split the main
calculation flow.

## 11. Maintenance rules

- Do not shorten the physical timestep inside chemistry to conceal convergence
  failure.
- Recovery changes the guess or conservation chart at the same timestep.
- Preserve the sequential explicit predictor unless the explicit-update option
  requests simultaneous updates.
- Do not pass previous-step history through the public chemistry input solely to
  improve Newton convergence.
- Keep reset recipes in `chemistry_newton_reset.c` and refer to them through the
  named enum.
- Keep LAPACK storage details, profiling, and diagnostic formatting out of the
  main control-flow function.
- Validate structural changes against saved evolution or batch output, not only
  by successful compilation.
