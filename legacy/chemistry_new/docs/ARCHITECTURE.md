# Portable chemistry architecture

Status: proposed design contract.  This document must be agreed before the
AREPO integration path is changed.

## 1. Goal and non-goals

The goal is a self-contained `chemistry_new/` directory which can be copied
unchanged under an AREPO source tree and can also be built by the standalone
standalone program. The same chemistry core, reaction network, rate functions,
cooling functions, tables, and default numerical algorithms must be used in
both hosts.

Copying the directory is not sufficient by itself: each host still adds a thin
build include and calls its own adapter.  No core source file may include an
AREPO or standalone header, inspect a hydrodynamic particle, perform an MPI
collective, or convert code units.

This refactor does not initially redesign the physical model, change a reaction
rate, change a tolerance, replace the existing implicit algorithm, or alter the
meaning/order of stored abundances.  Those are separate scientific changes and
must not be hidden inside the portability work.

## 2. Present boundary problems

The current tree has useful responsibility-based directories, but its actual
dependency direction is still host-first:

- core files include `allvars.h`, `proto.h`, and sometimes `tags.h`;
- runtime inputs are read from `All` (for example `T_CMB`, `Epsilon`, radiation
  rates, table directory, and unit factors);
- diagnostics use `ThisTask`, `FlagPrint`, and direct `printf` calls;
- some table readers allocate through `mymalloc` and contain MPI behavior;
- `chemistry_driver.c` combines the portable advance loop with the legacy
  `chemdata_in`/`chemdata_out` conversion;
- `chem_balance.c` owns hydrodynamic particle access and therefore is an
  adapter, despite currently living under `interface/`;
- public structure layout and species count are controlled by host compile
  definitions.

Putting these dependencies behind a new header which itself includes
`allvars.h` would only rename the coupling.  The values and services must cross
an explicit API boundary.

## 3. Dependency rule

Dependencies point inward only:

```text
AREPO hydrodynamics ──> AREPO adapter ─┐
                                      ├──> public API ──> chemistry core
standalone evolution ──> standalone adapter┘                │
                                                             ├── species/network/rates
                                                             ├── cooling/thermodynamics
                                                             ├── tables
                                                             └── solver implementations
```

Allowed dependencies are:

| Layer | May depend on | Must not depend on |
|---|---|---|
| public headers | C standard types, module-owned generated configuration | host headers, MPI, GSL headers |
| physics modules | public/internal chemistry headers, C math | particles, code units, global host state |
| solver implementations | RHS/Jacobian interfaces, optional solver-specific library | host state and host callbacks other than logging |
| adapters | public API and their host | another host's headers |
| standalone driver | standalone adapter/public API | AREPO |
| build integration | module source lists and host flags | copied or generated source code |

A source-level test will enforce that `src/`, `species/`, `reactions/`,
`rates/`, `cooling/`, `thermodynamics/`, and `solvers/` do not include
`allvars.h`, `proto.h`, or `tags.h`.

## 4. Data ownership and lifecycle

Four kinds of data are intentionally separate.

### `ChemistryConfig`

Immutable model and numerical policy chosen at initialization:

- enabled network/profile and cooling processes;
- subcycle safety factor and existing solver tolerances;
- abundance floors;
- table source and table-loading policy;
- default column/escape closure;
- error and logging callbacks.

Compile-time feature selection may remain during the compatibility phase, but
host preprocessor symbols are translated once by a module-owned
`chemistry_features.h`.  Physics source must not branch on an AREPO-only macro.
The long-term representation is a stable species enumeration plus an active
species/reaction mask; that change is postponed until after bitwise/tolerance
regression is established.

### `ChemistryUnits`

The public physics API uses CGS exclusively.  A units structure exists for
adapter conversion and validation, not for use in rate equations:

```c
typedef struct {
  double mass_g;
  double length_cm;
  double time_s;
  double energy_erg;
  double density_cgs;
  double specific_energy_cgs;
} ChemistryUnits;
```

AREPO and standalone adapters construct this structure from their host settings.
The core receives `rho_cgs`, seconds, kelvin, rates in s^-1, columns in cm^-2,
and specific thermal energy in erg g^-1.  Cosmological scale factors are
converted by the adapter; the core never reads `HubbleParam` or `cf_a3inv`.

### `ChemistryContext`

An opaque, initialized object which owns read-only shared data:

- validated configuration;
- reaction metadata and abundance floors;
- loaded/interpolated tables;
- selected solver operations;
- callback registrations.

After initialization it is read-only and can be shared by threads.  There is
no singleton context and no mutable `All` equivalent.  Initialization returns a
status instead of terminating the host.

### Step context, per-call state, and workspace

`ChemistryStepContext` is an immutable snapshot of values shared by all cells
in one chemistry sweep.  It is prepared once when the host's global time or
background changes and is then passed by pointer to every cell advance.  It
contains values such as CMB temperature, scale factor, uniform background
radiation rates, and derived coefficients such as the factor proportional to
`scale_factor^-4` used by Compton cooling.  It is caller-owned, requires no
allocation, and may be read concurrently.

This object is deliberately separate from `ChemistryContext`.  The latter is
unchanged for the lifetime of a run, while a step context changes with global
time.  Mutating a shared context once per step would hide temporal state and
would prevent two snapshots from being evaluated concurrently (for example in
a test, interpolation scheme, or task-parallel host).

`ChemistryState` is the state which advances: abundances, specific thermal
energy, and gamma.  Its abundance convention remains the current convention
during migration and is documented per species.

`ChemistryCellEnvironment` is an immutable snapshot of cell-specific external
conditions for one advance: physical density, metallicity, local radiation or
photo-rate additions, and any host-supplied column information.  It does not
repeat CMB temperature or scale factor.  It contains no unit-system fields and
no particle identifier required by the physics.  An optional opaque diagnostic
tag can be carried to logs.

`ChemistryWorkspace` owns all mutable scratch arrays, rates, continuum state,
and solver statistics.  The caller provides one workspace per simultaneous
advance (normally per thread), so advancing two cells cannot race on global
scratch storage.

`ChemistryResult` contains output timescales, cooling/heating diagnostics,
convergence information, and flags.  Diagnostics are output, not hidden
mutable state.

## 5. Public API

The public surface is small and versioned.  Detailed rate and cooling helpers
remain internal unless a demonstrated host use requires them.

```c
typedef enum {
  CHEMISTRY_OK = 0,
  CHEMISTRY_EINVAL,
  CHEMISTRY_EIO,
  CHEMISTRY_ENOMEM,
  CHEMISTRY_ENOCONVERGE,
  CHEMISTRY_EINTERNAL
} ChemistryStatus;

ChemistryStatus chemistry_context_create(
    const ChemistryConfig *config,
    ChemistryContext **context);

void chemistry_context_destroy(ChemistryContext *context);

size_t chemistry_workspace_size(const ChemistryContext *context);

ChemistryStatus chemistry_workspace_init(
    const ChemistryContext *context,
    void *storage,
    size_t storage_size,
    ChemistryWorkspace **workspace);

ChemistryStatus chemistry_step_prepare(
    const ChemistryContext *context,
    const ChemistryGlobalEnvironment *global_environment,
    ChemistryStepContext *step);

ChemistryStatus chemistry_advance(
    const ChemistryContext *context,
    const ChemistryStepContext *step,
    ChemistryWorkspace *workspace,
    const ChemistryCellEnvironment *cell,
    double duration_s,
    ChemistryState *state,
    ChemistryResult *result);
```

Caller-supplied workspace storage avoids allocation in a per-cell hot path.
An optional convenience allocator can be supplied, but is not the primitive
API.  `chemistry_step_prepare` is called once per global chemistry sweep, not
once per cell.  `chemistry_advance` receives pointers, so shared values are not
copied for each cell.  It has transactional failure semantics: on a hard error,
`state` is unchanged; a non-converged step reports the last accepted state only
when explicitly marked in `result`.

The legacy `compute_abundances_d()` signature is not this API.  It belongs in
the AREPO compatibility adapter so existing AREPO call sites can migrate
without simultaneously changing their data flow.

## 6. Host policy versus chemistry policy

The following choices belong outside the core:

- reading/writing `SphP`, `P`, or active-cell lists;
- conversion between entropy/conserved energy and CGS specific energy;
- cosmological density/time conversion;
- choosing a particle time step;
- MPI distribution, rank-specific output, and host memory allocation;
- parameter-file parsing;
- deciding whether a cell is eligible for chemistry.

Host values are classified by lifetime rather than copied wholesale from
`All`:

| Lifetime | Chemistry object | Current examples |
|---|---|---|
| run/initialization | `ChemistryConfig` / `ChemistryContext` | table root, abundance floors, solver tolerances, enabled model |
| global time/sweep | `ChemistryGlobalEnvironment` / `ChemistryStepContext` | `T_CMB`, `cf_atime`, uniform `J21`/`G0`, uniform photo-rates |
| cell | `ChemistryCellEnvironment` | physical density, metallicity, local radiation, columns |
| evolving cell | `ChemistryState` | abundances, thermal energy, gamma |
| temporary call | `ChemistryWorkspace` | rates, Jacobian, trial state, convergence scratch |

Unit factors remain in the adapter (or its adapter configuration), because the
core receives CGS values.  A uniform background and a local source contribution
are kept distinct in the two environment levels and combined by the core,
avoiding repeated storage without losing support for spatially varying fields.

The following belong inside the core because they change the local chemistry
solution:

- reaction and inverse-rate evaluation;
- cooling/heating evaluation;
- thermodynamic closure from the chemical state;
- positivity and element/charge conservation constraints;
- chemistry subcycling and the current explicit/implicit/Saha selection;
- convergence tests and accepted-step rollback.

Column density and escape probability are a deliberate seam.  The core accepts
host-provided columns when available, while a named `JEANS_LOCAL` policy
reproduces the current local Jeans-length closure.  This prevents a future
transport-derived column from requiring changes to reaction or cooling code.

## 7. ODE solver abstraction

The abstraction separates the chemical problem from the numerical method.  It
does not immediately replace the existing method with a third-party solver.

The network layer exposes an internal problem description:

```c
typedef struct {
  size_t variable_count;
  void *user_data;
  ChemistryStatus (*rhs)(double t, const double *y, double *dydt,
                         void *user_data);
  ChemistryStatus (*jacobian)(double t, const double *y, double *jac,
                              size_t leading_dimension, void *user_data);
  ChemistryStatus (*project)(double *y, void *user_data);
} ChemistryOdeProblem;

typedef struct {
  ChemistryStatus (*advance)(void *solver_state,
                             const ChemistryOdeProblem *problem,
                             double t0, double dt, double *y,
                             ChemistrySolverStats *stats);
  void (*reset)(void *solver_state);
} ChemistryOdeOps;
```

`rhs` and `jacobian` own no host state.  `project` applies the existing
positivity and conservation rules.  Energy-implicit coupling is represented as
a coupled problem/strategy, not as an AREPO special case.

The first implementation is `solvers/legacy/` and must preserve the present
iteration order, linear solve, tolerances, fallback, and rollback.  A GSL or
CVODE implementation may be added later behind the same operations and tested
against the legacy solver.  Solver selection is a context configuration, never
a global compile-time side effect in the host.

## 8. Table, allocation, logging, and errors

The core cannot assume a current working directory or an MPI rank.  Context
creation receives either a table root path or an I/O provider callback.  Table
files remain under `chemistry_new/tables/`, and the adapter resolves their
installed location.

All initialization allocation goes through an optional allocator in
`ChemistryConfig`; the default is `malloc/free`.  The hot advance path performs
no allocation.

Core code never calls `printf`, `abort`, `endrun`, or `terminate`.  It returns a
status and may emit a structured log event through a callback.  Rank filtering
and formatting are adapter responsibilities.  Assertions may remain only for
programmer invariants in test/debug builds.

## 9. Target directory tree

```text
chemistry_new/
  README.md
  include/chemistry/
    chemistry.h              public lifecycle and advance API
    chemistry_config.h       configuration and callbacks
    chemistry_species.h      species IDs and abundance convention
    chemistry_state.h        units/environment/state/result types
    chemistry_solver.h       public solver selection and statistics
    chemistry_version.h
  src/
    context.c
    advance.c
    internal.h
  species/
  reactions/
  rates/
  cooling/
  thermodynamics/
  tables/
  solvers/
    ode_internal.h
    legacy/
    saha/
  adapters/
    arepo/
      chemistry_arepo_adapter.c
      chemistry_arepo_adapter.h
    standalone/
      chemistry_standalone_adapter.c
      chemistry_standalone_adapter.h
  standalone/
    main.c
  mk/
    chemistry_new.mk          host-neutral source/header lists
    arepo.mk                  thin AREPO variable additions
    standalone.mk
  tests/
    unit/
    regression/
    integration/
    fixtures/
  docs/
    ARCHITECTURE.md
    MIGRATION.md
```

Physical files move only after include dependencies are removed.  Directory
moves and numerical edits are separate commits so regression failures can be
localized.

## 10. Adapter responsibilities

The standalone adapter:

1. reads the standalone cell and global run parameters;
2. constructs CGS units, environment, and state;
3. calls the public API;
4. writes accepted state and diagnostics back;
5. leaves collapse/density evolution in the standalone host.

The AREPO adapter performs the analogous conversion for one AREPO gas cell and
retains the legacy `compute_abundances_d()` wrapper during transition.  AREPO's
`chem_balance` loop is host integration and therefore lives in this adapter,
not in the core.  The adapter may use `allvars.h`, `proto.h`, `tags.h`, MPI, and
AREPO allocation facilities; no other module directory may do so.

## 11. Numerical-preservation migration

Each stage has a mandatory gate.  No stage combines a physics change with an
architecture change.

1. **Freeze reference outputs.** Record compiler/configuration, full-survey
   output, selected trajectories, table hashes, final states, diagnostics, and
   convergence flags.  Establish exact comparison where deterministic and
   explicit relative/absolute tolerances elsewhere.
2. **Introduce public types without changing calls.** Add module-owned headers
   and compile-time layout assertions against the legacy structures.
3. **Extract the standalone adapter.** Move only conversions and particle/global
   access; keep the same `compute_abundances_d()` implementation and call order.
4. **Pass runtime values explicitly.** Replace `All` reads one field at a time,
   with a regression run after each coherent group (units, radiation, solver
   policy, tables).
5. **Remove process services.** Replace direct allocation, output, termination,
   and MPI with context callbacks/status codes.
6. **Create context/workspace ownership.** Move tables and scratch state out of
   globals while retaining their values and initialization order.
7. **Wrap the legacy solver operations.** The solver abstraction first calls
   the unchanged algorithm; only its dependencies move.
8. **Move files into the target tree.** This is a mechanical step after the
   dependency test passes.
9. **Add the AREPO adapter and build include.** First compile/link smoke test,
   then a one-cell adapter equivalence test, then the existing AREPO problem
   regression.  AREPO call sites are changed only after equivalence succeeds.

Required tests are:

- species metadata, charge/element conservation, and reaction stoichiometry;
- rate/cooling golden values over boundary and representative grids;
- RHS/Jacobian finite-difference checks;
- explicit, implicit, Saha, rollback, and non-convergence solver tests;
- table path/loading failure tests;
- standalone adapter round-trip and full-survey regression;
- AREPO adapter unit-conversion and one-cell equivalence tests;
- a forbidden-host-include source check;
- a build matrix for the currently supported feature profiles.

## 12. Acceptance criteria

The portability refactor is complete only when:

- copying `chemistry_new/` plus adding the documented build include is enough
  to compile it in either host;
- public headers compile in a translation unit with no host headers;
- core files contain no host globals, particle access, MPI, or code-unit
  conversion;
- context sharing and workspace exclusivity are documented and tested;
- the standalone reference survey remains within its established regression
  criteria;
- the legacy and new adapter paths agree on identical cell inputs;
- failures return statuses and do not terminate the host;
- the old AREPO chemistry implementation can remain selectable until production
  comparison is complete.

## 13. Decisions to settle before implementation

Three decisions materially affect the implementation and should be made
explicitly:

1. **Species layout:** preserve the current feature-dependent array layout for
   the first portable release (lowest migration risk), or adopt a stable full
   species layout immediately (cleaner ABI, larger refactor).
2. **Workspace ownership:** caller-provided raw storage as proposed above, or an
   allocating `create/destroy` API only (simpler API, less control in AREPO hot
   paths).
3. **Column closure:** keep `JEANS_LOCAL` as the portable default for exact
   compatibility, while allowing the adapter to provide columns, or require
   every host to provide all columns from the start.

Recommended first-release choices are feature-dependent layout, caller-provided
workspace storage, and `JEANS_LOCAL` as the compatibility default.  They
minimize numerical and integration risk while leaving clean seams for the
long-term forms.
