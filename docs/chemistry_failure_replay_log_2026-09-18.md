# Reproducible chemistry-failure logging + log replay (2026-09-18)

## Motivation

The energy-implicit solver (`chemistry_noneq_energy_implicit` in
`chemistry_solver.c`) evaluates the implicit chemistry solve at many trial
temperatures while it searches for the energy root
(`give_g` -> `evaluate_g(..., report_chemistry_failure=0, ...)` ->
`chemistry_noneq_chem_implicit_trial`). That trial entry point has always
been called with `report_failure=0`, so a chemistry Newton solve that fails
to converge at one probed trial temperature produced **no output at all**
unless the *overall* energy solve also ultimately failed. `_with_guess`
(used by the continuation/guess-based recovery path) similarly passes
`report_failure=0`.

This is exactly the gap behind the earlier finding that AREPO particle
`1336657807`'s chemistry non-convergence in the 430000-470000 K trial range
left no trace in production logs: the energy loop simply moved on to another
trial temperature and, in most cases, still converged overall.

## What changed

`chemistry_solver.c` (byte-identical between ONEZONE `src/chemistry/solver/`
and AREPO `src/chemistry_ONECHEM/solver/`, per the 2026-09-18 merge):

- Added a shared helper, `chemistry_print_rerun_line()`, that emits one
  `/*CHEMISTRY_RERUN_V1*/key=value;...` replay line for a given environment,
  `n_nuclei`, `T_gas`, `dt`, and abundance state. This is now the single
  source of the replay-line format. The pre-existing energy-solver failure
  report (`chemistry_noneq_energy_implicit`, printed whenever
  `err > 1.0e-4 || work->diagnostics_enabled`) was refactored to call it;
  output is byte-identical to the previous inline block (verified below).
- `chemistry_noneq_chem_implicit_impl()` now **unconditionally** prints a
  `CHEM_NOT_CONVERGED ...` summary line plus one `chemistry_print_rerun_line()`
  replay line whenever the Newton solve does not converge, *before* the
  existing `if(allow_restarts && report_failure)` gate that controls the
  older, more verbose per-species dump
  (`chemistry_newton_diagnostics_report_failure`). This fires for all three
  entry points (`chemistry_noneq_chem_implicit`,
  `chemistry_noneq_chem_implicit_trial`,
  `chemistry_noneq_chem_implicit_with_guess`) regardless of
  `report_failure`/`allow_restarts`, so a rejected trial temperature during
  energy root-finding is now always reproducible, not just the terminal
  failure of the outer energy loop. Nothing about the Newton control flow,
  convergence criteria, or recovery policy was changed -- this is a pure
  logging addition.

`main.c` (ONEZONE standalone only):

- New CLI mode, `--chem-replay-log <path>`. It scans the given file for
  every line containing the `/*CHEMISTRY_RERUN_V1*/` marker (ignoring
  everything else -- human-readable prefixes, the older verbose per-species
  dump, unrelated program output), parses each one into its own record (so
  fields from one record can never leak into another), and replays each as
  an independent single implicit chemistry solve, printing one
  `CHEM_REPLAY_LOG_RESULT` line per record. A parse error on one line is
  reported and skipped; it does not abort the batch. This means a raw
  captured log -- AREPO stderr, or ONEZONE's own stdout -- can be fed back
  in directly with no manual line-extraction step, and a log with many
  failure events (e.g. many rejected trial temperatures from one energy
  solve) replays all of them in one pass.
- The existing `--chem-step` (and other modes) parser is untouched.

## Testing performed

All testing was done first in a scratch copy (`rsync`'d from the working
tree, including the pre-existing uncommitted `chemistry_newton_control.c`
change, which this work does not touch) before any file in either real
repository was modified.

- Compiled cleanly under three configurations: the high-ionization replay
  build (`-DHELIUM -DCARBON_OXYGEN_HIGH_IONIZATION ... -DCHEM_REPLAY_EXTERNAL_RADIATION`),
  a build matching AREPO's `energy_pivot_initial_states.fixture` conditions
  (`-DHD -DIONIZATION_PHOTON -DLW_RADIATION_LOCAL -DBACKGROUND_LW_RADIATION
  -DTHERMAL_FEEDBACK -DDUST_COOLING`), and a bare default build.
- Regression baseline (`make check-chemistry-initialization`): fails
  identically (`short-step solve failed` x3, exit 2) with and without this
  change, under the high-ionization build -- confirmed pre-existing and
  unrelated to this work, not a regression it introduces.
- AREPO's own `tests/energy_pivot_regression.c` (built via
  `tests/energy_pivot_regression.mk` against `energy_pivot_initial_states.fixture`):
  produces byte-identical stdout before and after the patch (all three
  fixture cases already converge, so no new `CHEM_NOT_CONVERGED` line fires
  on this fixture; this confirms the change is behavior-preserving on the
  happy path).
- Forced a genuine chemistry Newton failure (a disposable, scratch-only
  build with `MAX_ITER_CHEM` temporarily set to 1) to exercise the new
  unconditional log: confirmed `CHEM_NOT_CONVERGED` and the
  `/*CHEMISTRY_RERUN_V1*/` line both fire, and that the emitted replay line,
  fed straight back into `--chem-step` on a normal (non-crippled) build,
  reproduces the same physical result as the original input (round-trip
  fidelity).
- Built a synthetic ~180-line "raw log" (two independent
  `CHEMISTRY_RERUN_V1` records under different `diagnostic_id`s, embedded in
  realistic surrounding program output) and confirmed
  `--chem-replay-log` finds both, parses both, and replays each
  independently and correctly (no cross-contamination between records).
  Also confirmed graceful handling of a marker-free log (`records_found=0`)
  and a malformed record (specific parse-error message, record skipped,
  batch continues).

## Not done here

- The pre-existing uncommitted `chemistry_newton_control.c` fix
  (`constraint_gross`/`implied_gross` gross-scale normalization) is
  untouched and still uncommitted.
- Nothing has been committed in either repository.
- The original 430000-470000 K non-convergence investigation's remaining
  steps (finalize/commit the Newton fix, full regression suite, propagate
  and verify AREPO consistency) are still open.
