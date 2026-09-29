# Initialization contract

> Historical migration notes. Current C++ ownership and calling contracts are in
> [/Users/j/Desktop/codex/1ZONE/src/chemistry/README.md](/Users/j/Desktop/codex/1ZONE/src/chemistry/README.md).
> References below to the transitional C context or InitChem do not describe the current API.

## Decision

Chemistry core uses explicit eager initialization.  `chemistry_context_create`
is the only public operation allowed to load files, allocate run-lifetime
storage, construct tables, or validate the compiled network.  A successful
context is immutable and ready for every enabled chemistry path.

`chemistry_advance` and `chemistry_step_prepare` must not perform hidden
initialization.  In particular they must not open files, allocate a missing
table, or mutate a process-wide `initialized` flag.

This choice is important for AREPO because initialization can fail before the
time integration begins, occurs outside threaded cell loops, and has a clear
MPI error boundary.  Chemistry core itself performs no MPI collective and does
not terminate the process.  It returns a `ChemistryStatus`; the AREPO adapter
decides how ranks report and combine that failure.

## Lifecycle

```text
host parses configuration and establishes units
  -> chemistry_config_defaults
  -> host sets table_root and timescale_safety_factor
  -> chemistry_context_create       (once per process/rank)
  -> chemistry_workspace_init       (once per worker/concurrent call)

for each host-wide chemistry sweep:
  -> chemistry_step_prepare         (once for shared global conditions)
  -> chemistry_advance              (once per active cell, no I/O/allocation)

after every chemistry call has finished:
  -> chemistry_context_destroy
```

Context creation is transactional.  On any parse, validation, allocation, or
I/O failure it releases partial state, returns a precise status, and leaves the
output context pointer null.  It never publishes a partially initialized
context.

## What is initialized where

| Data or operation | Owner / phase |
|---|---|
| reaction network parse and stoichiometry validation | context creation |
| temperature grid and forward/inverse reaction-rate tables | context creation |
| molecular, continuum, metal, and dust cooling tables | context creation when compiled |
| LW opacity used by chemistry transport approximations | context creation when compiled |
| background spectrum and spectrum-to-rate lookup | sibling `RadiationRateContext` initialization |
| abundance floors and solver tolerances | compiled chemistry configuration |
| fixed-temperature derived quantities | context creation when compiled |
| CMB temperature and scale factor | step preparation |
| uniform LW intensity to derived-rate conversion | host, when the background changes |
| derived uniform LW rates | step preparation |
| density, metallicity, local radiation | cell environment |
| rate vectors, Jacobians, Newton scratch | workspace |

Compile-time-disabled modules are not initialized and allocate no storage.
This provides the useful memory saving without runtime lazy initialization.

## Table access

The public implementation uses `ChemistryConfig.table_root`. All required
files are opened synchronously during context creation. The path is not
retained, so the caller's string need not remain alive after
`chemistry_context_create` returns.

During the legacy migration,
`InitChem(const ChemistryInitConfig *, ChemRunContext *)` follows the same
usage rule and accepts the directory containing `tables/`. It currently
creates the context-owned reaction network; `DestroyChem` releases that
network.  It also allocates and eagerly fills the context-owned temperature,
reaction-rate, atomic-cooling, and molecular-line cooling tables. `DestroyChem`
also owns the eagerly loaded dust- and continuum-opacity table objects.
This removes `All.ChemistryTableDir` without placing evolving run or step
values in the initialization configuration.

`InitChem` deliberately does not initialize `RadiationRateContext`,
derive uniform LW rates, or write host abundance-floor arrays.  The host owns
that sibling context and places its derived values into the step context when
the host-wide environment changes. Abundance floors are part of the compiled
chemistry configuration. Thus neither is a side effect of loading tables.

Table-reading code returns a `ChemistryStatus` and does not read
`All.ChemistryTableDir` or call the host's termination routine. Some legacy
readers still print an error message to `stderr`; a logging callback is not
part of the current public API.

## Legacy compatibility

The standalone adapter may expose an `ensure configured` check because the old
`compute_abundances_d` signature cannot carry a context.  That check only
verifies that explicit host initialization already happened; it must not
trigger initialization from the first cell call.

The historical cooling-time fallbacks guarded by `flag_dust_init`,
`flag_tabulate`, `flag_read_opacities`, and `flag_read_SD93` have been removed.
The transitional `InitChem` eagerly initializes every enabled standalone cooling
table before any cell evaluation.  A source check prevents those process-wide
lazy flags from being reintroduced.

The table readers still emit some legacy `stderr` diagnostics when input is
malformed. Their storage is context-owned, failures return through
`ChemistryStatus`, and collective coordination remains outside the core.
