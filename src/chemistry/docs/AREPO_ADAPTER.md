# AREPO boundary

Production integration consists of two files:

- `/Users/j/Desktop/codex/1ZONE/src/chemistry/adapters/arepo/chem_balance.cc`
- `/Users/j/Desktop/codex/1ZONE/src/chemistry/adapters/arepo/chemistry_arepo.h`

The source owns AREPO lifecycle hooks, background radiation, particle import/export,
MPI communication, and direct calls to the public chemistry API. The header exposes
host lifecycle/radiation hooks and includes the public API. It defines no duplicate
cell/state/result types and does not expose internal core types.

Read the main flow as:

1. `chemistry_arepo_host_init`: create the Model and radiation resources.
2. `cooling_only`: prepare the shared step, select active particles and compute durations.
3. `particle2in_arepo`: convert host units directly into `ChemistryCellEnvironment`
   and `ChemistryState`, backed by a cell-local abundance buffer.
4. Construct a typed cell-local ChemWorkspace and call `Driver::advance` directly.
   The diagnostic retry uses the same workspace; the API resets it on entry.
5. `out2particle_arepo`: convert thermal energy to AREPO units and commit composition.
6. `chemistry_arepo_host_destroy`: release run-lifetime resources.

With CHEM_BALANCE enabled, the MPI payload stays host-owned. `chem_evaluate`
converts that payload directly into the same public API types and calls
`Driver::advance`; it does not use an intermediate adapter API.

Driver handles validation, internal-state conversion, subcycling, solver choice
and successful-state export. AREPO uses Model and Driver directly, with no core
C compatibility API. Only lifecycle/particle hooks called by the older C host
retain C linkage.

Boundary tests compile both particle and CHEM_BALANCE branches. Runtime tests
exercise background radiation, density/time/energy conversion, diagnostic retry,
and particle export for the local path and the imported-cell evaluation path.
These tests use mock host storage and a controlled public API evaluator; they do
not run MPI communication or a full AREPO simulation. Native C++ API tests separately
link the real chemistry implementation.
