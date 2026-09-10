# AREPO adapter and `chem_balance`

`chem_balance` remains an AREPO orchestration module.  It is not part of the
chemistry core and must not implement reaction, cooling, thermodynamic, or
abundance-normalization formulas.

## Ownership and lifetimes

- Each MPI rank initializes one `ChemArepoContext` after AREPO units and the
  table path are known.  Every rank loads/constructs its own rank-local tables;
  the context is not sent between ranks.
- `ChemArepoContext` is read-only while cells are advanced and is destroyed at
  run shutdown.
- `cooling_only()` constructs one `ChemGlobalStep` on its stack at the start of
  a chemistry sweep.  CMB temperature, scale factor, simulation time, and
  uniform radiation rates are written once, then shared by pointer.
- Each evaluation has a private `ChemWorkspace`.  No workspace is shared by
  simultaneous evaluations.

The current `chem_balance` evaluation loop remains intentionally serial while
the first AREPO integration is validated. Diagnostic enablement is now owned by
each `ChemWorkspace`, and a failed cell is retried with only that workspace's
diagnostics enabled. The former `FlagPrint` and `count_print` process globals
are gone. OpenMP can be reintroduced separately after the host counters and
diagnostic-output ordering are given explicit parallel semantics.

The existing no-argument `cooling_only()` hook means AREPO needs one
rank-lifetime owner for `ChemArepoContext`.  `chemistry_arepo_host.c` provides
that single typed object.  It is an adapter-owned handle, not hidden core state:
all core calls still receive it explicitly.  Its table members are not placed
back into `All` or exposed as separate globals.

`chemistry_arepo_host_prepare_step()` is the only place that maps shared AREPO
values into `ChemGlobalStep`. The radiation sibling caches uniform-background
rates and recomputes them only when J21 or effective temperature changes. Its
explicit HD rate is copied without reconstructing it from the H2 rate.

## Call flow

```text
cooling_only()
  prepare ChemGlobalStep once
  select active particles
  chem_balance()                         [when MPI exchange is enabled]
    pack AREPO wire input
    exchange wire input
    for each local/received cell:
      wire -> ChemArepoCell + ChemArepoState
      chemistry_arepo_advance(context, step, cell, dt, state, workspace, result)
      state + result -> wire output
    exchange wire output
    commit state to SphP
```

The non-`CHEM_BALANCE` path follows the same conversion and advance function,
but skips wire exchange.  It must not retain a second `compute_abundances_d()`
path.

## Boundary responsibilities

`chem_balance` owns:

- active-particle filtering and timestep conversion;
- physical-density calculation from AREPO's density storage;
- conversion of AREPO specific thermal energy to and from CGS, preserving the
  historical operation order of each path;
- MPI wire structures and exchange;
- copying `SphP[i].Abund[]` without renormalization;
- committing `Utherm`, `Energy`, `AbundMass`, gamma, and diagnostics.

`chemistry_arepo_advance` owns:

- construction of `ChemEnvironment` and `ChemState`;
- the single call into the chemistry core;
- conversion of the accepted state and workspace diagnostics back to adapter
  values.

The transitional standalone `compute_abundances_d()` function is not an AREPO
API. It contains standalone singleton and unit assumptions and must not be linked
as the implementation used by `chem_balance`.

## Required host header switch

Current AREPO includes `chemistry_SC/chemistry.h` directly from `allvars.h`.
That legacy header defines the species enum and cannot coexist with the
portable definition.  The required host-side change is a conditional include:

```c
#ifdef CHEMISTRY_PORTABLE
#include "chemistry/chemistry_core.h"
#else
#include "chemistry_SC/chemistry.h"
#endif
```

`chems.h` may continue to reserve `TGCHEM_NUM_ABUNDANCES` before this include.
`chemistry_core.h` now verifies at compile time that the reserved host capacity
equals chemistry's feature-dependent stored-species count.

## Minimal AREPO lifecycle changes

In addition to the header switch, AREPO must:

1. compile `CHEMISTRY_PORTABLE_CORE_SOURCES` and
   `CHEMISTRY_AREPO_ADAPTER_SOURCES` from `mk/arepo.mk`, while excluding the old
   `src/chemistry_SC` implementation sources;
2. replace the old `InitChem()` call in `begrun.c` with
   `chemistry_arepo_host_init()` after units and parameters are available, and
   treat a false return as initialization failure;
3. call `chemistry_arepo_host_destroy()` during orderly shutdown;
4. retain the existing `cooling_only(void)` and `compute_temp(int)` call sites;
   both symbols are supplied by the portable `chem_balance.c` boundary.

LW continuum opacity remains an AREPO radiation-transport responsibility. The
portable initializer does not call `tabulate_opacities_in_LW_band()`; the
selected `RayTrace_*` module must call its own implementation during radiation
initialization.

No AREPO caller receives table pointers, reaction objects, or per-step CMB and
scale-factor copies.
