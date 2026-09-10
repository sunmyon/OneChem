# Mutable global-state migration

The target is zero mutable global state in chemistry core.  Immutable physical
constants, compile-time coefficient arrays, and stateless function tables may
remain module-local `static const`; they have no initialization lifetime or
ownership ambiguity.

## Ownership map

| Legacy state | Target owner | Status |
|---|---|---|
| reaction array, reaction count, species reaction indices | run context (`ChemReactionNetwork`) | migrated |
| temperature grid and forward/inverse rate arrays | run context (`ChemistryRateTables`) | migrated |
| atomic cooling coefficients tabulated on the common grid | run context (`ChemistryRateTables`) | migrated |
| CO/OH/H2O molecular cooling lookup data | run context (`ChemistryCoolingTables`) | migrated |
| dust opacity lookup data and interpolation metadata | run context (`ChemistryDustTables`) | migrated |
| continuum opacity lookup data | run context (`ChemistryContinuumOpacityTables`) | migrated |
| LW transport opacity lookup data | host radiation-transport module | migrated out of core |
| blackbody FUV/LW derived lookup data | sibling `RadiationRateContext` | migrated out of core |
| BH power-law one-value caches | host radiation model | pending |
| table root | context creation input, passed explicitly to each loader | migrated |
| `FlagPrint` and calculation diagnostics | per-evaluation workspace | migrated |
| profiling traces and counters | workspace or host reduction service | pending |
| CMB, scale factor, time, uniform radiation | immutable step context | migrated |
| cell density, metallicity, local radiation | cell environment | migrated |

## First migrated component: reaction network

`ChemReactionNetwork` is an opaque heap-owned object created during eager
initialization.  It owns the parsed reactions and both species-to-reaction
indices.  Network evaluation functions require an explicit network pointer;
there is no process-wide reaction array or `N_RE` count.

The transitional `ChemRunContext` owns the network pointer.  `InitChem` creates
it and `DestroyChem` releases it.  The final `ChemistryContext` will absorb this
ownership without changing the network evaluation API.

## Initialization path and diagnostics

The table root is no longer cached in a module-static buffer. `InitChem` passes
the caller-owned root explicitly through every loader, and each file is opened
while initialization is in progress. The root does not become runtime state.

Diagnostic enablement is a field of `ChemWorkspace`. Adapters initialize it
for each cell, and solver/cooling helper calls receive it explicitly where they
need to emit legacy traces. A failed evaluation may be repeated with diagnostics
enabled without mutating another cell or thread. The former `FlagPrint` and
`count_print` symbols do not exist in portable chemistry.

## Rate and atomic-cooling tables

`ChemistryRateTables` owns the temperature grid, forward and inverse reaction
rates, atomic cooling coefficients, and formation energies used by Saha and
inverse-rate construction.  It is allocated once by `InitChem`, filled eagerly
by `chem_tabulate`, referenced read-only by rate/cooling/Saha evaluation, and
released by `DestroyChem`.  None of these arrays has external linkage.

## Molecular cooling tables

`ChemistryCoolingTables` currently owns the CO, OH, and H2O molecular line
cooling data.  The fixed grid metadata is installed during creation, the three
external files are read eagerly into that object, and all interpolation
functions receive a `const ChemistryCoolingTables *`.  The former mutable
module-static `CO_table`, `OH_table`, and `H2O_table` objects no longer exist.

The object name is intentionally broader than the first migrated component;
dust, continuum, and metal cooling subobjects will join this run-lifetime
owner as their readers are migrated.

## Dust opacity tables

`ChemistryDustTables` owns the temperature grid, Rosseland and Planck opacity
arrays, and the derived logarithmic lookup metadata.  Dust balance and opacity
interpolation require an explicit `const ChemistryDustTables *`; the former
`n_Tdust`, `Kappa_*`, and lookup-coordinate globals have been removed.

Each rank loads its own object during explicit context initialization. The core
reader contains no MPI calls; any alternative distribution policy belongs to
the host adapter.

## Continuum opacity tables

`ChemistryContinuumOpacityTables` owns the gas-temperature and density grids,
the Rosseland and Planck opacity matrices, and their interpolation bounds.
Continuum-opacity evaluation receives an explicit read-only table pointer from
the run context; the former `Kappa_*`, grid-count, and bound globals have been
removed.

As with dust tables, every rank constructs and owns its own table object. The
core contains no collective operation or task-number dependency.

## FUV/LW blackbody tables

The blackbody spectrum, effective-temperature grid, H-minus photodetachment
table, and G0 table are no longer chemistry state. They are owned by the
sibling module's explicit `RadiationRateContext`, one directory above
`src/chemistry/`. The host initializes that context once and calls
`radiation_rate_uniform_background` only when J21 or effective temperature
changes. Chemistry receives only the resulting H-minus, H2, and HD destruction
rates and G0 through its step context.

The radiation-rate module has no mutable file-scope state, performs no I/O, and has
no MPI dependency.  Every MPI task may construct its own small deterministic
model; a host that instead broadcasts precomputed data remains free to do so.
It also exposes the undiluted blackbody values and H-minus cross section needed
by point-source transport. The portable chemistry directory contains no
spectrum-to-rate conversion or historical `get_Hminus_PD_rate` compatibility
symbol.

The separate `BH_RADIATION` power-law model is still pending.  It should be
implemented alongside the host radiation model rather than returned to the
chemistry context.

## LW transport opacity

The former `cooling/continuum_opacity_LW.c` was an AREPO ray-tracing component,
not a chemistry cooling kernel: it constructed opacity state consumed by the
radiation solver and depended on host cell temperature/state. It has been
removed from portable chemistry, and `InitChem` no longer initializes it.
AREPO's selected `RayTrace_*` module owns both
`tabulate_opacities_in_LW_band()` and `calc_continuum_opacities_LW()` and must
initialize that state as part of radiation-transport initialization.

Each MPI rank owns one network object.  The pointer is shared read-only by the
threads on that rank; it is never broadcast as a pointer.  MPI distribution of
initialization input remains an adapter responsibility.

## Migration rule

Each component moves as one coherent unit:

1. define its owning context subobject;
2. make initialization transactional and add destruction;
3. pass a `const` pointer through every evaluation path;
4. remove the old file-scope storage;
5. run dependency checks, build, and full numerical regression.

Creating a single module-global `ChemistryContext` is explicitly not an
acceptable intermediate endpoint because it merely renames the old global bag.
