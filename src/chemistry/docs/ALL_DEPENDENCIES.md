# Migration of `All` dependencies

`All` is removed by lifetime and responsibility, not replaced with a second
global bag.  The initial inventory from `5d30197` is grouped as follows.

| Destination | Current fields |
|---|---|
| run context | loaded tables, subcycle factor, optional fixed temperature |
| step context | `T_CMB`, `cf_atime`, `Time`, `ComovingIntegrationOn`, uniform `J21`/`G0` and photo-rates |
| cell environment | physical density, metallicity, and local radiation after adapter conversion |
| adapter only | unit factors, `HubbleParam`, `cf_a3inv`, time-bin conversion, particle scheduling and conserved-energy commit |
| remove/replace service | rank output, snapshot side effects, allocation, termination |

The migration keeps host-global reads outside the chemistry operation and
passes two longer-lived objects by pointer through the transitional cell
environment:

- `ChemRunContext`, configured once per run;
- `ChemGlobalStep`, prepared once whenever global chemistry conditions change.

`ChemEnvironment` contains only pointers to those objects plus cell-specific
values.  It does not copy CMB temperature, scale factor, background rates,
subcycle safety, or unit conversion for every `compute_abundances_d()` call.
Direct `All` reads have been removed from cooling, reaction-rate evaluation,
solver decisions, and the legacy compute wrapper. Shielding columns are
derived inside the chemistry evaluation rather than supplied by the host.

The legacy lifecycle storage now lives only in the standalone adapter. The core
driver has no lifecycle setter and no unit factor.  The final public API passes
`ChemistryContext` and `ChemistryStepContext` explicitly and has no module
singleton.

Compile-time features are not represented as runtime switches in
`ChemRunContext`.  For example, `fixed_temperature` exists only under
`FIX_TEMPERATURE`; HD-only fields remain guarded by `HD`.  The context carries
values required by the selected build, not features excluded from that build.
The fixed H/He/D reference composition is similarly not a run input: it defines
the normalization of `Abund[]` and is exposed as immutable chemistry ABI
metadata.  The former external `XH`, `XHe`, and `XHD` symbols no longer exist.

## Migrated source boundary

The following implementation group no longer reaches `allvars.h`, `proto.h`,
or `tags.h`, directly or through its chemistry headers:

- core driver;
- eager initialization orchestration and table-root resolution;
- aggregate radiative/chemical cooling, H2 and molecular-line cooling kernels,
  plus LW continuum opacity;
- species registry;
- reaction-network parser, forward-rate table, rate evaluation, and shielding;
- explicit/implicit chemistry and energy solver implementation, including Saha;
- thermodynamic closure and partition functions;
- chemistry profiling support.

`make core-dependency-check` asks the compiler for each module's complete
header dependency graph and rejects a host-header regression.  Transitional
`InitChem` no longer reads `All`, and background-rate preparation is a step
responsibility rather than initialization. Abundance floors belong to the
compiled chemistry configuration, not host mutation from `InitChem`.

Remaining host dependencies are concentrated inside the dust and continuum
table/cooling implementations, the standalone adapter, and the old
hydrodynamics interface.  Background FUV spectrum-to-rate conversion now lives
in the standalone host module above this directory.
