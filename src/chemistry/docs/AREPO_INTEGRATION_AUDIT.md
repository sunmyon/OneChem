# AREPO integration audit

> Historical migration notes. Current C++ ownership and calling contracts are in
> [/Users/j/Desktop/codex/1ZONE/src/chemistry/README.md](/Users/j/Desktop/codex/1ZONE/src/chemistry/README.md).
> References below to the transitional C context or InitChem do not describe the current API.

This audit was performed against AREPO commit `984ccbfbf1` on branch
`sink-wind-feedback`.  The original AREPO checkout was kept unchanged; all
compile experiments used an archive copy under `/tmp`.

## Base chemistry connection

With radiation migration considered separately, portable chemistry needs only
the following AREPO touch points:

1. a `CHEMISTRY_PORTABLE` option in the run `Config.sh`;
2. a conditional include of `chemistry/mk/arepo.mk` from
   `makefiles/modules.make`, replacing the old `chemistry_SC` object list;
3. a conditional public-header selection in `src/allvars.h`;
4. replacement of `InitChem()` by `chemistry_arepo_host_init()` in
   `src/begrun.c`; and
5. `chemistry_arepo_host_destroy()` before `MPI_Finalize()` in `src/main.c`.

No new fields or parameters are required: this AREPO tree already contains
`All.ChemistryTableDir`, `All.Epsilon`, and their parameter registration.
The portable make fragment applies its generated configuration header and
private include directories only to portable chemistry objects.  It does not
modify the flags of unrelated AREPO objects.

The portable core, both AREPO adapter objects, `chem_balance`, and the modified
`allvars`, `begrun`, and `main` objects compile with the active AREPO generated
configuration.  A complete executable is not yet a useful integration gate:
the archived baseline configuration already fails in unrelated sink/metal
sources because fields referenced by those sources are absent from the active
particle structures.

## Radiation is a separate migration

The old `chemistry_SC/chemistry.h` is also serving as an accidental radiation
header.  All three `RayTrace_RSPH*` implementations include it directly.  The
currently selected `RayTrace_RSPH_START` variant has seven such source files.
They use it for species indices, numerical constants, `compute_temp()`, and the
H-minus spectrum interface; they do not require the old chemistry solver.

The symbols below must therefore not be silently retained in the portable
chemistry core:

- `Hminus_detachment_John1988()`;
- `make_rate_tables_LW_radiation()` and `get_Hminus_PD_rate()`;
- `update_background_J21()`;
- `init_J21_table()` and the time-dependent J21 table, when enabled; and
- the RayTrace-owned `tabulate_opacities_in_LW_band()` initialization call.

Their owner is the sibling `radiation/` module. It holds an explicit
rank-lifetime `RadiationRateContext` and produces the final rates
placed in `ChemGlobalStep` or `ChemArepoCell`.  Chemistry consumes those rates;
it does not know J21, effective-temperature spectra, or optical-depth transport.

That module is now implemented and regression-tested independently of AREPO.
It additionally exposes undiluted blackbody quantities and the H-minus cross
section required by RayTrace, without depending on chemistry headers.

## Two-stage migration

For the first integration build, one conditional legacy-header shim is the
smallest safe bridge.  Under `CHEMISTRY_PORTABLE`, the old
`chemistry_SC/chemistry.h` should expose only portable species definitions and
the constants still shared by RayTrace.  This allows the selected seven
RayTrace sources to compile unchanged.  The shim was tested in the temporary
tree and all seven objects compiled.

The next stage replaces those direct legacy includes with
`radiation/radiation_rates.h` plus the portable species header and removes
the shim. The generic spectrum-to-rate calculation is shared, while standalone
and AREPO retain their own context ownership and update policy. Keeping these
as two stages avoids combining chemistry validation with a broad RayTrace
source cleanup.

## Validation order

1. compile and test the portable core without AREPO;
2. compile all portable objects with the actual AREPO feature macros;
3. compile the AREPO lifecycle and `chem_balance` boundary objects;
4. compile the selected RayTrace variant through the temporary shim;
5. add and test the explicit radiation model and its initialization;
6. run a single-cell state/rate comparison against the old implementation;
7. only then attempt a complete AREPO executable and runtime comparison.
