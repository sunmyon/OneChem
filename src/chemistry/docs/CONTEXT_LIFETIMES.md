# Context lifetimes and update contract

The central rule is that a value is updated at the widest lifetime at which it
is constant.  `All` is therefore decomposed; it is not renamed.

| Object | Update frequency | Representative contents |
|---|---|---|
| `ChemistryContext` / transitional `ChemRunContext` | once per run | loaded tables, subcycle factor |
| `ChemistryStepContext` / transitional `ChemGlobalStep` | once per host-wide chemistry sweep | CMB temperature, scale factor, simulation time, uniform radiation |
| `ChemistryCellEnvironment` / transitional `ChemEnvironment` | once per active cell | density, metallicity, local radiation |
| `ChemistryState` | evolved for one cell | abundances, CGS thermal energy, gamma |
| `ChemistryWorkspace` | private to one concurrent call | rates, cooling terms, solver scratch |

For AREPO, one global integration step normally prepares one immutable step
context and every active cell receives the same pointer.  Assigning CMB
temperature or scale factor for every particle is neither required nor part of
the public API.  The transitional `ChemEnvironment` currently stores two
pointers (`run` and `global`) while old internal function signatures are being
migrated; it does not contain copies of the shared values.

The standalone compatibility adapter uses this lifecycle:

```c
ChemRunContext run = {0};
ChemStandaloneUnits units = {0};
RadiationRateContext radiation_context;

run.subcycle_safety_factor = All.Epsilon;
units.specific_energy_cgs = All.UnitEnergy_in_cgs / All.UnitMass_in_g;
chemistry_standalone_configure(&run, &units); /* once after set_units() */
radiation_rate_context_init(&radiation_context); /* once per run/task */

ChemGlobalStep step = {0};
step.cmb_temperature = All.T_CMB;
step.scale_factor = All.cf_atime;
step.simulation_time = All.Time;
step.comoving_integration = All.ComovingIntegrationOn;
RadiationBlackbodyInput background = {All.J21, Teff};
chemistry_standalone_prepare_step_with_background(
    &step, &radiation_context, &background);  /* once per shared step */

for(each_active_cell)
  compute_abundances_d(mode, &input[cell], &output[cell], dt[cell]);

radiation_rate_context_destroy(&radiation_context); /* once at shutdown */
```

The convenience background call fills every uniform radiation input consumed
by chemistry: H-minus, H2, and HD destruction rates plus G0.  In an AREPO
adapter the radiation solver may provide those final rates directly instead;
chemistry does not require J21 or Teff to reconstruct them.

`ChemStandaloneUnits` is adapter-owned because the core consumes CGS. It is not a
member of `ChemRunContext` or the public chemistry state API.

Compile-time selection remains compile-time selection.  `ChemRunContext` does
not contain booleans that revive excluded physics.  For example,
`fixed_temperature` exists only when `FIX_TEMPERATURE` is compiled, and HD-only
step fields exist only when `HD` is compiled.

The standalone adapter keeps a private singleton only to preserve the legacy
argument list of `compute_abundances_d`. AREPO owns an explicit public
`ChemistryContext` per MPI process and does not use that singleton.
