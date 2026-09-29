# Context lifetimes and update contract

The central rule is that a value is updated at the widest lifetime at which it
is constant.  `All` is therefore decomposed; it is not renamed.

| Object | Update frequency | Representative contents |
|---|---|---|
| `chemistry::Model` | once per run | loaded tables, subcycle factor |
| `ChemistryStepContext` / transitional `ChemGlobalStep` | once per host-wide chemistry sweep | CMB temperature, scale factor, simulation time, uniform radiation |
| `ChemistryCellEnvironment` / transitional `ChemEnvironment` | once per active cell | density, metallicity, local radiation |
| `ChemistryState` | evolved for one cell | abundances, CGS thermal energy, gamma |
| `ChemWorkspace` | private to one concurrent call | rates, cooling terms, solver scratch |

For AREPO, one global integration step normally prepares one immutable step
context and every active cell receives the same pointer.  Assigning CMB
temperature or scale factor for every particle is neither required nor part of
the public API.  The transitional `ChemEnvironment` currently stores two
pointers (`run` and `global`) while old internal function signatures are being
migrated; it does not contain copies of the shared values.

A `Model::Handle` owns a non-copyable, non-movable `chemistry::Model`.
The handle can move while resource addresses remain stable. Model directly owns
its resource objects; evaluations borrow const pointers through `ChemRunView`.
All borrowed views must expire before Model destruction. Resource factories return
owning handles, so failed construction automatically releases partial resources.

The standalone compatibility adapter uses this lifecycle (error handling omitted):

```cpp
ChemistryConfig config;
chemistry_config_defaults(&config);
config.table_root = All.ChemistryTableDir;
config.timescale_safety_factor = All.Epsilon;
auto created = chemistry::Model::create(config); // check status before using model
ChemStandaloneUnits units = {0};
RadiationRateContext radiation_context;

units.specific_energy_cgs = All.UnitEnergy_in_cgs / All.UnitMass_in_g;
chemistry_standalone_configure(&created.model->run(), &units);
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

chemistry_standalone_configure(NULL, NULL); /* drop borrowed pointers */
created.model.reset(); // normal scope exit also destroys it
radiation_rate_context_destroy(&radiation_context); /* once at shutdown */
```

The convenience background call fills every uniform radiation input consumed
by chemistry: H-minus, H2, and HD destruction rates plus G0.  In an AREPO
adapter the radiation solver may provide those final rates directly instead;
chemistry does not require J21 or Teff to reconstruct them.

`ChemStandaloneUnits` is adapter-owned because the core consumes CGS. It is not a
member of `ChemRunView` or the public chemistry state API.

Compile-time selection remains compile-time selection.  `ChemRunView` does
not contain booleans that revive excluded physics.  For example,
`fixed_temperature` exists only when `FIX_TEMPERATURE` is compiled, and HD-only
step fields exist only when `HD` is compiled.

The standalone adapter keeps a private singleton only to preserve the legacy
argument list of `compute_abundances_d`. AREPO owns an explicit
`Model::Handle` per MPI process and does not use that singleton.
