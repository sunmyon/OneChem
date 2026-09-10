# Host-independent chemistry API

External hydrodynamics codes include only `chemistry/chemistry.h`. They do not
include the internal solver headers or reproduce the chemistry data layout.

## Ownership and call order

```c
ChemistryConfig config;
ChemistryContext *context = NULL;

chemistry_config_defaults(&config);
/* Directory that contains the chemistry tables/ subdirectory. */
config.table_root = chemistry_table_directory;
config.subcycle_safety_factor = 0.01;
chemistry_context_create(&config, &context);

ChemistryGlobalEnvironment global = {0};
ChemistryStepContext step;
global.cmb_temperature = cmb_temperature_K;
global.scale_factor = scale_factor;
chemistry_step_prepare(context, &global, &step);

size_t bytes = chemistry_workspace_size(context);
/* Allocate `bytes` with normal C alignment once per concurrent caller. */
ChemistryWorkspace *workspace = NULL;
chemistry_workspace_init(context, storage, bytes, &workspace);

ChemistryCellEnvironment cell = {0};
cell.density_cgs = density_g_cm3;
cell.metallicity = metallicity_in_solar_units;
ChemistryState state = {abundances, chemistry_species_count(context),
                        specific_thermal_energy_cgs, gamma};
ChemistryResult result;
ChemistryAdvanceOptions options = {
    CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};

chemistry_advance(context, &step, workspace, &cell, duration_s, &options,
                  &state, &result);

chemistry_context_destroy(context);
```

The context owns its loaded network and table objects and is read-only after
creation. A step is prepared once for values shared by many cells. Each
simultaneous chemistry call requires a separate workspace. The host owns
abundance and workspace storage. A build with `CHEMISTRY_PROFILE` uses
process-wide profiling counters; concurrent profiled calls require external
synchronization.

## Data contract

- Density, time, specific thermal energy, rates, and heating terms use
  CGS units stated in the public structures.
- Abundances use `n_species / (rho / proton_mass)` and the feature-dependent
  ordering reported by `chemistry_species_name()`.
- `metallicity` is in solar units.
- `J21` describes the radiation field. Chemistry consumes the final `G0` and
  photo-rates supplied alongside it; it does not reconstruct local rates from
  local `J21`.
- The core computes its shielding columns from the local state; columns are not
  host inputs in the current public API.
- A non-OK solver return leaves `ChemistryState` unchanged. Diagnostic values
  may still be returned in `ChemistryResult`. Invalid arguments are rejected
  before either object is modified.

The HeI/HeII photoionization, He photoheating, and CO photodissociation fields
are reserved extension inputs. They default to zero and are intentionally not
connected to reactions until the host radiation model supplies them.

## Adapter boundary

An adapter performs only host-specific work: unit conversion, cell selection,
communication, and committing the accepted state. AREPO's `All`, `P`, `SphP`,
time-bin, and MPI logic therefore remain in the AREPO adapter. The numerical
chemistry, energy coupling, and subcycling are reached only through
`chemistry_advance()`.
