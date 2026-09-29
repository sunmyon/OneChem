# Native C++ interface

Include `/Users/j/Desktop/codex/1ZONE/src/chemistry/chemistry_model.hpp` for Model,
Driver and workspace declarations. The host and core must use matching species
and physics compile options. The implementation requires C++17.

```cpp
ChemistryConfig config{};
chemistry_config_defaults(&config);
config.table_root = table_root; // directory containing tables/
auto created = chemistry::Model::create(config);
if(created.status != CHEMISTRY_STATUS_OK) return created.status;

ChemistryGlobalEnvironment global{};
global.scale_factor = scale_factor;
global.cmb_temperature = cmb_temperature;
ChemistryStepContext step{};
auto status = chemistry::Driver::prepare(global, step);
if(status != CHEMISTRY_STATUS_OK) return status;

ChemWorkspace workspace{};
ChemistryAdvanceOptions options{CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};
ChemistryResult result{};
status = chemistry::Driver{}.advance(*created.model, step, workspace, cell,
                                    duration_s, options, state, result);
```

The host initializes cell and state before advance. State.abundances points to
an array of TGCHEM_NUM_ABUNDANCES elements. Config.table_root is used only during
creation. Model::Handle releases all resources on scope exit; no teardown API or
raw workspace allocation is needed.

## Contracts

- Density, duration, energy, heating and rates use CGS units stated in the value
  structures. Abundances use n_species / (rho / proton_mass).
- Species metadata is available through chemistry_species_definition; composition
  is immutable compile-time metadata. Metallicity is in solar units.
- A shared step supplies CMB temperature, scale factor, time and uniform radiation.
  Cell inputs supply local conditions and final photo-rates/heating rates.
- Chemistry computes shielding from its state. The host supplies radiation rates;
  chemistry does not reconstruct them from local J21.
- The host-value Driver overload validates/imports into temporary state. A solver
  failure leaves caller abundances, energy and gamma unchanged. Diagnostics may
  be returned even on failure. Invalid arguments are rejected before outputs change.
- Each simultaneous call needs a separate ChemWorkspace. Model data is read-only.
  Profiling mode uses process-wide counters and needs external synchronization.
- All native solver entry points return ChemistryStatus with OK meaning success.
  The internal-state Driver overload may mutate its state before a failure;
  callers needing transactional behavior use the host-value overload.

The old opaque ChemistryContext/ChemistryWorkspace, InitChem/DestroyChem and C
advance/evaluator wrappers have been removed. There is no compatibility flag.
The C-readable value headers do not imply that core functions have a C ABI.
Host C hooks and foreign libraries may still use C linkage at adapter boundaries.

HeI/HeII photo-rates/heating and CO dissociation fields remain subject to the
compiled reaction network; supplying a field does not enable absent physics.
