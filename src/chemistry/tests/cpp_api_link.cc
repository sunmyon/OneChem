// Include function-bearing headers independently before the umbrella header:
// linkage must follow the compatibility option regardless of include order.
#include <chemistry/chemistry_config.h>
#include <chemistry/chemistry_status.h>
#include <chemistry/chemistry.h>

#include "chemistry_model.hpp"
#include "species.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <type_traits>
#ifdef CHEMISTRY_STANDALONE_HOST
#include "chemistry_standalone_adapter.h"
#include "standalone_chemistry_compat.h"
#include "chemistry_standalone_types.h"
#endif
#ifdef CHEMISTRY_PROFILE
#include <thread>
#endif

extern "C" int chemistry_test_c_boundary(
    const ChemistryConfig *, const ChemistryGlobalEnvironment *);


static_assert(std::is_standard_layout<ChemistryConfig>::value, "C ABI");
static_assert(std::is_standard_layout<ChemistryState>::value, "C ABI");
static_assert(std::is_trivially_copyable<ChemistryState>::value, "C ABI");

#define CHECK(condition) do { \
  if (!(condition)) { \
    std::fprintf(stderr, "C/C++ API check failed: %s\n", #condition); \
    return 1; \
  } \
} while (0)

int main(int argc, char **argv)
{
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_set_failure(17);
  int worker_initial = -1, worker_final = -1;
  std::thread worker([&] {
    worker_initial = chemistry_profile_last_failure();
    chemistry_profile_set_failure(29);
    worker_final = chemistry_profile_last_failure();
  });
  worker.join();
  CHECK(worker_initial == 0 && worker_final == 29);
  CHECK(chemistry_profile_last_failure() == 17);
  chemistry_profile_set_failure(0);
#endif
  ChemistryConfig config{};
  chemistry_config_defaults(&config);
  config.table_root = argc > 1 ? argv[1] : "src/chemistry";
  ChemistryGlobalEnvironment global{};
  global.cmb_temperature = 2.7255;
  global.scale_factor = 1.;
  CHECK(chemistry_test_c_boundary(&config, &global));
  CHECK(std::strcmp(chemistry_status_string(CHEMISTRY_STATUS_IO_ERROR),
                    "I/O error") == 0);

  auto model = chemistry::Model::create(config);
  CHECK(model.status == CHEMISTRY_STATUS_OK && model.model);
#ifdef CHEMISTRY_STANDALONE_HOST
  ChemistryStandaloneAdapter first_adapter, second_adapter;
  ChemStandaloneUnits units{1.0};
  chemdata_in host_input{};
  chemdata_out host_output{};
  host_output.Utherm = 123.0;
  double full_guess[N_SP]{};
  double equilibrium[N_SP];
  for(double &value : equilibrium) value = -7.0;
  auto rejects_unprepared = [&](ChemistryStandaloneAdapter &adapter) {
    return adapter.advance(1, &host_input, &host_output, 1.0) == 1 &&
        adapter.solve_fixed(&host_input, &host_output, 100.0, 1.0) == 0 &&
        adapter.solve_fixed_guess(&host_input, &host_output, 100.0, 1.0,
                                  full_guess, 10) == 0 &&
        adapter.equilibrium_abundances(&host_input, &host_output, 100.0, 1.0,
                                      equilibrium) == 0 &&
        host_output.Utherm == 123.0 && equilibrium[0] == -7.0;
  };
  CHECK(rejects_unprepared(first_adapter));
  first_adapter.configure(&model.model->run(), &units);
  CHECK(first_adapter.set_timescale_safety_factor(0.25));
  CHECK(rejects_unprepared(first_adapter));
  ChemGlobalStep host_step{};
  first_adapter.prepare_step(&host_step);
  first_adapter.configure(&model.model->run(), &units);
  CHECK(rejects_unprepared(first_adapter));
  CHECK(!second_adapter.set_timescale_safety_factor(0.5));
  second_adapter.configure(&model.model->run(), &units);
  first_adapter.configure(nullptr, nullptr);
  CHECK(!first_adapter.set_timescale_safety_factor(0.25));
  CHECK(rejects_unprepared(first_adapter));
  CHECK(second_adapter.set_timescale_safety_factor(0.5));
#endif
  ChemWorkspace workspace{};
  ChemistryStepContext step{};
  CHECK(chemistry::Driver::prepare(global, step) == CHEMISTRY_STATUS_OK);
  ChemistryCellEnvironment cell{};
  cell.density_cgs = 1.e-21;
  double abundances[TGCHEM_NUM_ABUNDANCES]{};
  abundances[sHI] = XH;
  ChemistryState state{};
  state.abundances = abundances;
  state.abundance_count = TGCHEM_NUM_ABUNDANCES;
  state.specific_thermal_energy_cgs = 1.e12;
  state.gamma = 5. / 3.;
  ChemistryAdvanceOptions options{CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};
  ChemistryResult result{};
  const chemistry::Driver driver;
  CHECK(driver.advance(*model.model, step, workspace, cell, 0., options, state, result)
        == CHEMISTRY_STATUS_OK);
  CHECK(std::isfinite(state.specific_thermal_energy_cgs));
  cell.density_cgs = -1.;
  const double saved_energy = state.specific_thermal_energy_cgs;
  CHECK(driver.advance(*model.model, step, workspace, cell, 1., options, state, result)
        == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  CHECK(state.specific_thermal_energy_cgs == saved_energy);
  std::puts("PASS: native C++ Model / typed workspace / Driver interface");
  return 0;
}
