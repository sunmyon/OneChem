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
