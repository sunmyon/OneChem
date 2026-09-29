#include <chemistry/chemistry.h>
#include <string.h>

/* A real C caller checks value layouts; legacy-enabled builds also exercise
 * C linkage to status utilities implemented in C++. */
int chemistry_test_c_boundary(const ChemistryConfig *config,
                              const ChemistryGlobalEnvironment *global)
{
  return config->abi_version == CHEMISTRY_ABI_VERSION &&
      config->timescale_safety_factor > 0. &&
      global->cmb_temperature == 2.7255 && global->scale_factor == 1.
      ;
}
