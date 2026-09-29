#include <chemistry/chemistry.h>

static void use_public_types(void)
{
  ChemistryConfig config = {CHEMISTRY_ABI_VERSION, NULL, 0.1};
  ChemistryComposition composition = chemistry_default_composition();
  ChemistryGlobalEnvironment global = {0};
  ChemistryStepContext step = {0};
  ChemistryCellEnvironment cell = {0};
  ChemistryAdvanceOptions options = {
      CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};
  ChemistryState state = {0};
  ChemistryResult result = {0};

  (void)config;
  (void)composition;
  (void)global;
  (void)step;
  (void)cell;
  (void)options;
  (void)state;
  (void)result;
}

int main(void)
{
  use_public_types();
  return 0;
}
