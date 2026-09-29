#include "chemistry_arepo.h"

static void use_arepo_public_types(void)
{
  ChemistryConfig config = {0};
  ChemistryCellEnvironment cell = {0};
  ChemistryState state = {0};
  ChemistryResult result = {0};

  (void)config;
  (void)cell;
  (void)state;
  (void)result;
}

int main(void)
{
  use_arepo_public_types();
  return 0;
}
