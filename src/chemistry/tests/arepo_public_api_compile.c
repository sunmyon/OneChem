#include "chemistry_arepo.h"

static void use_arepo_public_types(void)
{
  ChemArepoConfig config = {0};
  ChemArepoContext context = {0};
  ChemArepoCell cell = {0};
  ChemArepoState state = {0};
  ChemArepoResult result = {0};

  (void)config;
  (void)context;
  (void)cell;
  (void)state;
  (void)result;
}

int main(void)
{
  use_arepo_public_types();
  return 0;
}
