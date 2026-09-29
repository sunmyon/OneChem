#include "chemistry_model.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "chemistry_build_config.h"
#include "forward_rates.h"
#include "reaction_network.h"
#include "species.h"

FILE *open_chemistry_table(const char *table_root, const char *filename)
{
  char path[1024];
  const int length = snprintf(path, sizeof(path), "%s/tables/%s", table_root,
                              filename);
  if(length < 0 || (size_t) length >= sizeof(path))
    return NULL;
  return fopen(path, "r");
}

int main(void)
{
  const ChemSpeciesDefinition *civ = chemistry_species_definition(sCIV);
  const ChemSpeciesDefinition *oiv = chemistry_species_definition(sOIV);
  assert(civ != NULL && civ->charge == 3 && civ->nuclei[CHEM_ELEMENT_C] == 1);
  assert(oiv != NULL && oiv->charge == 3 && oiv->nuclei[CHEM_ELEMENT_O] == 1);

  for(int rate = 141; rate <= 170; rate++) {
    if(rate >= 153 && rate <= 158)
      continue;
    if(rate == 147 || rate == 162 || rate == 168)
      continue;
    assert(forward_rate_table[rate] != NULL);
    const double value = forward_rate_table[rate](1.0e4);
    assert(isfinite(value) && value >= 0.0);
  }

  auto owner = chemistry::load_network(".");
  auto *network = owner.get();
  assert(network);

  double abundance[N_SP] = {0.0};
  double forward[CHEM_RATE_COUNT] = {0.0};
  double inverse[CHEM_RATE_COUNT] = {0.0};
  double rhs[N_SP] = {0.0};
  abundance[sCI] = 2.0e-4;
  forward[153] = 3.0;
  chemistry_network_rhs(network, CHEM_REACTION_EXPLICIT, 1.0, abundance,
                        forward, inverse, rhs);
  assert(fabs(rhs[sCI] + 6.0e-4) < 1.0e-16);
  assert(fabs(rhs[sCII] - 6.0e-4) < 1.0e-16);
  assert(fabs(rhs[sELEC] - 6.0e-4) < 1.0e-16);

  memset(rhs, 0, sizeof(rhs));
  memset(abundance, 0, sizeof(abundance));
  memset(forward, 0, sizeof(forward));
  abundance[sCO] = 3.0e-5;
  forward[171] = 4.0;
  chemistry_network_rhs(network, CHEM_REACTION_EXPLICIT, 1.0, abundance,
                        forward, inverse, rhs);
  assert(fabs(rhs[sCO] + 1.2e-4) < 1.0e-16);
  assert(fabs(rhs[sCI] - 1.2e-4) < 1.0e-16);
  assert(fabs(rhs[sOI] - 1.2e-4) < 1.0e-16);
  assert(rhs[sELEC] == 0.0);


  puts("metal ionization network check passed");
  return 0;
}
