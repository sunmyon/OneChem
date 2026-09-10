#include <string.h>

#include "allvars.h"
#include "species.h"

#define CHEM_DEFINE_SPECIES(symbol, name_, charge_, h_, he_, d_, c_, o_) \
  [s##symbol] = {name_, charge_, {h_, he_, d_, c_, o_}},

static const ChemSpeciesDefinition species_definitions[N_SP] = {
  CHEM_SPECIES_LIST(CHEM_DEFINE_SPECIES)
};

#undef CHEM_DEFINE_SPECIES

const ChemSpeciesDefinition *chemistry_species_definition(enum SPECIES species)
{
  if(species < 0 || species >= sEND)
    return NULL;
  return &species_definitions[species];
}

int chemistry_species_from_name(const char *name, enum SPECIES *species)
{
  for(int candidate = 0; candidate < sEND; candidate++) {
    if(strcmp(name, species_definitions[candidate].name) == 0) {
      *species = candidate;
      return 1;
    }
  }
  return 0;
}
