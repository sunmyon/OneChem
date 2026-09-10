#include <string.h>

#include "chemistry_constants.h"
#include "species.h"

#define CHEM_DEFINE_SPECIES(symbol, name_, charge_, h_, he_, d_, c_, o_) \
  [s##symbol] = {name_, charge_, {h_, he_, d_, c_, o_}},

static const ChemSpeciesDefinition species_definitions[N_SP] = {
  CHEM_SPECIES_LIST(CHEM_DEFINE_SPECIES)
};

#undef CHEM_DEFINE_SPECIES

#define CHEM_RESERVOIR_ENUM(symbol) s##symbol,
static const enum SPECIES hydrogen_reservoirs[] = {
  CHEM_HYDROGEN_RESERVOIRS(CHEM_RESERVOIR_ENUM)
};
static const enum SPECIES hydrogen_random_partition[] = {
  CHEM_HYDROGEN_RANDOM_PARTITION(CHEM_RESERVOIR_ENUM)
};
#ifdef HELIUM
static const enum SPECIES helium_reservoirs[] = {
  CHEM_HELIUM_RESERVOIRS(CHEM_RESERVOIR_ENUM)
};
static const enum SPECIES helium_random_partition[] = {
  CHEM_HELIUM_RANDOM_PARTITION(CHEM_RESERVOIR_ENUM)
};
#endif
#ifdef HD
static const enum SPECIES deuterium_reservoirs[] = {
  CHEM_DEUTERIUM_RESERVOIRS(CHEM_RESERVOIR_ENUM)
};
static const enum SPECIES deuterium_random_partition[] = {
  CHEM_DEUTERIUM_RANDOM_PARTITION(CHEM_RESERVOIR_ENUM)
};
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
static const enum SPECIES carbon_reservoirs[] = {
  CHEM_CARBON_RESERVOIRS(CHEM_RESERVOIR_ENUM)
};
static const enum SPECIES oxygen_reservoirs[] = {
  CHEM_OXYGEN_RESERVOIRS(CHEM_RESERVOIR_ENUM)
};
#endif
#undef CHEM_RESERVOIR_ENUM

#define CHEM_DEFINE_ELEMENT_CONSERVATION( \
    element_, total_, metallicity_, upper_factor_, reservoirs_, newton_count_, \
    random_species_, random_count_, normalization_, fallback_) \
  {element_, total_, metallicity_, upper_factor_, reservoirs_, \
   sizeof(reservoirs_) / sizeof(reservoirs_[0]), newton_count_, \
   random_species_, random_count_, \
   normalization_, fallback_},
const ChemElementConservationDefinition chemistry_element_conservation[] = {
  CHEM_ELEMENT_CONSERVATION_LIST(CHEM_DEFINE_ELEMENT_CONSERVATION)
};
#undef CHEM_DEFINE_ELEMENT_CONSERVATION

const size_t chemistry_element_conservation_count =
    sizeof(chemistry_element_conservation) /
    sizeof(chemistry_element_conservation[0]);

#define CHEM_DEFINE_SPECIES_INITIALIZATION( \
    symbol, floor_, mean_, upper_, ionization_, floor_guess_) \
  [s##symbol] = {floor_, mean_, upper_, ionization_, floor_guess_},
const ChemSpeciesInitialization chemistry_species_initialization[N_SP] = {
  CHEM_SPECIES_INITIALIZATION_LIST(CHEM_DEFINE_SPECIES_INITIALIZATION)
};
#undef CHEM_DEFINE_SPECIES_INITIALIZATION

#define CHEM_DEFINE_RANDOM_PARTITION(element_) \
  {CHEM_RANDOM_INITIALIZE_ELEMENT_PARTITION, element_, sEND, 0.0, 0.0},
#define CHEM_DEFINE_RANDOM_LOG(symbol, minimum_, maximum_) \
  {CHEM_RANDOM_INITIALIZE_LOG_SPECIES, CHEM_ELEMENT_COUNT, s##symbol, \
   minimum_, maximum_},
const ChemRandomInitialization chemistry_random_initialization[] = {
  CHEM_RANDOM_INITIALIZATION_LIST(CHEM_DEFINE_RANDOM_PARTITION,
                                  CHEM_DEFINE_RANDOM_LOG)
};
#undef CHEM_DEFINE_RANDOM_LOG
#undef CHEM_DEFINE_RANDOM_PARTITION

const size_t chemistry_random_initialization_count =
    sizeof(chemistry_random_initialization) /
    sizeof(chemistry_random_initialization[0]);

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
