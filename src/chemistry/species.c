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

double chemistry_net_charge(const double *abundance)
{
  double charge = 0.0;

  for(int species = 0; species < N_SP; species++)
    charge += species_definitions[species].charge * abundance[species];

  return charge;
}

double chemistry_neutral_electron_abundance(const double *abundance)
{
  double electron = 0.0;

  for(int species = 0; species < N_SP; species++)
    if(species != sELEC)
      electron += species_definitions[species].charge * abundance[species];

  return electron;
}

void chemistry_explicit_enforce_constraints(double metallicity,
                                            double *abundance)
{
  double element_total[CHEM_ELEMENT_COUNT] = {0.0};
  double element_content[CHEM_ELEMENT_COUNT] = {0.0};
  double element_scale[CHEM_ELEMENT_COUNT] = {1.0, 1.0, 1.0, 1.0, 1.0};

  for(size_t index = 0; index < chemistry_element_conservation_count; index++){
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[index];
    double total = XH * element->total_per_hydrogen;

    if(element->metallicity_scaled)
      total *= metallicity;

    element_total[element->element] = total;
  }

  /* Conservation is a property of the compiled network, not of the chosen
   * explicit update order. Include both evolved and algebraic intermediates. */
  for(int species = 0; species < N_SP; species++) {
    const ChemSpeciesDefinition *definition = &species_definitions[species];
    double species_abundance = abundance[species];

    if(definition->nuclei[CHEM_ELEMENT_H] > 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0 &&
       definition->nuclei[CHEM_ELEMENT_O] == 0)
      element_content[CHEM_ELEMENT_H] +=
          definition->nuclei[CHEM_ELEMENT_H] * species_abundance;

    element_content[CHEM_ELEMENT_HE] +=
        definition->nuclei[CHEM_ELEMENT_HE] * species_abundance;
    element_content[CHEM_ELEMENT_D] +=
        definition->nuclei[CHEM_ELEMENT_D] * species_abundance;
    element_content[CHEM_ELEMENT_C] +=
        definition->nuclei[CHEM_ELEMENT_C] * species_abundance;
  }

  for(int element = CHEM_ELEMENT_H; element <= CHEM_ELEMENT_C; element++)
    if(element_content[element] > 0.0)
      element_scale[element] =
          element_total[element] / element_content[element];

  for(int species = 0; species < N_SP; species++) {
    const ChemSpeciesDefinition *definition = &species_definitions[species];

    if(definition->nuclei[CHEM_ELEMENT_H] > 0 &&
       definition->nuclei[CHEM_ELEMENT_D] == 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0 &&
       definition->nuclei[CHEM_ELEMENT_O] == 0)
      abundance[species] *= element_scale[CHEM_ELEMENT_H];

    if(definition->nuclei[CHEM_ELEMENT_HE] > 0)
      abundance[species] *= element_scale[CHEM_ELEMENT_HE];

    if(definition->nuclei[CHEM_ELEMENT_D] > 0)
      abundance[species] *= element_scale[CHEM_ELEMENT_D];

    if(definition->nuclei[CHEM_ELEMENT_C] > 0)
      abundance[species] *= element_scale[CHEM_ELEMENT_C];
  }

  double oxygen_in_carbon = 0.0;
  double oxygen_only_content = 0.0;

  for(int species = 0; species < N_SP; species++) {
    const ChemSpeciesDefinition *definition = &species_definitions[species];
    double oxygen_content =
        definition->nuclei[CHEM_ELEMENT_O] * abundance[species];

    if(definition->nuclei[CHEM_ELEMENT_C] > 0)
      oxygen_in_carbon += oxygen_content;
    else
      oxygen_only_content += oxygen_content;
  }

  if(oxygen_only_content > 0.0 &&
     element_total[CHEM_ELEMENT_O] >= oxygen_in_carbon)
    element_scale[CHEM_ELEMENT_O] =
        (element_total[CHEM_ELEMENT_O] - oxygen_in_carbon) /
        oxygen_only_content;

  for(int species = 0; species < N_SP; species++) {
    const ChemSpeciesDefinition *definition = &species_definitions[species];

    if(definition->nuclei[CHEM_ELEMENT_O] > 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0)
      abundance[species] *= element_scale[CHEM_ELEMENT_O];
  }

  abundance[sELEC] = chemistry_neutral_electron_abundance(abundance);
}
