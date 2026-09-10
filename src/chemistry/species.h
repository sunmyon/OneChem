#pragma once

#include <stddef.h>

#include "chemistry_core.h"

typedef enum {
  CHEM_ELEMENT_H,
  CHEM_ELEMENT_HE,
  CHEM_ELEMENT_D,
  CHEM_ELEMENT_C,
  CHEM_ELEMENT_O,
  CHEM_ELEMENT_COUNT
} ChemElement;

typedef struct {
  const char *name;
  signed char charge;
  unsigned char nuclei[CHEM_ELEMENT_COUNT];
} ChemSpeciesDefinition;

typedef enum {
  CHEM_INITIAL_NORMALIZATION_NONE,
  CHEM_INITIAL_NORMALIZATION_BOTH_DIRECTIONS,
  CHEM_INITIAL_NORMALIZATION_EXCESS_ONLY
} ChemInitialNormalization;

typedef struct {
  ChemElement element;
  double total_per_hydrogen;
  int metallicity_scaled;
  double initialization_upper_bound_factor;
  const enum SPECIES *reservoirs;
  int reservoir_count;
  int newton_reservoir_count;
  const enum SPECIES *random_partition_species;
  int random_partition_species_count;
  ChemInitialNormalization initial_normalization;
  enum SPECIES initial_fallback_reservoir;
} ChemElementConservationDefinition;

typedef struct {
  double floor;
  double mean_per_baryon;
  double upper_bound_per_hydrogen;
  unsigned char scale_with_ionization;
  unsigned char use_floor_as_guess;
} ChemSpeciesInitialization;

typedef enum {
  CHEM_RANDOM_INITIALIZE_ELEMENT_PARTITION,
  CHEM_RANDOM_INITIALIZE_LOG_SPECIES
} ChemRandomInitializationKind;

typedef struct {
  ChemRandomInitializationKind kind;
  ChemElement element;
  enum SPECIES species;
  double minimum;
  double maximum;
} ChemRandomInitialization;

extern const ChemElementConservationDefinition
    chemistry_element_conservation[];
extern const size_t chemistry_element_conservation_count;
extern const ChemSpeciesInitialization chemistry_species_initialization[N_SP];
extern const ChemRandomInitialization chemistry_random_initialization[];
extern const size_t chemistry_random_initialization_count;

const ChemSpeciesDefinition *chemistry_species_definition(enum SPECIES species);
int chemistry_species_from_name(const char *name, enum SPECIES *species);
