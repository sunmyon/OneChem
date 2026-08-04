#pragma once

#include "chemistry.h"

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

const ChemSpeciesDefinition *chemistry_species_definition(enum SPECIES species);
int chemistry_species_from_name(const char *name, enum SPECIES *species);
