#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "chemistry/chemistry.h"
#include "species.h"
#include "chemistry_constants.h"

static int near(double a, double b)
{
  return fabs(a - b) <= 1.e-13 * fmax(fabs(a), fabs(b));
}

static double total(const double *a, ChemElement element)
{
  double sum = 0.0;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    sum += chemistry_species_definition(static_cast<SPECIES>(j))->nuclei[element] * a[j];
  return sum;
}

static void check_budgets(const double *a, double z)
{
  const double targets[] = {XH, He_number_fraction, D_number_fraction,
                            XH * A_C_sun * z, XH * A_O_sun * z};
  for(int element = CHEM_ELEMENT_HE; element < CHEM_ELEMENT_COUNT; element++) {
    int present = 0;
    for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
      present |= chemistry_species_definition(static_cast<SPECIES>(j))->nuclei[element] != 0;
    if(present) assert(near(total(a, static_cast<ChemElement>(element)), targets[element]));
  }
  double charge = 0.0;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    if(j != sELEC) charge += chemistry_species_definition(static_cast<SPECIES>(j))->charge * a[j];
  assert(near(a[sELEC], charge));
}

int main(void)
{
  double a[TGCHEM_NUM_ABUNDANCES], before[TGCHEM_NUM_ABUNDANCES];
  /* Exercise every compiled stored species, including neutral intermediates,
   * multiply charged ions, O2 stoichiometry, and negative ions. */
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) a[j] = 1.e-8 * (j + 1);
  a[sHI] = 0.5;
  a[sHII] = 0.1;
  a[sELEC] = -1.0; /* dependent input is replaced, never used as an ion */
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 0.2) == CHEMISTRY_STATUS_OK);
  check_budgets(a, 0.2);
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) {
    const ChemSpeciesDefinition *sp = chemistry_species_definition(static_cast<SPECIES>(j));
    if(j != sELEC && !sp->nuclei[CHEM_ELEMENT_HE] && !sp->nuclei[CHEM_ELEMENT_D] &&
       !sp->nuclei[CHEM_ELEMENT_C] && !sp->nuclei[CHEM_ELEMENT_O])
      assert(a[j] == before[j]); /* H-only species are untouched */
    for(int k = 0; k < j; k++) {
      const ChemSpeciesDefinition *other = chemistry_species_definition(static_cast<SPECIES>(k));
      for(int element = CHEM_ELEMENT_HE; element < CHEM_ELEMENT_COUNT; element++) {
        if(!sp->nuclei[element] || !other->nuclei[element]) continue;
        if(element == CHEM_ELEMENT_O && (sp->nuclei[CHEM_ELEMENT_C] || other->nuclei[CHEM_ELEMENT_C])) continue;
        assert(near(a[j] / a[k], before[j] / before[k]));
      }
    }
  }
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 0.2) == CHEMISTRY_STATUS_OK);
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) assert(near(a[j], before[j]));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 0.0) == CHEMISTRY_STATUS_OK);
  check_budgets(a, 0.0);
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) {
    const ChemSpeciesDefinition *sp = chemistry_species_definition(static_cast<SPECIES>(j));
    if(sp->nuclei[CHEM_ELEMENT_C] || sp->nuclei[CHEM_ELEMENT_O]) assert(a[j] == 0.0);
  }
#ifdef CARBON_OXYGEN_CHEMISTRY
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 1.0) == CHEMISTRY_STATUS_INCOMPATIBLE_STATE);
  assert(memcmp(a, before, sizeof(a)) == 0);
  a[sCO2] = 1.e-4;
  a[sOI] = 1.e-4;
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 1.0) == CHEMISTRY_STATUS_INCOMPATIBLE_STATE);
  assert(memcmp(a, before, sizeof(a)) == 0);
  a[sCO2] = 0;
  a[sOI] = 0;
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
  /* All C/O in high ions: the former hard-coded lists failed this state. */
  a[sCIII] = 1.e-4;
  a[sOIII] = 1.e-4;
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 1.0) == CHEMISTRY_STATUS_OK);
  check_budgets(a, 1.0);
#endif
#endif
  a[sHM] = 2.0;
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 0.0) == CHEMISTRY_STATUS_INCOMPATIBLE_STATE);
  assert(memcmp(a, before, sizeof(a)) == 0);
  a[sHI] = -1.0;
  memcpy(before, a, sizeof(a));
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, 0.0) == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  assert(memcmp(a, before, sizeof(a)) == 0);
  assert(chemistry_normalize_abundances(a, 0, 1.0) == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  assert(chemistry_normalize_abundances(a, TGCHEM_NUM_ABUNDANCES, NAN) == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  puts("Metadata-based elemental normalization passed");
  return 0;
}
