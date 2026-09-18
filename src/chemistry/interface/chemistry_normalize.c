/* Adapted from sink-wind-feedback 874ab26672. All family membership,
 * stoichiometric weights and charges come from the shared species metadata. */
#include <float.h>
#include <math.h>
#include <string.h>

#include "include/chemistry/chemistry.h"
#include "species.h"
#include "chemistry_constants.h"

static ChemistryStatus normalize_family(double *a, ChemElement element,
                                         int exclude_carbon, double target)
{
  double total = 0.0;
  int present = 0;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) {
    const ChemSpeciesDefinition *sp = chemistry_species_definition(j);
    if(!sp->nuclei[element] || (exclude_carbon && sp->nuclei[CHEM_ELEMENT_C]))
      continue;
    present = 1;
    total += sp->nuclei[element] * a[j];
  }
  if(!present)
    return CHEMISTRY_STATUS_OK;
  if(!isfinite(total) || !isfinite(target))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  if(total == 0.0 && target > 0.0)
    return CHEMISTRY_STATUS_INCOMPATIBLE_STATE;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) {
    const ChemSpeciesDefinition *sp = chemistry_species_definition(j);
    if(sp->nuclei[element] && !(exclude_carbon && sp->nuclei[CHEM_ELEMENT_C]))
      a[j] = total > 0.0 ? (a[j] / total) * target : 0.0;
  }
  return CHEMISTRY_STATUS_OK;
}

ChemistryStatus chemistry_normalize_abundances(double *abundance, size_t count,
                                               double metallicity)
{
  if(!abundance || count != TGCHEM_NUM_ABUNDANCES ||
     !isfinite(metallicity) || metallicity < 0.0)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  double corrected[TGCHEM_NUM_ABUNDANCES];
  memcpy(corrected, abundance, sizeof(corrected));
  /* Electrons are dependent; replace them from ionic charge below. */
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    if(j != sELEC && (!isfinite(corrected[j]) || corrected[j] < 0.0))
      return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  ChemistryStatus status = normalize_family(corrected, CHEM_ELEMENT_HE, 0, He_number_fraction);
  if(status != CHEMISTRY_STATUS_OK) return status;
  status = normalize_family(corrected, CHEM_ELEMENT_D, 0, D_number_fraction);
  if(status != CHEMISTRY_STATUS_OK) return status;
  status = normalize_family(corrected, CHEM_ELEMENT_C, 0, XH * A_C_sun * metallicity);
  if(status != CHEMISTRY_STATUS_OK) return status;

  /* Carbon-bearing species keep their common scale. Subtract ALL their
   * oxygen before scaling oxygen species without carbon. */
  double oxygen_in_carbon = 0.0;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++) {
    const ChemSpeciesDefinition *sp = chemistry_species_definition(j);
    if(sp->nuclei[CHEM_ELEMENT_C])
      oxygen_in_carbon += sp->nuclei[CHEM_ELEMENT_O] * corrected[j];
  }
  const double oxygen_target = XH * A_O_sun * metallicity;
  double oxygen_remaining = oxygen_target - oxygen_in_carbon;
  if(!isfinite(oxygen_remaining)) return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  if(oxygen_remaining < -32.0 * DBL_EPSILON * oxygen_target)
    return CHEMISTRY_STATUS_INCOMPATIBLE_STATE;
  /* Retain upstream's roundoff-only tolerance, not a species floor. */
  if(oxygen_remaining < 0.0) oxygen_remaining = 0.0;
  status = normalize_family(corrected, CHEM_ELEMENT_O, 1, oxygen_remaining);
  if(status != CHEMISTRY_STATUS_OK) return status;

  double net_charge = 0.0;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    if(j != sELEC)
      net_charge += chemistry_species_definition(j)->charge * corrected[j];
  if(!isfinite(net_charge) || net_charge < 0.0)
    return CHEMISTRY_STATUS_INCOMPATIBLE_STATE;
  corrected[sELEC] = net_charge;
  memcpy(abundance, corrected, sizeof(corrected));
  return CHEMISTRY_STATUS_OK;
}
