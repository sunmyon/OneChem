#ifndef CHEMISTRY_COMPOSITION_H
#define CHEMISTRY_COMPOSITION_H

/*
 * Abundances are particle counts per baryon mass unit:
 *
 *   abundance[i] = n_i / (rho / proton_mass).
 *
 * These values define the abundance ABI of this chemistry network.  They are
 * chemistry-owned constants, not host hydrodynamics inputs.
 */
#define CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON 0.76
#define CHEMISTRY_HELIUM_NUCLEI_PER_BARYON 0.0625
#define CHEMISTRY_DEUTERIUM_NUCLEI_PER_BARYON 2.25e-5
#define CHEMISTRY_DEUTERIUM_TO_HYDROGEN_RATIO 3.e-5

typedef struct ChemistryComposition {
  double hydrogen_nuclei_per_baryon;
  double helium_nuclei_per_baryon;
  double deuterium_nuclei_per_baryon;
} ChemistryComposition;

static inline ChemistryComposition chemistry_default_composition(void)
{
  const ChemistryComposition composition = {
      CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON,
      CHEMISTRY_HELIUM_NUCLEI_PER_BARYON,
      CHEMISTRY_DEUTERIUM_NUCLEI_PER_BARYON};
  return composition;
}

#endif
