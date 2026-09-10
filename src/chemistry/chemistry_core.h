#pragma once

#include "chemistry_build_config.h"
#include "species_definitions.h"
#include "include/chemistry/chemistry_composition.h"

#ifdef CHEMISTRY_SC

#define H_number_fraction CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON
#define He_number_fraction CHEMISTRY_HELIUM_NUCLEI_PER_BARYON
#define D_number_fraction CHEMISTRY_DEUTERIUM_NUCLEI_PER_BARYON

/* Transitional internal spellings.  These are compile-time aliases, not
 * external variables; new public code should use the descriptive names. */
#define XH CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON
#define XHe CHEMISTRY_HELIUM_NUCLEI_PER_BARYON

#ifndef CLIGHT
#define CLIGHT CLIGHT_REAL
#endif
#ifndef CLIGHT_REAL
#define CLIGHT_REAL 2.99792458e10
#endif

#define CHEM_DECLARE_SPECIES(symbol, name, charge, h, he, d, c, o) s##symbol,
enum SPECIES {
  CHEM_SPECIES_LIST(CHEM_DECLARE_SPECIES)
  sEND
};
#undef CHEM_DECLARE_SPECIES

#define CHEM_COUNT_STORED_SPECIES(symbol, name, charge, h, he, d, c, o) + 1
enum {
  sSTORED_END = 0 CHEM_SPECIES_STORED_LIST(CHEM_COUNT_STORED_SPECIES)
};
#undef CHEM_COUNT_STORED_SPECIES

#define N_SP sEND

/* AREPO's legacy chems.h may reserve the storage array before including this
 * header. Preserve that value long enough to verify the host allocation. */
#ifdef TGCHEM_NUM_ABUNDANCES
enum { CHEMISTRY_HOST_ABUNDANCE_CAPACITY = TGCHEM_NUM_ABUNDANCES };
#undef TGCHEM_NUM_ABUNDANCES
#define TGCHEM_NUM_ABUNDANCES sSTORED_END
typedef char chemistry_host_abundance_capacity_must_match[
    CHEMISTRY_HOST_ABUNDANCE_CAPACITY == TGCHEM_NUM_ABUNDANCES ? 1 : -1];
#else
#define TGCHEM_NUM_ABUNDANCES sSTORED_END
#endif

double compute_gamma_from_abund(const double *abund, double temperature);
double chemistry_gamma_from_abundances(const double *abund,
                                       double temperature);
double chemistry_molecular_weight(const double *abund);
double chemistry_temperature(const double *abund, double gamma,
                             double utherm_cgs);
double chemistry_jeans_radius(double rho_cgs, double temperature,
                              double thermal_mass_in_proton_masses);

#endif
