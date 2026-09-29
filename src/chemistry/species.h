#pragma once

#include <stddef.h>
#include "chemistry_build_config.h"
#include "include/chemistry/chemistry_composition.h"


#define CHEM_MAX_RANDOM_PARTITION_SPECIES 4

#define CHEM_SPECIES_ELECTRON(X) \
  X(ELEC, "ELEC", -1, 0, 0, 0, 0, 0)

#define CHEM_SPECIES_ATOMIC_HYDROGEN(X) \
  X(HI,   "HI",    0, 1, 0, 0, 0, 0) \
  X(HII,  "HII",   1, 1, 0, 0, 0, 0)

#define CHEM_SPECIES_BASE(X) \
  CHEM_SPECIES_ELECTRON(X) \
  CHEM_SPECIES_ATOMIC_HYDROGEN(X)

#define CHEM_SPECIES_MOLECULAR_HYDROGEN(X) \
  X(H2I, "H2I", 0, 2, 0, 0, 0, 0)

#define CHEM_SPECIES_NEGATIVE_HYDROGEN(X) \
  X(HM, "HM", -1, 1, 0, 0, 0, 0)

#ifdef HELIUM
#define CHEM_SPECIES_HELIUM(X) \
  X(HeI,   "HeI",   0, 0, 1, 0, 0, 0) \
  X(HeII,  "HeII",  1, 0, 1, 0, 0, 0) \
  X(HeIII, "HeIII", 2, 0, 1, 0, 0, 0)
#else
#define CHEM_SPECIES_HELIUM(X)
#endif

#ifdef H2_PLUS
#define CHEM_SPECIES_H2_PLUS(X) X(H2II, "H2II", 1, 2, 0, 0, 0, 0)
#else
#define CHEM_SPECIES_H2_PLUS(X)
#endif

#ifdef HD
#define CHEM_SPECIES_DEUTERIUM_BASE(X) \
  X(DI,  "DI",  0, 0, 0, 1, 0, 0) \
  X(DII, "DII", 1, 0, 0, 1, 0, 0) \
  X(HDI, "HDI", 0, 1, 0, 1, 0, 0)
#ifdef HD_PLUS
#define CHEM_SPECIES_HD_PLUS(X) X(HDII, "HDII", 1, 1, 0, 1, 0, 0)
#else
#define CHEM_SPECIES_HD_PLUS(X)
#endif
#ifdef DMINUS
#define CHEM_SPECIES_D_MINUS(X) X(DM, "DM", -1, 0, 0, 1, 0, 0)
#else
#define CHEM_SPECIES_D_MINUS(X)
#endif
#define CHEM_SPECIES_DEUTERIUM(X) \
  CHEM_SPECIES_DEUTERIUM_BASE(X) \
  CHEM_SPECIES_HD_PLUS(X) \
  CHEM_SPECIES_D_MINUS(X)
#else
#define CHEM_SPECIES_DEUTERIUM(X)
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
#define CHEM_SPECIES_CARBON_HIGH_IONS(X) \
  X(CIII, "CIII", 2, 0, 0, 0, 1, 0) \
  X(CIV,  "CIV",  3, 0, 0, 0, 1, 0)
#define CHEM_SPECIES_OXYGEN_HIGH_IONS(X) \
  X(OII,  "OII",  1, 0, 0, 0, 0, 1) \
  X(OIII, "OIII", 2, 0, 0, 0, 0, 1) \
  X(OIV,  "OIV",  3, 0, 0, 0, 0, 1)
#define CHEM_SPECIES_OII_INTERMEDIATE(X)
#else
#define CHEM_SPECIES_CARBON_HIGH_IONS(X)
#define CHEM_SPECIES_OXYGEN_HIGH_IONS(X)
#define CHEM_SPECIES_OII_INTERMEDIATE(X) \
  X(OII,   "OII",   1, 0, 0, 0, 0, 1)
#endif

#define CHEM_SPECIES_CARBON_OXYGEN_STORED(X) \
  X(CI,   "CI",    0, 0, 0, 0, 1, 0) \
  X(CII,  "CII",   1, 0, 0, 0, 1, 0) \
  CHEM_SPECIES_CARBON_HIGH_IONS(X) \
  X(CO,   "CO",    0, 0, 0, 0, 1, 1) \
  X(CO2,  "CO2",   0, 0, 0, 0, 1, 2) \
  X(OI,   "OI",    0, 0, 0, 0, 0, 1) \
  CHEM_SPECIES_OXYGEN_HIGH_IONS(X) \
  X(OH,   "OH",    0, 1, 0, 0, 0, 1) \
  X(H2O,  "H2O",   0, 2, 0, 0, 0, 1) \
  X(O2,   "O2",    0, 0, 0, 0, 0, 2)

#define CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X) \
  X(COII,  "COII",  1, 0, 0, 0, 1, 1) \
  X(CH,    "CH",    0, 1, 0, 0, 1, 0) \
  X(CH2,   "CH2",   0, 2, 0, 0, 1, 0) \
  CHEM_SPECIES_OII_INTERMEDIATE(X) \
  X(OHII,  "OHII",  1, 1, 0, 0, 0, 1) \
  X(H2OII, "H2OII", 1, 2, 0, 0, 0, 1) \
  X(H3OII, "H3OII", 1, 3, 0, 0, 0, 1) \
  X(O2II,  "O2II",  1, 0, 0, 0, 0, 2)
#else
#define CHEM_SPECIES_CARBON_HIGH_IONS(X)
#define CHEM_SPECIES_OXYGEN_HIGH_IONS(X)
#define CHEM_SPECIES_OII_INTERMEDIATE(X)
#define CHEM_SPECIES_CARBON_OXYGEN_STORED(X)
#define CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X)
#endif

/* Fully solved intermediates are part of the caller-owned state.  Only the
 * reduced equilibrium formulation keeps them as internal algebraic species. */
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
#define CHEM_SPECIES_FULLY_SOLVED_INTERMEDIATES(X) \
  CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X)
#define CHEM_SPECIES_INTERNAL_LIST(X)
#else
#define CHEM_SPECIES_FULLY_SOLVED_INTERMEDIATES(X)
#define CHEM_SPECIES_INTERNAL_LIST(X) \
  CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X)
#endif

#define CHEM_SPECIES_STORED_LIST(X) \
  CHEM_SPECIES_BASE(X) \
  CHEM_SPECIES_HELIUM(X) \
  CHEM_SPECIES_MOLECULAR_HYDROGEN(X) \
  CHEM_SPECIES_H2_PLUS(X) \
  CHEM_SPECIES_NEGATIVE_HYDROGEN(X) \
  CHEM_SPECIES_DEUTERIUM(X) \
  CHEM_SPECIES_CARBON_OXYGEN_STORED(X) \
  CHEM_SPECIES_FULLY_SOLVED_INTERMEDIATES(X)

#define CHEM_SPECIES_LIST(X) \
  CHEM_SPECIES_STORED_LIST(X) \
  CHEM_SPECIES_INTERNAL_LIST(X)

/* Conservation reservoirs. These are choices for the dependent abundance,
 * not duplicate lists of every species carrying an element. Element content
 * is always obtained from the nuclei metadata above. */
#define CHEM_HYDROGEN_RESERVOIRS(X) \
  X(HI) \
  X(HII) \
  X(H2I)

#ifdef H2_PLUS
#define CHEM_HYDROGEN_RANDOM_H2_PLUS(X) X(H2II)
#else
#define CHEM_HYDROGEN_RANDOM_H2_PLUS(X)
#endif

#define CHEM_HYDROGEN_RANDOM_PARTITION(X) \
  CHEM_HYDROGEN_RESERVOIRS(X) \
  CHEM_HYDROGEN_RANDOM_H2_PLUS(X)

#ifdef HELIUM
#define CHEM_HELIUM_RESERVOIRS(X) \
  X(HeI) \
  X(HeII)
#define CHEM_HELIUM_RANDOM_PARTITION(X) \
  X(HeI) \
  X(HeII) \
  X(HeIII)
#define CHEM_ELEMENT_CONSERVATION_HELIUM(X) \
  X(CHEM_ELEMENT_HE, XHe / XH, 0, 0.0, helium_reservoirs, 1, \
    helium_random_partition, \
    sizeof(helium_random_partition) / sizeof(helium_random_partition[0]), \
    CHEM_INITIAL_NORMALIZATION_EXCESS_ONLY, sHeII)
/* Preserve caller-provided helium fractions.  Electron-scaled seeds for
 * vanishing HeII or HeIII inject positive charge and create electron spikes. */
#define CHEM_SPECIES_INITIALIZATION_HELIUM(X) \
  X(HeI,   0.0, 0.0, 0.09, 0, 1) \
  X(HeII,  0.0, 0.0, 0.09, 0, 1) \
  X(HeIII, 0.0, 0.0, 0.09, 0, 1)
#else
#define CHEM_HELIUM_RESERVOIRS(X)
#define CHEM_HELIUM_RANDOM_PARTITION(X)
#define CHEM_ELEMENT_CONSERVATION_HELIUM(X)
#define CHEM_SPECIES_INITIALIZATION_HELIUM(X)
#endif

#ifdef HD
#define CHEM_DEUTERIUM_RESERVOIRS(X) \
  X(DI) \
  X(DII) \
  X(HDI)
#define CHEM_DEUTERIUM_RANDOM_PARTITION(X) \
  X(DI) \
  X(DII) \
  X(HDI)
#define CHEM_ELEMENT_CONSERVATION_DEUTERIUM(X) \
  X(CHEM_ELEMENT_D, D_number_fraction / XH, 0, 0.0, \
    deuterium_reservoirs, \
    sizeof(deuterium_reservoirs) / sizeof(deuterium_reservoirs[0]), \
    deuterium_random_partition, \
    sizeof(deuterium_random_partition) / \
        sizeof(deuterium_random_partition[0]), \
    CHEM_INITIAL_NORMALIZATION_EXCESS_ONLY, sHDI)
#define CHEM_SPECIES_INITIALIZATION_DEUTERIUM_BASE(X) \
  X(DI,  1.e-25, 4.e-5, D_number_fraction / XH * 1.001, 0, 1) \
  X(DII, 1.e-80, 6.e-9, D_number_fraction / XH * 1.001, 1, 1) \
  X(HDI, 1.e-30, 1.e-8, D_number_fraction / XH * 1.001, 0, 0)
#else
#define CHEM_DEUTERIUM_RESERVOIRS(X)
#define CHEM_DEUTERIUM_RANDOM_PARTITION(X)
#define CHEM_ELEMENT_CONSERVATION_DEUTERIUM(X)
#define CHEM_SPECIES_INITIALIZATION_DEUTERIUM_BASE(X)
#endif

#ifdef HD_PLUS
#define CHEM_SPECIES_INITIALIZATION_HD_PLUS(X) \
  X(HDII, 1.e-30, 2.e-18, 0.5 * D_number_fraction / XH, 1, 0)
#else
#define CHEM_SPECIES_INITIALIZATION_HD_PLUS(X)
#endif

#ifdef DMINUS
#define CHEM_SPECIES_INITIALIZATION_D_MINUS(X) \
  X(DM, 1.e-30, 5.e-16, D_number_fraction / XH, 1, 0)
#else
#define CHEM_SPECIES_INITIALIZATION_D_MINUS(X)
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
#define CHEM_CARBON_HIGH_ION_RESERVOIRS(X) X(CIII)
#define CHEM_OXYGEN_HIGH_ION_RESERVOIRS(X) X(OIV)
#else
#define CHEM_CARBON_HIGH_ION_RESERVOIRS(X)
#define CHEM_OXYGEN_HIGH_ION_RESERVOIRS(X)
#endif
#define CHEM_CARBON_RESERVOIRS(X) \
  X(CI) \
  X(CII) \
  CHEM_CARBON_HIGH_ION_RESERVOIRS(X) \
  X(CO) \
  X(CO2)
#define CHEM_OXYGEN_RESERVOIRS(X) \
  X(OI) \
  CHEM_OXYGEN_HIGH_ION_RESERVOIRS(X) \
  X(OH) \
  X(H2O) \
  X(O2)
#define CHEM_ELEMENT_CONSERVATION_CARBON_OXYGEN(X) \
  X(CHEM_ELEMENT_C, A_C_sun, 1, 1.001, carbon_reservoirs, \
    sizeof(carbon_reservoirs) / sizeof(carbon_reservoirs[0]), \
    NULL, 0, \
    CHEM_INITIAL_NORMALIZATION_NONE, sCI) \
  X(CHEM_ELEMENT_O, A_O_sun, 1, 1.001, oxygen_reservoirs, \
    sizeof(oxygen_reservoirs) / sizeof(oxygen_reservoirs[0]), \
    NULL, 0, \
    CHEM_INITIAL_NORMALIZATION_NONE, sOI)
#define CHEM_SPECIES_INITIALIZATION_CARBON_OXYGEN(X) \
  X(CI,  1.e-16, 0.0, 0.0, 0, 1) \
  X(CII, 1.e-16, 0.0, 0.0, 0, 1) \
  CHEM_SPECIES_INITIALIZATION_CARBON_HIGH_IONS(X) \
  X(OI,  1.e-16, 0.0, 0.0, 0, 1)
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
#define CHEM_SPECIES_INITIALIZATION_CARBON_HIGH_IONS(X) \
  X(CIII, 1.e-30, 0.0, 0.0, 0, 1) \
  X(CIV,  1.e-30, 0.0, 0.0, 0, 1) \
  X(OII,  1.e-30, 0.0, 0.0, 0, 1) \
  X(OIII, 1.e-30, 0.0, 0.0, 0, 1) \
  X(OIV,  1.e-30, 0.0, 0.0, 0, 1)
#else
#define CHEM_SPECIES_INITIALIZATION_CARBON_HIGH_IONS(X)
#endif
#else
#define CHEM_CARBON_RESERVOIRS(X)
#define CHEM_OXYGEN_RESERVOIRS(X)
#define CHEM_ELEMENT_CONSERVATION_CARBON_OXYGEN(X)
#define CHEM_SPECIES_INITIALIZATION_CARBON_OXYGEN(X)
#define CHEM_SPECIES_INITIALIZATION_CARBON_HIGH_IONS(X)
#endif

#ifdef H2_PLUS
#define CHEM_SPECIES_INITIALIZATION_H2_PLUS(X) \
  X(H2II, 1.e-55, 2.e-14, 0.5, 1, 0)
#else
#define CHEM_SPECIES_INITIALIZATION_H2_PLUS(X)
#endif

#define CHEM_ELEMENT_CONSERVATION_LIST(X) \
  X(CHEM_ELEMENT_H, 1.0, 0, 0.0, hydrogen_reservoirs, \
    sizeof(hydrogen_reservoirs) / sizeof(hydrogen_reservoirs[0]), \
    hydrogen_random_partition, \
    sizeof(hydrogen_random_partition) / \
        sizeof(hydrogen_random_partition[0]), \
    CHEM_INITIAL_NORMALIZATION_BOTH_DIRECTIONS, sH2I) \
  CHEM_ELEMENT_CONSERVATION_HELIUM(X) \
  CHEM_ELEMENT_CONSERVATION_DEUTERIUM(X) \
  CHEM_ELEMENT_CONSERVATION_CARBON_OXYGEN(X)

#define CHEM_SPECIES_INITIALIZATION_LIST(X) \
  X(ELEC, 1.e-80, 1.e-15, 2.0, 0, 0) \
  X(HI,   1.e-20, 0.76,  1.001, 0, 1) \
  X(HII,  1.e-20, 1.e-15, 1.001, 0, 1) \
  CHEM_SPECIES_INITIALIZATION_HELIUM(X) \
  X(HM,  1.e-55, 9.e-12, 1.0, 1, 0) \
  X(H2I, 1.e-50, 6.e-6, 0.50005, 0, 0) \
  CHEM_SPECIES_INITIALIZATION_H2_PLUS(X) \
  CHEM_SPECIES_INITIALIZATION_DEUTERIUM_BASE(X) \
  CHEM_SPECIES_INITIALIZATION_HD_PLUS(X) \
  CHEM_SPECIES_INITIALIZATION_D_MINUS(X) \
  CHEM_SPECIES_INITIALIZATION_CARBON_OXYGEN(X)

#ifdef HELIUM
#define CHEM_RANDOM_INITIALIZATION_HELIUM(PARTITION, LOG) \
  PARTITION(CHEM_ELEMENT_HE)
#else
#define CHEM_RANDOM_INITIALIZATION_HELIUM(PARTITION, LOG)
#endif

#ifdef HD_PLUS
#define CHEM_RANDOM_INITIALIZATION_HD_PLUS(PARTITION, LOG) \
  LOG(HDII, 1.e-30, 0.01 * D_number_fraction / XH)
#else
#define CHEM_RANDOM_INITIALIZATION_HD_PLUS(PARTITION, LOG)
#endif

#ifdef DMINUS
#define CHEM_RANDOM_INITIALIZATION_D_MINUS(PARTITION, LOG) \
  LOG(DM, 1.e-30, 0.01 * D_number_fraction / XH)
#else
#define CHEM_RANDOM_INITIALIZATION_D_MINUS(PARTITION, LOG)
#endif

#ifdef HD
#define CHEM_RANDOM_INITIALIZATION_DEUTERIUM(PARTITION, LOG) \
  PARTITION(CHEM_ELEMENT_D) \
  CHEM_RANDOM_INITIALIZATION_HD_PLUS(PARTITION, LOG) \
  CHEM_RANDOM_INITIALIZATION_D_MINUS(PARTITION, LOG)
#else
#define CHEM_RANDOM_INITIALIZATION_DEUTERIUM(PARTITION, LOG)
#endif

/* Ordering is part of the seeded random-reset behavior. */
#define CHEM_RANDOM_INITIALIZATION_LIST(PARTITION, LOG) \
  PARTITION(CHEM_ELEMENT_H) \
  LOG(HM, 1.e-55, 1.e-5) \
  CHEM_RANDOM_INITIALIZATION_HELIUM(PARTITION, LOG) \
  CHEM_RANDOM_INITIALIZATION_DEUTERIUM(PARTITION, LOG)

#define H_number_fraction CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON
#define He_number_fraction CHEMISTRY_HELIUM_NUCLEI_PER_BARYON
#define D_number_fraction CHEMISTRY_DEUTERIUM_NUCLEI_PER_BARYON

/* Transitional internal spellings.  These are compile-time aliases, not
 * external variables; new public code should use the descriptive names. */
#define XH CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON
#define XHe CHEMISTRY_HELIUM_NUCLEI_PER_BARYON

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
extern const ChemSpeciesInitialization *const chemistry_species_initialization;
extern const ChemRandomInitialization chemistry_random_initialization[];
extern const size_t chemistry_random_initialization_count;

const ChemSpeciesDefinition *chemistry_species_definition(enum SPECIES species);
int chemistry_species_from_name(const char *name, enum SPECIES *species);
double chemistry_net_charge(const double *abundance);
double chemistry_neutral_electron_abundance(const double *abundance);
void chemistry_explicit_enforce_constraints(double metallicity,
                                            double *abundance);
