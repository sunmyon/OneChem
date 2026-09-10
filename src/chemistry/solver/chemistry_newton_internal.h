#pragma once

#include "chemistry_core.h"
#include "chemistry_types.h"

#ifdef CARBON_OXYGEN_CHEMISTRY
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 6
#else
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 4
#endif

typedef struct {
  int count;
  double coefficient[CHEM_MAX_CONSERVATION_CONSTRAINTS][N_SP];
  double total[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  enum SPECIES dependent_species[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  double dependent_inverse[
      CHEM_MAX_CONSERVATION_CONSTRAINTS *
      CHEM_MAX_CONSERVATION_CONSTRAINTS];
  unsigned char is_dependent[N_SP];
} ChemConservationSystem;

/* State of one implicit Newton system.  Keeping the conservation chart and
 * its reduced species layout here ensures that they are always updated
 * together. */
typedef struct {
  double abundance[N_SP];
  double initial_abundance[N_SP];
  double cooling_initial_abundance[N_SP];
  double initialization_floor[N_SP];
  double newton_lower_bound[N_SP];
  double upper_bound[N_SP];
  double total_step[N_SP];
  enum SPECIES solver_species[N_SP];
  int solver_species_count;
  ChemConservationSystem conservation;
} ChemNewtonSystem;

int chemistry_invert_dependent_matrix(ChemConservationSystem *conservation);
int chemistry_build_conservation_system(
    ChemConservationSystem *conservation, double metallicity,
    const double *abundance);
int chemistry_build_coupled_conservation_system(
    ChemConservationSystem *conservation, double metallicity,
    const double *abundance);
int chemistry_species_is_algebraic(
    enum SPECIES species, const ChemConservationSystem *conservation);
int chemistry_update_conservation_dependent_species(
    double *abundance, const ChemConservationSystem *conservation);
int chemistry_update_implicit_algebraic_species(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation);
int chemistry_species_is_intermediate_equilibrium(enum SPECIES species);
int chemistry_update_intermediate_equilibrium_species(
    unsigned char mode, double nH, double *abundance, const ChemRates *rates);
int chemistry_reduce_jacobian_for_algebraic_species(
    double nH, const double *abundance, const ChemRates *rates,
    const enum SPECIES *solver_species, int n_solver_species,
    const ChemConservationSystem *conservation, double drfdy[N_SP][N_SP]);

#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
int chemistry_is_full_implicit_initialized_intermediate(
    enum SPECIES species);
#endif
