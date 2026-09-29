#pragma once

#include "chemistry_core.h"
#include "chemistry_types.h"

#ifdef CARBON_OXYGEN_CHEMISTRY
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 6
#else
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 4
#endif

typedef struct ChemConservationSystem {
  int count;
  double coefficient[CHEM_MAX_CONSERVATION_CONSTRAINTS][N_SP];
  double total[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  enum SPECIES dependent_species[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  double dependent_inverse[
      CHEM_MAX_CONSERVATION_CONSTRAINTS *
      CHEM_MAX_CONSERVATION_CONSTRAINTS];
  unsigned char is_dependent[N_SP];
#ifdef __cplusplus
  bool invert();
  bool build(double metallicity, const double *abundance);
  bool is_algebraic(enum SPECIES species) const;
  bool project(double *abundance) const;
  bool update_algebraic(double nH, double *abundance, const ChemRates *rates) const;
  bool reduce_jacobian(double nH, const double *abundance, const ChemRates *rates, const enum SPECIES *solver_species, int n_solver_species, double drfdy[N_SP][N_SP]) const;
#endif
} ChemConservationSystem;

/* Frozen physical initial state and bounds for a fixed-rate problem. */
typedef struct {
  double initial_abundance[N_SP];
  double cooling_initial_abundance[N_SP];
  double initialization_floor[N_SP];
  double newton_lower_bound[N_SP];
  double upper_bound[N_SP];
} ChemNewtonReference;

/* Mutable iterate/chart plus borrowed, read-only problem data. */
typedef struct {
  double abundance[N_SP];
  const double *initial_abundance;
  const double *cooling_initial_abundance;
  const double *initialization_floor;
  const double *newton_lower_bound;
  const double *upper_bound;
  double total_step[N_SP];
  enum SPECIES solver_species[N_SP];
  int solver_species_count;
  ChemConservationSystem conservation;
} ChemNewtonSystem;

int chemistry_species_is_intermediate_equilibrium(enum SPECIES species);
int chemistry_update_intermediate_equilibrium_species(
    unsigned char mode, double nH, double *abundance, const ChemRates *rates);

#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
int chemistry_is_full_implicit_initialized_intermediate(
    enum SPECIES species);
#endif

/* Shared Newton preparation and evaluation helpers. */

void chemistry_print_rerun_line(const ChemEnvironment *environment,
                                       double n_nuclei, double T_gas,
                                       double dt, const ChemState *state);

int chemistry_abundances_are_nonnegative(const double *abundance);

void chemistry_correct_full_implicit_intermediate_guess(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation);

int chemistry_prepare_implicit_newton_problem(
    const ChemEnvironment *environment, const ChemState *state,
    const ChemRates *rates, double nH, ChemNewtonReference *reference,
    ChemNewtonSystem *problem);

int chemistry_initialize_implicit_newton_guess(
    double nH, double dt, const double *newton_guess, const ChemRates *rates,
    ChemNewtonSystem *problem);

int chemistry_build_newton_jacobian(
    double nH, const double *abundance, const double *reaction_rhs,
    const ChemRates *rates, const ChemConservationSystem *conservation,
    const enum SPECIES *solver_species, int solver_species_count,
    int rebuild_newton_matrix, int *force_numerical_jacobian,
    double rhs_jacobian[N_SP][N_SP]);

void chemistry_update_implicit_timescale(
    ChemWorkspace *work, double dt, const double *abundance,
    const double *initial_abundance,
    const ChemConservationSystem *conservation);

int chemistry_export_implicit_solution(
    const ChemEnvironment *environment, const ChemStepData *step,
    ChemState *state, ChemWorkspace *work, const ChemState *input_state,
    double dt, double nH, double temperature,
    const double *initialization_floor, const double *cooling_initial_abundance,
    double *abundance, const ChemConservationSystem *conservation);

#ifdef __cplusplus
namespace chemistry { using ConservationSystem = ::ChemConservationSystem; }
#endif
