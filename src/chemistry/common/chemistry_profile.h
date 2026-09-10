#pragma once

#include "chemistry_types.h"

typedef enum {
  CHEM_PROFILE_TOTAL,
  CHEM_PROFILE_ENERGY_SOLVER,
  CHEM_PROFILE_ENERGY_RESIDUAL,
  CHEM_PROFILE_IMPLICIT_CHEMISTRY,
  CHEM_PROFILE_REACTION_RATES,
  CHEM_PROFILE_INTERMEDIATES,
  CHEM_PROFILE_INTERMEDIATE_ASSEMBLY,
  CHEM_PROFILE_INTERMEDIATE_LINEAR_SOLVE,
  CHEM_PROFILE_NEWTON_RHS,
  CHEM_PROFILE_NEWTON_JACOBIAN,
  CHEM_PROFILE_JACOBIAN_SOLVER_COLUMNS,
  CHEM_PROFILE_JACOBIAN_RESPONSE,
  CHEM_PROFILE_JACOBIAN_INTERMEDIATE_SOLVE,
  CHEM_PROFILE_JACOBIAN_SCHUR,
  CHEM_PROFILE_NEWTON_LINEAR_SOLVE,
  CHEM_PROFILE_NEWTON_STEP,
  CHEM_PROFILE_CHEMICAL_COOLING,
  CHEM_PROFILE_RADIATIVE_COOLING,
  CHEM_PROFILE_MOLECULAR_COOLING,
  CHEM_PROFILE_DUST_TEMPERATURE,
  CHEM_PROFILE_COUNT
} ChemistryProfileSection;

#ifdef CHEMISTRY_PROFILE

double profile_start(void);
void profile_end(ChemistryProfileSection section, double started_at);
void chemistry_profile_set_failure(int reason);
int chemistry_profile_last_failure(void);
void chemistry_profile_record_residual(const ChemEnvironment *environment,
                                       const ChemState *state,
                                       const ChemWorkspace *work,
                                       double temperature, double dt,
                                       double energy_old, double energy_new,
                                       double residual, int chemistry_status);
void chemistry_profile_report_energy_warning(void);
void chemistry_profile_record_equilibrium_initialization(double fraction,
                                                         double electron_before,
                                                         double electron_after);
void chemistry_profile_record_newton_solve(int iterations, int reset_count,
                                           int clipped_steps);
void chemistry_profile_record_clipped_species(int species);
void chemistry_profile_record_newton_limiter(int species, double nH);
void chemistry_profile_report(void);

#else

static inline double profile_start(void) { return 0.0; }

static inline void profile_end(ChemistryProfileSection section,
                               double started_at) {
  (void)section;
  (void)started_at;
}

static inline void chemistry_profile_set_failure(int reason) { (void)reason; }

static inline int chemistry_profile_last_failure(void) { return 0; }

static inline void chemistry_profile_record_residual(
    const ChemEnvironment *environment, const ChemState *state,
    const ChemWorkspace *work, double temperature, double dt, double energy_old,
    double energy_new, double residual, int chemistry_status) {
  (void)environment;
  (void)state;
  (void)work;
  (void)temperature;
  (void)dt;
  (void)energy_old;
  (void)energy_new;
  (void)residual;
  (void)chemistry_status;
}

static inline void chemistry_profile_report_energy_warning(void) {}

static inline void chemistry_profile_record_equilibrium_initialization(
    double fraction, double electron_before, double electron_after) {
  (void)fraction;
  (void)electron_before;
  (void)electron_after;
}

static inline void chemistry_profile_record_newton_solve(int iterations,
                                                         int reset_count,
                                                         int clipped_steps) {
  (void)iterations;
  (void)reset_count;
  (void)clipped_steps;
}

static inline void chemistry_profile_record_clipped_species(int species) {
  (void)species;
}

static inline void chemistry_profile_record_newton_limiter(int species,
                                                           double nH) {
  (void)species;
  (void)nH;
}

static inline void chemistry_profile_report(void) {}

#endif
