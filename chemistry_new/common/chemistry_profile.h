#pragma once

#ifdef CHEMISTRY_PROFILE

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

double chemistry_profile_time(void);
void chemistry_profile_add(ChemistryProfileSection section, double elapsed);
void chemistry_profile_set_failure(int reason);
int chemistry_profile_last_failure(void);
void chemistry_profile_record_residual(const ChemEnvironment *environment, const ChemState *state,
                                       const ChemWorkspace *work, double temperature, double dt,
                                       double energy_old, double energy_new, double residual,
                                       int chemistry_status);
void chemistry_profile_report_energy_warning(void);
void chemistry_profile_report(void);

#endif
