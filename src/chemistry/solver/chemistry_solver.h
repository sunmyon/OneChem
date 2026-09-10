#pragma once

#include "chemistry_types.h"

int chemistry_noneq_explicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
int chemistry_noneq_chem_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
int chemistry_noneq_chem_implicit_with_guess(const ChemEnvironment *environment, const ChemStepData *step,
                                             ChemState *state, ChemWorkspace *workspace, double dt,
                                             const double *newton_guess, int max_iterations);
int chemistry_update_intermediate_equilibrium_species(unsigned char mode, double nH,
                                                       double *abundance, const ChemRates *rates);
int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt);
#ifdef SAHA
int chemistry_saha(const ChemEnvironment *environment, ChemState *state, ChemWorkspace *workspace, double dt);
#endif
