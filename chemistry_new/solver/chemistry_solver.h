#pragma once

#include "chemistry_types.h"

int chemistry_noneq_exp_nometal(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
int chemistry_noneq_chem_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt);
#ifdef SAHA
int chemistry_saha(const ChemEnvironment *environment, ChemState *state, ChemWorkspace *workspace, double dt);
#endif
