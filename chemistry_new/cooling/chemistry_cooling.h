#pragma once

#include "chemistry_types.h"

int evaluate_radiative_cooling(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
int radiative_cooling_high_T(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
int chemical_cooling(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *workspace, double dt, double *abund_old);
