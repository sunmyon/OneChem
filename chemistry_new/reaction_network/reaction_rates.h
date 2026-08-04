#pragma once

#include "chemistry_types.h"

int evaluate_reaction_rates(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
int evaluate_chemistry_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *workspace);
