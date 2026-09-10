#pragma once

#include "chemistry_types.h"

ChemistryStatus chemistry_core_advance(int mode,
                                       const ChemEnvironment *environment,
                                       double duration_s,
                                       ChemState *state,
                                       ChemWorkspace *workspace);

/* Internal core helper used by host adapters exercising individual solvers. */
void chemistry_prepare_step_data(const ChemEnvironment *environment,
                                 const ChemState *initial_state,
                                 double initial_temperature,
                                 double duration_s,
                                 ChemStepData *step);
