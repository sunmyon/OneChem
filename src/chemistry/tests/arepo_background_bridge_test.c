#include <math.h>
#include <string.h>

#include "allvars.h"
#include "chemistry_arepo_host.h"

struct global_data_all_processes All;

ChemistryStatus chemistry_arepo_context_init(const ChemArepoConfig *config,
                                             ChemArepoContext *context)
{
  if(config == NULL || context == NULL)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  memset(context, 0, sizeof(*context));
  context->chemistry = (ChemistryContext *)context;
  return CHEMISTRY_STATUS_OK;
}

void chemistry_arepo_context_destroy(ChemArepoContext *context)
{
  if(context != NULL)
    memset(context, 0, sizeof(*context));
}

ChemistryStatus chemistry_step_prepare(
    const ChemistryContext *context,
    const ChemistryGlobalEnvironment *global,
    ChemistryStepContext *step)
{
  if(context == NULL || global == NULL || step == NULL)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  memset(step, 0, sizeof(*step));
  step->global = *global;
  return CHEMISTRY_STATUS_OK;
}

static int nearly_equal(double actual, double expected)
{
  return fabs(actual - expected) <= 1.e-13 * fabs(expected);
}

int main(void)
{
  ChemistryStepContext step;

  All.Epsilon = 0.01;
  All.J21 = 1.0;
  All.T_CMB = 30.0;
  All.cf_atime = 0.1;
  All.Time = 42.0;

  if(!chemistry_arepo_host_init() ||
     !chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_H2_diss, 1.39e-12) ||
     !nearly_equal(step.global.background_rate_HD_diss, 1.55e-12))
    return 1;

  if(!chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_HD_diss, 1.55e-12))
    return 1;

  All.J21 = 2.0;
  if(!chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_H2_diss, 2.78e-12) ||
     !nearly_equal(step.global.background_rate_HD_diss, 3.10e-12))
    return 1;

  chemistry_arepo_host_destroy();
  return 0;
}
