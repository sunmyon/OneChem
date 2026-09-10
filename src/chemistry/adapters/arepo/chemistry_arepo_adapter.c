#include <math.h>
#include <stddef.h>
#include <string.h>

#include "chemistry_arepo_adapter.h"

ChemistryStatus chemistry_arepo_context_init(const ChemArepoConfig *config,
                                             ChemArepoContext *context)
{
  if(config == NULL || context == NULL || config->table_root == NULL ||
     config->table_root[0] == '\0' ||
     !(config->subcycle_safety_factor > 0.0))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  memset(context, 0, sizeof(*context));
  ChemistryConfig chemistry_config;
  chemistry_config_defaults(&chemistry_config);
  chemistry_config.table_root = config->table_root;
  chemistry_config.subcycle_safety_factor = config->subcycle_safety_factor;
  return chemistry_context_create(&chemistry_config, &context->chemistry);
}

void chemistry_arepo_context_destroy(ChemArepoContext *context)
{
  if(context == NULL)
    return;

  chemistry_context_destroy(context->chemistry);
  memset(context, 0, sizeof(*context));
}

static void import_environment(const ChemArepoCell *cell,
                               ChemistryCellEnvironment *environment)
{
  memset(environment, 0, sizeof(*environment));
  environment->density_cgs = cell->rho_cgs;
  environment->metallicity = cell->metallicity;
  environment->local_J21 = cell->J21;
  environment->local_rate_Hminus_pd = cell->rate_Hminus_pd;
  environment->local_rate_H2_diss = cell->rate_H2_diss;
#ifdef HD
  environment->local_rate_HD_diss = cell->rate_HD_diss;
#endif
  environment->local_G0 = cell->G0;
  environment->stellar_radiation_flux = cell->stellar_radiation_flux;
  environment->photoionization_rate_HI = cell->photoionization_rate_HI;
  environment->photoheating_rate_HI = cell->photoheating_rate_HI;
  environment->photoionization_rate_HeI = cell->photoionization_rate_HeI;
  environment->photoheating_rate_HeI = cell->photoheating_rate_HeI;
  environment->photoionization_rate_HeII = cell->photoionization_rate_HeII;
  environment->photoheating_rate_HeII = cell->photoheating_rate_HeII;
  environment->local_rate_CO_diss = cell->rate_CO_diss;
  environment->diagnostics_enabled = cell->diagnostics_enabled;
  environment->diagnostic_id = (uint64_t)cell->diagnostic_id;
}

static void export_result(const ChemistryResult *chemistry_result,
                          ChemArepoResult *result)
{
  result->cooling_time_s = chemistry_result->cooling_time_s;
  result->chemistry_time_s = chemistry_result->chemistry_time_s;
  result->dust_temperature = chemistry_result->dust_temperature;
  result->dust_kappa_planck = chemistry_result->dust_kappa_planck;
  result->implicit_flag = chemistry_result->implicit_flag;
  result->convergence_flag = chemistry_result->convergence_flag;
}

ChemistryStatus chemistry_arepo_advance(const ChemArepoContext *context,
                                        const ChemistryStepContext *step,
                                        const ChemArepoCell *cell,
                                        double duration_s,
                                        ChemArepoState *state,
                                        ChemArepoResult *result)
{
  if(context == NULL || context->chemistry == NULL || step == NULL ||
     cell == NULL || state == NULL || result == NULL ||
     !(cell->rho_cgs > 0.0) || !(duration_s >= 0.0))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  ChemistryCellEnvironment environment;
  ChemistryState chemistry_state = {
      state->abund,
      TGCHEM_NUM_ABUNDANCES,
      state->utherm_cgs,
      state->gamma};
  ChemistryResult chemistry_result;
  const ChemistryAdvanceOptions options = {
      CHEMISTRY_MODE_DEFAULT,
      CHEMISTRY_SOLVER_AUTO};

  const size_t workspace_size = chemistry_workspace_size(context->chemistry);
  const size_t workspace_words =
      (workspace_size + sizeof(max_align_t) - 1) / sizeof(max_align_t);
  max_align_t workspace_storage[workspace_words];
  ChemistryWorkspace *workspace = NULL;

  memset(result, 0, sizeof(*result));
  import_environment(cell, &environment);

  ChemistryStatus status = chemistry_workspace_init(
      context->chemistry, workspace_storage, sizeof(workspace_storage),
      &workspace);
  if(status != CHEMISTRY_STATUS_OK)
    return status;

  status = chemistry_advance(
      context->chemistry, step, workspace, &environment, duration_s, &options,
      &chemistry_state, &chemistry_result);
  state->utherm_cgs = chemistry_state.specific_thermal_energy_cgs;
  state->gamma = chemistry_state.gamma;
  export_result(&chemistry_result, result);
  return status;
}
