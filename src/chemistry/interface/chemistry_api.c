#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "include/chemistry/chemistry.h"

#include "chemistry_core.h"
#include "chemistry_driver.h"
#include "chemistry_init.h"
#include "chemistry_types.h"
#include "species.h"

struct ChemistryContext {
  ChemRunContext run;
};

struct ChemistryWorkspace {
  ChemWorkspace core;
};

void chemistry_config_defaults(ChemistryConfig *config)
{
  if(config == NULL)
    return;

  memset(config, 0, sizeof(*config));
  config->abi_version = CHEMISTRY_ABI_VERSION;
  config->timescale_safety_factor = 0.1;
}

ChemistryStatus chemistry_context_create(const ChemistryConfig *config,
                                         ChemistryContext **context)
{
  if(context == NULL)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  *context = NULL;
  if(config == NULL || config->abi_version != CHEMISTRY_ABI_VERSION ||
     config->table_root == NULL || config->table_root[0] == '\0' ||
     !(config->timescale_safety_factor > 0.0) ||
     !isfinite(config->timescale_safety_factor))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  ChemistryContext *created = malloc(sizeof(*created));
  if(created == NULL)
    return CHEMISTRY_STATUS_OUT_OF_MEMORY;

  memset(created, 0, sizeof(*created));

  const ChemistryInitConfig init = {config->table_root};
  ChemistryStatus status = InitChem(&init, &created->run);
  if(status != CHEMISTRY_STATUS_OK) {
    free(created);
    return status;
  }

  created->run.timescale_safety_factor = config->timescale_safety_factor;
  *context = created;
  return CHEMISTRY_STATUS_OK;
}

void chemistry_context_destroy(ChemistryContext *context)
{
  if(context == NULL)
    return;

  DestroyChem(&context->run);
  free(context);
}

size_t chemistry_species_count(const ChemistryContext *context)
{
  return context != NULL ? TGCHEM_NUM_ABUNDANCES : 0;
}

const char *chemistry_species_name(const ChemistryContext *context,
                                   size_t species_index)
{
  if(context == NULL || species_index >= TGCHEM_NUM_ABUNDANCES)
    return NULL;

  const ChemSpeciesDefinition *definition =
      chemistry_species_definition((enum SPECIES)species_index);
  return definition != NULL ? definition->name : NULL;
}

size_t chemistry_workspace_size(const ChemistryContext *context)
{
  return context != NULL ? sizeof(ChemistryWorkspace) : 0;
}

ChemistryStatus chemistry_workspace_init(const ChemistryContext *context,
                                          void *storage,
                                          size_t storage_size,
                                          ChemistryWorkspace **workspace)
{
  if(workspace == NULL)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  *workspace = NULL;
  if(context == NULL || storage == NULL ||
     storage_size < sizeof(ChemistryWorkspace) ||
     (uintptr_t)storage % _Alignof(ChemistryWorkspace) != 0)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  memset(storage, 0, sizeof(ChemistryWorkspace));
  *workspace = storage;
  return CHEMISTRY_STATUS_OK;
}

void chemistry_workspace_reset(ChemistryWorkspace *workspace)
{
  if(workspace != NULL)
    memset(workspace, 0, sizeof(*workspace));
}

ChemistryStatus chemistry_step_prepare(
    const ChemistryContext *context,
    const ChemistryGlobalEnvironment *global_environment,
    ChemistryStepContext *step)
{
  if(context == NULL || global_environment == NULL || step == NULL ||
     !(global_environment->scale_factor > 0.0) ||
     !isfinite(global_environment->scale_factor) ||
     !(global_environment->cmb_temperature >= 0.0) ||
     !isfinite(global_environment->cmb_temperature))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  memset(step, 0, sizeof(*step));
  step->global = *global_environment;
  return CHEMISTRY_STATUS_OK;
}

static void import_global_environment(const ChemistryStepContext *step,
                                      ChemGlobalStep *global)
{
  memset(global, 0, sizeof(*global));
  global->cmb_temperature = step->global.cmb_temperature;
  global->scale_factor = step->global.scale_factor;
  global->simulation_time = step->global.simulation_time;
  global->comoving_integration = step->global.comoving_integration;
  global->background_J21 = step->global.background_J21;
  global->background_G0 = step->global.background_G0;
  global->background_rate_Hminus_pd =
      step->global.background_rate_Hminus_pd;
  global->background_rate_H2_diss =
      step->global.background_rate_H2_diss;
#ifdef HD
  global->background_rate_HD_diss =
      step->global.background_rate_HD_diss;
#endif
}

static void import_cell_environment(const ChemistryContext *context,
                                    const ChemGlobalStep *global,
                                    const ChemistryCellEnvironment *cell,
                                    ChemEnvironment *environment)
{
  memset(environment, 0, sizeof(*environment));
  environment->run = &context->run;
  environment->global = global;
  environment->rho_cgs = cell->density_cgs;
  environment->metallicity = cell->metallicity;
  environment->J21 = cell->local_J21;
  environment->G0 = cell->local_G0;
  environment->rate_Hminus_pd = cell->local_rate_Hminus_pd;
  environment->rate_H2_diss = cell->local_rate_H2_diss;
#ifdef HD
  environment->rate_HD_diss = cell->local_rate_HD_diss;
#endif
  environment->photoionization_rate_HI = cell->photoionization_rate_HI;
  environment->photoheating_rate_HI = cell->photoheating_rate_HI;
#ifdef HELIUM
  environment->photoionization_rate_HeI = cell->photoionization_rate_HeI;
  environment->photoheating_rate_HeI = cell->photoheating_rate_HeI;
  environment->photoionization_rate_HeII = cell->photoionization_rate_HeII;
  environment->photoheating_rate_HeII = cell->photoheating_rate_HeII;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  environment->rate_CO_diss = cell->local_rate_CO_diss;
#endif
  environment->stellar_radiation_flux = cell->stellar_radiation_flux;
  environment->ID = (int)cell->diagnostic_id;
}

ChemistryStatus chemistry_advance(
    const ChemistryContext *context,
    const ChemistryStepContext *step,
    ChemistryWorkspace *workspace,
    const ChemistryCellEnvironment *cell,
    double duration_s,
    const ChemistryAdvanceOptions *options,
    ChemistryState *state,
    ChemistryResult *result)
{
  if(context == NULL || step == NULL || workspace == NULL || cell == NULL ||
     options == NULL || state == NULL || result == NULL ||
     state->abundances == NULL ||
     state->abundance_count != TGCHEM_NUM_ABUNDANCES ||
     !(cell->density_cgs > 0.0) || !isfinite(cell->density_cgs) ||
     !(duration_s >= 0.0) || !isfinite(duration_s) ||
     options->mode != CHEMISTRY_MODE_DEFAULT ||
     options->solver != CHEMISTRY_SOLVER_AUTO)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  chemistry_workspace_reset(workspace);
  memset(result, 0, sizeof(*result));

  ChemGlobalStep global;
  ChemEnvironment environment;
  ChemState core_state;
  import_global_environment(step, &global);
  import_cell_environment(context, &global, cell, &environment);

  memset(&core_state, 0, sizeof(core_state));
  memcpy(core_state.abund, state->abundances,
         TGCHEM_NUM_ABUNDANCES * sizeof(*state->abundances));
  core_state.utherm_cgs = state->specific_thermal_energy_cgs;
  core_state.gamma = state->gamma;
  workspace->core.diagnostics_enabled = cell->diagnostics_enabled;

  const double initial_thermal_energy = core_state.utherm_cgs;
  ChemistryStatus status = chemistry_core_advance(
      (int)options->mode, &environment, duration_s, &core_state,
      &workspace->core);

  if(status == CHEMISTRY_STATUS_OK) {
    memcpy(state->abundances, core_state.abund,
           TGCHEM_NUM_ABUNDANCES * sizeof(*state->abundances));
    state->specific_thermal_energy_cgs = core_state.utherm_cgs;
    state->gamma = core_state.gamma;
  }

  const double cooling_radiative =
      fabs(workspace->core.evaluation.t_cool_rad);
  const double cooling_chemical =
      fabs(workspace->core.evaluation.t_cool_chem);
  result->cooling_time_s = cooling_radiative < cooling_chemical ?
      cooling_radiative : cooling_chemical;
  result->chemistry_time_s = fabs(workspace->core.evaluation.t_chem);
  result->thermal_energy_change_cgs = status == CHEMISTRY_STATUS_OK ?
      core_state.utherm_cgs - initial_thermal_energy : 0.0;
  result->dust_temperature = workspace->core.continuum.dust_temperature;
  result->dust_kappa_planck = workspace->core.continuum.dust_kappa_planck;
  result->continuum_optical_depth = workspace->core.continuum.tau;
  result->continuum_reduction_factor =
      workspace->core.continuum.reduction_factor;
  result->implicit_flag = workspace->core.evaluation.flag_implicit;
  result->convergence_flag = workspace->core.evaluation.flag_convergence;
  result->used_implicit_energy = result->implicit_flag == 1;
  result->used_implicit_chemistry = result->implicit_flag != 0;
  result->converged = status == CHEMISTRY_STATUS_OK &&
                      result->convergence_flag == 0;
  return status;
}
