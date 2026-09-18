#include <math.h>
#include <string.h>

/* This file is deliberately host-facing.  These headers must not leak back
 * through chemistry_driver.c or any other core module. */
#include "allvars.h"

#include "chemistry_constants.h"
#include "chemistry_driver.h"
#include "chemistry_models.h"
#include "chemistry_standalone_adapter.h"
#include "chemistry_standalone_types.h"
#include "chemistry_profile.h"
#include "chemistry_solver.h"
#include "reaction_network.h"
#include "reaction_rates.h"
#include "standalone_math.h"

#ifdef CHEMISTRY_SC

static ChemRunContext standalone_run;
static ChemGlobalStep standalone_step;
static ChemStandaloneUnits standalone_units;
static int standalone_run_is_configured;
static int standalone_step_is_prepared;

void chemistry_standalone_configure(const ChemRunContext *run,
                                    const ChemStandaloneUnits *units)
{
  if(run == NULL || units == NULL || !(units->specific_energy_cgs > 0.0)) {
    standalone_run_is_configured = 0;
    return;
  }

  standalone_run = *run;
  standalone_units = *units;
  standalone_run_is_configured = 1;
}

int chemistry_standalone_set_timescale_safety_factor(double value)
{
  if(!standalone_run_is_configured || !isfinite(value) || value <= 0.0)
    return 0;
  standalone_run.timescale_safety_factor = value;
  return 1;
}

void chemistry_standalone_prepare_step(const ChemGlobalStep *step)
{
  if(step == NULL) {
    standalone_step_is_prepared = 0;
    return;
  }

  standalone_step = *step;
  standalone_step_is_prepared = 1;
}

static int adapter_is_ready(void)
{
  return standalone_run_is_configured && standalone_step_is_prepared;
}

static void import_environment(const chemdata_in *input,
                               ChemEnvironment *environment)
{
  memset(environment, 0, sizeof(*environment));
  environment->run = &standalone_run;
  environment->global = &standalone_step;
  environment->rho_cgs = input->rho_cgs;
#ifdef METALLICITY
  environment->metallicity = input->Metallicity;
#endif
#ifdef LW_RADIATION_LOCAL
  environment->J21 = input->J21;
  environment->rate_Hminus_pd = input->rate_Hminus_pd;
  environment->rate_H2_diss = input->rate_H2_diss;
#ifdef HD
  environment->rate_HD_diss = input->rate_HD_diss;
#endif
#ifdef DUST_COOLING
  environment->G0 = input->G0;
#endif
#endif
#ifdef THERMAL_FEEDBACK
  environment->stellar_radiation_flux = input->stellar_rad_flux;
#endif
#ifdef IONIZATION_PHOTON
  environment->photoionization_rate_HI = input->k22;
  environment->photoheating_rate_HI = input->piHI;
#endif
#ifdef HELIUM
  environment->photoionization_rate_HeI = input->photoionization_rate_HeI;
  environment->photoheating_rate_HeI = input->photoheating_rate_HeI;
  environment->photoionization_rate_HeII = input->photoionization_rate_HeII;
  environment->photoheating_rate_HeII = input->photoheating_rate_HeII;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  for(int stage = 0; stage < 3; stage++) {
    environment->photoionization_rate_C[stage] = input->photoionization_rate_C[stage];
    environment->photoheating_rate_C[stage] = input->photoheating_rate_C[stage];
    environment->photoionization_rate_O[stage] = input->photoionization_rate_O[stage];
    environment->photoheating_rate_O[stage] = input->photoheating_rate_O[stage];
  }
  environment->rate_CO_diss = input->rate_CO_diss;
#endif
  environment->ID = input->ID;
}

static void import_state(const chemdata_out *output, ChemState *state)
{
  memset(state->abund, 0, sizeof(state->abund));
  memcpy(state->abund, output->Abund, sizeof(output->Abund));
  state->utherm_cgs =
      output->Utherm * standalone_units.specific_energy_cgs;
  state->gamma = output->Gamma;
}

static void export_state(const ChemState *state, chemdata_out *output)
{
  memcpy(output->Abund, state->abund, sizeof(output->Abund));
  output->Utherm =
      state->utherm_cgs / standalone_units.specific_energy_cgs;
  output->Gamma = state->gamma;
}

static void export_workspace(const ChemWorkspace *workspace,
                             chemdata_out *output)
{
  output->t_cool = standalone_min(fabs(workspace->evaluation.t_cool_rad),
                                  fabs(workspace->evaluation.t_cool_chem));
  output->t_chem = fabs(workspace->evaluation.t_chem);
  output->Tdust = workspace->continuum.dust_temperature;
  output->kappa_dust = workspace->continuum.dust_kappa_planck;
  output->tau_cont = workspace->continuum.tau;
  output->continuum_reduction_factor = workspace->continuum.reduction_factor;
  memcpy(output->forward_rates, workspace->rates.forward,
         sizeof(output->forward_rates));
  memcpy(output->inverse_rates, workspace->rates.inverse,
         sizeof(output->inverse_rates));
  output->diagnostics = workspace->diagnostics;
}

int compute_abundances_d(int mode, const chemdata_in *input,
                         chemdata_out *output, double dtime)
{
  if(!adapter_is_ready())
    return 1;
#ifdef CHEMISTRY_PROFILE
  const double profile_timer = profile_start();
#endif
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};

  import_environment(input, &environment);
  import_state(output, &state);
  workspace.diagnostics_enabled = input->flag_print;

  const ChemistryStatus status = chemistry_core_advance(
      mode, &environment, dtime, &state, &workspace);
  export_state(&state, output);
  export_workspace(&workspace, output);

#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_TOTAL, profile_timer);
#endif
  /* Preserve the historical standalone ABI: zero succeeds, one fails. */
  return status == CHEMISTRY_STATUS_OK ? 0 : 1;
}

int compute_chemistry_implicit_d(const chemdata_in *input,
                                 chemdata_out *output,
                                 double temperature, double dtime)
{
  if(!adapter_is_ready())
    return 1;
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;

  import_environment(input, &environment);
  import_state(output, &state);
  workspace.diagnostics_enabled = input->flag_print;
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight =
      chemistry_molecular_weight(state.abund);
  workspace.evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace.evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace.evaluation.t_chem = MAX_REAL_NUMBER;
  chemistry_prepare_step_data(&environment, &state, temperature, dtime, &step);

  const int status = chemistry_noneq_chem_implicit(
      &environment, &step, &state, &workspace, dtime);
  export_state(&state, output);
  export_workspace(&workspace, output);
  return status;
}

int compute_equilibrium_abundances_d(const chemdata_in *input,
                                     const chemdata_out *output,
                                     double temperature, double dtime,
                                     double *full_abundance)
{
  if(!adapter_is_ready())
    return 1;
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;

  import_environment(input, &environment);
  import_state(output, &state);
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight =
      chemistry_molecular_weight(state.abund);
  chemistry_prepare_step_data(&environment, &state, temperature, dtime, &step);
  evaluate_reaction_rates(&environment, &step, &state, &workspace);

  for(int species = 0; species < N_SP; species++)
    full_abundance[species] = state.abund[species] / XH;
  return chemistry_update_intermediate_equilibrium_species(
      CHEM_REACTION_IMPLICIT, XH * environment.rho_cgs / PROTONMASS,
      full_abundance, &workspace.rates);
}

int compute_chemistry_implicit_guess_d(const chemdata_in *input,
                                       chemdata_out *output,
                                       double temperature, double dtime,
                                       const double *newton_guess,
                                       int max_iterations)
{
  if(!adapter_is_ready())
    return 1;
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;
  double full_guess[N_SP] = {0};

  import_environment(input, &environment);
  import_state(output, &state);
  workspace.diagnostics_enabled = input->flag_print;
  memcpy(full_guess, newton_guess, sizeof(output->Abund));
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight =
      chemistry_molecular_weight(state.abund);
  workspace.evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace.evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace.evaluation.t_chem = MAX_REAL_NUMBER;
  chemistry_prepare_step_data(&environment, &state, temperature, dtime, &step);

  const int status = chemistry_noneq_chem_implicit_with_guess(
      &environment, &step, &state, &workspace, dtime, full_guess,
      max_iterations);
  export_state(&state, output);
  export_workspace(&workspace, output);
  return status;
}

#endif
