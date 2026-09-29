// Driver orchestrates a cell step; Model and solvers remain independently owned.
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "chemistry_model.hpp"
#include "chemistry_cooling.h"
#include "cooling_func.h"
#include "reaction_rates.h"
#include "shielding.h"
#include "chemistry_solver.h"
#include "fixed_rate_solver.hpp"
#include "thermodynamics.h"


#define MAX_ITERATION 1000
#define HIGH_TEMPERATURE_THRESHOLD 1.0e5


void chemistry::Driver::prepare_cell(const ChemEnvironment &env, const ChemState &initial,
                                     double initial_temperature, double duration_s,
                                     ChemStepData &output) {
  const auto *environment = &env;
  const auto *initial_state = &initial;
  auto *step = &output;

  const double n_nuclei = environment->rho_cgs / PROTONMASS;

  memset(step, 0, sizeof(*step));
  step->jeans_length = get_core_radius(environment->rho_cgs, initial_temperature);
  step->column_HI = step->jeans_length * n_nuclei * initial_state->abund[sHI];
  step->column_H2 = step->jeans_length * n_nuclei * initial_state->abund[sH2I];
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  step->column_HD = step->jeans_length * n_nuclei * initial_state->abund[sHDI];
#endif

#ifdef CHEMISTRY_LW_SELF_SHIELDING
  step->shield_HI = EachShield(step->column_HI, 1.49e-1, 1.62e+0, 2.85e23);
  step->shield_H2 = EachShield(step->column_H2, 5.20e-3, 2.38e-1, 2.34e19);
  step->shield_H2self = SelfShield(step->column_H2, initial_temperature,
                                   9.12 * sqrt(initial_temperature / 1.e4));
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  step->shield_HDself = SelfShield(step->column_HD, initial_temperature,
                                   7.45 * sqrt(initial_temperature / 1.e4));
#endif
#endif

  step->G0 = environment->G0 + environment->global->background_G0;
  step->k22 = environment->photoionization_rate_HI;
  step->piHI = environment->photoheating_rate_HI;
#ifdef CHEMISTRY_SPECIES_HELIUM
  step->photoionization_rate_HeI = environment->photoionization_rate_HeI;
  step->photoheating_rate_HeI = environment->photoheating_rate_HeI;
  step->photoionization_rate_HeII = environment->photoionization_rate_HeII;
  step->photoheating_rate_HeII = environment->photoheating_rate_HeII;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  for(int stage = 0; stage < 3; stage++) {
    step->photoionization_rate_C[stage] =
        environment->photoionization_rate_C[stage];
    step->photoheating_rate_C[stage] = environment->photoheating_rate_C[stage];
    step->photoionization_rate_O[stage] =
        environment->photoionization_rate_O[stage];
    step->photoheating_rate_O[stage] = environment->photoheating_rate_O[stage];
  }
#endif

#ifdef IONIZING_PHOTON
  const double recombination_rate = 2.6e-13;
  const double neutral_fraction_limit = 1.e-3;
  const double neutral_fraction_eq = recombination_rate * XH * n_nuclei /
                                     (step->k22 + 1.e-300);

  if(step->k22 * duration_s > 1.0 && neutral_fraction_eq < neutral_fraction_limit) {
    const double factor = neutral_fraction_eq / neutral_fraction_limit;
    step->k22 *= factor;
    step->piHI *= factor;
  }
#endif
}

void chemistry::Driver::set_high_temperature_hydrogen_state(ChemState &state)
{
  state.abund[sHI] = 1.e-10 * XH;
  state.abund[sHII] = XH;
  state.abund[sH2I] = 1.e-10 * XH;
  state.abund[sHM] = 1.e-12 * XH;
  state.abund[sELEC] = XH;
  state.gamma = 5.0 / 3.0;
}

void chemistry::Driver::evaluate_high_temperature_cooling(const ChemEnvironment &env, const ChemStepData &cell_step, const ChemState &input, double temperature, ChemWorkspace &output)
{
  const auto *environment = &env;
  const auto *step = &cell_step;
  auto *state = &input;
  auto *workspace = &output;

  const double n_nuclei = environment->rho_cgs / PROTONMASS;
  double energy_rate = -8.418e-26 * pow(temperature, 0.11) * n_nuclei * n_nuclei;

  energy_rate -= 1.43e-27 * sqrt(temperature) *
                 (1.0 + 4.0 * XHe / XH) * n_nuclei * n_nuclei;

#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA
  energy_rate -= continuum_cooling_rate_thin(
      temperature, environment->rho_cgs, (double *) state->abund,
      workspace->diagnostics_enabled);
#endif
#if CHEMISTRY_METAL_LINE_COOLING_SD93
  energy_rate -= get_cooling_rate_SD93(environment->run->continuum_opacity_tables, temperature) * environment->metallicity *
                 n_nuclei * n_nuclei;
#endif
  energy_rate += step->piHI * state->abund[sHI] * n_nuclei;
#ifdef CHEMISTRY_SPECIES_HELIUM
  energy_rate += step->photoheating_rate_HeI * state->abund[sHeI] * n_nuclei;
  energy_rate += step->photoheating_rate_HeII * state->abund[sHeII] * n_nuclei;
#endif

  workspace->evaluation.rhoU_dot_rad = energy_rate;
  if(fabs(energy_rate) < 1.e-40) {
    workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  } else {
    const double energy_density = state->utherm_cgs * environment->rho_cgs;
    workspace->evaluation.t_cool_rad = energy_density / energy_rate;
  }
}

int chemistry::Driver::evolve_high_temperature(const ChemEnvironment &env, const ChemStepData &cell_step, ChemState &input, double &remaining, ChemWorkspace &output)
{
  const auto *environment = &env;
  const auto *step = &cell_step;
  auto *state = &input;
  auto *workspace = &output;
  auto *remaining_time_s = &remaining;

  set_high_temperature_hydrogen_state(*state);

  while(*remaining_time_s > 0.0) {
    update_thermodynamic_state(state, &workspace->thermodynamics);
    const double temperature = workspace->thermodynamics.temperature;
    if(temperature < HIGH_TEMPERATURE_THRESHOLD)
      break;

    evaluate_high_temperature_cooling(*environment, *step, *state, temperature, *workspace);

    double dt = *remaining_time_s;
    if(workspace->evaluation.t_cool_rad != MAX_REAL_NUMBER)
      dt = dmin(dt, 0.1 * fabs(workspace->evaluation.t_cool_rad));

    if(!(dt > 0.0) || !isfinite(dt))
      return 1;

    const double delta_utherm = workspace->evaluation.rhoU_dot_rad /
                                environment->rho_cgs * dt;
    if(state->utherm_cgs + delta_utherm <= 0.0)
      return 1;

    state->utherm_cgs += delta_utherm;
    *remaining_time_s -= dt;
  }

  return 0;
}

ChemSolverPlan chemistry::Driver::select_solvers(double n_nuclei, double time_cooling, double time_chemistry, double time_remaining)
{
  ChemSolverPlan plan = {CHEMISTRY_SOLVER_EXPLICIT, ENERGY_SOLVER_EXPLICIT,
                         dmin(time_cooling, dmin(time_chemistry, time_remaining))};

#if defined(CHEMISTRY_ENERGY_IMPLICIT) && !defined(CHEMISTRY_HYBRID_IMPLICIT)
  plan.energy = ENERGY_SOLVER_IMPLICIT;
  plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  plan.timestep = time_remaining;
#endif

#if defined(CHEMISTRY_IMPLICIT) && !defined(CHEMISTRY_HYBRID_IMPLICIT)
  plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  plan.timestep = dmin(time_cooling, time_remaining);
#endif

#ifdef CHEMISTRY_HYBRID_IMPLICIT
  plan.timestep = time_remaining;
  if(time_remaining > time_cooling || n_nuclei > N_EQCHEM) {
    plan.energy = ENERGY_SOLVER_IMPLICIT;
    plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  } else if(time_remaining > time_chemistry) {
    plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  } else if(n_nuclei > 1.e9) {
    if(time_chemistry > time_cooling) {
      plan.energy = ENERGY_SOLVER_IMPLICIT;
      plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
    } else {
      plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
    }
  }
#endif

#ifdef CHEMISTRY_SAHA
  if(n_nuclei > N_EQCHEM && plan.chemistry == CHEMISTRY_SOLVER_IMPLICIT)
    plan.chemistry = CHEMISTRY_SOLVER_SAHA;
#endif

  return plan;
}


ChemistryStatus chemistry::Driver::advance(const ChemEnvironment &env, double duration_s,
                                          ChemState &input, ChemWorkspace &output, int mode) const {
  const auto *environment = &env;
  auto *state = &input;
  auto *workspace = &output;

  if(environment->run == NULL || environment->global == NULL ||
     !(duration_s >= 0.0) || !isfinite(duration_s) ||
     (mode != 0 && mode != 1 && mode != 2))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  ChemistryStatus ifunc;
  ChemistryStatus status = CHEMISTRY_STATUS_OK;
  workspace->status = CHEMISTRY_STATUS_OK;
  /* Reconcile input composition once, before thermodynamics, rates and the
   * old state used for chemical heating. This also covers direct core calls. */
#ifdef CHEMISTRY_NORMALIZE_ABUNDANCES
  status = chemistry_normalize_abundances(
      state->abund, TGCHEM_NUM_ABUNDANCES, environment->metallicity);
  if(status != CHEMISTRY_STATUS_OK)
    return status;
#else
  /* Keep the existing charge-only reconstruction when normalization is off. */
  const double neutral_electron =
      chemistry_neutral_electron_abundance(state->abund);
  if(!isfinite(neutral_electron) || neutral_electron < 0.0)
    return CHEMISTRY_STATUS_INCOMPATIBLE_STATE;
  state->abund[sELEC] = neutral_electron;
#endif
  ChemStepData step;
  update_thermodynamic_state(state, &workspace->thermodynamics);
  double T_gas = workspace->thermodynamics.temperature;
  double dtime = duration_s;
  const double rho_cgs = environment->rho_cgs;
  const double n_nuclei = rho_cgs / PROTONMASS;

  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;

  prepare_cell(*environment, *state, T_gas, duration_s, step);
  if(workspace->diagnostics_enabled){
    printf("shield_HI=%g H2=%g H2self=%g column_HI=%g H2=%g T_gas=%g\n"
	   , step.shield_HI, step.shield_H2, step.shield_H2self, step.column_HI, step.column_H2, T_gas);
  }

  if(T_gas > HIGH_TEMPERATURE_THRESHOLD) {
    if(evolve_high_temperature(*environment, step, *state, dtime, *workspace))
      return CHEMISTRY_STATUS_NOT_CONVERGED;
    if(dtime <= 0.0)
      return CHEMISTRY_STATUS_OK;

    update_thermodynamic_state(state, &workspace->thermodynamics);
    T_gas = workspace->thermodynamics.temperature;
  }

  // Reuse the cell step data throughout low-temperature subcycling.
  status = evaluate_timesteps(*environment, step, *state, *workspace);
  if(status != CHEMISTRY_STATUS_OK)
    return status;
  if(mode==0 || mode ==2)
    return CHEMISTRY_STATUS_OK;

  if(workspace->diagnostics_enabled)
    printf("[ChemID%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , environment->ID, __LINE__, n_nuclei, T_gas
	   , state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I], state->abund[sHM], dtime);

  double t_now = 0.;    /*    set current time */

  for(int iter=0;iter<MAX_ITERATION;iter++){
    /*---------  Set the internal time step for subcycling */
    double time1 = environment->run->timescale_safety_factor * dmin(fabs(workspace->evaluation.t_cool_rad),
                                                                    fabs(workspace->evaluation.t_cool_chem));
    double time2 = environment->run->timescale_safety_factor * fabs(workspace->evaluation.t_chem);
    double time_elapsed = t_now;
    double time3 = dtime - time_elapsed;

    if(time3 <= 0.)
      goto LOOP_OUT;

    /* take the minimum */
    double dt = dmin(time1, dmin(time2,time3));

    if(iter >= MAX_ITERATION - 10){
      printf("[iter%d]time=%g %g %g => dt=%g\n",iter, time1, time2, time3, dt);
      printf("t_now=%g,dtime=%g\n",t_now, dtime);
      printf("n_nuclei=%g,Tgas=%g\n", n_nuclei, T_gas);
    }

    ChemSolverPlan solvers = select_solvers(n_nuclei, time1, time2, time3);
    dt = solvers.timestep;

    workspace->evaluation.flag_implicit = solvers.energy == ENERGY_SOLVER_IMPLICIT ? 1 :
                                          solvers.chemistry != CHEMISTRY_SOLVER_EXPLICIT ? 2 : 0;

    if(workspace->diagnostics_enabled)
      printf("[ChemID%d] line=%d iter=%d energy_solver=%d chemistry_solver=%d\n",
             environment->ID, __LINE__, iter, solvers.energy, solvers.chemistry);

    double T_gas_in = T_gas;

#if defined(CHEMISTRY_ENERGY_IMPLICIT)
    if(solvers.energy == ENERGY_SOLVER_IMPLICIT){
      workspace->evaluation.flag_convergence = EnergySolver{}.solve(*environment, step, *state, *workspace, solvers.chemistry, dt) != CHEMISTRY_STATUS_OK;
      if(workspace->status != CHEMISTRY_STATUS_OK)
        return workspace->status;
      if(workspace->evaluation.flag_convergence)
        status = CHEMISTRY_STATUS_NOT_CONVERGED;
    }else
#endif
      {
        ChemState state_old = *state;
	int isnan_flag=0;

#ifdef CHEMISTRY_SAHA
	if(solvers.chemistry == CHEMISTRY_SOLVER_SAHA){
	  ifunc = SahaSolver{}.solve(*environment, *state, *workspace, dt);
	  if(ifunc != CHEMISTRY_STATUS_OK)
	    status = CHEMISTRY_STATUS_NOT_CONVERGED;
	}else
#endif
	if(solvers.chemistry == CHEMISTRY_SOLVER_IMPLICIT){
	  ifunc = FixedRateSolver{}.advance(*environment, step, *state, *workspace, dt, {});
	  if(ifunc != CHEMISTRY_STATUS_OK)
	    status = CHEMISTRY_STATUS_NOT_CONVERGED;
	}else
	  {
	    ifunc = ExplicitSolver{}.solve(*environment, step, *state, *workspace, dt);
	    if(ifunc != CHEMISTRY_STATUS_OK){
	      isnan_flag = 1;
	      status = CHEMISTRY_STATUS_NOT_CONVERGED;
	      workspace->diagnostics_enabled = 1;
	    }
	  }

	state->gamma = chemistry_gamma_from_abundances(state->abund, T_gas_in);

	/* Update temperature */
	ChemistryStatus cooling_status =
	    chemistry::CoolingEvaluator(*environment, *state, *workspace).radiative(step);
	if(cooling_status != CHEMISTRY_STATUS_OK)
	  return cooling_status;

	double t_cool = dmin(fabs(workspace->evaluation.t_cool_rad), fabs(workspace->evaluation.t_cool_chem));

	double del_rhoU_rad  = workspace->evaluation.rhoU_dot_rad * dt;
	double del_rhoU_chem = workspace->evaluation.rhoU_dot_chem * dt; /* erg/cm3 */
	double del_rhoU_tot  = del_rhoU_rad + del_rhoU_chem;
	double del_U_tot     = del_rhoU_tot / rho_cgs; /* delta U in erg/gram */

	state->utherm_cgs += del_U_tot;
	if(workspace->diagnostics_enabled)
	  printf("(line%d) t_cool=%g t_chem=%g\n", __LINE__, t_cool, workspace->evaluation.t_chem);

#if defined(CHEMISTRY_HYBRID_IMPLICIT) && defined(CHEMISTRY_ENERGY_IMPLICIT)
	/************* DEBUGGING... [CHON200417] ************/
	if(0.1 * state_old.utherm_cgs < fabs(del_U_tot) || isnan_flag){
	  solvers.energy = ENERGY_SOLVER_IMPLICIT;
	  solvers.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
#ifdef CHEMISTRY_SAHA
	  if(n_nuclei > N_EQCHEM)
	    solvers.chemistry = CHEMISTRY_SOLVER_SAHA;
#endif
	  *state = state_old;

	  workspace->evaluation.flag_convergence = EnergySolver{}.solve(*environment, step, *state, *workspace, solvers.chemistry, dt) != CHEMISTRY_STATUS_OK;
	  if(workspace->status != CHEMISTRY_STATUS_OK)
	    return workspace->status;
	  if(workspace->evaluation.flag_convergence)
	    status = CHEMISTRY_STATUS_NOT_CONVERGED;
	  workspace->evaluation.flag_implicit = 1;
	}
#endif
      }

    t_now += dt;

    update_thermodynamic_state(state, &workspace->thermodynamics);
    T_gas = workspace->thermodynamics.temperature;

#ifdef CHEMISTRY_TEMPERATURE_MAX
    if(T_gas > CHEMISTRY_TEMPERATURE_MAX)
      state->utherm_cgs *= CHEMISTRY_TEMPERATURE_MAX / T_gas;
#endif

#ifdef CHEMISTRY_FIXED_TEMPERATURE
    state->utherm_cgs *= environment->run->fixed_temperature / T_gas;
    workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
    workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#endif

    /* stop if time = a_end */
    if(dtime - t_now < dtime/MAX_ITERATION)
      goto LOOP_OUT;
  }  /* iter */

 LOOP_OUT:;

  if(workspace->diagnostics_enabled)
    printf("[ChemID%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , environment->ID, __LINE__, n_nuclei, T_gas
	   , state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I], state->abund[sHM], dtime);

#ifdef CHEMISTRY_IMPLICIT
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;
#else
  workspace->evaluation.t_chem = fabs(workspace->evaluation.t_chem);
#endif

#ifdef CHEMISTRY_ENERGY_IMPLICIT
  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#else
  workspace->evaluation.t_cool_rad = fabs(workspace->evaluation.t_cool_rad);
  workspace->evaluation.t_cool_chem = fabs(workspace->evaluation.t_cool_chem);
#endif

  return status;
}

ChemistryStatus chemistry::Driver::evaluate_timesteps(const ChemEnvironment &env,
                                          const ChemStepData &cell_step,
                                          const ChemState &input,
                                          ChemWorkspace &output)
{
  const auto *environment = &env;
  const auto *step = &cell_step;
  auto *state = &input;
  auto *workspace = &output;

  update_thermodynamic_state(state, &workspace->thermodynamics);

  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;

#if !defined(CHEMISTRY_ENERGY_IMPLICIT) || defined(CHEMISTRY_HYBRID_IMPLICIT)
  chemistry::ReactionRateEvaluator(*environment, *state, *workspace).evaluate(*step);

#if !defined(CHEMISTRY_IMPLICIT) || defined(CHEMISTRY_HYBRID_IMPLICIT)
  chemistry::ReactionRateEvaluator(*environment, *state, *workspace).evaluate_timescale();
#endif

  ChemistryStatus cooling_status =
      chemistry::CoolingEvaluator(*environment, *state, *workspace).radiative(*step);
  if(cooling_status != CHEMISTRY_STATUS_OK)
    return cooling_status;
#endif

#ifdef CHEMISTRY_FIXED_TEMPERATURE
  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#endif

  return CHEMISTRY_STATUS_OK;
}


ChemistryStatus chemistry::Driver::prepare(const ChemistryGlobalEnvironment &global_input,
                                           ChemistryStepContext &step_output) {
  const auto *global_environment = &global_input;
  auto *step = &step_output;

  if(!(global_environment->scale_factor > 0.0) ||
     !isfinite(global_environment->scale_factor) ||
     !(global_environment->cmb_temperature >= 0.0) ||
     !isfinite(global_environment->cmb_temperature))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  memset(step, 0, sizeof(*step));
  step->global = *global_environment;
  return CHEMISTRY_STATUS_OK;
}

void chemistry::Driver::import_global_environment(const ChemistryStepContext &input,
                                                    ChemGlobalStep &output)
{
  const auto *step = &input;
  auto *global = &output;
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
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  global->background_rate_HD_diss =
      step->global.background_rate_HD_diss;
#endif
}

void chemistry::Driver::import_cell_environment(const Model &model,
    const ChemGlobalStep &global_input, const ChemistryCellEnvironment &cell_input,
    ChemEnvironment &output)
{
  const auto *global = &global_input;
  const auto *cell = &cell_input;
  auto *environment = &output;
  memset(environment, 0, sizeof(*environment));
  environment->run = &model.run();
  environment->global = global;
  environment->rho_cgs = cell->density_cgs;
  environment->metallicity = cell->metallicity;
  environment->J21 = cell->local_J21;
  environment->G0 = cell->local_G0;
  environment->rate_Hminus_pd = cell->local_rate_Hminus_pd;
  environment->rate_H2_diss = cell->local_rate_H2_diss;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  environment->rate_HD_diss = cell->local_rate_HD_diss;
#endif
  environment->photoionization_rate_HI = cell->photoionization_rate_HI;
  environment->photoheating_rate_HI = cell->photoheating_rate_HI;
#ifdef CHEMISTRY_SPECIES_HELIUM
  environment->photoionization_rate_HeI = cell->photoionization_rate_HeI;
  environment->photoheating_rate_HeI = cell->photoheating_rate_HeI;
  environment->photoionization_rate_HeII = cell->photoionization_rate_HeII;
  environment->photoheating_rate_HeII = cell->photoheating_rate_HeII;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  for(int stage = 0; stage < 3; stage++) {
    environment->photoionization_rate_C[stage] =
        cell->photoionization_rate_C[stage];
    environment->photoheating_rate_C[stage] = cell->photoheating_rate_C[stage];
    environment->photoionization_rate_O[stage] =
        cell->photoionization_rate_O[stage];
    environment->photoheating_rate_O[stage] = cell->photoheating_rate_O[stage];
  }
  environment->rate_CO_diss = cell->local_rate_CO_diss;
#endif
  environment->stellar_radiation_flux = cell->stellar_radiation_flux;
  environment->ID = (int)cell->diagnostic_id;
}

ChemistryStatus chemistry::Driver::advance(const Model &model, const ChemistryStepContext &step_input,
    ChemWorkspace &core_workspace, const ChemistryCellEnvironment &cell_input, double duration_s,
    const ChemistryAdvanceOptions &options_input, ChemistryState &state_input,
    ChemistryResult &result_output) const {
  const auto *step = &step_input;
  const auto *cell = &cell_input;
  const auto *options = &options_input;
  auto *state = &state_input;
  auto *result = &result_output;

  if(state->abundances == NULL ||
     state->abundance_count != TGCHEM_NUM_ABUNDANCES ||
     !(cell->density_cgs > 0.0) || !isfinite(cell->density_cgs) ||
     !(duration_s >= 0.0) || !isfinite(duration_s) ||
     options->mode != CHEMISTRY_MODE_DEFAULT ||
     options->solver != CHEMISTRY_SOLVER_AUTO)
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  core_workspace = {};
  memset(result, 0, sizeof(*result));

  ChemGlobalStep global;
  ChemEnvironment environment;
  ChemState core_state;
  import_global_environment(*step, global);
  import_cell_environment(model, global, *cell, environment);

  memset(&core_state, 0, sizeof(core_state));
  memcpy(core_state.abund, state->abundances,
         TGCHEM_NUM_ABUNDANCES * sizeof(*state->abundances));
  core_state.utherm_cgs = state->specific_thermal_energy_cgs;
  core_state.gamma = state->gamma;
  core_workspace.diagnostics_enabled = cell->diagnostics_enabled;

  const double initial_thermal_energy = core_state.utherm_cgs;
  ChemistryStatus status = advance(environment, duration_s, core_state, core_workspace, (int)options->mode);

  if(status == CHEMISTRY_STATUS_OK) {
    memcpy(state->abundances, core_state.abund,
           TGCHEM_NUM_ABUNDANCES * sizeof(*state->abundances));
    state->specific_thermal_energy_cgs = core_state.utherm_cgs;
    state->gamma = core_state.gamma;
  }

  const double cooling_radiative =
    fabs(core_workspace.evaluation.t_cool_rad);
  const double cooling_chemical =
    fabs(core_workspace.evaluation.t_cool_chem);
  result->cooling_time_s = cooling_radiative < cooling_chemical ?
    cooling_radiative : cooling_chemical;
  result->chemistry_time_s = fabs(core_workspace.evaluation.t_chem);
  result->thermal_energy_change_cgs = status == CHEMISTRY_STATUS_OK ?
    core_state.utherm_cgs - initial_thermal_energy : 0.0;
  result->dust_temperature = core_workspace.continuum.dust_temperature;
  result->dust_kappa_planck = core_workspace.continuum.dust_kappa_planck;
  result->continuum_optical_depth = core_workspace.continuum.tau;
  result->continuum_reduction_factor =
    core_workspace.continuum.reduction_factor;
  result->implicit_flag = core_workspace.evaluation.flag_implicit;
  result->convergence_flag = core_workspace.evaluation.flag_convergence;
  result->used_implicit_energy = result->implicit_flag == 1;
  result->used_implicit_chemistry = result->implicit_flag != 0;
  result->converged = status == CHEMISTRY_STATUS_OK &&
    result->convergence_flag == 0;
  return status;
}
