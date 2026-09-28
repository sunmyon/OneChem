/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "chemistry_core.h"
#include "include/chemistry/chemistry.h"
#include "species.h"
#include "chemistry_constants.h"
#include "chemistry_cooling.h"
#include "chemistry_diagnostics.h"
#include "chemistry_models.h"
#include "chemistry_profile.h"
#include "chemistry_types.h"
#include "cooling_func.h"
#include "chemistry_driver.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "reaction_rates.h"
#include "shielding.h"
#include "chemistry_solver.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

#define MAX_ITERATION 1000
#define HIGH_TEMPERATURE_THRESHOLD 1.0e5

static ChemistryStatus evaluate_timesteps(const ChemEnvironment *environment,
                                          const ChemStepData *step,
                                          const ChemState *state,
                                          ChemWorkspace *workspace);

void chemistry_prepare_step_data(const ChemEnvironment *environment, const ChemState *initial_state, double initial_temperature, double duration_s, ChemStepData *step)
{
  const double n_nuclei = environment->rho_cgs / PROTONMASS;

  memset(step, 0, sizeof(*step));
  step->jeans_length = get_core_radius(environment->rho_cgs, initial_temperature);
  step->column_HI = step->jeans_length * n_nuclei * initial_state->abund[sHI];
  step->column_H2 = step->jeans_length * n_nuclei * initial_state->abund[sH2I];
#ifdef HD
  step->column_HD = step->jeans_length * n_nuclei * initial_state->abund[sHDI];
#endif

#ifdef SHIELD_LW
  step->shield_HI = EachShield(step->column_HI, 1.49e-1, 1.62e+0, 2.85e23);
  step->shield_H2 = EachShield(step->column_H2, 5.20e-3, 2.38e-1, 2.34e19);
  step->shield_H2self = SelfShield(step->column_H2, initial_temperature,
                                   9.12 * sqrt(initial_temperature / 1.e4));
#ifdef HD
  step->shield_HDself = SelfShield(step->column_HD, initial_temperature,
                                   7.45 * sqrt(initial_temperature / 1.e4));
#endif
#endif

  step->G0 = environment->G0 + environment->global->background_G0;
  step->k22 = environment->photoionization_rate_HI;
  step->piHI = environment->photoheating_rate_HI;
#ifdef HELIUM
  step->photoionization_rate_HeI = environment->photoionization_rate_HeI;
  step->photoheating_rate_HeI = environment->photoheating_rate_HeI;
  step->photoionization_rate_HeII = environment->photoionization_rate_HeII;
  step->photoheating_rate_HeII = environment->photoheating_rate_HeII;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  for(int stage = 0; stage < 3; stage++) {
    step->photoionization_rate_C[stage] =
        environment->photoionization_rate_C[stage];
    step->photoheating_rate_C[stage] = environment->photoheating_rate_C[stage];
    step->photoionization_rate_O[stage] =
        environment->photoionization_rate_O[stage];
    step->photoheating_rate_O[stage] = environment->photoheating_rate_O[stage];
  }
#endif

#ifdef IONIZATION_PHOTON
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

static void set_high_temperature_hydrogen_state(ChemState *state)
{
  state->abund[sHI] = 1.e-10 * XH;
  state->abund[sHII] = XH;
  state->abund[sH2I] = 1.e-10 * XH;
  state->abund[sHM] = 1.e-12 * XH;
  state->abund[sELEC] = XH;
  state->gamma = 5.0 / 3.0;
}

static void evaluate_high_temperature_cooling(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, double temperature, ChemWorkspace *workspace)
{
  const double n_nuclei = environment->rho_cgs / PROTONMASS;
  double energy_rate = -8.418e-26 * pow(temperature, 0.11) * n_nuclei * n_nuclei;

  energy_rate -= 1.43e-27 * sqrt(temperature) *
                 (1.0 + 4.0 * XHe / XH) * n_nuclei * n_nuclei;

#ifdef CONTINUUM_COOLING_MATSUKOBA
  energy_rate -= continuum_cooling_rate_thin(
      temperature, environment->rho_cgs, (double *) state->abund,
      workspace->diagnostics_enabled);
#endif
#ifdef METAL_LINE_COOLING_SD93
  energy_rate -= get_cooling_rate_SD93(temperature) * environment->metallicity *
                 n_nuclei * n_nuclei;
#endif
  energy_rate += step->piHI * state->abund[sHI] * n_nuclei;
#ifdef HELIUM
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

static int evolve_high_temperature(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, double *remaining_time_s, ChemWorkspace *workspace)
{
  set_high_temperature_hydrogen_state(state);

  while(*remaining_time_s > 0.0) {
    update_thermodynamic_state(state, &workspace->thermodynamics);
    const double temperature = workspace->thermodynamics.temperature;
    if(temperature < HIGH_TEMPERATURE_THRESHOLD)
      break;

    evaluate_high_temperature_cooling(environment, step, state, temperature, workspace);

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

static ChemSolverPlan select_solvers(double n_nuclei, double time_cooling, double time_chemistry, double time_remaining)
{
  ChemSolverPlan plan = {CHEMISTRY_SOLVER_EXPLICIT, ENERGY_SOLVER_EXPLICIT,
                         dmin(time_cooling, dmin(time_chemistry, time_remaining))};

#if defined(ENERGY_IMPLICIT) && !defined(HYBRID_IMPLICIT)
  plan.energy = ENERGY_SOLVER_IMPLICIT;
  plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  plan.timestep = time_remaining;
#endif

#if defined(CHEMISTRY_IMPLICIT) && !defined(HYBRID_IMPLICIT)
  plan.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
  plan.timestep = dmin(time_cooling, time_remaining);
#endif

#ifdef HYBRID_IMPLICIT
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

#ifdef SAHA
  if(n_nuclei > N_EQCHEM && plan.chemistry == CHEMISTRY_SOLVER_IMPLICIT)
    plan.chemistry = CHEMISTRY_SOLVER_SAHA;
#endif

  return plan;
}


ChemistryStatus chemistry_core_advance(int mode,
                                       const ChemEnvironment *environment,
                                       double duration_s,
                                       ChemState *state,
                                       ChemWorkspace *workspace)
{
  if(environment == NULL || environment->run == NULL ||
     environment->global == NULL || state == NULL || workspace == NULL ||
     !(duration_s >= 0.0) || !isfinite(duration_s) ||
     (mode != 0 && mode != 1 && mode != 2))
    return CHEMISTRY_STATUS_INVALID_ARGUMENT;

  int ifunc;
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

  chemistry_prepare_step_data(environment, state, T_gas, duration_s, &step);
  if(workspace->diagnostics_enabled){
    printf("shield_HI=%g H2=%g H2self=%g column_HI=%g H2=%g T_gas=%g\n"
	   , step.shield_HI, step.shield_H2, step.shield_H2self, step.column_HI, step.column_H2, T_gas);
  }

  if(T_gas > HIGH_TEMPERATURE_THRESHOLD) {
    if(evolve_high_temperature(environment, &step, state, &dtime, workspace))
      return CHEMISTRY_STATUS_NOT_CONVERGED;
    if(dtime <= 0.0)
      return CHEMISTRY_STATUS_OK;

    update_thermodynamic_state(state, &workspace->thermodynamics);
    T_gas = workspace->thermodynamics.temperature;
  }

  /* The low-temperature path below is intentionally left for the next
   * migration step.  It should consume environment, step, state and
   * workspace without rebuilding the fixed step data. */
  status = evaluate_timesteps(environment, &step, state, workspace);
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

#if defined(ENERGY_IMPLICIT)
    if(solvers.energy == ENERGY_SOLVER_IMPLICIT){
      workspace->evaluation.flag_convergence = chemistry_noneq_energy_implicit(environment, &step, state, workspace, solvers.chemistry, dt);
      if(workspace->status != CHEMISTRY_STATUS_OK)
        return workspace->status;
      if(workspace->evaluation.flag_convergence)
        status = CHEMISTRY_STATUS_NOT_CONVERGED;
    }else
#endif
      {
        ChemState state_old = *state;
	int isnan_flag=0;

#ifdef SAHA
	if(solvers.chemistry == CHEMISTRY_SOLVER_SAHA){
	  ifunc = chemistry_saha(environment, state, workspace, dt);
	  if(ifunc != 1)
	    status = CHEMISTRY_STATUS_NOT_CONVERGED;
	}else
#endif
	if(solvers.chemistry == CHEMISTRY_SOLVER_IMPLICIT){
	  ifunc = chemistry_noneq_chem_implicit(environment, &step, state, workspace, dt);
	  if(ifunc != 1)
	    status = CHEMISTRY_STATUS_NOT_CONVERGED;
	}else
	  {
	    ifunc = chemistry_noneq_explicit(environment, &step, state, workspace, dt);
	    if(ifunc == 2){
	      isnan_flag = 1;
	      status = CHEMISTRY_STATUS_NOT_CONVERGED;
	      workspace->diagnostics_enabled = 1;
	    }
	  }

	state->gamma = chemistry_gamma_from_abundances(state->abund, T_gas_in);

	/* Update temperature */
	ChemistryStatus cooling_status =
	    evaluate_radiative_cooling(environment, &step, state, workspace);
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

#if defined(HYBRID_IMPLICIT) && defined(ENERGY_IMPLICIT)
	/************* DEBUGGING... [CHON200417] ************/
	if(0.1 * state_old.utherm_cgs < fabs(del_U_tot) || isnan_flag){
	  solvers.energy = ENERGY_SOLVER_IMPLICIT;
	  solvers.chemistry = CHEMISTRY_SOLVER_IMPLICIT;
#ifdef SAHA
	  if(n_nuclei > N_EQCHEM)
	    solvers.chemistry = CHEMISTRY_SOLVER_SAHA;
#endif
	  *state = state_old;

	  workspace->evaluation.flag_convergence = chemistry_noneq_energy_implicit(environment, &step, state, workspace, solvers.chemistry, dt);
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

#ifdef SET_TMAX
    if(T_gas > SET_TMAX)
      state->utherm_cgs *= SET_TMAX / T_gas;
#endif

#ifdef FIX_TEMPERATURE
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

#ifdef ENERGY_IMPLICIT
  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#else
  workspace->evaluation.t_cool_rad = fabs(workspace->evaluation.t_cool_rad);
  workspace->evaluation.t_cool_chem = fabs(workspace->evaluation.t_cool_chem);
#endif

  return status;
}

static ChemistryStatus evaluate_timesteps(const ChemEnvironment *environment,
                                          const ChemStepData *step,
                                          const ChemState *state,
                                          ChemWorkspace *workspace)
{
  update_thermodynamic_state(state, &workspace->thermodynamics);

  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;

#if !defined(ENERGY_IMPLICIT) || defined(HYBRID_IMPLICIT)
  if(evaluate_reaction_rates(environment, step, state, workspace))
    return CHEMISTRY_STATUS_INTERNAL_ERROR;

#if !defined(CHEMISTRY_IMPLICIT) || defined(HYBRID_IMPLICIT)
  if(evaluate_chemistry_timescale(environment, state, workspace))
    return CHEMISTRY_STATUS_INTERNAL_ERROR;
#endif

  ChemistryStatus cooling_status =
      evaluate_radiative_cooling(environment, step, state, workspace);
  if(cooling_status != CHEMISTRY_STATUS_OK)
    return cooling_status;
#endif

#ifdef FIX_TEMPERATURE
  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#endif

  return CHEMISTRY_STATUS_OK;
}
#endif /* CHEMISTRY_SC */
