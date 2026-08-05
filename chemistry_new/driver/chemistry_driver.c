/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"

#include "chemistry.h"
#include "chemistry_constants.h"
#include "chemistry_cooling.h"
#include "chemistry_diagnostics.h"
#include "chemistry_models.h"
#include "chemistry_profile.h"
#include "chemistry_types.h"
#include "cooling_func.h"
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

static int evaluate_timesteps(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
int chemistry_advance(int mode, const ChemEnvironment *environment, double duration_s, ChemState *state, ChemWorkspace *workspace);

static void initialize_step_data(const ChemEnvironment *environment, const ChemState *initial_state, double initial_temperature, double duration_s, ChemStepData *step)
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

  step->G0 = environment->G0;
  step->k22 = environment->photoionization_rate_HI;
  step->piHI = environment->photoheating_rate_HI;

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
      temperature, environment->rho_cgs, (double *) state->abund);
#endif
#ifdef METAL_LINE_COOLING_SD93
  energy_rate -= get_cooling_rate_SD93(temperature) * environment->metallicity *
                 n_nuclei * n_nuclei;
#endif
  energy_rate += step->piHI * state->abund[sHI] * n_nuclei;

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

static void import_environment(const chemdata_in *input, ChemEnvironment *environment)
{
  memset(environment, 0, sizeof(*environment));
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
#ifdef BACKGROUND_LW_RADIATION
  environment->G0 += All.G0;
#endif
#ifdef IONIZATION_PHOTON
  environment->photoionization_rate_HI = input->k22;
  environment->photoheating_rate_HI = input->piHI;
#endif

  environment->ID = input->ID;
  environment->cmb_temperature = All.T_CMB;
  environment->scale_factor = All.cf_atime;
}

static void import_state(const chemdata_out *output, ChemState *state)
{
  memset(state->abund, 0, sizeof(state->abund));
  memcpy(state->abund, output->Abund, sizeof(output->Abund));
  state->utherm_cgs = output->Utherm *
                      (All.UnitEnergy_in_cgs / All.UnitMass_in_g);
  state->gamma = output->Gamma;
}

static void export_state(const ChemState *state, chemdata_out *output)
{
  memcpy(output->Abund, state->abund, sizeof(output->Abund));
  output->Utherm = state->utherm_cgs /
                   (All.UnitEnergy_in_cgs / All.UnitMass_in_g);
  output->Gamma = state->gamma;
}

static void export_workspace(const ChemWorkspace *workspace, chemdata_out *output)
{
  output->t_cool = dmin(fabs(workspace->evaluation.t_cool_rad),
                        fabs(workspace->evaluation.t_cool_chem));
  output->t_chem = fabs(workspace->evaluation.t_chem);
  output->Tdust = workspace->continuum.dust_temperature;
  output->kappa_dust = workspace->continuum.dust_kappa_planck;
#ifdef ONEZONE
  output->tau_cont = workspace->continuum.tau;
  output->continuum_reduction_factor = workspace->continuum.reduction_factor;
  memcpy(output->forward_rates, workspace->rates.forward, sizeof(output->forward_rates));
  memcpy(output->inverse_rates, workspace->rates.inverse, sizeof(output->inverse_rates));
  output->diagnostics = workspace->diagnostics;
#endif
}

int compute_abundances_d(int mode, const chemdata_in *input, chemdata_out *output, double dtime)
{
#ifdef CHEMISTRY_PROFILE
  double profile_start = chemistry_profile_time();
#endif
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};

  FlagPrint = input->flag_print;

  import_environment(input, &environment);
  import_state(output, &state);

  /* Transitional boundary: move this conversion into chem_balance later. */
  int status = chemistry_advance(mode, &environment, dtime, &state, &workspace);

  export_state(&state, output);
  export_workspace(&workspace, output);

#ifdef CHEMISTRY_PROFILE
  chemistry_profile_add(CHEM_PROFILE_TOTAL, chemistry_profile_time() - profile_start);
#endif
  return status;
}

int compute_chemistry_implicit_d(const chemdata_in *input, chemdata_out *output,
                                 double temperature, double dtime)
{
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;

  FlagPrint = input->flag_print;
  import_environment(input, &environment);
  import_state(output, &state);
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight = chemistry_molecular_weight(state.abund);
  workspace.evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace.evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace.evaluation.t_chem = MAX_REAL_NUMBER;
  initialize_step_data(&environment, &state, temperature, dtime, &step);

  int status = chemistry_noneq_chem_implicit(&environment, &step, &state, &workspace, dtime);
  export_state(&state, output);
  export_workspace(&workspace, output);
  return status;
}

int compute_equilibrium_abundances_d(const chemdata_in *input, const chemdata_out *output,
                                     double temperature, double dtime, double *full_abundance)
{
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;

  import_environment(input, &environment);
  import_state(output, &state);
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight = chemistry_molecular_weight(state.abund);
  initialize_step_data(&environment, &state, temperature, dtime, &step);
  evaluate_reaction_rates(&environment, &step, &state, &workspace);

  for(int species = 0; species < N_SP; species++)
    full_abundance[species] = state.abund[species] / XH;
  return chemistry_update_intermediate_equilibrium_species(
      CHEM_REACTION_IMPLICIT, XH * environment.rho_cgs / PROTONMASS,
      full_abundance, &workspace.rates);
}

int compute_chemistry_implicit_guess_d(const chemdata_in *input, chemdata_out *output,
                                       double temperature, double dtime,
                                       const double *newton_guess, int max_iterations)
{
  ChemEnvironment environment = {0};
  ChemState state;
  ChemWorkspace workspace = {0};
  ChemStepData step;
  double full_guess[N_SP] = {0};

  FlagPrint = input->flag_print;
  import_environment(input, &environment);
  import_state(output, &state);
  memcpy(full_guess, newton_guess, sizeof(output->Abund));
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight = chemistry_molecular_weight(state.abund);
  workspace.evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace.evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace.evaluation.t_chem = MAX_REAL_NUMBER;
  initialize_step_data(&environment, &state, temperature, dtime, &step);

  int status = chemistry_noneq_chem_implicit_with_guess(
      &environment, &step, &state, &workspace, dtime, full_guess, max_iterations);
  export_state(&state, output);
  export_workspace(&workspace, output);
  return status;
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


int chemistry_advance(int mode, const ChemEnvironment *environment, double duration_s, ChemState *state, ChemWorkspace *workspace)
{
  int ifunc, status = 0;
  ChemStepData step;
  update_thermodynamic_state(state, &workspace->thermodynamics);
  double T_gas = workspace->thermodynamics.temperature;
  double dtime = duration_s;
  const double rho_cgs = environment->rho_cgs;
  const double n_nuclei = rho_cgs / PROTONMASS;

  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;

  initialize_step_data(environment, state, T_gas, duration_s, &step);

  if(T_gas > HIGH_TEMPERATURE_THRESHOLD) {
    status = evolve_high_temperature(environment, &step, state, &dtime, workspace);
    if(status != 0 || dtime <= 0.0)
      return status;

    update_thermodynamic_state(state, &workspace->thermodynamics);
    T_gas = workspace->thermodynamics.temperature;
  }

  /* The low-temperature path below is intentionally left for the next
   * migration step.  It should consume environment, step, state and
   * workspace without rebuilding the fixed step data. */
  status |= evaluate_timesteps(environment, &step, state, workspace);
  if(mode==0 || mode ==2)
    return status;

  if(FlagPrint)
    printf("[Task%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , ThisTask, __LINE__, n_nuclei, T_gas
	   , state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I], state->abund[sHM], dtime);

  double t_now = 0.;    /*    set current time */

  for(int iter=0;iter<MAX_ITERATION;iter++){
    /*---------  Set the internal time step for subcycling */
    double time1 = All.Epsilon * dmin(fabs(workspace->evaluation.t_cool_rad),
                                      fabs(workspace->evaluation.t_cool_chem));
    double time2 = All.Epsilon * fabs(workspace->evaluation.t_chem);
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

    if(FlagPrint)
      printf("[Task%d] line=%d iter=%d energy_solver=%d chemistry_solver=%d\n",
             ThisTask, __LINE__, iter, solvers.energy, solvers.chemistry);

    double T_gas_in = T_gas;

#if defined(ENERGY_IMPLICIT)
    if(solvers.energy == ENERGY_SOLVER_IMPLICIT){
      workspace->evaluation.flag_convergence = chemistry_noneq_energy_implicit(environment, &step, state, workspace, solvers.chemistry, dt);
      if(workspace->evaluation.flag_convergence)
        status = 1;
    }else
#endif
      {
        ChemState state_old = *state;
	int isnan_flag=0;

#ifdef SAHA
	if(solvers.chemistry == CHEMISTRY_SOLVER_SAHA){
	  ifunc = chemistry_saha(environment, state, workspace, dt);
	  if(ifunc != 1)
	    status = 1;
	}else
#endif
	if(solvers.chemistry == CHEMISTRY_SOLVER_IMPLICIT){
	  ifunc = chemistry_noneq_chem_implicit(environment, &step, state, workspace, dt);
	  if(ifunc != 1)
	    status = 1;
	}else
	  {
	    ifunc = chemistry_noneq_explicit(environment, &step, state, workspace, dt);
	    if(ifunc == 2){
	      isnan_flag = 1;
	      status = 1;
	      FlagPrint = 1;  //we allow to display debug info for this erroneous case
	    }
	  }

	state->gamma = compute_gamma_from_abund(state->abund, T_gas_in);

	/* Update temperature */
	int cooling_status = evaluate_radiative_cooling(environment, &step, state, workspace);
	if(cooling_status)
	  status = 1;

	double t_cool = dmin(fabs(workspace->evaluation.t_cool_rad), fabs(workspace->evaluation.t_cool_chem));

	double del_rhoU_rad  = workspace->evaluation.rhoU_dot_rad * dt;
	double del_rhoU_chem = workspace->evaluation.rhoU_dot_chem * dt; /* erg/cm3 */
	double del_rhoU_tot  = del_rhoU_rad + del_rhoU_chem;
	double del_U_tot     = del_rhoU_tot / rho_cgs; /* delta U in erg/gram */

	state->utherm_cgs += del_U_tot;
	if(FlagPrint)
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
	  if(workspace->evaluation.flag_convergence)
	    status = 1;
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
    state->utherm_cgs *= All.InitTemperature / T_gas;
    workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
    workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#endif

    /* stop if time = a_end */
    if(dtime - t_now < dtime/MAX_ITERATION)
      goto LOOP_OUT;
  }  /* iter */

 LOOP_OUT:;

  if(FlagPrint)
    printf("[Task%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , ThisTask, __LINE__, n_nuclei, T_gas
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

static int evaluate_timesteps(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace)
{
  int status = 0;
  update_thermodynamic_state(state, &workspace->thermodynamics);

  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace->evaluation.t_chem = MAX_REAL_NUMBER;

#if !defined(ENERGY_IMPLICIT) || defined(HYBRID_IMPLICIT)
  status |= evaluate_reaction_rates(environment, step, state, workspace);

#if !defined(CHEMISTRY_IMPLICIT) || defined(HYBRID_IMPLICIT)
  status |= evaluate_chemistry_timescale(environment, state, workspace);
#endif

  status |= evaluate_radiative_cooling(environment, step, state, workspace);
#endif

#ifdef FIX_TEMPERATURE
  workspace->evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace->evaluation.t_cool_chem = MAX_REAL_NUMBER;
#endif

  return status;
}
#endif /* CHEMISTRY_SC */
