#include "fixed_rate_solver.hpp"
/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"
#include "chemistry_cooling.h"
#include "chemistry_types.h"

#include "cooling_func.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "reaction_rates.h"
#include "species.h"
#include "chemistry_newton_control.h"
#include "chemistry_newton_diagnostics.h"
#include "chemistry_newton_internal.h"
#include "chemistry_solver.h"
#include "saha.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

#define CHEMISTRY_NEWTON_TRACE_ABUNDANCE 1.e-12

static void repair_initial_abundances(double *abundance);

extern "C" int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
/* Fixed-rate Newton factorization now lives in solver/linear_algebra.cc.
 * The dgesv declaration above remains for the other chemistry paths.
 * Those paths still own their matrix assembly and solve calls.
 * Fixed-rate factor reuse is controlled by the C++ NewtonSession;
 * no factorization state is shared with the energy solver.
 */



/* Emit one reproducible "CHEMISTRY_RERUN_V1" replay line for the given
 * physical conditions and abundance state.  This is the single source of
 * the replay-line format: both the energy-solver failure report and the
 * chemistry Newton-solver failure report call this function, so a replay
 * line is always available to feed back into --chem-step / --energy-step,
 * regardless of which level detected the non-convergence. */
void chemistry_print_rerun_line(const ChemEnvironment *environment,
                                       double n_nuclei, double T_gas,
                                       double dt, const ChemState *state)
{
  printf("/*CHEMISTRY_RERUN_V1*/diagnostic_id=%d;comoving_integration=%d;"
         "simulation_time=%.17g;scale_factor=%.17g;n_nuclei=%.17g;"
         "T_gas=%.17g;cmb_temperature=%.17g;metallicity=%.17g;"
         "dt=%.17g;gamma=%.17g;utherm_cgs=%.17g;"
         "timescale_safety_factor=%.17g;",
         environment->ID, environment->global->comoving_integration,
         environment->global->simulation_time,
         environment->global->scale_factor, n_nuclei,
         T_gas, environment->global->cmb_temperature,
         environment->metallicity, dt, state->gamma,
         state->utherm_cgs,
         environment->run->timescale_safety_factor);

  for(int diagnostic_species = 0; diagnostic_species < sSTORED_END;
      diagnostic_species++)
    printf("yi[s%s]=%.17g;",
           chemistry_species_definition(static_cast<SPECIES>(diagnostic_species))->name,
           state->abund[diagnostic_species]);

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  printf("local_J21=%.17g;local_rate_Hminus_pd=%.17g;"
         "local_rate_H2_diss=%.17g;",
         environment->J21, environment->rate_Hminus_pd,
         environment->rate_H2_diss);
#ifdef HD
  printf("local_rate_HD_diss=%.17g;", environment->rate_HD_diss);
#endif
#ifdef DUST_COOLING
  printf("local_G0=%.17g;", environment->G0);
#endif
#endif
#ifdef BACKGROUND_LW_RADIATION
  printf("background_J21=%.17g;background_rate_Hminus_pd=%.17g;"
         "background_rate_H2_diss=%.17g;",
         environment->global->background_J21,
         environment->global->background_rate_Hminus_pd,
         environment->global->background_rate_H2_diss);
#ifdef HD
  printf("background_rate_HD_diss=%.17g;",
         environment->global->background_rate_HD_diss);
#endif
#ifdef METALLICITY
  printf("background_G0=%.17g;", environment->global->background_G0);
#endif
#endif
#endif

#ifdef IONIZATION_PHOTON
  printf("photoionization_rate_HI=%.17g;photoheating_rate_HI=%.17g;",
         environment->photoionization_rate_HI,
         environment->photoheating_rate_HI);
#endif
#ifdef HELIUM
  printf("photoionization_rate_HeI=%.17g;photoheating_rate_HeI=%.17g;"
         "photoionization_rate_HeII=%.17g;photoheating_rate_HeII=%.17g;",
         environment->photoionization_rate_HeI,
         environment->photoheating_rate_HeI,
         environment->photoionization_rate_HeII,
         environment->photoheating_rate_HeII);
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  const char *carbon_stage[3] = {"CI", "CII", "CIII"};
  const char *oxygen_stage[3] = {"OI", "OII", "OIII"};
  for(int stage = 0; stage < 3; stage++) {
    printf("photoionization_rate_%s=%.17g;photoheating_rate_%s=%.17g;",
           carbon_stage[stage], environment->photoionization_rate_C[stage],
           carbon_stage[stage], environment->photoheating_rate_C[stage]);
    printf("photoionization_rate_%s=%.17g;photoheating_rate_%s=%.17g;",
           oxygen_stage[stage], environment->photoionization_rate_O[stage],
           oxygen_stage[stage], environment->photoheating_rate_O[stage]);
  }
  printf("local_rate_CO_diss=%.17g;", environment->rate_CO_diss);
#endif
#ifdef THERMAL_FEEDBACK
  printf("stellar_radiation_flux=%.17g;",
         environment->stellar_radiation_flux);
#endif
  printf("\n");
}

#ifdef ENERGY_IMPLICIT
#define MAX_ITER 500
#define BETA_ENERGY 1.0

void chemistry::EnergySolver::restore_workspace(ChemWorkspace *destination,
                              const ChemWorkspace *source)
{
  ChemistryStatus status = destination->status;
  *destination = *source;
  if(status != CHEMISTRY_STATUS_OK)
    destination->status = status;
}
/* This is not the normal energy convergence tolerance.  It only identifies
 * a discontinuous chemistry-branch jump that leaves a residual many orders
 * of magnitude larger than the thermal energy. */
#define ENERGY_GROSS_RESIDUAL_FACTOR 1.e6

ChemistryStatus chemistry::EnergySolver::solve(const ChemEnvironment &env, const ChemStepData &cell_step,
                                   ChemState &input, ChemWorkspace &output,
                                   ChemistrySolver chemistry_solver, double dt) const {
  const EnergyInput problem{env, cell_step, chemistry_solver, dt};
  const auto *environment = &env;
  const auto *step = &cell_step;
  auto *state = &input;
  auto *work = &output;

#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  const double rho_gas = environment->rho_cgs;
  const double T_gas = work->thermodynamics.temperature;

  int iter = 0, ITER_reset = 30;
  double denom=0., err, err_chem;
  double final_chemistry_error = 1.0;
  int final_evaluation_available = 0;
  double e_new;
  double T_gas_init, T_old, T_new, dT, d_T;
  const double T_max=1.e6, T_min=0.1;
  double T_large, T_small;
  double g = 0., g0;

  ChemWorkspace work_old = *work;
  ChemState state_old = *state;
  T_gas_init =  T_gas;


  T_large = 5000, T_small = 1.5;

  T_old = T_gas_init, dT = 1.0e-3 * T_old;
  T_new = T_old;

#ifdef IONIZATION_PHOTON
  double nH = XH * rho_gas / PROTONMASS;
  double k2 = 2.6e-13;
  double x_neut_eq = MAX_REAL_NUMBER;
  if(step->k22 > 0.)
    x_neut_eq = k2 * nH / step->k22;
  double T_temp = step->piHI * dt/ BOLTZMANN;

  if(x_neut_eq < 0.1){
    if(step->k22 > 0. && step->k22 * dt > 0.2)
      T_temp = dmin(T_temp, step->piHI / step->k22 / BOLTZMANN);
  }else{
    if(T_temp > 2.e4){
      double Lambda = continuum_cooling_rate_thin_Lambda(2.e4);
      double fneut_crit = step->piHI / 2./ Lambda/ nH;
      if(fneut_crit > 1.e-6)
	T_temp = dmin(T_temp, 1.e5);

      if(fneut_crit > 1.e-4)
	T_temp = dmin(T_temp, 5.e4);

      if(fneut_crit > 1.e-2)
	T_temp = dmin(T_temp, 2.5e4);
    }
  }

  if(T_temp > 2.e5)
    T_temp = 2.e5;

  if(work->diagnostics_enabled)
    printf("[2]piHI = %g, %g T_temp=%g k22=%g x_neut_eq=%g dt=%g\n", step->piHI, step->piHI, T_temp, step->k22,x_neut_eq, dt);

  if(T_temp > T_gas){
    T_old = T_temp;
    dT = 0.1 * T_temp;
    T_new = T_old;
  }
#endif

  /***************************************************************/

#ifdef FLAG_BRENT
  double g_plus = 1.0e30, g_minus = -1.0e30;
  double T_plus, T_minus;
  int flag_brent = 0, sign;
#endif
  int flag_converged = 0;

  while(1){
    err_chem = trial_residual(problem, state, work, T_new, &e_new, &g0);
    if(work->status != CHEMISTRY_STATUS_OK) {
      err = NAN;
      goto END;
    }

    if(err_chem < 1.e-2)
      break;

    if(isnan(dT) || isinf(dT)){
      dT = 0.1 * T_new;
      printf("(LINE=%d) dT=nan -> %g g=%g g0=%g denom=%g ID=%d\n",
	     __LINE__, dT, g, g0, denom, environment->ID);
    }

    if(T_new + dT > T_max)
      T_new = dmax(T_large,T_min) , T_large /= 1.1;
    /*  else if(T_new + dT > 30000) T_new /= 1.1;*/
    else if(T_new + dT < 1.5)
      T_new = dmin(T_small,T_max), T_small *= 1.1;
    else T_new += dT, T_new *=1.001 ;

    restore_workspace(work, &work_old);
    *state = state_old;

    if(isnan(T_new)){
      printf("[ChemID%d](iter=%d) T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	     environment->ID, iter, state->utherm_cgs, state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I]);
      printf("(LINE%d)T_gas=%g T_new=%g T_l=%g T_s=%g dT=%g\n",
	     __LINE__, T_gas_init, T_new, T_large, T_small, dT);
      fflush(stdout);
    }

    if(iter++ > 20){
      printf("didn't converged chem_loop\n"), err = 1.0;
      goto END;
    }
  }


#ifdef FLAG_BRENT

  if(g0 > 0) g_plus = g0, T_plus = T_new, sign = 1;
  else g_minus = g0, T_minus = T_new, sign = -1;
#endif

  restore_workspace(work, &work_old);
  *state = state_old;

  T_old = T_new;
  T_new = T_old + dT;


  int j;
  for(j=iter;j<MAX_ITER;j++){
#ifdef FLAG_BRENT
    if( flag_brent == 1){
      restore_workspace(work, &work_old);
      *state = state_old;
      err = solve_bisection(problem, state, work, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

      if(err < 1.e-5)
	flag_converged = 1;

      break;
    }

    if(j==ITER_reset) break;
#endif


    /************************************************************/
    /**** When the chemistry network does not be converged, *****/
    /**** we iteratively search for the temperature,        *****/
    /**** at which the chemical reaction be converged.      *****/
    /************************************************************/
    while(1){
      double T_new_temp;
      err_chem = trial_residual(problem, state, work, T_new, &e_new, &g);
      if(work->status != CHEMISTRY_STATUS_OK) {
        err = NAN;
        goto END;
      }
#ifdef ID_CHEM_DEBUG
      if(work->diagnostics_enabled)
	printf("[%d] T=%g dT=%g err_chem=%g\n", j, T_new, dT, err_chem);
#endif

      denom = g - g0;
      if(fabs(denom) < 1.e-10){
	if(denom > 0)
	  denom = 0.01;
	else
	  denom = -0.01;
      }

      dT = BETA_ENERGY*g*(T_old - T_new)/denom;

      if(err_chem < 1e-4)
	break;

      if(isnan(dT) || isinf(dT)){
	dT = 0.1 * T_new;
	printf("(LINE=%d) dT=nan -> %g g=%g g0=%g denom=%g ID=%d\n",
	       __LINE__, dT, g, g0, denom, environment->ID);
      }

      T_new_temp = T_new;

      if(T_new + dT > T_max)
	T_new = dmax(T_large,T_min), T_large /= 1.1;
      else if(T_new + dT < 1.5){
	if(err_chem < 1.e-4 && T_new < 1.1 * T_small)
	  T_new = 0.5 * (T_new + T_old);
	else
	  T_new = dmin(T_small,T_max), T_small *= 1.1;
      }else{
	T_new = T_new_temp + dT;
      }

      if(j++ > MAX_ITER){
	printf("didn't converged T_new=%g dT=%g.\n"
	       , T_new, dT), err = 1.0;
	goto END;
      }

      restore_workspace(work, &work_old);
      *state = state_old;
    }


    /************************************************************/
    /************************************************************/
    /************************************************************/

#ifdef FLAG_BRENT
    if(g > 0){
      if(g < g_plus) g_plus = g, T_plus = T_new;
    }else{
      if(g > g_minus) g_minus = g, T_minus = T_new;
    }

    if( g*sign < 0 && j>10)  flag_brent = 1;
#endif


    g0 = g;
    T_old = T_new;

    if(T_new + dT > T_max){
      T_new = dmax(T_large,T_min), T_large /= 1.1;
      err = 1.0;
    }else if(T_new + dT < T_small){
      T_new = dmin(T_small,T_max), T_small *= 1.1;
      err = 1.0;
    }else{
      err   = fabs(dT/T_new);
      T_new += dT;
    }


    if(err < 1.0e-8 || (err < 1.0e-5 && j > 20)){
      flag_converged = 1;
      break;
    }

    restore_workspace(work, &work_old);
    *state = state_old;
  }


#ifdef FLAG_BRENT
  if(flag_converged == 0){
    T_new = T_gas_init;

    for(;j<MAX_ITER;j++){
      restore_workspace(work, &work_old);
      *state = state_old;

      d_T = dmax(0.2 * T_new, 2000);
      if(sign > 0)
	T_new = dmax(T_new-d_T,T_new/10.);
      else
	T_new = dmin(T_new + d_T,T_new * 10.);

      err_chem = trial_residual(problem, state, work, T_new, &e_new, &g);
      if(work->status != CHEMISTRY_STATUS_OK) {
        err = NAN;
        goto END;
      }

      if(err_chem < 1.e-4 && sign * g < 0.){
	if(sign > 0)
	  g_minus = g, T_minus = T_new;
	else
	  g_plus = g, T_plus = T_new;

	flag_brent = 1;
	break;
      }

      if(isnan(T_new)){
	printf("[ChemID%d](iter=%d) T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	       environment->ID, j, state->utherm_cgs, state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I]);
	printf("(LINE%d)T_gas=%g T_new=%g T_l=%g T_s=%g dT=%g\n"
	       , __LINE__, T_gas_init, T_new, T_large, T_small, dT);
	fflush(stdout);
      }

      if(T_new < T_min || T_new > T_max){
	printf("didn't converged or T_new(%g) is too low.\n", T_new);
	err = 1., flag_converged = 0;

	break;
      }
    }

    if(flag_brent == 1){
      restore_workspace(work, &work_old);
      *state = state_old;
      err = solve_bisection(problem, state, work, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

      if(err < 1.e-5) flag_converged = 1;
    }
  }
#endif

  final_chemistry_error = commit_residual(problem, state, work, T_new, &e_new, &g);
  final_evaluation_available = 1;
  state->utherm_cgs = e_new;

 END:;
  double normalized_energy_residual = final_evaluation_available ?
      fabs(g) / dmax(dmax(fabs(e_new), fabs(state_old.utherm_cgs)), 1.0) :
      MAX_REAL_NUMBER;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  if(chemistry_solver == CHEMISTRY_SOLVER_IMPLICIT &&
     (err > 1.e-4 || !final_evaluation_available ||
      final_chemistry_error >= 1.e-4 ||
      normalized_energy_residual >= ENERGY_GROSS_RESIDUAL_FACTOR)) {
    ChemState normal_state = *state;
    ChemWorkspace normal_workspace = *work;
    double normal_temperature = T_new;
    double normal_energy = e_new;
    double normal_residual = g;
    double normal_error = err;
    ChemState continuation_state;
    ChemWorkspace continuation_workspace;
    double continuation_temperature, continuation_energy,
           continuation_residual;
    if(solve_continuation(problem, &state_old, &work_old, T_gas_init, T_min, T_max, &continuation_state, &continuation_workspace, &continuation_temperature, &continuation_energy, &continuation_residual)) {
      *state = continuation_state;
      restore_workspace(work, &continuation_workspace);
      state->utherm_cgs = continuation_energy;
      T_new = continuation_temperature;
      e_new = continuation_energy;
      g = continuation_residual;
      err = fabs(g) /
          dmax(dmax(fabs(e_new), fabs(state_old.utherm_cgs)), 1.0);
      final_chemistry_error = work->evaluation.err_chem;
      final_evaluation_available = 1;
    } else {
      *state = normal_state;
      restore_workspace(work, &normal_workspace);
      T_new = normal_temperature;
      e_new = normal_energy;
      g = normal_residual;
      err = normal_error;
    }
  }
#endif
  if(!isfinite(err) || !isfinite(T_new)) {
    restore_workspace(work, &work_old);
    *state = state_old;
    work->status = CHEMISTRY_STATUS_INTERNAL_ERROR;
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_ENERGY_SOLVER, profile_timer);
#endif
    return CHEMISTRY_STATUS_NOT_CONVERGED;
  }

  if(err > 1.0e-4 || work->diagnostics_enabled){
    printf("%s err=%g id=%d T_old=%g T_new=%g\n",
	   err > 1.0e-4 ? "not converged(in energy)." : "converged(in energy diagnostic).",
	   err, environment->ID, T_gas_init, T_new);
    chemistry_print_rerun_line(environment, rho_gas / PROTONMASS, T_gas_init,
                               dt, &state_old);
  }

  if(err > 1.0e-4){
    if(chemistry_solver == CHEMISTRY_SOLVER_IMPLICIT &&
       (!final_evaluation_available || final_chemistry_error >= 1.e-4)) {
      ChemState diagnostic_state = state_old;
      ChemWorkspace diagnostic_workspace = work_old;
      double diagnostic_energy, diagnostic_residual;
      commit_residual_with_report(problem, &diagnostic_state, &diagnostic_workspace, T_new, &diagnostic_energy, &diagnostic_residual);
    }
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_report_energy_warning();
#endif
    restore_workspace(work, &work_old);
    *state = state_old;
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_ENERGY_SOLVER, profile_timer);
#endif
    return CHEMISTRY_STATUS_NOT_CONVERGED;
  }

#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_ENERGY_SOLVER, profile_timer);
#endif
  return CHEMISTRY_STATUS_OK;
}
double chemistry::EnergySolver::evaluate_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g, ResidualMode mode, const double *newton_guess, int max_chemistry_iterations)
{
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  const bool commit = mode != ResidualMode::trial;
  const bool report_chemistry_failure = mode == ResidualMode::commit_with_report;
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  ChemWorkspace work_trial = *work;
  ChemState state_trial = *state;
  double e_old = state_trial.utherm_cgs;
  ChemistryStatus chemistry_status = CHEMISTRY_STATUS_OK;

  work_trial.thermodynamics.temperature = T_gas_new;
#ifdef SAHA
  if(chemistry_solver == CHEMISTRY_SOLVER_SAHA)
    chemistry_status = chemistry::SahaSolver{}.solve(*environment, state_trial, work_trial, dt);
  else
#endif
    {
#ifdef CHEMISTRY_PROFILE
      double chemistry_timer = profile_start();
#endif
      // Temperature trials may reject chemistry without emitting failure reports.
      const chemistry::FixedRateOptions options = newton_guess
          ? chemistry::FixedRateOptions{newton_guess, max_chemistry_iterations, 0, 0}
          : chemistry::FixedRateOptions{nullptr, chemistry::default_chemistry_iterations, 1,
                                       report_chemistry_failure};
      chemistry_status = chemistry::FixedRateSolver{}.advance(
          *environment, *step, state_trial, work_trial, dt, options);
#ifdef CHEMISTRY_PROFILE
      profile_end(CHEM_PROFILE_IMPLICIT_CHEMISTRY, chemistry_timer);
#endif
    }
  if(chemistry_status != CHEMISTRY_STATUS_OK) {
    /* A failed chemical trial has no updated chemical energy rate. Never
     * evaluate an energy residual using the copied rate from the old state. */
    *e_new = NAN;
    *g = NAN;
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_ENERGY_RESIDUAL, profile_timer);
#endif
    return 1.0;
  }
  ChemistryStatus cooling_status =
      chemistry::CoolingEvaluator(*environment, state_trial, work_trial).radiative(*step);
  if(cooling_status != CHEMISTRY_STATUS_OK) {
    /* Residual probes use copies of the workspace. Preserve a fatal cooling
     * evaluation error explicitly so it cannot be mistaken for an ordinary
     * rejected root-finding point. */
    work->status = cooling_status;
    *e_new = NAN;
    *g = NAN;
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_ENERGY_RESIDUAL, profile_timer);
#endif
    return 1.0;
  }

  double T_gamma = chemistry_temperature(state_trial.abund, state_trial.gamma, state_trial.utherm_cgs);
  state_trial.gamma =
      chemistry_gamma_from_abundances(state_trial.abund, T_gamma);

  double gamma_new = state_trial.gamma;
  double xmu_new = chemistry_molecular_weight(state_trial.abund);

  *e_new = BOLTZMANN * T_gas_new/( (gamma_new - 1.0)*xmu_new*PROTONMASS );
  double de  = (work_trial.evaluation.rhoU_dot_chem + work_trial.evaluation.rhoU_dot_rad)/environment->rho_cgs * dt;

#ifdef IONIZATION_PHOTON
  double e_max = BOLTZMANN * 1.e6 /( (gamma_new - 1.0)*xmu_new*PROTONMASS );
  if(e_old + de > e_max && de > 0.){
    de =  dmax(e_max - e_old, 0.);
  }
#endif

  *g  = *e_new - e_old - de;

#ifdef CHEMISTRY_PROFILE
  chemistry_profile_record_residual(environment, &state_trial, &work_trial, T_gas_new, dt,
                                    e_old, *e_new, *g, chemistry_status == CHEMISTRY_STATUS_OK);
  profile_end(CHEM_PROFILE_ENERGY_RESIDUAL, profile_timer);
#endif

  if(work->diagnostics_enabled)
    printf("resid(e_new-e_old-de)=%g e_new=%g old=%g de=%g Tgas=%g rhoU_dot_chem=%g cool=%g rho_gas=%g dt=%g xmu=%g gamma=%g\n"
	   , *g, *e_new, e_old, de, T_gas_new, work_trial.evaluation.rhoU_dot_chem, work_trial.evaluation.rhoU_dot_rad, environment->rho_cgs, dt, xmu_new, gamma_new);

  if(commit){
    *state = state_trial;
    *work = work_trial;
  }
  return work_trial.evaluation.err_chem;
}

double chemistry::EnergySolver::trial_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g){
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  return evaluate_residual(problem, state, work, T_gas_new, e_new, g, ResidualMode::trial, nullptr, 0);
}

double chemistry::EnergySolver::residual_with_guess(const EnergyInput &problem, ChemState *state,
                                ChemWorkspace *work, double T_gas_new, double *e_new,
                                double *g, const double *newton_guess,
                                int max_chemistry_iterations){
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  return evaluate_residual(problem, state, work, T_gas_new, e_new, g, ResidualMode::commit, newton_guess, max_chemistry_iterations);
}

double chemistry::EnergySolver::commit_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g){
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  return evaluate_residual(problem, state, work, T_gas_new, e_new, g, ResidualMode::commit, nullptr, 0);
}

double chemistry::EnergySolver::commit_residual_with_report(const EnergyInput &problem, ChemState *state,
                              ChemWorkspace *work,
                              double T_gas_new, double *e_new, double *g){
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  return evaluate_residual(problem, state, work, T_gas_new, e_new, g, ResidualMode::commit_with_report, nullptr, 0);
}



int chemistry::EnergySolver::evaluate_energy_continuation_point(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace, double temperature,
    const double *newton_guess, EnergyContinuationPoint *point)
{
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  point->state = *initial_state;
  point->workspace = *initial_workspace;
  point->temperature = temperature;
  double chemistry_error;
  if(newton_guess)
    chemistry_error = residual_with_guess(problem, &point->state, &point->workspace, temperature, &point->energy, &point->residual, newton_guess, 1000);
  else
    chemistry_error = commit_residual(problem, &point->state, &point->workspace, temperature, &point->energy, &point->residual);
  if(chemistry_error >= 1.e-4 || !isfinite(point->energy) ||
     !isfinite(point->residual))
    return 0;
  point->state.utherm_cgs = point->energy;
  return 1;
}

int chemistry::EnergySolver::refine_energy_continuation_bracket(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace,
    EnergyContinuationPoint lower, EnergyContinuationPoint upper,
    EnergyContinuationPoint *root)
{
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  for(int iteration = 0; iteration < 80; iteration++) {
    double temperature = 0.5 * (lower.temperature + upper.temperature);
    const EnergyContinuationPoint *nearest =
        temperature - lower.temperature < upper.temperature - temperature ?
        &lower : &upper;
    EnergyContinuationPoint middle;
    if(!evaluate_energy_continuation_point(problem, initial_state, initial_workspace, temperature, nearest->state.abund, &middle)) {
      const EnergyContinuationPoint *other = nearest == &lower ? &upper : &lower;
      if(!evaluate_energy_continuation_point(problem, initial_state, initial_workspace, temperature, other->state.abund, &middle))
        return 0;
    }
    if(lower.residual * middle.residual <= 0.0)
      upper = middle;
    else
      lower = middle;

    EnergyContinuationPoint candidate =
        fabs(lower.residual) < fabs(upper.residual) ? lower : upper;
    double residual_scale = dmax(
        dmax(fabs(candidate.energy), fabs(initial_state->utherm_cgs)), 1.0);
    if(fabs(candidate.residual) / residual_scale < 1.e-6 &&
       fabs(upper.temperature - lower.temperature) /
           dmax(candidate.temperature, 1.0) < 1.e-8) {
      *root = candidate;
      return 1;
    }
  }
  return 0;
}

int chemistry::EnergySolver::solve_continuation(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace, double temperature_initial,
    double temperature_min, double temperature_max, ChemState *state,
    ChemWorkspace *workspace, double *temperature, double *energy,
    double *residual)
{
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  EnergyContinuationPoint start;
  if(!evaluate_energy_continuation_point(problem, initial_state, initial_workspace, temperature_initial, NULL, &start))
    return 0;

  double start_scale = dmax(
      dmax(fabs(start.energy), fabs(initial_state->utherm_cgs)), 1.0);
  if(fabs(start.residual) / start_scale < 1.e-6) {
    *state = start.state;
    *workspace = start.workspace;
    *temperature = start.temperature;
    *energy = start.energy;
    *residual = start.residual;
    return 1;
  }

  int primary_direction = start.residual < 0.0 ? 1 : -1;
  for(int direction_index = 0; direction_index < 2; direction_index++) {
    int direction = direction_index == 0 ? primary_direction : -primary_direction;
    EnergyContinuationPoint previous = start;
    for(int iteration = 0; iteration < 20000; iteration++) {
      double nominal_step = dmax(1.e-4, 1.e-3 * previous.temperature);
      double step_size = nominal_step;
      EnergyContinuationPoint trial;
      int valid = 0;
      for(int retry = 0; retry < 24; retry++) {
        double trial_temperature = previous.temperature + direction * step_size;
        if(trial_temperature <= temperature_min ||
           trial_temperature >= temperature_max)
          break;
        if(evaluate_energy_continuation_point(problem, initial_state, initial_workspace, trial_temperature, previous.state.abund, &trial)) {
          valid = 1;
          break;
        }
        step_size *= 0.5;
      }
      if(!valid)
        break;
      if(previous.residual * trial.residual <= 0.0) {
        EnergyContinuationPoint lower = previous.temperature < trial.temperature ?
            previous : trial;
        EnergyContinuationPoint upper = previous.temperature < trial.temperature ?
            trial : previous;
        EnergyContinuationPoint root;
        if(refine_energy_continuation_bracket(problem, initial_state, initial_workspace, lower, upper, &root)) {
          *state = root.state;
          *workspace = root.workspace;
          *temperature = root.temperature;
          *energy = root.energy;
          *residual = root.residual;
          return 1;
        }
      }
      previous = trial;
    }
  }
  return 0;
}


double chemistry::EnergySolver::solve_bisection(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_min, double g_min, double T_max, double g_max, double *T_new0, double *e_new0, int max_iter){
  const auto *environment = &problem.environment;
  const auto *step = &problem.step;
  const auto chemistry_solver = problem.solver;
  const double dt = problem.duration;

  int iter;
  double err = 1.0;
  double T_a, T_b;
  double T_cand;
  double val_a;
  double val_new;
  double e_new;

  if(g_min * g_max > 0){
    printf("something wrong in solve_bisection."
	   " g_min and g_max must have opposite sign.\n");
    return 1.e30;
  }

  T_a = T_min, val_a = g_min;
  T_b = T_max;
  T_cand = 0.5 * (T_a + T_b);
  e_new = *e_new0;

  iter = 0;
  while(iter < max_iter && err > 1.0e-6){
    T_cand = 0.5 * (T_a + T_b);

    err = trial_residual(problem, state, work, T_cand, &e_new, &val_new);
    if(work->status != CHEMISTRY_STATUS_OK)
      return 1.e30;

    /* trial_residual returns a non-finite residual when the chemistry trial at T_cand
     * failed.  Such a point carries no sign information: comparing it with
     * `<` is always false, so it would be taken as the same side as val_a,
     * val_a would become non-finite, and every later midpoint would move
     * T_a toward T_b, ending with a "converged" bracket that is not a root.
     * Fail instead and let the caller's recovery path handle it. */
    if(!isfinite(val_new)){
      printf("[ChemID%d] solve_bisection: non-finite residual at T=%g "
             "(chemistry trial failed); bracket [%g, %g] cannot be refined.\n",
             environment->ID, T_cand, T_a, T_b);
      return 1.e30;
    }

    if(val_a * val_new < 0)
      T_b = T_cand;
    else T_a = T_cand, val_a = val_new;

    err = fabs( (T_a - T_b)/T_b );

    if( err < 1.0e-5) break;

    iter++;
  }

  *T_new0 = T_cand;
  *e_new0 = e_new;

  return err;
}
#endif


#ifdef CHEM_INTERMEDIATE_EQUILIBRIUM
static const enum SPECIES intermediate_equilibrium_species[] = {
  sHM,
#ifdef CARBON_OXYGEN_CHEMISTRY
  sCOII, sCH, sCH2, sOII, sOHII, sH2OII, sH3OII, sO2II,
#endif
};
static const size_t intermediate_equilibrium_species_count = sizeof(intermediate_equilibrium_species) / sizeof(intermediate_equilibrium_species[0]);
#else
static const enum SPECIES intermediate_equilibrium_species[1] = {sEND};
static const size_t intermediate_equilibrium_species_count = 0;
#endif

int chemistry_species_is_intermediate_equilibrium(enum SPECIES species)
{
  for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
    if(intermediate_equilibrium_species[i] == species)
      return 1;
  return 0;
}

static int chemistry_solve_intermediate_equilibrium_list(unsigned char mode,
                                                         double nH,
                                                         double *abundance,
                                                         const ChemRates *rates,
                                                         const enum SPECIES *species_list,
                                                         size_t species_count)
{
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#define INTERMEDIATE_RETURN(value) do { profile_end(CHEM_PROFILE_INTERMEDIATES, profile_timer); return (value); } while(0)
#else
#define INTERMEDIATE_RETURN(value) return (value)
#endif
  double updated[N_SP];

  if(species_count == 0)
    INTERMEDIATE_RETURN(1);

  if(species_count == 1) {
    enum SPECIES species = species_list[0];
    if(!chemistry_network_intermediate_equilibrium_abundance(rates->network,
        species, mode, nH, abundance, rates->forward, rates->inverse,
        &updated[0], NULL))
      INTERMEDIATE_RETURN(0);
    abundance[species] = updated[0];
    INTERMEDIATE_RETURN(1);
  }

  int n_intermediate = (int)species_count;
  double trial[N_SP];
  double rhs[N_SP];
  double matrix[N_SP * N_SP];
  int pivots[N_SP];
  int one = 1;
  int info;

  memcpy(trial, abundance, sizeof(trial));
  for(int i = 0; i < n_intermediate; i++)
    trial[species_list[i]] = 0.0;

#ifdef CHEMISTRY_PROFILE
  double assembly_timer = profile_start();
#endif
  chemistry_network_intermediate_system(rates->network, mode, nH, trial,
                                        rates->forward, rates->inverse,
                                        species_list, n_intermediate, rhs, matrix);

  for(int row = 0; row < n_intermediate; row++)
    updated[row] = -rhs[row];
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_INTERMEDIATE_ASSEMBLY, assembly_timer);
  double intermediate_linear_timer = profile_start();
#endif

  dgesv_(&n_intermediate, &one, matrix, &n_intermediate, pivots, updated, &n_intermediate, &info);
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_INTERMEDIATE_LINEAR_SOLVE,
              intermediate_linear_timer);
#endif
  if(info != 0)
    INTERMEDIATE_RETURN(0);

  double negative_tolerance = 1.e-30;
#ifdef CARBON_OXYGEN_CHEMISTRY
  double element_scale = 0.0;
  const enum SPECIES stored_metal_species[] = {sCI, sCII, sCO, sCO2, sOI, sOH, sH2O, sO2};
  for(size_t i = 0; i < sizeof(stored_metal_species) / sizeof(stored_metal_species[0]); i++)
    element_scale = dmax(element_scale, fabs(abundance[stored_metal_species[i]]));
  negative_tolerance = dmax(negative_tolerance, 1.e-12 * element_scale);
#endif

  for(int i = 0; i < n_intermediate; i++) {
    if(!isfinite(updated[i]) || updated[i] < -negative_tolerance)
      INTERMEDIATE_RETURN(0);
    abundance[species_list[i]] = dmax(updated[i], 0.0);
  }
  INTERMEDIATE_RETURN(1);
#undef INTERMEDIATE_RETURN
}

int chemistry_update_intermediate_equilibrium_species(unsigned char mode, double nH, double *abundance, const ChemRates *rates)
{
  return chemistry_solve_intermediate_equilibrium_list(
      mode, nH, abundance, rates, intermediate_equilibrium_species,
      intermediate_equilibrium_species_count);
}

static void chemistry_add_conservation_constraint(ChemConservationSystem *conservation, enum SPECIES dependent_species, double total)
{
  int constraint = conservation->count++;
  conservation->dependent_species[constraint] = dependent_species;
  conservation->total[constraint] = total;
  conservation->is_dependent[dependent_species] = 1;
}

bool ChemConservationSystem::invert() {
  auto *conservation = this;
  const int count = conservation->count;
  double augmented[CHEM_MAX_CONSERVATION_CONSTRAINTS][2 * CHEM_MAX_CONSERVATION_CONSTRAINTS];

  for(int constraint = 0; constraint < count; constraint++)
    for(int column = 0; column < 2 * count; column++)
      augmented[constraint][column] = column < count ?
          conservation->coefficient[constraint][conservation->dependent_species[column]] :
          (column - count == constraint ? 1.0 : 0.0);

  for(int pivot = 0; pivot < count; pivot++) {
    int pivot_row = pivot;
    for(int row = pivot + 1; row < count; row++)
      if(fabs(augmented[row][pivot]) > fabs(augmented[pivot_row][pivot]))
        pivot_row = row;
    if(fabs(augmented[pivot_row][pivot]) < 1.e-30)
      return 0;
    if(pivot_row != pivot)
      for(int column = 0; column < 2 * count; column++) {
        double temporary = augmented[pivot][column];
        augmented[pivot][column] = augmented[pivot_row][column];
        augmented[pivot_row][column] = temporary;
      }

    double inverse_pivot = 1.0 / augmented[pivot][pivot];
    for(int column = 0; column < 2 * count; column++)
      augmented[pivot][column] *= inverse_pivot;
    for(int row = 0; row < count; row++) {
      if(row == pivot)
        continue;
      double factor = augmented[row][pivot];
      for(int column = 0; column < 2 * count; column++)
        augmented[row][column] -= factor * augmented[pivot][column];
    }
  }

  for(int dependent = 0; dependent < count; dependent++)
    for(int constraint = 0; constraint < count; constraint++)
      conservation->dependent_inverse[constraint * count + dependent] = augmented[dependent][count + constraint];
  return 1;
}

static enum SPECIES chemistry_largest_reservoir(
    const double *abundance, ChemElement element,
    const enum SPECIES *reservoirs, int reservoir_count)
{
  enum SPECIES largest = reservoirs[0];
  double largest_content =
      chemistry_species_definition(largest)->nuclei[element] *
      abundance[largest];

  for(int i = 1; i < reservoir_count; i++) {
    enum SPECIES species = reservoirs[i];
    double content =
        chemistry_species_definition(static_cast<SPECIES>(species))->nuclei[element] *
        abundance[species];

    if(content > largest_content) {
      largest = species;
      largest_content = content;
    }
  }

  return largest;
}

bool ChemConservationSystem::build(double metallicity, const double *abundance) {
  auto *conservation = this;
  memset(conservation, 0, sizeof(*conservation));

  for(size_t element_index = 0; element_index < chemistry_element_conservation_count; element_index++)
    {
      const ChemElementConservationDefinition *element = &chemistry_element_conservation[element_index];

#ifdef CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION
      if(element->element == CHEM_ELEMENT_C || element->element == CHEM_ELEMENT_O)
        continue;
#endif

      enum SPECIES dependent =
          chemistry_largest_reservoir(abundance, element->element, element->reservoirs, element->newton_reservoir_count);

      int constraint = conservation->count;
      double total   = element->total_per_hydrogen;

      if(element->metallicity_scaled)
        total *= metallicity;

      chemistry_add_conservation_constraint(conservation, dependent, total);

      for(int species = 0; species < N_SP; species++)
        conservation->coefficient[constraint][species] = chemistry_species_definition(static_cast<SPECIES>(species))->nuclei[element->element];
    }

  enum SPECIES charge_dependent  = sELEC;
  double largest_negative_charge = -1.0;
  for(int species = 0; species < N_SP; species++)
    {
      const ChemSpeciesDefinition *definition = chemistry_species_definition(static_cast<SPECIES>(species));
      if(definition->charge >= 0 || chemistry_species_is_intermediate_equilibrium(static_cast<SPECIES>(species)) || conservation->is_dependent[species])
        continue;
      double charge_content = fabs((double)definition->charge) * dmax(abundance[species], 0.0);
      if(charge_content > largest_negative_charge)
        {
          largest_negative_charge = charge_content;
          charge_dependent        = static_cast<SPECIES>(species);
        }
    }
  int charge_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, charge_dependent, 0.0);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[charge_constraint][species] = chemistry_species_definition(static_cast<SPECIES>(species))->charge;

#if defined(CARBON_OXYGEN_CHEMISTRY) && !defined(CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION)
  /* Select C/O columns jointly. A nonsingular matrix is insufficient: the
   * reconstructed abundances must also lie in the nonnegative domain. */
  const enum SPECIES candidates[] = {
    sCI, sCII, sCO, sCO2, sOI, sOH, sH2O, sO2,
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
    sCIII, sOIV,
#endif
  };
  int cr = -1, orow = -1;
  for(int row = 0; row < conservation->count; row++)
    {
      if(conservation->coefficient[row][sCI] == 1.0)
        cr = row;
      if(conservation->coefficient[row][sOI] == 1.0)
        orow = row;
    }
  if(cr >= 0 && orow >= 0)
    {
      ChemConservationSystem best;
      double best_score = -1.0;
      size_t count      = sizeof(candidates) / sizeof(candidates[0]);
      for(size_t c = 0; c < count; c++)
        {
          if(!chemistry_species_definition(candidates[c])->nuclei[CHEM_ELEMENT_C])
            continue;
          for(size_t o = 0; o < count; o++)
            {
              if(c == o || !chemistry_species_definition(candidates[o])->nuclei[CHEM_ELEMENT_O])
                continue;
              ChemConservationSystem trial  = *conservation;
              trial.dependent_species[cr]   = candidates[c];
              trial.dependent_species[orow] = candidates[o];
              memset(trial.is_dependent, 0, sizeof(trial.is_dependent));
              for(int row = 0; row < trial.count; row++)
                trial.is_dependent[trial.dependent_species[row]] = 1;
              if(!(&trial)->invert())
                continue;
              double projected[N_SP];
              memcpy(projected, abundance, sizeof(projected));
              if(!(&trial)->project(projected))
                continue;
              int feasible = 1;
              for(int species = 0; species < N_SP; species++)
                if(!isfinite(projected[species]) || projected[species] < 0.0)
                  feasible = 0;
              if(!feasible)
                continue;
              /* Maximize the smaller normalized reservoir margin. The inverse
               * weights account for amplification of elemental residuals. */
              double score     = INFINITY;
              const int rows[] = {cr, orow};
              for(int index = 0; index < 2; index++)
                {
                  int row      = rows[index];
                  double scale = 0.0;
                  for(int constraint = 0; constraint < trial.count; constraint++)
                    scale += fabs(trial.dependent_inverse[constraint * trial.count + row]) * fabs(trial.total[constraint]);
                  score = dmin(score, scale > 0.0 ? projected[trial.dependent_species[row]] / scale : 0.0);
                }
              if(score > best_score)
                {
                  best       = trial;
                  best_score = score;
                }
            }
        }
      if(best_score < 0.0)
        return 0;
      *conservation = best;
      return 1;
    }
#endif
  return (conservation)->invert();
}





bool ChemConservationSystem::is_algebraic(enum SPECIES species) const {
  const auto *conservation = this;
  return chemistry_species_is_intermediate_equilibrium(static_cast<SPECIES>(species)) ||
         conservation->is_dependent[species];
}

bool ChemConservationSystem::project(double *abundance) const {
  const auto *conservation = this;
  double remaining[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  int valid = 1;
  for(int constraint = 0; constraint < conservation->count; constraint++) {
    remaining[constraint] = conservation->total[constraint];
    for(int species = 0; species < N_SP; species++)
      if(!conservation->is_dependent[species])
        remaining[constraint] -= conservation->coefficient[constraint][species] * abundance[species];
  }

  for(int dependent = 0; dependent < conservation->count; dependent++) {
    double value = 0.0;
    for(int constraint = 0; constraint < conservation->count; constraint++)
      value += conservation->dependent_inverse[constraint * conservation->count + dependent] * remaining[constraint];
    abundance[conservation->dependent_species[dependent]] = value;
    if(!isfinite(value) || value < 0.0)
      valid = 0;
  }

  return valid;
}

int chemistry_abundances_are_nonnegative(const double *abundance)
{
  for(int species = 0; species < N_SP; species++)
    if(!isfinite(abundance[species]) || abundance[species] < 0.0)
      return 0;
  return 1;
}

#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
static const enum SPECIES full_implicit_equilibrium_initialization_species[] = {
  sHM,
#ifdef CARBON_OXYGEN_CHEMISTRY
  sCOII, sCH, sCH2, sOII, sOHII, sH2OII, sH3OII, sO2II,
#endif
};

int chemistry_is_full_implicit_initialized_intermediate(
    enum SPECIES species)
{
  for(size_t i = 0;
      i < sizeof(full_implicit_equilibrium_initialization_species) /
          sizeof(full_implicit_equilibrium_initialization_species[0]);
      i++)
    if(full_implicit_equilibrium_initialization_species[i] == species)
      return 1;
  return 0;
}

static int chemistry_initialize_full_implicit_intermediates(
    double nH, double input_electron, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation)
{
  const int count = (int)(sizeof(full_implicit_equilibrium_initialization_species) /
                          sizeof(full_implicit_equilibrium_initialization_species[0]));
  double base[N_SP], equilibrium[N_SP];
  memcpy(base, abundance, sizeof(base));

  /* The external intermediate history is deliberately ignored.  Return its
   * nuclei and charge to the locally selected conservation reservoirs first. */
  for(int i = 0; i < count; i++)
    base[full_implicit_equilibrium_initialization_species[i]] = 0.0;
  if(!(conservation)->project(base))
    return 0;
  for(int dependent = 0; dependent < conservation->count; dependent++)
    if(base[conservation->dependent_species[dependent]] < 0.0)
      return 0;

  memcpy(equilibrium, base, sizeof(equilibrium));
  /* The input electron abundance contains the charge carried by the
   * intermediates before their history is discarded above.  Use that
   * physical electron abundance for the equilibrium reconstruction; the
   * subsequent conservation projection recomputes the exactly neutral value
   * from the newly reconstructed intermediates. */
  equilibrium[sELEC] = input_electron;
  int equilibrium_inside_conservation_domain =
      chemistry_solve_intermediate_equilibrium_list(
          CHEM_REACTION_IMPLICIT, nH, equilibrium, rates,
          full_implicit_equilibrium_initialization_species, (size_t)count) &&
      (conservation)->project(equilibrium);
  for(int dependent = 0;
      equilibrium_inside_conservation_domain &&
      dependent < conservation->count; dependent++) {
    enum SPECIES species = conservation->dependent_species[dependent];
    if(equilibrium[species] < 0.0)
      equilibrium_inside_conservation_domain = 0;
  }

  /* Equilibrium values are only a Newton starting candidate.  If that
   * candidate is outside the elemental or charge-conservation domain, start
   * from the already projected non-equilibrium base state instead. */
  memcpy(abundance,
         equilibrium_inside_conservation_domain ? equilibrium : base,
         sizeof(equilibrium));
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_record_equilibrium_initialization(
      equilibrium_inside_conservation_domain ? 1.0 : 0.0,
      input_electron, abundance[sELEC]);
#endif
  return 1;
}

#endif

bool ChemConservationSystem::update_algebraic(double nH, double *abundance, const ChemRates *rates) const {
  const auto *conservation = this;
  if(intermediate_equilibrium_species_count == 0) {
    if(!(conservation)->project(abundance))
      return 0;
    return chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates) &&
           chemistry_abundances_are_nonnegative(abundance);
  }

  /* Conservation-dependent abundances, especially the electron abundance,
   * also depend on charged equilibrium intermediates.  Iterate the two
   * algebraic systems together; half damping stabilizes the common
   * ion-recombination fixed point y_ion proportional to 1 / y_e. */
  for(int iteration = 0; iteration < 128; iteration++) {
    double previous_dependent[CHEM_MAX_CONSERVATION_CONSTRAINTS];
    for(int dependent = 0; dependent < conservation->count; dependent++)
      previous_dependent[dependent] = abundance[conservation->dependent_species[dependent]];

    if(!(conservation)->project(abundance))
      return 0;

    double maximum_change = 0.0;
    for(int dependent = 0; dependent < conservation->count; dependent++) {
      enum SPECIES species = conservation->dependent_species[dependent];
      double target = abundance[species];
      if(target < 0.0)
        return 0;
      double scale = dmax(dmax(fabs(previous_dependent[dependent]), fabs(target)), 1.e-30);
      maximum_change = dmax(maximum_change,
                            fabs(target - previous_dependent[dependent]) / scale);
      abundance[species] = 0.5 * (previous_dependent[dependent] + target);
    }

    if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates))
      return 0;

    if(maximum_change < 1.e-10) {
      if(!(conservation)->project(abundance))
        return 0;
      if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates))
        return 0;
      return chemistry_abundances_are_nonnegative(abundance);
    }
  }

  return 0;
}

bool ChemConservationSystem::reduce_jacobian(double nH, const double *abundance, const ChemRates *rates, const enum SPECIES *solver_species, int n_solver_species, double drfdy[N_SP][N_SP]) const {
  const auto *conservation = this;
  int n_intermediate = (int)intermediate_equilibrium_species_count;
  int n_algebraic = n_intermediate + conservation->count;
  if(n_algebraic == 0)
    return 1;

  enum SPECIES algebraic_species[N_SP];
  for(int intermediate = 0; intermediate < n_intermediate; intermediate++)
    algebraic_species[intermediate] = intermediate_equilibrium_species[intermediate];
  for(int dependent = 0; dependent < conservation->count; dependent++)
    algebraic_species[n_intermediate + dependent] = conservation->dependent_species[dependent];

  double algebraic_response[N_SP * N_SP];
  memset(algebraic_response, 0, sizeof(algebraic_response));

#ifdef CHEMISTRY_PROFILE
  double response_timer = profile_start();
#endif
  for(int solver_j = 0; solver_j < n_solver_species; solver_j++) {
    enum SPECIES column_species = solver_species[solver_j];
    for(int dependent = 0; dependent < conservation->count; dependent++)
      for(int constraint = 0; constraint < conservation->count; constraint++)
        algebraic_response[solver_j * n_algebraic + n_intermediate + dependent] +=
            conservation->dependent_inverse[constraint * conservation->count + dependent] *
            conservation->coefficient[constraint][column_species];
  }
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_JACOBIAN_RESPONSE, response_timer);
#endif

  if(n_intermediate > 0) {
    double intermediate_jacobian[N_SP * N_SP];
    double intermediate_response[N_SP * N_SP];
    double conservation_intermediate_response[N_SP * CHEM_MAX_CONSERVATION_CONSTRAINTS];
    int intermediate_pivots[N_SP];

    memset(conservation_intermediate_response, 0,
           sizeof(conservation_intermediate_response));
    for(int intermediate = 0; intermediate < n_intermediate; intermediate++) {
      enum SPECIES intermediate_species = intermediate_equilibrium_species[intermediate];
      for(int dependent = 0; dependent < conservation->count; dependent++)
        for(int constraint = 0; constraint < conservation->count; constraint++)
          conservation_intermediate_response[intermediate * conservation->count + dependent] +=
              conservation->dependent_inverse[constraint * conservation->count + dependent] *
              conservation->coefficient[constraint][intermediate_species];
    }

    for(int intermediate_j = 0; intermediate_j < n_intermediate; intermediate_j++) {
      enum SPECIES column_species = intermediate_equilibrium_species[intermediate_j];
      for(int intermediate_i = 0; intermediate_i < n_intermediate; intermediate_i++) {
        enum SPECIES row_species = intermediate_equilibrium_species[intermediate_i];
        double reduced_derivative = drfdy[column_species][row_species];
        for(int dependent = 0; dependent < conservation->count; dependent++) {
          enum SPECIES dependent_species = conservation->dependent_species[dependent];
          reduced_derivative -= drfdy[dependent_species][row_species] *
              conservation_intermediate_response[intermediate_j * conservation->count + dependent];
        }
        intermediate_jacobian[intermediate_j * n_intermediate + intermediate_i] = reduced_derivative;
      }
    }

    for(int solver_j = 0; solver_j < n_solver_species; solver_j++) {
      enum SPECIES column_species = solver_species[solver_j];
      for(int intermediate_i = 0; intermediate_i < n_intermediate; intermediate_i++) {
        enum SPECIES row_species = intermediate_equilibrium_species[intermediate_i];
        double response = drfdy[column_species][row_species];
        for(int dependent = 0; dependent < conservation->count; dependent++) {
          enum SPECIES dependent_species = conservation->dependent_species[dependent];
          response -= drfdy[dependent_species][row_species] *
              algebraic_response[solver_j * n_algebraic + n_intermediate + dependent];
        }
        intermediate_response[solver_j * n_intermediate + intermediate_i] = response;
      }
    }

#ifdef CHEMISTRY_PROFILE
    double intermediate_solve_timer = profile_start();
#endif
    int info;
    dgesv_(&n_intermediate, &n_solver_species, intermediate_jacobian, &n_intermediate, intermediate_pivots, intermediate_response, &n_intermediate, &info);
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_JACOBIAN_INTERMEDIATE_SOLVE,
                intermediate_solve_timer);
#endif
    if(info != 0)
      return 0;

    for(int solver_j = 0; solver_j < n_solver_species; solver_j++)
      for(int intermediate = 0; intermediate < n_intermediate; intermediate++)
        algebraic_response[solver_j * n_algebraic + intermediate] = intermediate_response[solver_j * n_intermediate + intermediate];

    for(int solver_j = 0; solver_j < n_solver_species; solver_j++)
      for(int dependent = 0; dependent < conservation->count; dependent++)
        for(int intermediate = 0; intermediate < n_intermediate; intermediate++)
          algebraic_response[solver_j * n_algebraic + n_intermediate + dependent] -=
              conservation_intermediate_response[intermediate * conservation->count + dependent] *
              algebraic_response[solver_j * n_algebraic + intermediate];
  }

#ifdef CHEMISTRY_PROFILE
  double schur_timer = profile_start();
#endif
  unsigned char active_algebraic_count[N_SP];
  unsigned char active_algebraic[N_SP][N_SP];
  for(int solver_i = 0; solver_i < n_solver_species; solver_i++) {
    enum SPECIES row_species = solver_species[solver_i];
    active_algebraic_count[solver_i] = 0;
    for(int algebraic_k = 0; algebraic_k < n_algebraic; algebraic_k++) {
      if(drfdy[algebraic_species[algebraic_k]][row_species] != 0.0)
        active_algebraic[solver_i][active_algebraic_count[solver_i]++] = algebraic_k;
    }
  }

  for(int solver_j = 0; solver_j < n_solver_species; solver_j++) {
    enum SPECIES column_species = solver_species[solver_j];
    for(int solver_i = 0; solver_i < n_solver_species; solver_i++) {
      enum SPECIES row_species = solver_species[solver_i];
      for(int active = 0; active < active_algebraic_count[solver_i]; active++) {
        int algebraic_k = active_algebraic[solver_i][active];
        drfdy[column_species][row_species] -= drfdy[algebraic_species[algebraic_k]][row_species] * algebraic_response[solver_j * n_algebraic + algebraic_k];
      }
    }
  }
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_JACOBIAN_SCHUR, schur_timer);
#endif

  return 1;
}

#ifdef CHEM_EXPLICIT_SEQUENTIAL
void chemistry::ExplicitSolver::update_species(enum SPECIES species, double nH, double dt, const double *initial_abundance, double *abundance, const ChemRates *rates)
{
  double production;
  double destruction;
  chemistry_network_species_production_destruction(rates->network, species,
      nH, abundance, rates->forward, rates->inverse, &production,
      &destruction);
  abundance[species] = (initial_abundance[species] + production * dt) /
                       (1.0 + destruction * dt);
  if(abundance[species] < 0.0)
    abundance[species] = 0.0;
}
#else

int chemistry::ExplicitSolver::update_algebraic_species(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation)
{
  for(int iteration = 0; iteration < 128; iteration++) {
    double previous_dependent[CHEM_MAX_CONSERVATION_CONSTRAINTS];
    for(int dependent = 0; dependent < conservation->count; dependent++)
      previous_dependent[dependent] =
          abundance[conservation->dependent_species[dependent]];

    if(!(conservation)->project(abundance))
      return 0;

    double maximum_change = 0.0;
    for(int dependent = 0; dependent < conservation->count; dependent++) {
      enum SPECIES species = conservation->dependent_species[dependent];
      double target = abundance[species];
      if(target < -1.e-14 * dmax(conservation->total[dependent], 1.e-30))
        return 0;
      target = dmax(target, 0.0);
      double scale = dmax(dmax(fabs(previous_dependent[dependent]),
                               fabs(target)), 1.e-30);
      maximum_change = dmax(maximum_change,
                            fabs(target - previous_dependent[dependent]) /
                                scale);
      abundance[species] =
          0.5 * (previous_dependent[dependent] + target);
    }

    if(!chemistry_update_intermediate_equilibrium_species(
            CHEM_REACTION_EXPLICIT, nH, abundance, rates))
      return 0;

    if(maximum_change < 1.e-10) {
      if(!(conservation)->project(abundance))
        return 0;
      for(int dependent = 0; dependent < conservation->count; dependent++) {
        enum SPECIES species = conservation->dependent_species[dependent];
        if(abundance[species] < -1.e-14 *
               dmax(conservation->total[dependent], 1.e-30))
          return 0;
        abundance[species] = dmax(abundance[species], 0.0);
      }
      if(!chemistry_update_intermediate_equilibrium_species(
             CHEM_REACTION_EXPLICIT, nH, abundance, rates))
        return 0;
      return (conservation)->project(abundance);
    }
  }

  return 0;
}

/* One linearly implicit Euler step.  All reaction rates and derivatives are
 * evaluated at the same old state, then the coupled species correction is
 * obtained from one dense solve:
 *
 *   (I - dt dR/dy) delta_y = dt R(y_old).
 *
 * Conservation reservoirs and equilibrium intermediates are algebraic, so
 * they are eliminated from the matrix and reconstructed after the solve. */
int chemistry::ExplicitSolver::update_simultaneously(
    double nH, double metallicity, double dt, double *abundance,
    const ChemRates *rates)
{
  ChemConservationSystem conservation;
  if(!(&conservation)->build(metallicity, abundance) ||
     !update_algebraic_species(
         nH, abundance, rates, &conservation))
    return 0;

  enum SPECIES solver_species[N_SP];
  int solver_species_count = 0;
  for(int species = 0; species < N_SP; species++)
    if(!(&conservation)->is_algebraic(static_cast<SPECIES>(species)))
      solver_species[solver_species_count++] = static_cast<SPECIES>(species);

  double reaction_rhs[N_SP];
  double rhs_jacobian[N_SP][N_SP];
  chemistry_network_rhs(rates->network, CHEM_REACTION_EXPLICIT, nH,
                        abundance, rates->forward, rates->inverse,
                        reaction_rhs);
  chemistry_network_jacobian(rates->network, CHEM_REACTION_EXPLICIT, nH,
                             abundance, rates->forward, rates->inverse,
                             rhs_jacobian);
  if(!(&conservation)->reduce_jacobian(nH, abundance, rates, solver_species, solver_species_count, rhs_jacobian))
    return 0;

  double matrix[N_SP * N_SP];
  double correction[N_SP];
  int pivots[N_SP];
  int one = 1;
  int info;
  for(int column = 0; column < solver_species_count; column++) {
    enum SPECIES column_species = solver_species[column];
    for(int row = 0; row < solver_species_count; row++) {
      enum SPECIES row_species = solver_species[row];
      matrix[column * solver_species_count + row] =
          (column == row ? 1.0 : 0.0) -
          dt * rhs_jacobian[column_species][row_species];
    }
  }
  for(int row = 0; row < solver_species_count; row++)
    correction[row] = dt * reaction_rhs[solver_species[row]];

  dgesv_(&solver_species_count, &one, matrix, &solver_species_count, pivots,
         correction, &solver_species_count, &info);
  if(info != 0)
    return 0;

  /* Do not let an under-resolved trace species reduce the complete network
   * step to machine zero. Keep an already negligible positive abundance when
   * its linear correction crosses the positivity boundary; finite reservoirs
   * still participate in the common line search below. */
  for(int row = 0; row < solver_species_count; row++) {
    enum SPECIES species = solver_species[row];
    if(abundance[species] < 1.e-30 &&
       correction[row] < -abundance[species])
      correction[row] = 0.0;
  }

  /* A single coupled solve can overshoot a positivity boundary.  Shorten the
   * complete vector step, preserving its direction and simultaneous nature,
   * until both evolved and reconstructed algebraic abundances are valid. */
  for(int line_search = 0; line_search < 48; line_search++) {
    double fraction = ldexp(1.0, -line_search);
    double trial[N_SP];
    memcpy(trial, abundance, sizeof(trial));
    int valid = 1;
    for(int row = 0; row < solver_species_count; row++) {
      enum SPECIES species = solver_species[row];
      trial[species] += fraction * correction[row];
      if(!isfinite(trial[species]) || trial[species] < 0.0) {
        valid = 0;
        break;
      }
    }
    if(valid && update_algebraic_species(
                    nH, trial, rates, &conservation)) {
      memcpy(abundance, trial, sizeof(trial));
      return 1;
    }
  }

  return 0;
}
#endif

ChemistryStatus chemistry::ExplicitSolver::solve(const ChemEnvironment &env, const ChemStepData &cell_step,
    ChemState &input, ChemWorkspace &output, double dt) const
{
  const auto *environment = &env;
  const auto *step = &cell_step;
  auto *state = &input;
  auto *work = &output;
  const double nH = XH * environment->rho_cgs / PROTONMASS;
  const double temperature = work->thermodynamics.temperature;
  double abund_old[N_SP];
#ifdef CHEM_EXPLICIT_SEQUENTIAL
  double initial_abundance[N_SP];
#endif
  double abundance[N_SP];

  chemistry::ReactionRateEvaluator(*environment, *state, *work).evaluate(*step);

  for(int species = 0; species < N_SP; species++) {
    abund_old[species] = state->abund[species];
#ifdef CHEM_EXPLICIT_SEQUENTIAL
    initial_abundance[species] = state->abund[species] / XH;
#endif
    abundance[species] = state->abund[species] / XH;
  }

  if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
    return CHEMISTRY_STATUS_NOT_CONVERGED;

#ifdef CHEM_EXPLICIT_SEQUENTIAL
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_intermediate_equilibrium(static_cast<SPECIES>(species))) {
      if(!chemistry_update_intermediate_equilibrium_species(
              CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
        return CHEMISTRY_STATUS_NOT_CONVERGED;
    } else {
      update_species(
          static_cast<SPECIES>(species), nH, dt, initial_abundance, abundance, &work->rates);
    }
  }
#else
  if(!update_simultaneously(
         nH, environment->metallicity, dt, abundance, &work->rates))
    return CHEMISTRY_STATUS_NOT_CONVERGED;
#endif

#ifdef CHEM_EXPLICIT_SEQUENTIAL
  if(!chemistry_update_intermediate_equilibrium_species(
         CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
    return CHEMISTRY_STATUS_NOT_CONVERGED;
#endif

  for(int species = 0; species < N_SP; species++)
    state->abund[species] = abundance[species] * XH;

#ifdef CHEM_EXPLICIT_SEQUENTIAL
  for(int iteration = 0; iteration < 16; iteration++) {
    double previous_intermediate_equilibrium[sizeof(intermediate_equilibrium_species) / sizeof(intermediate_equilibrium_species[0])];
    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
      previous_intermediate_equilibrium[i] = state->abund[intermediate_equilibrium_species[i]];

    chemistry_explicit_enforce_constraints(environment->metallicity,
                                           state->abund);
    for(int species = 0; species < N_SP; species++)
      abundance[species] = state->abund[species] / XH;
    if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
      return CHEMISTRY_STATUS_NOT_CONVERGED;
    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
      state->abund[intermediate_equilibrium_species[i]] = abundance[intermediate_equilibrium_species[i]] * XH;

    double max_change = 0.0;
    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++) {
      double value = state->abund[intermediate_equilibrium_species[i]];
      double change = fabs(value - previous_intermediate_equilibrium[i]) / (fabs(value) + 1.e-80);
      if(change > max_change)
        max_change = change;
    }
    if(max_change < 1.e-12)
      break;
  }
#endif

  double charge_scale = 0.0;
  for(int species = 0; species < N_SP; species++) {
    if(!isfinite(state->abund[species]) || state->abund[species] < 0.0)
      return CHEMISTRY_STATUS_NOT_CONVERGED;
    charge_scale += fabs((double)chemistry_species_definition(static_cast<SPECIES>(species))->charge *
                         state->abund[species]);
  }
  if(fabs(chemistry_net_charge(state->abund)) >
     256.0 * DBL_EPSILON * dmax(charge_scale, 1.e-30))
    return CHEMISTRY_STATUS_NOT_CONVERGED;

  if(isnan(state->abund[sHI])) {
    printf("[ChemID%d]HI is nan (dt=%g, rho_gas=%g, T_gas=%g)\n",
           environment->ID, dt, environment->rho_cgs, temperature);
    return CHEMISTRY_STATUS_NOT_CONVERGED;
  }
  if(isnan(temperature)) {
    printf("[ChemID%d]line=%d T_gas is nan (dt=%g, rho_gas=%g, T_gas=%g)\n",
           environment->ID, __LINE__, dt, environment->rho_cgs, temperature);
    return CHEMISTRY_STATUS_NOT_CONVERGED;
  }

  chemistry::ReactionRateEvaluator(*environment, *state, *work).evaluate_timescale();
  ChemistryStatus cooling_status =
      chemistry::CoolingEvaluator(*environment, *state, *work).chemical(dt, abund_old);
  if(cooling_status != CHEMISTRY_STATUS_OK)
    return cooling_status;
  return CHEMISTRY_STATUS_OK;
}

#define ITER_RESET 40

#ifdef SAHA
ChemistryStatus chemistry::SahaSolver::solve(const ChemEnvironment &env, ChemState &input,
    ChemWorkspace &output, double dt) const
{
  const auto *environment = &env;
  auto *state = &input;
  auto *workspace = &output;
  const double n_nuclei = environment->rho_cgs / PROTONMASS;
  const double nH = XH * n_nuclei;
  const double T_gas = workspace->thermodynamics.temperature;
  const double y_H2_init = state->abund[sH2I] / XH;
  const double y_Hp_init = state->abund[sHII] / XH;
  const double y_Hm_init = state->abund[sHM] / XH;
  double y_input[N_SP];

  y_input[sELEC] = state->abund[sELEC] / XH;
  y_input[sHI] = state->abund[sHI] / XH;
  y_input[sHII] = state->abund[sHII] / XH;
  y_input[sH2I] = state->abund[sH2I] / XH;
  y_input[sHM] = state->abund[sHM] / XH;

  workspace->evaluation.err_chem = chem_saha_5sp(
      environment->run->rate_tables, nH, T_gas, y_input);

  state->abund[sHI] = y_input[sHI] * XH;
  state->abund[sHII] = y_input[sHII] * XH;
  state->abund[sH2I] = y_input[sH2I] * XH;
  state->abund[sHM] = y_input[sHM] * XH;
  state->abund[sELEC] = y_input[sELEC] * XH;

  double edot_H2form = 0.;
  double edot_ionize = 0.;
  double edot_HMform = 0.;
  if(dt >= 1.e-30){
    edot_H2form = 4.48 * eV_to_erg * (y_input[sH2I] - y_H2_init) * nH / dt;
    edot_ionize = -13.598 * eV_to_erg * (y_input[sHII] - y_Hp_init) * nH / dt;
    edot_HMform = 0.78 * eV_to_erg * (y_input[sHM] - y_Hm_init) * nH / dt;
  }

  workspace->evaluation.rhoU_dot_chem = edot_H2form + edot_ionize + edot_HMform;
#ifdef DEBUG_CHEMISTRY_SC
  workspace->diagnostics.edot_ciHI = edot_ionize;
  workspace->diagnostics.edot_H2form_gas = edot_H2form;
#endif

  return CHEMISTRY_STATUS_OK;
}
#endif

void chemistry_correct_full_implicit_intermediate_guess(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation){
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
  double uncorrected_abundance[N_SP];
  memcpy(uncorrected_abundance, abundance, sizeof(uncorrected_abundance));

  int valid = chemistry_solve_intermediate_equilibrium_list(
      CHEM_REACTION_IMPLICIT, nH, abundance, rates,
      full_implicit_equilibrium_initialization_species,
      sizeof(full_implicit_equilibrium_initialization_species) /
          sizeof(full_implicit_equilibrium_initialization_species[0]));
  if(valid)
    valid = (conservation)->project(abundance);
  for(int dependent = 0; valid && dependent < conservation->count; dependent++)
    if(abundance[conservation->dependent_species[dependent]] < 0.0)
      valid = 0;

  if(!valid)
    memcpy(abundance, uncorrected_abundance, sizeof(uncorrected_abundance));
#else
  (void)nH;
  (void)abundance;
  (void)rates;
  (void)conservation;
#endif
}

int chemistry_prepare_implicit_newton_problem(
    const ChemEnvironment *environment, const ChemState *state,
    const ChemRates *rates, double nH, ChemNewtonReference *reference,
    ChemNewtonSystem *problem){
  problem->initial_abundance = reference->initial_abundance;
  problem->cooling_initial_abundance = reference->cooling_initial_abundance;
  problem->initialization_floor = reference->initialization_floor;
  problem->newton_lower_bound = reference->newton_lower_bound;
  problem->upper_bound = reference->upper_bound;
  problem->solver_species_count = 0;

  /* Freeze the physical old state before constructing any Newton guess.
   * A valid zero (or sub-floor abundance) is not missing initial data. */
  for(int species = 0; species < N_SP; species++){
    if(!isfinite(state->abund[species]) || state->abund[species] < 0.0){
      chemistry_profile_set_failure(1);
      return 0;
    }
    problem->abundance[species] = state->abund[species] / XH;
    reference->initial_abundance[species] = problem->abundance[species];
    reference->cooling_initial_abundance[species] = state->abund[species];
  }

  if(!chemistry_update_intermediate_equilibrium_species(
          CHEM_REACTION_IMPLICIT, nH, problem->abundance, rates)){
    chemistry_profile_set_failure(1);
    return 0;
  }

  /* Start from bounds that are valid for every compiled species. */
  for(int species = 0; species < N_SP; species++){
    reference->initialization_floor[species] = 0.0;

    /* Export floors do not define the physical old state or a Newton bound.
     * Newton may resolve any nonnegative abundance, including exact zero. */
    reference->newton_lower_bound[species] = 0.0;
    reference->upper_bound[species] = 2.0;
  }

  for(int species = 0; species < N_SP; species++){
    const ChemSpeciesInitialization *parameters =
        &chemistry_species_initialization[species];

    reference->initialization_floor[species] = parameters->floor;

    if(parameters->upper_bound_per_hydrogen > 0.0)
      reference->upper_bound[species] =
          parameters->upper_bound_per_hydrogen;
  }

  for(size_t element_index = 0;
      element_index < chemistry_element_conservation_count; element_index++){
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[element_index];

    if(element->initialization_upper_bound_factor == 0.0)
      continue;

    double upper_bound = element->total_per_hydrogen;

    if(element->metallicity_scaled)
      upper_bound *= environment->metallicity;

    upper_bound *= element->initialization_upper_bound_factor;
    upper_bound = dmax(upper_bound, 1.e-80);

    for(int species = 0; species < N_SP; species++){
      if(chemistry_species_definition(static_cast<SPECIES>(species))->nuclei[element->element] > 0)
        reference->upper_bound[species] =
            dmin(reference->upper_bound[species], upper_bound);
    }
  }

  /* Equilibrium reconstruction seeds the guess, never the physical old state. */
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
  const double equilibrium_initialization_electron = problem->abundance[sELEC];
  for(size_t index = 0;
       index < sizeof(full_implicit_equilibrium_initialization_species) /
                   sizeof(full_implicit_equilibrium_initialization_species[0]);
       index++)
    problem
        ->abundance[full_implicit_equilibrium_initialization_species[index]] =
        0.0;
#endif

  /* Select algebraic reservoirs and make the first state consistent. */
  repair_initial_abundances(problem->abundance);

  if(!(&problem->conservation)->build(environment->metallicity, problem->abundance)
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
      || !chemistry_initialize_full_implicit_intermediates(
             nH, equilibrium_initialization_electron, problem->abundance, rates,
             &problem->conservation)
#else
      || !(&problem->conservation)->update_algebraic(nH, problem->abundance, rates)
#endif
  ){
    chemistry_profile_set_failure(2);
    return 0;
  }

  /* Freeze only the reduced layout; the physical old state is already saved. */
  for(int species = 0; species < N_SP; species++){
    if(!(&problem->conservation)->is_algebraic(static_cast<SPECIES>(species)))
      problem->solver_species[problem->solver_species_count++] = static_cast<SPECIES>(species);
  }

  return 1;
}

int chemistry_initialize_implicit_newton_guess(
    double nH, double dt, const double *newton_guess, const ChemRates *rates,
    ChemNewtonSystem *problem){
  double predictor_seed[N_SP];
  memcpy(predictor_seed, problem->abundance, sizeof(predictor_seed));
  if(!newton_guess && !chemistry_newton_build_explicit_guess(
                           nH, dt, predictor_seed,
                           problem->abundance, problem->newton_lower_bound,
                           problem->upper_bound, rates, &problem->conservation))
    return 0;

  if(newton_guess){
    for(int species = 0; species < N_SP; species++)
      if(!(&problem->conservation)->is_algebraic(static_cast<SPECIES>(species)))
        problem->abundance[species] = newton_guess[species] / XH;
    if(!(&problem->conservation)->update_algebraic(nH, problem->abundance, rates))
      return 0;
  }

  chemistry_correct_full_implicit_intermediate_guess(
      nH, problem->abundance, rates, &problem->conservation);
  for(int species = 0; species < N_SP; species++)
    problem->total_step[species] =
        problem->abundance[species] - problem->initial_abundance[species];

  return 1;
}

static int chemistry_build_numerical_newton_jacobian(
    double nH, const double *abundance, const double *reaction_rhs,
    const ChemRates *rates, const ChemConservationSystem *conservation,
    const enum SPECIES *solver_species, int solver_species_count,
    double rhs_jacobian[N_SP][N_SP]){
  const double relative_perturbation = 1.0e-4;
  const double zero_abundance_perturbation = 1.0e-10;
  const int maximum_perturbation_attempts = 20;

  for(int solver_column = 0; solver_column < solver_species_count;
       solver_column++){
    const int perturbed_species = solver_species[solver_column];
    double perturbation =
        fabs(abundance[perturbed_species]) < 1.0e-100
            ? zero_abundance_perturbation
            : relative_perturbation * abundance[perturbed_species];

    double forward_rhs[N_SP];
    double backward_rhs[N_SP];
    int forward_valid = 0;
    int backward_valid = 0;

    for(int attempt = 0; attempt < maximum_perturbation_attempts &&
                          !forward_valid && !backward_valid;
         attempt++){
      double perturbed_abundance[N_SP];

      memcpy(perturbed_abundance, abundance, sizeof(perturbed_abundance));
      perturbed_abundance[perturbed_species] += 0.5 * perturbation;
      forward_valid = (conservation)->update_algebraic(nH, perturbed_abundance, rates);
      if(forward_valid){
        chemistry_network_rhs(rates->network, CHEM_REACTION_IMPLICIT, nH,
                              perturbed_abundance, rates->forward,
                              rates->inverse, forward_rhs);
      }

      memcpy(perturbed_abundance, abundance, sizeof(perturbed_abundance));
      perturbed_abundance[perturbed_species] -= 0.5 * perturbation;
      backward_valid = perturbed_abundance[perturbed_species] >= 0.0 &&
                       (conservation)->update_algebraic(nH, perturbed_abundance, rates);
      if(backward_valid){
        chemistry_network_rhs(rates->network, CHEM_REACTION_IMPLICIT, nH,
                              perturbed_abundance, rates->forward,
                              rates->inverse, backward_rhs);
      }

      if(!forward_valid && !backward_valid)
        perturbation *= 0.5;
    }

    if(!forward_valid && !backward_valid)
      return 0;

    for(int species = 0; species < N_SP; species++){
      double rhs_change;
      double rhs_scale;

      if(forward_valid && backward_valid){
        rhs_jacobian[perturbed_species][species] =
            (forward_rhs[species] - backward_rhs[species]) / perturbation;
        rhs_scale =
            dmax(fabs(forward_rhs[species]), fabs(backward_rhs[species]));
        rhs_change = fabs(forward_rhs[species] - backward_rhs[species]);
      }else if(forward_valid){
        rhs_jacobian[perturbed_species][species] =
            2.0 * (forward_rhs[species] - reaction_rhs[species]) / perturbation;
        rhs_scale =
            dmax(fabs(forward_rhs[species]), fabs(reaction_rhs[species]));
        rhs_change = fabs(forward_rhs[species] - reaction_rhs[species]);
      }else {
        rhs_jacobian[perturbed_species][species] =
            2.0 * (reaction_rhs[species] - backward_rhs[species]) /
            perturbation;
        rhs_scale =
            dmax(fabs(reaction_rhs[species]), fabs(backward_rhs[species]));
        rhs_change = fabs(reaction_rhs[species] - backward_rhs[species]);
      }

      if(rhs_scale != 0.0 && rhs_change / rhs_scale < 1.0e-15)
        rhs_jacobian[perturbed_species][species] = 0.0;
    }
  }

  return 1;
}

int chemistry_build_newton_jacobian(
    double nH, const double *abundance, const double *reaction_rhs,
    const ChemRates *rates, const ChemConservationSystem *conservation,
    const enum SPECIES *solver_species, int solver_species_count,
    int rebuild_newton_matrix, int *force_numerical_jacobian,
    double rhs_jacobian[N_SP][N_SP]){
  memset(rhs_jacobian, 0, sizeof(double) * N_SP * N_SP);

  double jacobian_timer = profile_start();

  if(rebuild_newton_matrix){
    int used_analytic_jacobian = 0;

#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
    used_analytic_jacobian = !*force_numerical_jacobian;
    if(used_analytic_jacobian){
      double solver_columns_timer = profile_start();
      chemistry_network_jacobian(rates->network, CHEM_REACTION_IMPLICIT, nH,
                                 abundance, rates->forward, rates->inverse,
                                 rhs_jacobian);
      profile_end(CHEM_PROFILE_JACOBIAN_SOLVER_COLUMNS,
                  solver_columns_timer);

      if(!(conservation)->reduce_jacobian(nH, abundance, rates, solver_species, solver_species_count, rhs_jacobian)){
        *force_numerical_jacobian = 1;
        used_analytic_jacobian = 0;
      }
    }
#else
    (void)force_numerical_jacobian;
#endif

    if(!used_analytic_jacobian &&
        !chemistry_build_numerical_newton_jacobian(
            nH, abundance, reaction_rhs, rates, conservation, solver_species,
            solver_species_count, rhs_jacobian)){
      return 0;
    }
  }

  profile_end(CHEM_PROFILE_NEWTON_JACOBIAN, jacobian_timer);
  return 1;
}

void chemistry_update_implicit_timescale(
    ChemWorkspace *work, double dt, const double *abundance,
    const double *initial_abundance,
    const ChemConservationSystem *conservation){
  double minimum_timescale = 1.0e20;

  for(int species = 0; species < N_SP; species++){
    if((conservation)->is_algebraic(static_cast<SPECIES>(species)))
      continue;
    if(abundance[species] != initial_abundance[species]){
      double species_timescale =
          fabs((abundance[species] + initial_abundance[species]) / 2.0 /
               (abundance[species] - initial_abundance[species])) *
          dt;
      if(species_timescale < minimum_timescale)
        minimum_timescale = species_timescale;
    }

    work->evaluation.t_chem = minimum_timescale;
  }
}

int chemistry_export_implicit_solution(
    const ChemEnvironment *environment, const ChemStepData *step,
    ChemState *state, ChemWorkspace *work, const ChemState *input_state,
    double dt, double nH, double temperature,
    const double *initialization_floor, const double *cooling_initial_abundance,
    double *abundance, const ChemConservationSystem *conservation){
  double converged_abundance[N_SP];
  memcpy(converged_abundance, abundance, sizeof(converged_abundance));
  for(int species = 0; species < N_SP; species++){
    if((conservation)->is_algebraic(static_cast<SPECIES>(species)))
      continue;
    if(abundance[species] < initialization_floor[species])
      abundance[species] = initialization_floor[species];
  }

  /* Floors apply only to independent species.  Reconstruct algebraic species
   * afterward so conservation and equilibrium remain internally consistent. */
  if(!(conservation)->update_algebraic(nH, abundance, &work->rates)){
    // Export floors are heuristic. They may exceed an elemental budget (for
    // example, a positive metal floor at zero metallicity). Retain the valid
    // converged solution when the floored candidate cannot conserve nuclei.
    memcpy(abundance, converged_abundance, sizeof(converged_abundance));
    if(!(conservation)->update_algebraic(nH, abundance, &work->rates)){
      *state = *input_state;
      return 0;
    }
  }

  if(!chemistry_abundances_are_nonnegative(abundance)){
    *state = *input_state;
    return 0;
  }

  for(int species = 0; species < N_SP; species++)
    state->abund[species] = abundance[species] * XH;

  ChemistryStatus cooling_status =
      chemistry::CoolingEvaluator(*environment, *state, *work).chemical(dt, cooling_initial_abundance);
  if(cooling_status != CHEMISTRY_STATUS_OK)
    return 0;

#ifdef CALC_CONTINUUM_OPACITY
  double number_density[N_SP];
  for(int species = 0; species < N_SP; species++)
    number_density[species] = abundance[species] * nH;

#ifdef UPDATE_L_CORE_IN_SOLVER
  double core_length = get_core_radius(environment->rho_cgs, temperature);
#else
  double core_length = step->jeans_length;
#endif
  double freefall_time = compute_freefall_time(environment->rho_cgs);
  solve_level_population_H_new(temperature, number_density, core_length,
                               freefall_time);
#else
  (void)step;
  (void)temperature;
#endif

  return 1;
}


static void repair_initial_abundances(double *yi){
  const double soft_normalization_tolerance = 1.e-3;

  /* Repair pathological trace-species guesses before measuring totals. */
  if(yi[sHM] > 1.e-3)
    yi[sHM] = 1.e-10;

  /* Soft normalization repairs grossly inconsistent input totals without
   * selecting the final conservation reservoir. Hydrogen is normalized for
   * either sign of the error; He and D are only reduced when overabundant. */
  for(size_t group = 0;
      group < chemistry_element_conservation_count; group++){
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[group];

    if(element->initial_normalization == CHEM_INITIAL_NORMALIZATION_NONE)
      continue;

    double target = element->total_per_hydrogen;
    double content = 0.0;

    for(int species = 0; species < N_SP; species++){
      const ChemSpeciesDefinition *definition =
          chemistry_species_definition(static_cast<SPECIES>(species));
      int belongs_to_group =
          definition->nuclei[element->element] > 0;

      if(element->element == CHEM_ELEMENT_H)
        belongs_to_group = belongs_to_group &&
            definition->nuclei[CHEM_ELEMENT_HE] == 0 &&
            definition->nuclei[CHEM_ELEMENT_D] == 0 &&
            definition->nuclei[CHEM_ELEMENT_C] == 0 &&
            definition->nuclei[CHEM_ELEMENT_O] == 0;

      if(belongs_to_group)
        content += definition->nuclei[element->element] * yi[species];
    }

    int normalize =
        element->initial_normalization ==
            CHEM_INITIAL_NORMALIZATION_EXCESS_ONLY
            ? content > target
            : fabs(target - content) > soft_normalization_tolerance;

    if(!normalize)
      continue;

    double scale = target / content;

    for(int species = 0; species < N_SP; species++){
      const ChemSpeciesDefinition *definition =
          chemistry_species_definition(static_cast<SPECIES>(species));
      int belongs_to_group =
          definition->nuclei[element->element] > 0;

      if(element->element == CHEM_ELEMENT_H)
        belongs_to_group = belongs_to_group &&
            definition->nuclei[CHEM_ELEMENT_HE] == 0 &&
            definition->nuclei[CHEM_ELEMENT_D] == 0 &&
            definition->nuclei[CHEM_ELEMENT_C] == 0 &&
            definition->nuclei[CHEM_ELEMENT_O] == 0;

      if(belongs_to_group)
        yi[species] *= scale;
    }
  }

  /* Repair charge imbalance. Positive excess is absorbed by electrons.
   * Negative excess first removes electrons, then ionizes hydrogen. */
  double negative_charge = 0.0;
  double positive_charge = 0.0;

  for(int species = 0; species < N_SP; species++){
    int charge = chemistry_species_definition(static_cast<SPECIES>(species))->charge;

    if(charge < 0)
      negative_charge -= charge * yi[species];
    else if(charge > 0)
      positive_charge += charge * yi[species];
  }

  if(negative_charge > 1.001 * positive_charge){
    double negative_excess = negative_charge - positive_charge;

    if(negative_excess < 0.9 * yi[sELEC]){
      yi[sELEC] -= negative_excess;
    }else{
      negative_excess -= 0.9 * yi[sELEC];
      yi[sELEC] *= 0.1;

      if(yi[sHI] > negative_excess){
        yi[sHI] -= negative_excess;
        yi[sHII] += negative_excess;
      }else{
        yi[sH2I] -= 0.5 * negative_excess;
        yi[sHII] += negative_excess;
      }
    }
  }else if(positive_charge > 1.00001 * negative_charge){
    double positive_excess = positive_charge - negative_charge;
    yi[sELEC] += positive_excess;
  }

  /* Reconstruct each element in a reservoir large enough to absorb the
   * residual without subtracting nearly equal trace abundances. */
  for(size_t group = 0;
      group < chemistry_element_conservation_count; group++){
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[group];

    if(element->initial_normalization == CHEM_INITIAL_NORMALIZATION_NONE)
      continue;

    double target = element->total_per_hydrogen;
    double reservoir_threshold = 0.1 * target;
    enum SPECIES reservoir = element->initial_fallback_reservoir;

    if(yi[element->reservoirs[0]] > reservoir_threshold)
      reservoir = element->reservoirs[0];
    else if(element->reservoir_count > 1 &&
            yi[element->reservoirs[1]] > reservoir_threshold)
      reservoir = element->reservoirs[1];

    double remaining = target;

    for(int species = 0; species < N_SP; species++){
      const ChemSpeciesDefinition *definition =
          chemistry_species_definition(static_cast<SPECIES>(species));
      int belongs_to_group =
          definition->nuclei[element->element] > 0;

      if(element->element == CHEM_ELEMENT_H)
        belongs_to_group = belongs_to_group &&
            definition->nuclei[CHEM_ELEMENT_HE] == 0 &&
            definition->nuclei[CHEM_ELEMENT_D] == 0 &&
            definition->nuclei[CHEM_ELEMENT_C] == 0 &&
            definition->nuclei[CHEM_ELEMENT_O] == 0;

      if(belongs_to_group && species != reservoir)
        remaining -= definition->nuclei[element->element] * yi[species];
    }

    int reservoir_nuclei =
        chemistry_species_definition(reservoir)
            ->nuclei[element->element];
    yi[reservoir] = remaining / reservoir_nuclei;
  }
}



#endif /* CHEMISTRY_SC */
