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
#include "chemistry_diagnostics.h"
#include "chemistry_models.h"
#include "chemistry_profile.h"
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
#include "chemistry_newton_reset.h"
#include "chemistry_solver.h"
#include "saha.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

#define CHEMISTRY_NEWTON_TRACE_ABUNDANCE 1.e-12

static void repair_initial_abundances(double *abundance);
static int chemistry_noneq_chem_implicit_trial(const ChemEnvironment *environment,
                                               const ChemStepData *step,
                                               ChemState *state,
                                               ChemWorkspace *work,
                                               double dt);
int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
#ifdef CHEM_NEWTON_LU_REUSE
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int dgetrs_(char *trans, int *n, int *nrhs, const double *a, int *lda,
            const int *ipiv, double *b, int *ldb, int *info);
#endif

#ifdef ENERGY_IMPLICIT
double give_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature, double *energy, double *residual);
double solve_bisection(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature_min, double residual_min, double temperature_max, double residual_max, double *temperature, double *energy, int max_iter);
static double commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature, double *energy, double *residual);
static double report_commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature, double *energy, double *residual);
static int solve_energy_by_chemistry_continuation(const ChemEnvironment *environment,
                                                  const ChemStepData *step,
                                                  const ChemState *initial_state,
                                                  const ChemWorkspace *initial_workspace,
                                                  ChemistrySolver chemistry_solver,
                                                  double dt, double temperature_initial,
                                                  double temperature_min,
                                                  double temperature_max,
                                                  ChemState *state,
                                                  ChemWorkspace *workspace,
                                                  double *temperature,
                                                  double *energy,
                                                  double *residual);
#endif

#ifdef ENERGY_IMPLICIT
#define MAX_ITER 500
#define BETA_ENERGY 1.0

static void restore_workspace(ChemWorkspace *destination,
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

int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt){
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
  const double T_max=3.e6, T_min=0.1;
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

  while(1){
    err_chem = give_g(environment, step, state, work, chemistry_solver, dt, T_new, &e_new, &g0);
    if(work->status != CHEMISTRY_STATUS_OK) {
      err = NAN;
      goto END;
    }

    if(err_chem < 1.e-2)
      break;

    if(isnan(dT) || isinf(dT)){
      dT = 0.1 * T_new;
      printf("(LINE=%d) dT=nan -> %g g=%g g0=%g denom=%g\n",
	     __LINE__, dT, g, g0, denom);
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
  double g_plus = 1.0e30, g_minus = -1.0e30;
  double T_plus, T_minus;
  int flag_brent = 0, sign;

  if(g0 > 0) g_plus = g0, T_plus = T_new, sign = 1;
  else g_minus = g0, T_minus = T_new, sign = -1;
#endif

  restore_workspace(work, &work_old);
  *state = state_old;

  T_old = T_new;
  T_new = T_old + dT;

  int flag_converged = 0;

  int j;
  for(j=iter;j<MAX_ITER;j++){
#ifdef FLAG_BRENT
    if( flag_brent == 1){
      restore_workspace(work, &work_old);
      *state = state_old;
      err = solve_bisection(environment, step, state, work, chemistry_solver, dt, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

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
      err_chem = give_g(environment, step, state, work, chemistry_solver, dt, T_new, &e_new, &g);
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
	printf("(LINE=%d) dT=nan -> %g g=%g g0=%g denom=%g\n",
	       __LINE__, dT, g, g0, denom);
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

      err_chem = give_g(environment, step, state, work, chemistry_solver, dt, T_new, &e_new, &g);
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
      err = solve_bisection(environment, step, state, work, chemistry_solver, dt, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

      if(err < 1.e-5) flag_converged = 1;
    }
  }
#endif

  final_chemistry_error = commit_g(environment, step, state, work,
                                   chemistry_solver, dt, T_new, &e_new, &g);
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
    if(solve_energy_by_chemistry_continuation(
           environment, step, &state_old, &work_old, chemistry_solver, dt,
           T_gas_init, T_min, T_max, &continuation_state,
           &continuation_workspace, &continuation_temperature,
           &continuation_energy, &continuation_residual)) {
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
    return 1;
  }

  if(err > 1.0e-4 || work->diagnostics_enabled){
    printf("%s err=%g id=%d T_old=%g T_new=%g\n",
	   err > 1.0e-4 ? "not converged(in energy)." : "converged(in energy diagnostic).",
	   err, environment->ID, T_gas_init, T_new);
    printf("/*CHEMISTRY_RERUN_V1*/diagnostic_id=%d;comoving_integration=%d;"
           "simulation_time=%.17g;scale_factor=%.17g;n_nuclei=%.17g;"
           "T_gas=%.17g;cmb_temperature=%.17g;metallicity=%.17g;"
           "dt=%.17g;gamma=%.17g;utherm_cgs=%.17g;",
           environment->ID, environment->global->comoving_integration,
           environment->global->simulation_time,
           environment->global->scale_factor, rho_gas / PROTONMASS,
           T_gas_init, environment->global->cmb_temperature,
           environment->metallicity, dt, state_old.gamma,
           state_old.utherm_cgs);

    for(int diagnostic_species = 0; diagnostic_species < sSTORED_END;
        diagnostic_species++)
      printf("yi[s%s]=%.17g;",
             chemistry_species_definition(diagnostic_species)->name,
             state_old.abund[diagnostic_species]);

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
#ifdef THERMAL_FEEDBACK
    printf("stellar_radiation_flux=%.17g;",
           environment->stellar_radiation_flux);
#endif
    printf("\n");
  }

  if(err > 1.0e-4){
    if(chemistry_solver == CHEMISTRY_SOLVER_IMPLICIT &&
       (!final_evaluation_available || final_chemistry_error >= 1.e-4)) {
      ChemState diagnostic_state = state_old;
      ChemWorkspace diagnostic_workspace = work_old;
      double diagnostic_energy, diagnostic_residual;
      report_commit_g(environment, step, &diagnostic_state,
                      &diagnostic_workspace, chemistry_solver, dt, T_new,
                      &diagnostic_energy, &diagnostic_residual);
    }
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_report_energy_warning();
#endif
    restore_workspace(work, &work_old);
    *state = state_old;
#ifdef CHEMISTRY_PROFILE
    profile_end(CHEM_PROFILE_ENERGY_SOLVER, profile_timer);
#endif
    return 1;
  }

#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_ENERGY_SOLVER, profile_timer);
#endif
  return 0;
}
static double evaluate_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g, int commit, int report_chemistry_failure, const double *newton_guess, int max_chemistry_iterations)
{
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  ChemWorkspace work_trial = *work;
  ChemState state_trial = *state;
  double e_old = state_trial.utherm_cgs;
  int chemistry_status = 1;

  work_trial.thermodynamics.temperature = T_gas_new;
#ifdef SAHA
  if(chemistry_solver == CHEMISTRY_SOLVER_SAHA)
    chemistry_saha(environment, &state_trial, &work_trial, dt);
  else
#endif
    {
#ifdef CHEMISTRY_PROFILE
      double chemistry_timer = profile_start();
#endif
      if(newton_guess)
        chemistry_status = chemistry_noneq_chem_implicit_with_guess(
            environment, step, &state_trial, &work_trial, dt,
            newton_guess, max_chemistry_iterations);
      else if(report_chemistry_failure)
        chemistry_status = chemistry_noneq_chem_implicit(
            environment, step, &state_trial, &work_trial, dt);
      else
        /* Energy root finding deliberately probes temperatures whose
         * chemistry solve may fail.  Preserve that failure as a rejected
         * residual point, but defer the detailed report until all energy
         * recovery paths have failed. */
        chemistry_status = chemistry_noneq_chem_implicit_trial(
            environment, step, &state_trial, &work_trial, dt);
#ifdef CHEMISTRY_PROFILE
      if(!chemistry_status)
        chemistry_status = -chemistry_profile_last_failure();
      profile_end(CHEM_PROFILE_IMPLICIT_CHEMISTRY, chemistry_timer);
#endif
    }
  ChemistryStatus cooling_status =
      evaluate_radiative_cooling(environment, step, &state_trial, &work_trial);
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
                                    e_old, *e_new, *g, chemistry_status);
  profile_end(CHEM_PROFILE_ENERGY_RESIDUAL, profile_timer);
#endif

  if(work->diagnostics_enabled)
    printf("resid(e_new-e_old-de)=%g e_new=%g old=%g de=%g Tgas=%g rhoU_dot_chem=%g cool=%g rho_gas=%g dt=%g xmu=%g gamma=%g\n"
	   , *g, *e_new, e_old, de, T_gas_new, work_trial.evaluation.rhoU_dot_chem, work_trial.evaluation.rhoU_dot_rad, environment->rho_cgs, dt, xmu_new, gamma_new);

  if(commit){
    *state = state_trial;
    *work = work_trial;
  }
  if((commit || newton_guess) && !chemistry_status)
    return 1.0;
  return work_trial.evaluation.err_chem;
}

double give_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt,
                    T_gas_new, e_new, g, 0, 0, NULL, 0);
}

static double give_g_with_guess(const ChemEnvironment *environment,
                                const ChemStepData *step, ChemState *state,
                                ChemWorkspace *work,
                                ChemistrySolver chemistry_solver,
                                double dt, double T_gas_new, double *e_new,
                                double *g, const double *newton_guess,
                                int max_chemistry_iterations){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt,
                    T_gas_new, e_new, g, 1, 0, newton_guess,
                    max_chemistry_iterations);
}

static double commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt,
                    T_gas_new, e_new, g, 1, 0, NULL, 0);
}

static double report_commit_g(const ChemEnvironment *environment,
                              const ChemStepData *step, ChemState *state,
                              ChemWorkspace *work,
                              ChemistrySolver chemistry_solver, double dt,
                              double T_gas_new, double *e_new, double *g){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt,
                    T_gas_new, e_new, g, 1, 1, NULL, 0);
}

typedef struct {
  double temperature;
  double energy;
  double residual;
  ChemState state;
  ChemWorkspace workspace;
} EnergyContinuationPoint;

static int evaluate_energy_continuation_point(
    const ChemEnvironment *environment, const ChemStepData *step,
    const ChemState *initial_state, const ChemWorkspace *initial_workspace,
    ChemistrySolver chemistry_solver, double dt, double temperature,
    const double *newton_guess, EnergyContinuationPoint *point)
{
  point->state = *initial_state;
  point->workspace = *initial_workspace;
  point->temperature = temperature;
  double chemistry_error;
  if(newton_guess)
    chemistry_error = give_g_with_guess(
        environment, step, &point->state, &point->workspace,
        chemistry_solver, dt, temperature, &point->energy, &point->residual,
        newton_guess, 1000);
  else
    chemistry_error = commit_g(
        environment, step, &point->state, &point->workspace,
        chemistry_solver, dt, temperature, &point->energy, &point->residual);
  if(chemistry_error >= 1.e-4 || !isfinite(point->energy) ||
     !isfinite(point->residual))
    return 0;
  point->state.utherm_cgs = point->energy;
  return 1;
}

static int refine_energy_continuation_bracket(
    const ChemEnvironment *environment, const ChemStepData *step,
    const ChemState *initial_state, const ChemWorkspace *initial_workspace,
    ChemistrySolver chemistry_solver, double dt,
    EnergyContinuationPoint lower, EnergyContinuationPoint upper,
    EnergyContinuationPoint *root)
{
  for(int iteration = 0; iteration < 80; iteration++) {
    double temperature = 0.5 * (lower.temperature + upper.temperature);
    const EnergyContinuationPoint *nearest =
        temperature - lower.temperature < upper.temperature - temperature ?
        &lower : &upper;
    EnergyContinuationPoint middle;
    if(!evaluate_energy_continuation_point(
           environment, step, initial_state, initial_workspace,
           chemistry_solver, dt, temperature, nearest->state.abund, &middle)) {
      const EnergyContinuationPoint *other = nearest == &lower ? &upper : &lower;
      if(!evaluate_energy_continuation_point(
             environment, step, initial_state, initial_workspace,
             chemistry_solver, dt, temperature, other->state.abund, &middle))
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

static int solve_energy_by_chemistry_continuation(
    const ChemEnvironment *environment, const ChemStepData *step,
    const ChemState *initial_state, const ChemWorkspace *initial_workspace,
    ChemistrySolver chemistry_solver, double dt, double temperature_initial,
    double temperature_min, double temperature_max, ChemState *state,
    ChemWorkspace *workspace, double *temperature, double *energy,
    double *residual)
{
  EnergyContinuationPoint start;
  if(!evaluate_energy_continuation_point(
         environment, step, initial_state, initial_workspace,
         chemistry_solver, dt, temperature_initial, NULL, &start))
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
        if(evaluate_energy_continuation_point(
               environment, step, initial_state, initial_workspace,
               chemistry_solver, dt, trial_temperature,
               previous.state.abund, &trial)) {
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
        if(refine_energy_continuation_bracket(
               environment, step, initial_state, initial_workspace,
               chemistry_solver, dt, lower, upper, &root)) {
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


double solve_bisection(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_min, double g_min, double T_max, double g_max, double *T_new0, double *e_new0, int max_iter){
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

    err = give_g(environment, step, state, work, chemistry_solver, dt, T_cand, &e_new, &val_new);
    if(work->status != CHEMISTRY_STATUS_OK)
      return 1.e30;

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

int chemistry_invert_dependent_matrix(ChemConservationSystem *conservation)
{
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
        chemistry_species_definition(species)->nuclei[element] *
        abundance[species];

    if(content > largest_content) {
      largest = species;
      largest_content = content;
    }
  }

  return largest;
}

#ifdef CARBON_OXYGEN_CHEMISTRY
static enum SPECIES chemistry_coupled_oxygen_reservoir(
    const double *abundance, enum SPECIES carbon_reservoir,
    enum SPECIES oxygen_reservoir)
{
  const enum SPECIES candidates[] = {
    sCI, sCII, sCO, sCO2, sOI, sOH, sH2O, sO2
  };
  const ChemSpeciesDefinition *carbon =
      chemistry_species_definition(carbon_reservoir);
  enum SPECIES best = oxygen_reservoir;
  double best_content = -1.0;

  for(size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); i++) {
    enum SPECIES candidate = candidates[i];
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(candidate);
    int determinant =
        carbon->nuclei[CHEM_ELEMENT_C] *
            definition->nuclei[CHEM_ELEMENT_O] -
        definition->nuclei[CHEM_ELEMENT_C] *
            carbon->nuclei[CHEM_ELEMENT_O];
    if(candidate == carbon_reservoir || determinant == 0)
      continue;
    double content =
        definition->nuclei[CHEM_ELEMENT_O] * abundance[candidate];
    if(content > best_content) {
      best_content = content;
      best = candidate;
    }
  }
  return best_content >= 0.0 ? best : oxygen_reservoir;
}
#endif

static int chemistry_build_conservation_system_mode(
    ChemConservationSystem *conservation, double metallicity,
    const double *abundance, int coupled_oxygen)
{
  memset(conservation, 0, sizeof(*conservation));

#if defined(CARBON_OXYGEN_CHEMISTRY) && \
    !defined(CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION)
  enum SPECIES carbon_dependent = sCI;
#else
  (void)coupled_oxygen;
#endif

  for(size_t element_index = 0;
      element_index < chemistry_element_conservation_count;
      element_index++) {
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[element_index];

#ifdef CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION
    if(element->element == CHEM_ELEMENT_C ||
       element->element == CHEM_ELEMENT_O)
      continue;
#endif

    enum SPECIES dependent = chemistry_largest_reservoir(
        abundance, element->element, element->reservoirs,
        element->newton_reservoir_count);

#if defined(CARBON_OXYGEN_CHEMISTRY) && \
    !defined(CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION)
    if(element->element == CHEM_ELEMENT_C)
      carbon_dependent = dependent;
    else if(element->element == CHEM_ELEMENT_O && coupled_oxygen)
      dependent = chemistry_coupled_oxygen_reservoir(
          abundance, carbon_dependent, dependent);
#endif

    int constraint = conservation->count;
    double total = element->total_per_hydrogen;

    if(element->metallicity_scaled)
      total *= metallicity;

    chemistry_add_conservation_constraint(
        conservation, dependent, total);

    for(int species = 0; species < N_SP; species++)
      conservation->coefficient[constraint][species] =
          chemistry_species_definition(species)->nuclei[element->element];
  }

  enum SPECIES charge_dependent = sELEC;
  double largest_negative_charge = -1.0;
  for(int species = 0; species < N_SP; species++) {
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(species);
    if(definition->charge >= 0 ||
       chemistry_species_is_intermediate_equilibrium(species) ||
       conservation->is_dependent[species])
      continue;
    double charge_content =
        fabs((double)definition->charge) * dmax(abundance[species], 0.0);
    if(charge_content > largest_negative_charge) {
      largest_negative_charge = charge_content;
      charge_dependent = species;
    }
  }
  int charge_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, charge_dependent, 0.0);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[charge_constraint][species] =
        chemistry_species_definition(species)->charge;

  return chemistry_invert_dependent_matrix(conservation);
}

int chemistry_build_conservation_system(ChemConservationSystem *conservation,
                                        double metallicity,
                                        const double *abundance)
{
  return chemistry_build_conservation_system_mode(
      conservation, metallicity, abundance, 0);
}

int chemistry_build_coupled_conservation_system(
    ChemConservationSystem *conservation, double metallicity,
    const double *abundance)
{
  return chemistry_build_conservation_system_mode(
      conservation, metallicity, abundance, 1);
}

int chemistry_species_is_algebraic(enum SPECIES species, const ChemConservationSystem *conservation)
{
  return chemistry_species_is_intermediate_equilibrium(species) ||
         conservation->is_dependent[species];
}

int chemistry_update_conservation_dependent_species(double *abundance, const ChemConservationSystem *conservation)
{
  double remaining[CHEM_MAX_CONSERVATION_CONSTRAINTS];
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
    if(!isfinite(value))
      return 0;
  }

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
  if(!chemistry_update_conservation_dependent_species(base, conservation))
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
  if(!chemistry_solve_intermediate_equilibrium_list(
         CHEM_REACTION_IMPLICIT, nH, equilibrium, rates,
         full_implicit_equilibrium_initialization_species, (size_t)count) ||
     !chemistry_update_conservation_dependent_species(equilibrium, conservation))
    return 0;

  /* This is only an initial guess for the full Newton system.  If its C/O
   * intermediate part lies outside the conservation domain, keep the HM
   * equilibrium guess and replace only the C/O intermediates by trace
   * abundances. */
  int equilibrium_inside_conservation_domain = 1;
  for(int dependent = 0; dependent < conservation->count; dependent++) {
    enum SPECIES species = conservation->dependent_species[dependent];
    if(equilibrium[species] < 0.0)
      equilibrium_inside_conservation_domain = 0;
  }
  memcpy(abundance, equilibrium, sizeof(equilibrium));
  if(!equilibrium_inside_conservation_domain) {
#ifdef CARBON_OXYGEN_CHEMISTRY
    const enum SPECIES carbon_oxygen_intermediates[] = {
      sCOII, sCH, sCH2, sOII, sOHII, sH2OII, sH3OII, sO2II
    };
    double carbon_total = 0.0;
    double oxygen_total = 0.0;
    for(int species = 0; species < N_SP; species++) {
      const ChemSpeciesDefinition *definition =
          chemistry_species_definition(species);
      carbon_total +=
          definition->nuclei[CHEM_ELEMENT_C] * base[species];
      oxygen_total +=
          definition->nuclei[CHEM_ELEMENT_O] * base[species];
    }
    const double floor_fraction = 1.e-12;
    for(size_t i = 0;
        i < sizeof(carbon_oxygen_intermediates) /
            sizeof(carbon_oxygen_intermediates[0]); i++) {
      enum SPECIES species = carbon_oxygen_intermediates[i];
      const ChemSpeciesDefinition *definition =
          chemistry_species_definition(species);
      double species_floor = INFINITY;
      if(definition->nuclei[CHEM_ELEMENT_C] > 0)
        species_floor = dmin(
            species_floor,
            floor_fraction * carbon_total /
                definition->nuclei[CHEM_ELEMENT_C]);
      if(definition->nuclei[CHEM_ELEMENT_O] > 0)
        species_floor = dmin(
            species_floor,
            floor_fraction * oxygen_total /
                definition->nuclei[CHEM_ELEMENT_O]);
      abundance[species] =
          isfinite(species_floor) ? species_floor : 0.0;
    }
#endif
  }
  if(!chemistry_update_conservation_dependent_species(abundance, conservation))
    return 0;
  for(int dependent = 0; dependent < conservation->count; dependent++)
    if(abundance[conservation->dependent_species[dependent]] < 0.0)
      return 0;
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_record_equilibrium_initialization(
      equilibrium_inside_conservation_domain ? 1.0 : 0.0,
      input_electron, abundance[sELEC]);
#endif
  return 1;
}

#endif

int chemistry_update_implicit_algebraic_species(double nH, double *abundance, const ChemRates *rates, const ChemConservationSystem *conservation)
{
  if(intermediate_equilibrium_species_count == 0) {
    if(!chemistry_update_conservation_dependent_species(abundance, conservation))
      return 0;
    return chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates);
  }

  /* Conservation-dependent abundances, especially the electron abundance,
   * also depend on charged equilibrium intermediates.  Iterate the two
   * algebraic systems together; half damping stabilizes the common
   * ion-recombination fixed point y_ion proportional to 1 / y_e. */
  for(int iteration = 0; iteration < 128; iteration++) {
    double previous_dependent[CHEM_MAX_CONSERVATION_CONSTRAINTS];
    for(int dependent = 0; dependent < conservation->count; dependent++)
      previous_dependent[dependent] = abundance[conservation->dependent_species[dependent]];

    if(!chemistry_update_conservation_dependent_species(abundance, conservation))
      return 0;

    double maximum_change = 0.0;
    for(int dependent = 0; dependent < conservation->count; dependent++) {
      enum SPECIES species = conservation->dependent_species[dependent];
      double target = abundance[species];
      double roundoff_tolerance = 0.0;
      for(int constraint = 0; constraint < conservation->count; constraint++)
        if(conservation->total[constraint] > 0.0 &&
           conservation->coefficient[constraint][species] != 0.0)
          roundoff_tolerance = dmax(roundoff_tolerance,
              64.0 * DBL_EPSILON * conservation->total[constraint]);
      if(target < 0.0 && target >= -roundoff_tolerance)
        target = abundance[species] = 0.0;
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
      if(!chemistry_update_conservation_dependent_species(abundance, conservation))
        return 0;
      if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates))
        return 0;
      return 1;
    }
  }

  return 0;
}

int chemistry_reduce_jacobian_for_algebraic_species(double nH, const double *abundance, const ChemRates *rates, const enum SPECIES *solver_species, int n_solver_species, const ChemConservationSystem *conservation, double drfdy[N_SP][N_SP])
{
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

static void chemistry_explicit_update_species(enum SPECIES species, double nH, double dt, const double *initial_abundance, double *abundance, const ChemRates *rates)
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

#define CHEM_EXPLICIT_SPECIES_ENUM(symbol, name, charge, h, he, d, c, o) \
  s##symbol,
static const enum SPECIES explicit_species_order[] = {
  CHEM_SPECIES_EXPLICIT_SEQUENCE(CHEM_EXPLICIT_SPECIES_ENUM)
};
#undef CHEM_EXPLICIT_SPECIES_ENUM

static void chemistry_explicit_enforce_constraints(const ChemEnvironment *environment, ChemState *state)
{
  double element_total[CHEM_ELEMENT_COUNT] = {0.0};
  double element_content[CHEM_ELEMENT_COUNT] = {0.0};
  double element_scale[CHEM_ELEMENT_COUNT] = {1.0, 1.0, 1.0, 1.0, 1.0};

  const size_t explicit_species_count =
      sizeof(explicit_species_order) / sizeof(explicit_species_order[0]);

  for(size_t index = 0; index < chemistry_element_conservation_count; index++){
    const ChemElementConservationDefinition *element =
        &chemistry_element_conservation[index];
    double total = XH * element->total_per_hydrogen;

    if(element->metallicity_scaled)
      total *= environment->metallicity;

    element_total[element->element] = total;
  }

  /* Accumulate the same normalization groups used by the original explicit
   * scheme. Hydrogen bound to D contributes to the H total but is scaled only
   * with the D family. Hydrogen bound to C/O is left to the metal families. */
  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(species);
    double abundance = state->abund[species];

    if(definition->nuclei[CHEM_ELEMENT_H] > 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0 &&
       definition->nuclei[CHEM_ELEMENT_O] == 0)
      element_content[CHEM_ELEMENT_H] +=
          definition->nuclei[CHEM_ELEMENT_H] * abundance;

    element_content[CHEM_ELEMENT_HE] +=
        definition->nuclei[CHEM_ELEMENT_HE] * abundance;
    element_content[CHEM_ELEMENT_D] +=
        definition->nuclei[CHEM_ELEMENT_D] * abundance;
    element_content[CHEM_ELEMENT_C] +=
        definition->nuclei[CHEM_ELEMENT_C] * abundance;
  }

  for(int element = CHEM_ELEMENT_H; element <= CHEM_ELEMENT_C; element++) {
    if(element_content[element] > 0.0)
      element_scale[element] =
          element_total[element] / element_content[element];
  }

  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(species);

    if(definition->nuclei[CHEM_ELEMENT_H] > 0 &&
       definition->nuclei[CHEM_ELEMENT_D] == 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0 &&
       definition->nuclei[CHEM_ELEMENT_O] == 0)
      state->abund[species] *= element_scale[CHEM_ELEMENT_H];

    if(definition->nuclei[CHEM_ELEMENT_HE] > 0)
      state->abund[species] *= element_scale[CHEM_ELEMENT_HE];

    if(definition->nuclei[CHEM_ELEMENT_D] > 0)
      state->abund[species] *= element_scale[CHEM_ELEMENT_D];

    if(definition->nuclei[CHEM_ELEMENT_C] > 0)
      state->abund[species] *= element_scale[CHEM_ELEMENT_C];
  }

  /* Carbon is normalized first. Oxygen carried by the resulting CO and CO2
   * is then held fixed while the oxygen-only reservoir is normalized. */
  double oxygen_in_carbon = 0.0;
  double oxygen_only_content = 0.0;

  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(species);
    double oxygen_content =
        definition->nuclei[CHEM_ELEMENT_O] * state->abund[species];

    if(definition->nuclei[CHEM_ELEMENT_C] > 0)
      oxygen_in_carbon += oxygen_content;
    else
      oxygen_only_content += oxygen_content;
  }

  if(oxygen_only_content > 0.0 &&
     element_total[CHEM_ELEMENT_O] >= oxygen_in_carbon) {
    element_scale[CHEM_ELEMENT_O] =
        (element_total[CHEM_ELEMENT_O] - oxygen_in_carbon) /
        oxygen_only_content;
  }

  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];
    const ChemSpeciesDefinition *definition =
        chemistry_species_definition(species);

    if(definition->nuclei[CHEM_ELEMENT_O] > 0 &&
       definition->nuclei[CHEM_ELEMENT_C] == 0)
      state->abund[species] *= element_scale[CHEM_ELEMENT_O];
  }

  /* Charge neutrality is reconstructed from the stored species updated by
   * this explicit scheme. Equilibrium-only intermediates are internal. */
  double net_non_electron_charge = 0.0;

  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];

    if(species != sELEC)
      net_non_electron_charge +=
          chemistry_species_definition(species)->charge *
          state->abund[species];
  }

  state->abund[sELEC] = net_non_electron_charge;
}

int chemistry_noneq_explicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt)
{
  const double nH = XH * environment->rho_cgs / PROTONMASS;
  const double temperature = work->thermodynamics.temperature;
  double abund_old[N_SP];
  double initial_abundance[N_SP];
  double abundance[N_SP];

  evaluate_reaction_rates(environment, step, state, work);

  for(int species = 0; species < N_SP; species++) {
    abund_old[species] = state->abund[species];
    initial_abundance[species] = state->abund[species] / XH;
    abundance[species] = initial_abundance[species];
  }

  if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
    return 2;

#ifdef CHEM_EXPLICIT_SIMULTANEOUS
  double production[N_SP];
  double destruction[N_SP];
  chemistry_network_production_destruction(work->rates.network, nH,
      abundance, work->rates.forward, work->rates.inverse, production,
      destruction);
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_intermediate_equilibrium(species))
      continue;
    abundance[species] = (initial_abundance[species] + production[species] * dt) /
                         (1.0 + destruction[species] * dt);
    if(abundance[species] < 0.0)
      abundance[species] = 0.0;
  }
#else
  const size_t explicit_species_count =
      sizeof(explicit_species_order) / sizeof(explicit_species_order[0]);

  for(size_t index = 0; index < explicit_species_count; index++) {
    enum SPECIES species = explicit_species_order[index];

    if(chemistry_species_is_intermediate_equilibrium(species)) {
      if(!chemistry_update_intermediate_equilibrium_species(
              CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
        return 2;
    } else {
      chemistry_explicit_update_species(
          species, nH, dt, initial_abundance, abundance, &work->rates);
    }
  }
#endif

  if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
    return 2;

  for(int species = 0; species < N_SP; species++)
    state->abund[species] = abundance[species] * XH;

  for(int iteration = 0; iteration < 16; iteration++) {
    double previous_intermediate_equilibrium[sizeof(intermediate_equilibrium_species) / sizeof(intermediate_equilibrium_species[0])];
    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
      previous_intermediate_equilibrium[i] = state->abund[intermediate_equilibrium_species[i]];

    chemistry_explicit_enforce_constraints(environment, state);
    for(int species = 0; species < N_SP; species++)
      abundance[species] = state->abund[species] / XH;
    if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
      return 2;
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

  if(isnan(state->abund[sHI])) {
    printf("[ChemID%d]HI is nan (dt=%g, rho_gas=%g, T_gas=%g)\n",
           environment->ID, dt, environment->rho_cgs, temperature);
    return 2;
  }
  if(isnan(temperature)) {
    printf("[ChemID%d]line=%d T_gas is nan (dt=%g, rho_gas=%g, T_gas=%g)\n",
           environment->ID, __LINE__, dt, environment->rho_cgs, temperature);
    return 2;
  }

  evaluate_chemistry_timescale(environment, state, work);
  ChemistryStatus cooling_status =
      chemical_cooling(environment, state, work, dt, abund_old);
  if(cooling_status != CHEMISTRY_STATUS_OK)
    return 0;
  return 1;
}

#ifndef MAX_ITER_CHEM
#define MAX_ITER_CHEM 3000
#endif
#define ITER_RESET 40

#ifdef SAHA
int chemistry_saha(const ChemEnvironment *environment, ChemState *state, ChemWorkspace *workspace, double dt)
{
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

  return 1;
}
#endif

typedef enum {
  NEWTON_RECOVERY_NONE,
  NEWTON_RECOVERY_INVALID_BASIN,
  NEWTON_RECOVERY_INVALID_STEP
} NewtonRecoveryReason;

typedef enum {
  NEWTON_BASIN_RESTART,
  NEWTON_BASIN_CONTINUE,
  NEWTON_BASIN_FAILED
} NewtonBasinOutcome;

typedef enum {
  NEWTON_STEP_CONTINUE,
  NEWTON_STEP_RECOVER_BASIN
} NewtonStepOutcome;

enum {
  NEWTON_REQUIRED_CONFIRMATIONS = 2,
  NEWTON_CLIPPING_HOLD_ITERATIONS = 2,
  NEWTON_BOUNDARY_STEPS_BEFORE_RESELECTION = 2,
  NEWTON_BASIN_ITERATION_LIMIT = 70
};

typedef struct {
  int iteration;
  int basin_start;

  int convergence_count;
  int clipping_hold;
  int boundary_streak;
  int chart_reselection_requested;

  int lu_valid;
  int force_numerical_jacobian;

  ChemNewtonResetIndex reset;
  int deterministic_reset_count;
  int random_restart_count;
  int clipped_step_count;
  int max_iterations;
  int allow_restarts;

  uint64_t random_state;
  double step_error;
  NewtonRecoveryReason recovery;
} NewtonControl;

typedef struct {
  double residual_error;
  ChemNewtonStepResult result;
  double correction[N_SP];
  int limiting_species;
} NewtonIterationData;

static void reset_for_newton_basin(NewtonControl *newton){
  newton->basin_start = newton->iteration;
  newton->convergence_count = 0;
  newton->clipping_hold = 0;
  newton->boundary_streak = 0;
  newton->chart_reselection_requested = 0;
  newton->lu_valid = 0;
  newton->force_numerical_jacobian = 0;
  newton->step_error = 100.0;
  newton->recovery = NEWTON_RECOVERY_NONE;
}

static void reset_after_chart_change(NewtonControl *newton){
  newton->convergence_count = 0;
  newton->boundary_streak = 0;
  newton->chart_reselection_requested = 0;
  newton->lu_valid = 0;
}

static void reset_for_newton_recovery(NewtonControl *newton){
  newton->convergence_count = 0;
  newton->step_error = 100.0;
}

static void chemistry_correct_full_implicit_intermediate_guess(
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
    valid = chemistry_update_conservation_dependent_species(abundance,
                                                            conservation);
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

static int chemistry_prepare_implicit_newton_problem(
    const ChemEnvironment *environment, const ChemState *state,
    const ChemRates *rates, double nH, ChemNewtonSystem *problem){
  double initialization_mean[N_SP];
  unsigned char initialize_at_floor[N_SP] = {0};

  problem->solver_species_count = 0;

  /* Normalize the input composition. */
  for(int species = 0; species < N_SP; species++)
    problem->abundance[species] = state->abund[species] / XH;

  if(!chemistry_update_intermediate_equilibrium_species(
          CHEM_REACTION_IMPLICIT, nH, problem->abundance, rates)){
    chemistry_profile_set_failure(1);
    return 0;
  }

  /* Start from bounds that are valid for every compiled species. */
  for(int species = 0; species < N_SP; species++){
    problem->initialization_floor[species] = 0.0;

    /* The initialization floor only provides a usable starting guess.
     * Newton itself may resolve any nonnegative abundance below it. */
    problem->newton_lower_bound[species] = 0.0;
    initialization_mean[species] = dmax(state->abund[species], 1.e-80);
    problem->upper_bound[species] = 2.0;
  }

  const ChemSpeciesInitialization *electron_initialization =
      &chemistry_species_initialization[sELEC];
  double ionization_scale =
      problem->abundance[sELEC] > electron_initialization->floor
          ? problem->abundance[sELEC] /
                electron_initialization->mean_per_baryon
          : 1.0;

  for(int species = 0; species < N_SP; species++){
    const ChemSpeciesInitialization *parameters =
        &chemistry_species_initialization[species];

    problem->initialization_floor[species] = parameters->floor;

    if(parameters->mean_per_baryon > 0.0 ||
       parameters->use_floor_as_guess)
      initialization_mean[species] =
          parameters->mean_per_baryon *
          (parameters->scale_with_ionization ? ionization_scale : 1.0);

    if(parameters->upper_bound_per_hydrogen > 0.0)
      problem->upper_bound[species] =
          parameters->upper_bound_per_hydrogen;

    initialize_at_floor[species] = parameters->use_floor_as_guess;
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
      if(chemistry_species_definition(species)->nuclei[element->element] > 0)
        problem->upper_bound[species] =
            dmin(problem->upper_bound[species], upper_bound);
    }
  }

  /* Replace unusable input values only for construction of the first guess. */
  for(int species = 0; species < N_SP; species++){
    initialization_mean[species] /= XH;

    if(problem->abundance[species] > problem->initialization_floor[species])
      continue;

    problem->abundance[species] = initialize_at_floor[species]
                                      ? problem->initialization_floor[species]
                                      : initialization_mean[species];
  }

  /* The caller's intermediate abundances are not part of the initial state
   * in the full-implicit initialization experiment. */
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

  if(!chemistry_build_conservation_system(
          &problem->conservation, environment->metallicity, problem->abundance)
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
      || !chemistry_initialize_full_implicit_intermediates(
             nH, equilibrium_initialization_electron, problem->abundance, rates,
             &problem->conservation)
#else
      || !chemistry_update_implicit_algebraic_species(
             nH, problem->abundance, rates, &problem->conservation)
#endif
  ){
    chemistry_profile_set_failure(2);
    return 0;
  }

  /* Freeze the chemistry-call initial state and reduced solver layout. */
  for(int species = 0; species < N_SP; species++){
    if(!chemistry_species_is_algebraic(species, &problem->conservation))
      problem->solver_species[problem->solver_species_count++] = species;

    problem->initial_abundance[species] = problem->abundance[species];
    problem->cooling_initial_abundance[species] =
        problem->initial_abundance[species] * XH;
  }

  return 1;
}

static int chemistry_initialize_implicit_newton_guess(
    double nH, double dt, const double *newton_guess, const ChemRates *rates,
    ChemNewtonSystem *problem){
  if(!newton_guess && !chemistry_newton_build_explicit_guess(
                           nH, dt, problem->initial_abundance,
                           problem->abundance, problem->newton_lower_bound,
                           problem->upper_bound, rates, &problem->conservation))
    return 0;

  if(newton_guess){
    for(int species = 0; species < N_SP; species++)
      if(!chemistry_species_is_algebraic(species, &problem->conservation))
        problem->abundance[species] = newton_guess[species] / XH;
    if(!chemistry_update_implicit_algebraic_species(
            nH, problem->abundance, rates, &problem->conservation))
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
      forward_valid = chemistry_update_implicit_algebraic_species(
          nH, perturbed_abundance, rates, conservation);
      if(forward_valid){
        chemistry_network_rhs(rates->network, CHEM_REACTION_IMPLICIT, nH,
                              perturbed_abundance, rates->forward,
                              rates->inverse, forward_rhs);
      }

      memcpy(perturbed_abundance, abundance, sizeof(perturbed_abundance));
      perturbed_abundance[perturbed_species] -= 0.5 * perturbation;
      backward_valid = perturbed_abundance[perturbed_species] >= 0.0 &&
                       chemistry_update_implicit_algebraic_species(
                           nH, perturbed_abundance, rates, conservation);
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

static int chemistry_build_newton_jacobian(
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

      if(!chemistry_reduce_jacobian_for_algebraic_species(
              nH, abundance, rates, solver_species, solver_species_count,
              conservation, rhs_jacobian)){
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

static int chemistry_solve_reduced_newton_system(
    int solver_species_count, int rebuild_newton_matrix, int *cached_lu_valid,
    int *lapack_pivots, double *newton_matrix, double *linear_rhs){
  int right_hand_side_count = 1;
  int linear_solver_status = 1;

  double linear_solve_timer = profile_start();

#ifdef CHEM_NEWTON_LU_REUSE
  if(rebuild_newton_matrix){
    dgetrf_(&solver_species_count, &solver_species_count, newton_matrix,
            &solver_species_count, lapack_pivots, &linear_solver_status);
    *cached_lu_valid = linear_solver_status == 0;
  }

  if(*cached_lu_valid){
    char transpose = 'N';
    dgetrs_(&transpose, &solver_species_count, &right_hand_side_count,
            newton_matrix, &solver_species_count, lapack_pivots, linear_rhs,
            &solver_species_count, &linear_solver_status);
  }
#else
  (void)rebuild_newton_matrix;
  (void)cached_lu_valid;
  dgesv_(&solver_species_count, &right_hand_side_count, newton_matrix,
         &solver_species_count, lapack_pivots, linear_rhs,
         &solver_species_count, &linear_solver_status);
#endif

  profile_end(CHEM_PROFILE_NEWTON_LINEAR_SOLVE, linear_solve_timer);
  return linear_solver_status;
}

static int chemistry_step_has_relevant_clip(
    const double *abundance, const double *newton_step,
    const int *clipped_species, double convergence_abundance_threshold,
    unsigned char *relevant_clipped_species){
  int relevant_step_was_clipped = 0;

  memset(relevant_clipped_species, 0,
         N_SP * sizeof(*relevant_clipped_species));

  for(int species = 0; species < N_SP; species++){
    const double abundance_before_step =
        abundance[species] - newton_step[species];
    int relevant_clip =
        clipped_species[species] &&
        dmax(fabs(abundance_before_step), fabs(abundance[species])) >=
            convergence_abundance_threshold;

#if defined(CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION) &&                  \
    !defined(CHEM_EQUILIBRIUM_SPECIES_CONVERGENCE)
    if(chemistry_is_full_implicit_initialized_intermediate(species))
      relevant_clip = 0;
#endif

    if(relevant_clip){
      relevant_step_was_clipped = 1;
      relevant_clipped_species[species] = 1;
    }
  }

  return relevant_step_was_clipped;
}

static void chemistry_update_implicit_timescale(
    ChemWorkspace *work, double dt, const double *abundance,
    const double *initial_abundance,
    const ChemConservationSystem *conservation){
  double minimum_timescale = 1.0e20;

  for(int species = 0; species < N_SP; species++){
    if(chemistry_species_is_algebraic(species, conservation))
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

static int chemistry_export_implicit_solution(
    const ChemEnvironment *environment, const ChemStepData *step,
    ChemState *state, ChemWorkspace *work, const ChemState *input_state,
    double dt, double nH, double temperature,
    const double *initialization_floor, double *cooling_initial_abundance,
    double *abundance, const ChemConservationSystem *conservation){
  for(int species = 0; species < N_SP; species++){
    if(chemistry_species_is_algebraic(species, conservation))
      continue;
    if(abundance[species] < initialization_floor[species])
      abundance[species] = initialization_floor[species];
  }

  /* Floors apply only to independent species.  Reconstruct algebraic species
   * afterward so conservation and equilibrium remain internally consistent. */
  if(!chemistry_update_implicit_algebraic_species(nH, abundance, &work->rates,
                                                   conservation)){
    *state = *input_state;
    return 0;
  }

  for(int species = 0; species < N_SP; species++){
    if(conservation->is_dependent[species] &&
        abundance[species] < 0.0 &&
        abundance[species] >= -CHEMISTRY_NEWTON_TRACE_ABUNDANCE){
      abundance[species] = 0.0;
    }
  }

  for(int species = 0; species < N_SP; species++)
    state->abund[species] = abundance[species] * XH;

  ChemistryStatus cooling_status =
      chemical_cooling(environment, state, work, dt, cooling_initial_abundance);
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

static int prepare_newton_basin(
    double nH, double metallicity, const ChemRates *rates,
    ChemNewtonSystem *system, NewtonControl *newton){
  reset_for_newton_basin(newton);

  if(newton->reset != CHEM_NEWTON_RESET_NONE){
    chemistry_newton_apply_deterministic_reset(
        newton->reset, metallicity, system->initial_abundance,
        system->abundance);

    chemistry_newton_reselect_conservation_chart(metallicity, system);

    chemistry_correct_full_implicit_intermediate_guess(
        nH, system->abundance, rates, &system->conservation);

    for(int species = 0; species < N_SP; species++){
      system->total_step[species] =
          system->abundance[species] - system->initial_abundance[species];
    }
  }

  int basin_is_valid = chemistry_newton_normalize_current_state(
      system->abundance, system->upper_bound);

  if(basin_is_valid){
    basin_is_valid = chemistry_update_implicit_algebraic_species(
        nH, system->abundance, rates, &system->conservation);
  }

  if(!basin_is_valid){
    newton->recovery = NEWTON_RECOVERY_INVALID_BASIN;
    return 0;
  }

  for(int species = 0; species < N_SP; species++){
    if(chemistry_species_is_algebraic(species, &system->conservation))
      system->total_step[species] = 0.0;
  }

  return 1;
}

static NewtonStepOutcome resolve_invalid_newton_step(
    double metallicity, ChemNewtonSystem *system, NewtonControl *newton){
  if(chemistry_newton_reselect_conservation_chart(metallicity, system)){
    reset_after_chart_change(newton);
    newton->iteration++;

    return NEWTON_STEP_CONTINUE;
  }

  newton->recovery = NEWTON_RECOVERY_INVALID_STEP;
  return NEWTON_STEP_RECOVER_BASIN;
}

static NewtonStepOutcome finish_newton_step(
    double nH, ChemNewtonSystem *system, NewtonControl *newton,
    NewtonIterationData *iteration_data){
  const int step_was_clipped = iteration_data->result.step_was_clipped;

  if(step_was_clipped)
    newton->lu_valid = 0;

  newton->step_error = chemistry_newton_material_step_norm(
      system->abundance, iteration_data->correction, &system->conservation,
      &iteration_data->limiting_species);

  for(int species = 0; species < N_SP; species++){
    system->total_step[species] =
        system->abundance[species] - system->initial_abundance[species];

    if(chemistry_species_is_algebraic(species, &system->conservation))
      system->total_step[species] = 0.0;
  }

  unsigned char relevant_clipped_species[N_SP];
  int relevant_step_was_clipped = chemistry_step_has_relevant_clip(
      system->abundance, iteration_data->correction,
      iteration_data->result.clipped_species,
      CHEMISTRY_NEWTON_TRACE_ABUNDANCE, relevant_clipped_species);

  if(relevant_step_was_clipped){
    newton->clipping_hold = NEWTON_CLIPPING_HOLD_ITERATIONS;
    newton->clipped_step_count++;
  }

  newton->boundary_streak =
      step_was_clipped ? newton->boundary_streak + 1 : 0;

  if(newton->boundary_streak >= NEWTON_BOUNDARY_STEPS_BEFORE_RESELECTION)
    newton->chart_reselection_requested = 1;

  if(iteration_data->residual_error > 1.0 ||
      relevant_step_was_clipped){
    newton->convergence_count = 0;
  }else if(newton->clipping_hold > 0){
    newton->clipping_hold--;
    newton->convergence_count = 0;
  }else if(newton->step_error <= 1.0){
    newton->convergence_count++;
  }else{
    newton->convergence_count = 0;
  }

  /* Profiling only; these calls do not change the Newton state. */
  if(newton->step_error > 0.0)
    chemistry_profile_record_newton_limiter(
        iteration_data->limiting_species, nH);

  for(int species = 0; species < N_SP; species++){
    if(relevant_clipped_species[species])
      chemistry_profile_record_clipped_species(species);
  }

  if(newton->allow_restarts && newton->convergence_count == 0 &&
     newton->iteration - newton->basin_start + 1 >=
         NEWTON_BASIN_ITERATION_LIMIT &&
     newton->reset <= newton->deterministic_reset_count){
    newton->recovery = NEWTON_RECOVERY_NONE;
    return NEWTON_STEP_RECOVER_BASIN;
  }

  newton->iteration++;
  return NEWTON_STEP_CONTINUE;
}

static NewtonBasinOutcome recover_newton_basin(
    double nH, double metallicity, const ChemRates *rates,
    ChemNewtonSystem *system, NewtonControl *newton){
  reset_for_newton_recovery(newton);

  if(newton->allow_restarts &&
     newton->reset < newton->deterministic_reset_count){
    newton->reset++;
    return NEWTON_BASIN_RESTART;
  }

  if(newton->recovery == NEWTON_RECOVERY_INVALID_STEP)
    chemistry_profile_set_failure(6);

  if(newton->allow_restarts && newton->iteration < newton->max_iterations){
    int found_random_guess = chemistry_find_random_newton_guess(
        nH, metallicity, system->initial_abundance, system->abundance,
        rates, &system->conservation, &newton->random_state,
        &newton->iteration, newton->max_iterations,
        &newton->random_restart_count);

    if(found_random_guess){
      newton->reset++;
      newton->force_numerical_jacobian = 0;
      return NEWTON_BASIN_RESTART;
    }
  }

  if(newton->recovery != NEWTON_RECOVERY_NONE){
    if(newton->recovery == NEWTON_RECOVERY_INVALID_BASIN)
      chemistry_profile_set_failure(3);

    return NEWTON_BASIN_FAILED;
  }

  newton->iteration++;
  return NEWTON_BASIN_CONTINUE;
}

static int chemistry_noneq_chem_implicit_impl(
    const ChemEnvironment *environment, const ChemStepData *step,
    ChemState *state, ChemWorkspace *work, double dt,
    const double *newton_guess, int max_iterations, int allow_restarts,
    int report_failure){
  /* Phase 1: construct the normalized implicit problem and its initial guess. */
  chemistry_profile_set_failure(0);

  const ChemState state_old = *state;

  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas / PROTONMASS;
  const double nH = XH * n_nuclei;
  const double T_gas = work->thermodynamics.temperature;

  evaluate_reaction_rates(environment, step, state, work);

  ChemNewtonSystem system;
  if(!chemistry_prepare_implicit_newton_problem(environment, state,
                                                 &work->rates, nH, &system) ||
      !chemistry_initialize_implicit_newton_guess(nH, dt, newton_guess,
                                                  &work->rates, &system)){
    *state = state_old;
    return 0;
  }

  /* Abundance state and physical bounds. */
  double *abundance = system.abundance;
  const double *initial_abundance = system.initial_abundance;
  double *cooling_initial_abundance = system.cooling_initial_abundance;

  const double *initialization_floor = system.initialization_floor;
  const double *lower_bound = system.newton_lower_bound;
  const double *upper_bound = system.upper_bound;

  /* Reduced-system layout. */
  double *total_step = system.total_step;
  enum SPECIES *solver_species = system.solver_species;
  ChemConservationSystem *conservation = &system.conservation;

  /* Dense work arrays for one Newton iteration. */
  int lapack_pivots[N_SP];

  double rhs_jacobian[N_SP][N_SP];
  double reaction_rhs[N_SP];
  double species_production[N_SP];
  double species_loss[N_SP];

  double newton_matrix[N_SP * N_SP];
  double linear_rhs[N_SP];

  NewtonIterationData iteration_data = {0};

  /* Phase 2: iterate Newton basins. Recovery changes the guess at the same
   * physical timestep. */
  NewtonControl newton = {
      .iteration = 1,
      .basin_start = 1,
      .reset = CHEM_NEWTON_RESET_NONE,
      .deterministic_reset_count = CHEM_NEWTON_RESET_LAST,
      .max_iterations = max_iterations,
      .allow_restarts = allow_restarts,
      .random_state = chemistry_newton_random_seed(nH, T_gas, dt),
      .step_error = 100.0,
      .recovery = NEWTON_RECOVERY_NONE,
  };

  ChemNewtonDiagnostics diagnostics;
  NewtonBasinOutcome basin_outcome;

restart_newton_basin:
  basin_outcome = NEWTON_BASIN_CONTINUE;

  if(!prepare_newton_basin(
          nH, environment->metallicity, &work->rates, &system, &newton)){
    basin_outcome = recover_newton_basin(
        nH, environment->metallicity, &work->rates, &system, &newton);
  }

  if(basin_outcome == NEWTON_BASIN_CONTINUE){
    /* Diagnostics only: print and snapshot the accepted basin entry. */
    print_newton_basin_entry(work, T_gas, abundance);

#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
    memcpy(diagnostics.basin_entry, abundance, sizeof(diagnostics.basin_entry));
#endif
  }

  /* Main Newton loop.
   * Each pass evaluates the residual, solves the reduced linear system,
   * applies a bounded correction, and updates convergence state. */
  while(basin_outcome == NEWTON_BASIN_CONTINUE &&
        newton.iteration <= max_iterations &&
        newton.convergence_count < NEWTON_REQUIRED_CONFIRMATIONS){
    /* Re-select a conservation reservoir after repeated boundary contact. */
    if(newton.chart_reselection_requested){
      if(chemistry_newton_reselect_conservation_chart(
              environment->metallicity, &system))
        reset_after_chart_change(&newton);
    }

    /* Reaction RHS and residual. */
    double rhs_timer = profile_start();

    chemistry_network_rhs_gross(
        work->rates.network, CHEM_REACTION_IMPLICIT, nH, abundance,
        work->rates.forward, work->rates.inverse, reaction_rhs,
        species_production, species_loss);

    profile_end(CHEM_PROFILE_NEWTON_RHS, rhs_timer);

    iteration_data.residual_error = chemistry_newton_evaluate_residual(
        dt, abundance, initial_abundance, upper_bound, reaction_rhs,
        species_production, species_loss, &work->rates, conservation);

    /* Reduced Jacobian and linear solve. */
    int rebuild_newton_matrix = 1;

#ifdef CHEM_NEWTON_LU_REUSE
    /* Reuse is limited to the mandatory second convergence confirmation after
     * an already tiny correction. This guard preserves the thermochemical
     * branch selected by a freshly rebuilt Jacobian. */
    rebuild_newton_matrix =
        !newton.lu_valid || newton.convergence_count != 1 ||
        newton.step_error >= 1.0e-12;
#endif

    if(!chemistry_build_newton_jacobian(
            nH, abundance, reaction_rhs, &work->rates, conservation,
            solver_species, system.solver_species_count,
            rebuild_newton_matrix, &newton.force_numerical_jacobian,
            rhs_jacobian)){
      chemistry_profile_set_failure(5);
      *state = state_old;

      return 0;
    }

    /* Diagnostics only. */
    print_newton_iteration(work, newton.iteration, T_gas, dt,
                           newton.step_error,
                           iteration_data.limiting_species, abundance,
                           initial_abundance, reaction_rhs);
    print_newton_jacobian(work, newton.iteration, rhs_jacobian);

    /* Backward Euler Newton system:
     *   (I - dt dR/dy) delta_y = dt R(y) - (y - y_initial). */
    if(rebuild_newton_matrix){
      for(int column = 0; column < system.solver_species_count; column++){
        const int column_species = solver_species[column];

        for(int row = 0; row < system.solver_species_count; row++){
          const int row_species = solver_species[row];
          const int matrix_index =
              column * system.solver_species_count + row;

          newton_matrix[matrix_index] =
              column == row
                  ? 1.0 - dt * rhs_jacobian[column_species][row_species]
                  : -dt * rhs_jacobian[column_species][row_species];
        }
      }
    }

    for(int row = 0; row < system.solver_species_count; row++){
      const int species = solver_species[row];
      linear_rhs[row] = fma(reaction_rhs[species], dt, -total_step[species]);
    }

    int linear_solver_status = chemistry_solve_reduced_newton_system(
        system.solver_species_count, rebuild_newton_matrix, &newton.lu_valid,
        lapack_pivots, newton_matrix, linear_rhs);

    memset(iteration_data.correction, 0,
           sizeof(iteration_data.correction));

    for(int row = 0; row < system.solver_species_count; row++)
      iteration_data.correction[solver_species[row]] = linear_rhs[row];

    /* Apply the bounded correction. */
    newton.step_error = 0.0;

    double newton_step_timer = profile_start();

    int step_is_valid =
        linear_solver_status == 0 &&
        chemistry_newton_apply_step(
            nH, abundance, iteration_data.correction, lower_bound, upper_bound,
            &work->rates, conservation, &iteration_data.result);

    profile_end(CHEM_PROFILE_NEWTON_STEP, newton_step_timer);

    NewtonStepOutcome step_outcome;

    if(step_is_valid){
      step_outcome = finish_newton_step(
          nH, &system, &newton, &iteration_data);
    }else{
      step_outcome = resolve_invalid_newton_step(
          environment->metallicity, &system, &newton);
    }

    if(step_outcome == NEWTON_STEP_RECOVER_BASIN){
      basin_outcome = recover_newton_basin(
          nH, environment->metallicity, &work->rates, &system, &newton);
    }
  }

  if(basin_outcome == NEWTON_BASIN_RESTART)
    goto restart_newton_basin;

  if(basin_outcome == NEWTON_BASIN_FAILED){
    *state = state_old;
    return 0;
  }

  /* Phase 3: certify the final iterate and report a failed solve. */
  if(!chemistry_update_implicit_algebraic_species(
          nH, abundance, &work->rates, conservation)){
    chemistry_profile_set_failure(7);
    *state = state_old;

    return 0;
  }

  int converged =
      newton.convergence_count >= NEWTON_REQUIRED_CONFIRMATIONS;

  if(converged){
    chemistry_profile_record_newton_solve(
        newton.iteration - 1, newton.reset, newton.clipped_step_count);
  }

  print_random_restart_success(
      &diagnostics, converged, newton.iteration, newton.reset,
      newton.random_restart_count);

  if(!converged){
    work->evaluation.err_chem = 1.0;

    if(allow_restarts && report_failure){
      chemistry_newton_diagnostics_report_failure(
          environment, work, &state_old, dt, n_nuclei, nH, T_gas,
          newton.step_error, iteration_data.limiting_species,
          newton.iteration, newton.reset, newton.random_restart_count,
          allow_restarts, abundance, initial_abundance,
          iteration_data.correction, CHEMISTRY_NEWTON_TRACE_ABUNDANCE,
          conservation);
    }

    if(allow_restarts)
      chemistry_profile_set_failure(8);

    *state = state_old;

    return 0;
  }

  chemistry_update_implicit_timescale(work, dt, abundance, initial_abundance,
                                      conservation);

  /* Phase 4: export the converged state and evaluate chemical cooling. */
  work->evaluation.err_chem = 0.0;

  return chemistry_export_implicit_solution(
      environment, step, state, work, &state_old, dt, nH, T_gas,
      initialization_floor, cooling_initial_abundance, abundance,
      conservation);
}

int chemistry_noneq_chem_implicit(const ChemEnvironment *environment,
                                  const ChemStepData *step, ChemState *state,
                                  ChemWorkspace *work, double dt){
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            NULL, MAX_ITER_CHEM, 1, 1);
}

static int
chemistry_noneq_chem_implicit_trial(const ChemEnvironment *environment,
                                    const ChemStepData *step, ChemState *state,
                                    ChemWorkspace *work, double dt){
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            NULL, MAX_ITER_CHEM, 1, 0);
}

int chemistry_noneq_chem_implicit_with_guess(const ChemEnvironment *environment,
                                             const ChemStepData *step,
                                             ChemState *state,
                                             ChemWorkspace *work, double dt,
                                             const double *newton_guess,
                                             int max_iterations){
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            newton_guess, max_iterations, 0, 0);
}

#undef MAX_ITER_CHEM

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
          chemistry_species_definition(species);
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
          chemistry_species_definition(species);
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
    int charge = chemistry_species_definition(species)->charge;

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
          chemistry_species_definition(species);
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
