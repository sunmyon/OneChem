/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <signal.h>
#include <stdint.h>

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
#include "species.h"
#include "chemistry_solver.h"
#include "saha.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

#define CHEMISTRY_NEWTON_TRACE_ABUNDANCE 1.e-12
#define CHEMISTRY_HYDROGEN_IONIZATION_FLOOR 1.e-20
#define CHEMISTRY_CARBON_OXYGEN_FLOOR 1.e-16

static void correct_initial_abundances(double *yi);
static void give_initial_guess(double temperature, double nH, double dt, double *yi, double *initial_guess, double *rates);
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
/* This is not the normal energy convergence tolerance.  It only identifies
 * a discontinuous chemistry-branch jump that leaves a residual many orders
 * of magnitude larger than the thermal energy. */
#define ENERGY_GROSS_RESIDUAL_FACTOR 1.e6

int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt){
#ifdef CHEMISTRY_PROFILE
  double profile_start = chemistry_profile_time();
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
  double g = 0., g0, array[1000];

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

  if(FlagPrint)
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

    *work = work_old;
    *state = state_old;

    if(isnan(T_new)){
      printf("[task%03d](iter=%d) T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	     ThisTask, iter, state->utherm_cgs,  state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I]);
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

  *work = work_old;
  *state = state_old;

  T_old = T_new;
  T_new = T_old + dT;

  int flag_converged = 0;

  int j;
  for(j=iter;j<MAX_ITER;j++){
#ifdef FLAG_BRENT
    if( flag_brent == 1){
      *work = work_old;
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
#ifdef ID_CHEM_DEBUG
      if(FlagPrint)
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

      *work = work_old;
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


    array[j] = g, g0 = g;
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

    *work = work_old;
    *state = state_old;
  }


#ifdef FLAG_BRENT
  if(flag_converged == 0){
    T_new = T_gas_init;

    for(;j<MAX_ITER;j++){
      *work = work_old;
      *state = state_old;

      d_T = dmax(0.2 * T_new, 2000);
      if(sign > 0)
	T_new = dmax(T_new-d_T,T_new/10.);
      else
	T_new = dmin(T_new + d_T,T_new * 10.);

      err_chem = give_g(environment, step, state, work, chemistry_solver, dt, T_new, &e_new, &g);

      if(err_chem < 1.e-4 && sign * g < 0.){
	if(sign > 0)
	  g_minus = g, T_minus = T_new;
	else
	  g_plus = g, T_plus = T_new;

	flag_brent = 1;
	break;
      }

      if(isnan(T_new)){
	printf("[task%03d](iter=%d) T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	       ThisTask, j, state->utherm_cgs,  state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I]);
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
      *work = work_old;
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
      *work = continuation_workspace;
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
      *work = normal_workspace;
      T_new = normal_temperature;
      e_new = normal_energy;
      g = normal_residual;
      err = normal_error;
    }
  }
#endif
  if(isnan(err) || isnan(T_new))
    raise(11);

  if(err > 1.0e-4 || FlagPrint == 1){
    printf("%s err=%g id=%d T_old=%g T_new=%g Task=%d\n"
	   "All.ComovingIntegrationOn=%d, All.Time=%g, dt=%g, n_nuclei=%g, T_gas=%g, gamma=%g;\n"
	   "yi[0]=%g, yi[1]=%g, yi[2]=%g;\n"
	   "yi[%d]=%g, yi[%d]=%g;\n",
	   err > 1.0e-4 ? "not converged(in energy)." : "converged(in energy diagnostic).",
	   err, environment->ID, T_gas_init, T_new, ThisTask,
	   All.ComovingIntegrationOn, All.Time, dt, rho_gas/PROTONMASS, T_gas_init, state_old.gamma,
	   state_old.abund[sELEC], state_old.abund[sHI], state_old.abund[sHII],
	   sH2I, state_old.abund[sH2I], sHM, state_old.abund[sHM]);
#ifdef H2_PLUS
    printf("yi[%d]=%g;\n", sH2II, state_old.abund[sH2II]);
#endif
#ifdef HELIUM
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n"
	   , sHeI, state_old.abund[sHeI], sHeII, state_old.abund[sHeII], sHeIII, state_old.abund[sHeIII]);
#endif
#ifdef HD
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n",
	   sDI,   state_old.abund[sDI],
	   sDII,  state_old.abund[sDII],
	   sHDI,  state_old.abund[sHDI]);
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
    printf("environment->rate_Hminus_pd=%g, environment->rate_H2_diss=%g;\n"
	   ,environment->rate_Hminus_pd, environment->rate_H2_diss);
#ifdef DUST_COOLING
    printf("environment->G0=%g;\n", environment->G0);
#endif
#endif
#ifdef BACKGROUND_LW_RADIATION
    printf("All.rate_Hminus_pd=%g, All.rate_H2_diss=%g, All.J21=%g;\n"
	   ,All.rate_Hminus_pd, All.rate_H2_diss, All.J21);
#ifdef METALLICITY
    printf("All.G0=%g;\n", All.G0);
#endif
#endif
#endif

#ifdef IONIZATION_PHOTON
    printf("step->k22=%g, step->piHI=%g;\n"
           , step->k22, step->piHI);
#endif

#ifdef METALLICITY
    printf("environment->Metallicity=%g;\n", environment->metallicity);
#endif

    printf("yf[0]=%g, yf[1]=%g, yf[2]=%g, yf[3]=%g, yf[4]=%g;\n",
	   state->abund[sELEC], state->abund[sHI], state->abund[sHII], state->abund[sH2I], state->abund[sHM]);
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
    *work = work_old;
    *state = state_old;
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_ENERGY_SOLVER, chemistry_profile_time() - profile_start);
#endif
    return 1;
  }

#ifdef CHEMISTRY_PROFILE
  chemistry_profile_add(CHEM_PROFILE_ENERGY_SOLVER, chemistry_profile_time() - profile_start);
#endif
  return 0;
}
#pragma opt


static double evaluate_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g, int commit, int report_chemistry_failure, const double *newton_guess, int max_chemistry_iterations)
{
#ifdef CHEMISTRY_PROFILE
  double profile_start = chemistry_profile_time();
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
      double chemistry_start = chemistry_profile_time();
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
      chemistry_profile_add(CHEM_PROFILE_IMPLICIT_CHEMISTRY, chemistry_profile_time() - chemistry_start);
#endif
    }
  evaluate_radiative_cooling(environment, step, &state_trial, &work_trial);

  double T_gamma = chemistry_temperature(state_trial.abund, state_trial.gamma, state_trial.utherm_cgs);
  state_trial.gamma = compute_gamma_from_abund(state_trial.abund, T_gamma);

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
  chemistry_profile_add(CHEM_PROFILE_ENERGY_RESIDUAL, chemistry_profile_time() - profile_start);
#endif

  if(FlagPrint)
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
  double val_a, val_b;
  double val_new;
  double e_new;

  if(g_min * g_max > 0){
    printf("something wrong in solve_bisection."
	   " g_min and g_max must have opposite sign.\n");
    return 1.e30;
  }

  T_a = T_min, val_a = g_min;
  T_b = T_max, val_b = g_max;
  T_cand = 0.5 * (T_a + T_b);
  e_new = *e_new0;

  iter = 0;
  while(iter < max_iter && err > 1.0e-6){
    T_cand = 0.5 * (T_a + T_b);

    err = give_g(environment, step, state, work, chemistry_solver, dt, T_cand, &e_new, &val_new);

    if(val_a * val_new < 0) T_b = T_cand, val_b = val_new;
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

static int chemistry_species_is_intermediate_equilibrium(enum SPECIES species)
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
  double profile_start = chemistry_profile_time();
#define INTERMEDIATE_RETURN(value) do { chemistry_profile_add(CHEM_PROFILE_INTERMEDIATES, chemistry_profile_time() - profile_start); return (value); } while(0)
#else
#define INTERMEDIATE_RETURN(value) return (value)
#endif
  double updated[N_SP];

  if(species_count == 0)
    INTERMEDIATE_RETURN(1);

  if(species_count == 1) {
    enum SPECIES species = species_list[0];
    if(!chemistry_network_intermediate_equilibrium_abundance(species, mode, nH, abundance, rates->forward, rates->inverse, &updated[0], NULL))
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
  double assembly_start = chemistry_profile_time();
#endif
  chemistry_network_intermediate_system(mode, nH, trial, rates->forward, rates->inverse,
                                        species_list, n_intermediate, rhs, matrix);

  for(int row = 0; row < n_intermediate; row++)
    updated[row] = -rhs[row];
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_add(CHEM_PROFILE_INTERMEDIATE_ASSEMBLY, chemistry_profile_time() - assembly_start);
  double intermediate_linear_start = chemistry_profile_time();
#endif

  dgesv_(&n_intermediate, &one, matrix, &n_intermediate, pivots, updated, &n_intermediate, &info);
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_add(CHEM_PROFILE_INTERMEDIATE_LINEAR_SOLVE, chemistry_profile_time() - intermediate_linear_start);
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

#ifdef CARBON_OXYGEN_CHEMISTRY
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 6
#else
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 4
#endif

typedef struct {
  int count;
  double coefficient[CHEM_MAX_CONSERVATION_CONSTRAINTS][N_SP];
  double total[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  enum SPECIES dependent_species[CHEM_MAX_CONSERVATION_CONSTRAINTS];
  double dependent_inverse[CHEM_MAX_CONSERVATION_CONSTRAINTS * CHEM_MAX_CONSERVATION_CONSTRAINTS];
  unsigned char is_dependent[N_SP];
} ChemConservationSystem;

static void chemistry_add_conservation_constraint(ChemConservationSystem *conservation, enum SPECIES dependent_species, double total)
{
  int constraint = conservation->count++;
  conservation->dependent_species[constraint] = dependent_species;
  conservation->total[constraint] = total;
  conservation->is_dependent[dependent_species] = 1;
}

static int chemistry_invert_dependent_matrix(ChemConservationSystem *conservation)
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

static enum SPECIES chemistry_largest_reservoir(const double *abundance,
                                                 const enum SPECIES *species,
                                                 const double *multiplicity,
                                                 int count)
{
  enum SPECIES largest = species[0];
  double largest_content = multiplicity[0] * abundance[largest];
  for(int i = 1; i < count; i++) {
    double content = multiplicity[i] * abundance[species[i]];
    if(content > largest_content) {
      largest = species[i];
      largest_content = content;
    }
  }
  return largest;
}

static int chemistry_build_conservation_system(ChemConservationSystem *conservation,
                                                double metallicity,
                                                const double *abundance)
{
  memset(conservation, 0, sizeof(*conservation));

  const enum SPECIES hydrogen_reservoirs[] = {
    sHI, sHII, sH2I,
#ifdef H2_PLUS
    sH2II,
#endif
  };
  const double hydrogen_multiplicity[] = {
    1.0, 1.0, 2.0,
#ifdef H2_PLUS
    2.0,
#endif
  };
  enum SPECIES hydrogen_dependent = chemistry_largest_reservoir(
      abundance, hydrogen_reservoirs, hydrogen_multiplicity,
      (int)(sizeof(hydrogen_reservoirs) / sizeof(hydrogen_reservoirs[0])));
  int hydrogen_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, hydrogen_dependent, 1.0);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[hydrogen_constraint][species] =
        chemistry_species_definition(species)->nuclei[CHEM_ELEMENT_H];

#ifdef HELIUM
  int helium_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, sHeI, XHe / XH);
  conservation->coefficient[helium_constraint][sHeI] = 1.0;
  conservation->coefficient[helium_constraint][sHeII] = 1.0;
  conservation->coefficient[helium_constraint][sHeIII] = 1.0;
#endif

#ifdef HD
  const enum SPECIES deuterium_reservoirs[] = {
    sDI, sDII, sHDI,
#ifdef HD_PLUS
    sHDII,
#endif
#ifdef DMINUS
    sDM,
#endif
  };
  const double deuterium_multiplicity[] = {
    1.0, 1.0, 1.0,
#ifdef HD_PLUS
    1.0,
#endif
#ifdef DMINUS
    1.0,
#endif
  };
  enum SPECIES deuterium_dependent = chemistry_largest_reservoir(
      abundance, deuterium_reservoirs, deuterium_multiplicity,
      (int)(sizeof(deuterium_reservoirs) / sizeof(deuterium_reservoirs[0])));
  int deuterium_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, deuterium_dependent,
                                        D_number_fraction / XH);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[deuterium_constraint][species] =
        chemistry_species_definition(species)->nuclei[CHEM_ELEMENT_D];
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
#ifndef CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION
  const enum SPECIES carbon_reservoirs[] = {sCI, sCII, sCO, sCO2};
  const double carbon_multiplicity[] = {1.0, 1.0, 1.0, 1.0};
  const enum SPECIES oxygen_reservoirs[] = {sOI, sOH, sH2O, sO2};
  const double oxygen_multiplicity[] = {1.0, 1.0, 1.0, 2.0};
  enum SPECIES carbon_dependent = chemistry_largest_reservoir(
      abundance, carbon_reservoirs, carbon_multiplicity, 4);
  enum SPECIES oxygen_dependent = chemistry_largest_reservoir(
      abundance, oxygen_reservoirs, oxygen_multiplicity, 4);
  int carbon_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, carbon_dependent, A_C_sun * metallicity);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[carbon_constraint][species] =
        chemistry_species_definition(species)->nuclei[CHEM_ELEMENT_C];

  int oxygen_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, oxygen_dependent, A_O_sun * metallicity);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[oxygen_constraint][species] =
        chemistry_species_definition(species)->nuclei[CHEM_ELEMENT_O];
#endif
#else
  (void)metallicity;
#endif

  int charge_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, sELEC, 0.0);
  for(int species = 0; species < N_SP; species++)
    conservation->coefficient[charge_constraint][species] =
        chemistry_species_definition(species)->charge;

  return chemistry_invert_dependent_matrix(conservation);
}

static int chemistry_species_is_conservation_dependent(enum SPECIES species, const ChemConservationSystem *conservation)
{
  return conservation->is_dependent[species];
}

static int chemistry_species_is_algebraic(enum SPECIES species, const ChemConservationSystem *conservation)
{
  return chemistry_species_is_intermediate_equilibrium(species) || chemistry_species_is_conservation_dependent(species, conservation);
}

static int chemistry_update_conservation_dependent_species(double *abundance, const ChemConservationSystem *conservation)
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

static int chemistry_is_full_implicit_initialized_intermediate(
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

  /* If the formal equilibrium intermediate population would consume more of
   * a conserved reservoir than is available, reduce all initialized
   * intermediates together.  This remains history-free and preserves the
   * equilibrium proportions used to seed the full implicit system. */
  double fraction = 1.0;
  for(int dependent = 0; dependent < conservation->count; dependent++) {
    enum SPECIES species = conservation->dependent_species[dependent];
    if(equilibrium[species] < 0.0) {
      double denominator = base[species] - equilibrium[species];
      if(denominator <= 0.0)
        return 0;
      fraction = dmin(fraction, base[species] / denominator);
    }
  }
  if(fraction < 1.0)
    fraction = dmax(0.0, fraction * (1.0 - 64.0 * DBL_EPSILON));
  memcpy(abundance, base, sizeof(base));
  for(int i = 0; i < count; i++) {
    enum SPECIES species = full_implicit_equilibrium_initialization_species[i];
    abundance[species] = fraction * equilibrium[species];
  }
  if(!chemistry_update_conservation_dependent_species(abundance, conservation))
    return 0;
  for(int dependent = 0; dependent < conservation->count; dependent++)
    if(abundance[conservation->dependent_species[dependent]] < 0.0)
      return 0;
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_record_equilibrium_initialization(
      fraction, input_electron, abundance[sELEC]);
#endif
  return 1;
}

#endif

static int chemistry_update_implicit_algebraic_species(double nH, double *abundance, const ChemRates *rates, const ChemConservationSystem *conservation)
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

static int chemistry_apply_implicit_newton_step(double nH, double *abundance, double *delta, const double *lower_bound, const double *upper_bound, const ChemRates *rates, const ChemConservationSystem *conservation, int allow_tiny_negative_dependent, int *invalid_species, int *step_was_clipped, int *clipped_species)
{
  double initial[N_SP];
  memcpy(initial, abundance, sizeof(initial));

  *invalid_species = -1;
  *step_was_clipped = 0;
  memset(clipped_species, 0, N_SP * sizeof(*clipped_species));
  double trial[N_SP];
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_algebraic(species, conservation)) {
      trial[species] = initial[species];
      continue;
    }
    if(!isfinite(delta[species])) {
      *invalid_species = species;
      return 0;
    }
    double limited_lower_bound = dmax(0.1 * initial[species], lower_bound[species]);
    if(initial[species] + delta[species] < limited_lower_bound) {
      delta[species] = limited_lower_bound - initial[species];
      *step_was_clipped = 1;
      clipped_species[species] = 1;
    } else if(initial[species] + delta[species] > upper_bound[species]) {
      delta[species] = dmin(0.9 * initial[species], 0.9 * (upper_bound[species] - initial[species]));
      *step_was_clipped = 1;
      clipped_species[species] = 1;
    }
    trial[species] = initial[species] + delta[species];
  }

  int valid = 0;
  double step_factor = 1.0;
  for(int attempt = 0; attempt < 48; attempt++) {
    for(int species = 0; species < N_SP; species++)
      trial[species] = chemistry_species_is_algebraic(species, conservation) ?
          initial[species] : initial[species] + step_factor * delta[species];
    valid = chemistry_update_implicit_algebraic_species(nH, trial, rates, conservation);
    if(valid) {
      for(int species = 0; species < N_SP; species++)
        if(!isfinite(trial[species]) ||
           (trial[species] < 0.0 &&
            (!allow_tiny_negative_dependent ||
             !chemistry_species_is_conservation_dependent(species, conservation) ||
             trial[species] < -CHEMISTRY_NEWTON_TRACE_ABUNDANCE)) ||
           (!chemistry_species_is_intermediate_equilibrium(species) &&
            trial[species] > upper_bound[species])) {
          valid = 0;
          *invalid_species = species;
          break;
        }
    }
    if(valid)
      break;
    step_factor *= 0.5;
  }
#ifdef CHEMISTRY_PROFILE
  if(!valid)
    chemistry_profile_set_failure(600);
#endif

  if(!valid)
    return 0;
  if(step_factor < 1.0)
    *step_was_clipped = 1;
  for(int species = 0; species < N_SP; species++) {
    delta[species] = trial[species] - initial[species];
    abundance[species] = trial[species];
  }
  return 1;
}

static uint64_t chemistry_random_u64(uint64_t *state)
{
  uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static double chemistry_random_unit(uint64_t *state)
{
  return (chemistry_random_u64(state) >> 11) * 0x1.0p-53;
}

static double chemistry_random_weight(uint64_t *state)
{
  return 1.e-12 + chemistry_random_unit(state);
}

static void chemistry_random_partition(uint64_t *random_state, int count, const double *multiplicity, double total, double *values)
{
  double weights[4];
  double normalization = 0.0;
  for(int i = 0; i < count; i++) {
    weights[i] = chemistry_random_weight(random_state);
    normalization += multiplicity[i] * weights[i];
  }
  for(int i = 0; i < count; i++)
    values[i] = total * weights[i] / normalization;
}

static double chemistry_random_fraction(uint64_t *random_state)
{
  double first = chemistry_random_weight(random_state);
  double second = chemistry_random_weight(random_state);
  return first / (first + second);
}

static double chemistry_random_log_value(uint64_t *random_state, double minimum, double maximum)
{
  return exp(log(minimum) + chemistry_random_unit(random_state) * log(maximum / minimum));
}

static int chemistry_randomize_newton_guess(double nH, double metallicity, double *abundance, const ChemRates *rates, const ChemConservationSystem *conservation, uint64_t *random_state)
{
  const double unit_multiplicity[4] = {1.0, 1.0, 1.0, 1.0};
  const double hydrogen_multiplicity[4] = {1.0, 1.0, 2.0, 2.0};
  const double oxygen_multiplicity[4] = {1.0, 1.0, 1.0, 2.0};
  double values[4];

#ifdef H2_PLUS
  chemistry_random_partition(random_state, 4, hydrogen_multiplicity, 1.0, values);
  abundance[sHI] = values[0];
  abundance[sHII] = values[1];
  abundance[sH2I] = values[2];
  abundance[sH2II] = values[3];
#else
  chemistry_random_partition(random_state, 3, hydrogen_multiplicity, 1.0, values);
  abundance[sHI] = values[0];
  abundance[sHII] = values[1];
  abundance[sH2I] = values[2];
#endif
  abundance[sHM] = chemistry_random_log_value(random_state, 1.e-55, 1.e-5);

#ifdef HELIUM
  chemistry_random_partition(random_state, 3, unit_multiplicity, XHe / XH, values);
  abundance[sHeI] = values[0];
  abundance[sHeII] = values[1];
  abundance[sHeIII] = values[2];
#endif

#ifdef HD
  chemistry_random_partition(random_state, 3, unit_multiplicity, D_number_fraction / XH, values);
  abundance[sDI] = values[0];
  abundance[sDII] = values[1];
  abundance[sHDI] = values[2];
#ifdef HD_PLUS
  abundance[sHDII] = chemistry_random_log_value(random_state, 1.e-30, 0.01 * D_number_fraction / XH);
#endif
#ifdef DMINUS
  abundance[sDM] = chemistry_random_log_value(random_state, 1.e-30, 0.01 * D_number_fraction / XH);
#endif
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  double carbon_left = A_C_sun * metallicity;
  double oxygen_left = A_O_sun * metallicity;
  abundance[sCO2] = dmin(carbon_left, 0.5 * oxygen_left) * chemistry_random_fraction(random_state);
  carbon_left -= abundance[sCO2];
  oxygen_left -= 2.0 * abundance[sCO2];
  abundance[sCO] = dmin(carbon_left, oxygen_left) * chemistry_random_fraction(random_state);
  carbon_left -= abundance[sCO];
  oxygen_left -= abundance[sCO];
  chemistry_random_partition(random_state, 2, unit_multiplicity, carbon_left, values);
  abundance[sCI] = values[0];
  abundance[sCII] = values[1];
  chemistry_random_partition(random_state, 4, oxygen_multiplicity, oxygen_left, values);
  abundance[sOI] = values[0];
  abundance[sOH] = values[1];
  abundance[sH2O] = values[2];
  abundance[sO2] = values[3];
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
  abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
  abundance[sH3OII] = abundance[sO2II] = 0.0;
#endif
#else
  (void)metallicity;
#endif

  return chemistry_update_implicit_algebraic_species(nH, abundance, rates, conservation);
}

static int chemistry_find_random_newton_guess(double nH, double metallicity,
                                               const double *initial_abundance,
                                               double *abundance,
                                               const ChemRates *rates,
                                               const ChemConservationSystem *conservation,
                                               uint64_t *random_state,
                                               int *iteration, int max_iterations,
                                               int *random_attempt_count)
{
  while(*iteration < max_iterations) {
    memcpy(abundance, initial_abundance, N_SP * sizeof(*abundance));
    (*iteration)++;
    (*random_attempt_count)++;
    if(chemistry_randomize_newton_guess(nH, metallicity, abundance, rates,
                                        conservation, random_state))
      return 1;
  }
  return 0;
}

static int chemistry_reduce_jacobian_for_algebraic_species(double nH, const double *abundance, const ChemRates *rates, const enum SPECIES *solver_species, int n_solver_species, const ChemConservationSystem *conservation, double drfdy[N_SP][N_SP])
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
  double response_start = chemistry_profile_time();
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
  chemistry_profile_add(CHEM_PROFILE_JACOBIAN_RESPONSE, chemistry_profile_time() - response_start);
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
    double intermediate_solve_start = chemistry_profile_time();
#endif
    int info;
    dgesv_(&n_intermediate, &n_solver_species, intermediate_jacobian, &n_intermediate, intermediate_pivots, intermediate_response, &n_intermediate, &info);
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_JACOBIAN_INTERMEDIATE_SOLVE, chemistry_profile_time() - intermediate_solve_start);
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
  double schur_start = chemistry_profile_time();
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
  chemistry_profile_add(CHEM_PROFILE_JACOBIAN_SCHUR, chemistry_profile_time() - schur_start);
#endif

  return 1;
}

static void chemistry_explicit_update_species(enum SPECIES species, double nH, double dt, const double *initial_abundance, double *abundance, const ChemRates *rates)
{
  double production;
  double destruction;
  chemistry_network_species_production_destruction(species, nH, abundance, rates->forward, rates->inverse, &production, &destruction);
  abundance[species] = (initial_abundance[species] + production * dt) /
                       (1.0 + destruction * dt);
  if(abundance[species] < 0.0)
    abundance[species] = 0.0;
}

static const enum SPECIES explicit_species_order[] = {
  sHI, sHII,
#ifdef HELIUM
  sHeI, sHeII, sHeIII,
#endif
  sELEC, sHM,
#ifdef H2_PLUS
  sH2II,
#endif
  sH2I,
#ifdef HD
  sDI, sDII, sHDI,
#ifdef HD_PLUS
  sHDII,
#endif
#ifdef DMINUS
  sDM,
#endif
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  sCI, sCII, sCO, sCO2, sOI, sOH, sH2O, sO2,
#endif
};

static void chemistry_explicit_enforce_constraints(const ChemEnvironment *environment, ChemState *state)
{
  double correction = H_number_fraction /
      (state->abund[sHI] + state->abund[sHII] + state->abund[sHM] + 2.0 * state->abund[sH2I]
#ifdef H2_PLUS
       + 2.0 * state->abund[sH2II]
#endif
#ifdef HD
       + state->abund[sHDI]
#ifdef HD_PLUS
       + state->abund[sHDII]
#endif
#endif
      );

  state->abund[sHI] *= correction;
  state->abund[sHII] *= correction;
  state->abund[sHM] *= correction;
  state->abund[sH2I] *= correction;
#ifdef H2_PLUS
  state->abund[sH2II] *= correction;
#endif

#ifdef HD
  correction = D_number_fraction /
      (state->abund[sDI] + state->abund[sDII] + state->abund[sHDI]
#ifdef HD_PLUS
       + state->abund[sHDII]
#endif
#ifdef DMINUS
       + state->abund[sDM]
#endif
      );
  state->abund[sDI] *= correction;
  state->abund[sDII] *= correction;
  state->abund[sHDI] *= correction;
#ifdef HD_PLUS
  state->abund[sHDII] *= correction;
#endif
#ifdef DMINUS
  state->abund[sDM] *= correction;
#endif
#endif

#ifdef HELIUM
  correction = He_number_fraction /
      (state->abund[sHeI] + state->abund[sHeII] + state->abund[sHeIII]);
  state->abund[sHeI] *= correction;
  state->abund[sHeII] *= correction;
  state->abund[sHeIII] *= correction;
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  double carbon_total = XH * A_C_sun * environment->metallicity;
  double carbon_sum = state->abund[sCI] + state->abund[sCII] + state->abund[sCO] + state->abund[sCO2];
  if(carbon_sum > 0.0) {
    correction = carbon_total / carbon_sum;
    state->abund[sCI] *= correction;
    state->abund[sCII] *= correction;
    state->abund[sCO] *= correction;
    state->abund[sCO2] *= correction;
  }

  double oxygen_total = XH * A_O_sun * environment->metallicity;
  double oxygen_in_carbon = state->abund[sCO] + 2.0 * state->abund[sCO2];
  double oxygen_sum = state->abund[sOI] + state->abund[sOH] + state->abund[sH2O] + 2.0 * state->abund[sO2];
  if(oxygen_sum > 0.0 && oxygen_total >= oxygen_in_carbon) {
    correction = (oxygen_total - oxygen_in_carbon) / oxygen_sum;
    state->abund[sOI] *= correction;
    state->abund[sOH] *= correction;
    state->abund[sH2O] *= correction;
    state->abund[sO2] *= correction;
  }
#else
  (void)environment;
#endif

  state->abund[sELEC] = state->abund[sHII] - state->abund[sHM];
#ifdef H2_PLUS
  state->abund[sELEC] += state->abund[sH2II];
#endif
#ifdef HD
  state->abund[sELEC] += state->abund[sDII];
#ifdef HD_PLUS
  state->abund[sELEC] += state->abund[sHDII];
#endif
#ifdef DMINUS
  state->abund[sELEC] -= state->abund[sDM];
#endif
#endif
#ifdef HELIUM
  state->abund[sELEC] += state->abund[sHeII] + 2.0 * state->abund[sHeIII];
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  state->abund[sELEC] += state->abund[sCII];
#endif
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
  chemistry_network_production_destruction(nH, abundance, work->rates.forward, work->rates.inverse, production, destruction);
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_intermediate_equilibrium(species))
      continue;
    abundance[species] = (initial_abundance[species] + production[species] * dt) /
                         (1.0 + destruction[species] * dt);
    if(abundance[species] < 0.0)
      abundance[species] = 0.0;
  }
#else
  for(size_t order_index = 0; order_index < sizeof(explicit_species_order) / sizeof(explicit_species_order[0]); order_index++) {
    enum SPECIES species = explicit_species_order[order_index];
    if(chemistry_species_is_intermediate_equilibrium(species)) {
      if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_EXPLICIT, nH, abundance, &work->rates))
        return 2;
    } else {
      chemistry_explicit_update_species(species, nH, dt, initial_abundance, abundance, &work->rates);
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
    printf("[task%03d]HI is nan (ID=%d, dt=%g, rho_gas=%g, T_gas=%g)\n",
           ThisTask, environment->ID, dt, environment->rho_cgs, temperature);
    return 2;
  }
  if(isnan(temperature)) {
    printf("[task%03d]line=%d T_gas is nan (dt=%g, rho_gas=%g, T_gas=%g)\n",
           ThisTask, __LINE__, dt, environment->rho_cgs, temperature);
    return 2;
  }

  evaluate_chemistry_timescale(environment, state, work);
  chemical_cooling(environment, state, work, dt, abund_old);
  return 1;
}

#ifndef MAX_ITER_CHEM
#define MAX_ITER_CHEM 1000
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

  workspace->evaluation.err_chem = chem_saha_5sp(nH, T_gas, y_input);

  state->abund[sHI] = y_input[sHI] * XH;
  state->abund[sHII] = y_input[sHII] * XH;
  state->abund[sH2I] = y_input[sH2I] * XH;
  state->abund[sHM] = y_input[sHM] * XH;
  state->abund[sELEC] = y_input[sELEC] * XH;

  double edot_H2form = 0.;
  double edot_ionize = 0.;
  double edot_HMform = 0.;
  if(dt >= 1.e-30) {
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

static int chemistry_noneq_chem_implicit_impl(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt, const double *newton_guess, int max_iterations, int allow_restarts, int report_failure)
{
#ifdef CHEMISTRY_PROFILE
  chemistry_profile_set_failure(0);
#endif
  enum SPECIES species;
  const ChemState state_old = *state;

  int n_sp = N_SP;
  int inc  = 1, info;

  const double rho_gas = environment->rho_cgs;
  const double n_nuclei  = rho_gas/PROTONMASS;
  const double nH = XH * n_nuclei;
  const double T_gas = work->thermodynamics.temperature;
  uint64_t random_state = UINT64_C(0x6a09e667f3bcc909);
  uint64_t value_bits;
  memcpy(&value_bits, &nH, sizeof(value_bits));
  random_state ^= value_bits;
  memcpy(&value_bits, &T_gas, sizeof(value_bits));
  random_state ^= value_bits + UINT64_C(0x9e3779b97f4a7c15);
  memcpy(&value_bits, &dt, sizeof(value_bits));
  random_state ^= value_bits + UINT64_C(0x3c6ef372fe94f82b);

  evaluate_reaction_rates(environment, step, state, work);

  double yi[n_sp];
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    yi[i_sp] = state->abund[i_sp] / XH;

  if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, yi, &work->rates)) {
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_set_failure(1);
#endif
    *state = state_old;
    return 0;
  }

  int err_flag = 0, ipiv[n_sp];
  double drfdy[n_sp][n_sp],dy[n_sp];
  double y_init[n_sp], abund_init[n_sp];
  double y_mean[n_sp],y_max[n_sp],y_min[n_sp];
  double yi_ion;
  double rf_tot[n_sp];
  double AA[n_sp*n_sp],ddy[n_sp], solver_rhs[n_sp];
  enum SPECIES solver_species[n_sp];
  int n_solver_species = 0;
  ChemConservationSystem conservation;

  for(int i_sp = 0; i_sp < n_sp; i_sp++) {
    y_min[i_sp] = 0.0;
    y_mean[i_sp] = dmax(state->abund[i_sp], 1.e-80);
    y_max[i_sp] = 2.0;
  }

  y_min[sELEC] = 1.e-80,  y_mean[sELEC]  = 1.e-15   ,  y_max[sELEC]  = 2.0;

  if(yi[sELEC] > y_min[sELEC])
    yi_ion = yi[sELEC]/y_mean[sELEC];
  else
    yi_ion = 1.;

  y_min[sHI]   = CHEMISTRY_HYDROGEN_IONIZATION_FLOOR, y_mean[sHI]  = 0.76,   y_max[sHI]  = 1.001;
  y_min[sHII]  = CHEMISTRY_HYDROGEN_IONIZATION_FLOOR, y_mean[sHII] = 1.e-15, y_max[sHII] = 1.001;
#ifdef HELIUM
  y_min[sHeI]  = 1.e-100,  y_mean[sHeI]   = 0.06          ,  y_max[sHeI]   = 0.09;
  y_min[sHeII] = 1.e-100,  y_mean[sHeII]  = 4.e-12 *yi_ion,  y_max[sHeII]  = 0.09;
  y_min[sHeIII]= 1.e-100,  y_mean[sHeIII] = 4.e-13 *yi_ion,  y_max[sHeIII] = 0.09;
#endif
  y_min[sHM]   = 1.e-55,  y_mean[sHM]    = 9.e-12 *yi_ion,  y_max[sHM]    = 1.0;
  y_min[sH2I]  = 1.e-50,  y_mean[sH2I]   = 6.e-06        ,  y_max[sH2I]   = 0.50005;
#ifdef H2_PLUS
  y_min[sH2II] = 1.e-55,  y_mean[sH2II]  = 2.e-14 *yi_ion,  y_max[sH2II]  = 0.5;
#endif
#ifdef HD
  y_min[sDI]   = 1.e-25,  y_mean[sDI]    = 4.e-05        ,  y_max[sDI]    = D_number_fraction/XH * 1.001;
  y_min[sDII]  = 1.e-80,  y_mean[sDII]   = 6.e-09 *yi_ion,  y_max[sDII]   = D_number_fraction/XH * 1.001;
  y_min[sHDI]  = 1.e-30,  y_mean[sHDI]   = 1.e-08        ,  y_max[sHDI]   = D_number_fraction/XH * 1.001;
#ifdef HD_PLUS
  y_min[sHDII] = 1.e-30,  y_mean[sHDII]  = 2.e-18 *yi_ion,  y_max[sHDII]  = 0.5*D_number_fraction/XH;
#endif
#ifdef DMINUS
  y_min[sDM]   = 1.e-30,  y_mean[sDM]    = 5.e-16 *yi_ion,  y_max[sDM]    = D_number_fraction/XH;
#endif
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  y_min[sCI] = y_min[sCII] = y_min[sOI] = CHEMISTRY_CARBON_OXYGEN_FLOOR;
  double carbon_max = dmax(A_C_sun * environment->metallicity * 1.001, 1.e-80);
  double oxygen_max = dmax(A_O_sun * environment->metallicity * 1.001, 1.e-80);
  y_max[sCI] = y_max[sCII] = y_max[sCO] = y_max[sCO2] = carbon_max;
  y_max[sOI] = y_max[sOH] = y_max[sH2O] = y_max[sO2] = oxygen_max;
  y_max[sCOII] = y_max[sCH] = y_max[sCH2] = carbon_max;
  y_max[sOII] = y_max[sOHII] = y_max[sH2OII] = y_max[sH3OII] = y_max[sO2II] = oxygen_max;
#endif

  for(int i_sp=0;i_sp<n_sp;i_sp++) {
    y_mean[i_sp] /= XH;
    y_init[i_sp] = yi[i_sp];

    if(yi[i_sp] <= y_min[i_sp]/* || yi[i_sp] > y_max[i_sp]*/){
      int initialize_at_floor = i_sp == sHI || i_sp == sHII;
#ifdef HD
      initialize_at_floor = initialize_at_floor || i_sp == sDI || i_sp == sDII;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
      initialize_at_floor = initialize_at_floor || i_sp == sCI || i_sp == sCII || i_sp == sOI;
#endif
      yi[i_sp] = initialize_at_floor ? y_min[i_sp] : y_mean[i_sp];
    }
  }

  /* In the initialization experiment the caller's intermediate abundances
   * are not part of this chemistry step's initial state. */
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
  const double equilibrium_initialization_electron = yi[sELEC];
  for(size_t i = 0;
      i < sizeof(full_implicit_equilibrium_initialization_species) /
          sizeof(full_implicit_equilibrium_initialization_species[0]);
      i++)
    yi[full_implicit_equilibrium_initialization_species[i]] = 0.0;
#endif

  correct_initial_abundances(yi);
  if(!chemistry_build_conservation_system(&conservation, environment->metallicity, yi)
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
     || !chemistry_initialize_full_implicit_intermediates(
            nH, equilibrium_initialization_electron, yi, &work->rates,
            &conservation)
#else
     || !chemistry_update_implicit_algebraic_species(
            nH, yi, &work->rates, &conservation)
#endif
     ) {
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_set_failure(2);
#endif
    *state = state_old;
    return 0;
  }
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    if(!chemistry_species_is_algebraic(i_sp, &conservation))
      solver_species[n_solver_species++] = i_sp;
  for(int i_sp=0;i_sp<n_sp;i_sp++){
    y_init[i_sp] = yi[i_sp];
    abund_init[i_sp] = y_init[i_sp] * XH;
  }

  if(newton_guess) {
    for(int i_sp = 0; i_sp < n_sp; i_sp++)
      if(!chemistry_species_is_algebraic(i_sp, &conservation))
        yi[i_sp] = newton_guess[i_sp] / XH;
    if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
      *state = state_old;
      return 0;
    }
  }

#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
  {
    double guess_before_equilibrium_correction[N_SP];
    memcpy(guess_before_equilibrium_correction, yi,
           sizeof(guess_before_equilibrium_correction));
    int correction_valid = chemistry_solve_intermediate_equilibrium_list(
        CHEM_REACTION_IMPLICIT, nH, yi, &work->rates,
        full_implicit_equilibrium_initialization_species,
        sizeof(full_implicit_equilibrium_initialization_species) /
            sizeof(full_implicit_equilibrium_initialization_species[0]));
    if(correction_valid)
      correction_valid = chemistry_update_conservation_dependent_species(
          yi, &conservation);
    for(int dependent = 0;
        correction_valid && dependent < conservation.count; dependent++)
      if(yi[conservation.dependent_species[dependent]] < 0.0)
        correction_valid = 0;
    if(!correction_valid)
      memcpy(yi, guess_before_equilibrium_correction, sizeof(yi));
  }
#endif

  //give_initial_guess(T_gas, nH, dt, yi, y_init, k);
  for(int i_sp=0;i_sp<n_sp;i_sp++){
    dy[i_sp] = yi[i_sp] - y_init[i_sp];
  }

  double err_max= 100.0;
  int flag_conv = 0, flag_reset = 0, itr = 1;
  int basin_iteration_start = itr;
  int species_clip_streak[N_SP] = {0};
  int convergence_hold_after_clipping = 0;
  int random_restart_count = 0;
#ifdef CHEM_NEWTON_LU_REUSE
  int cached_newton_lu_valid = 0;
#endif
#ifdef CHEMISTRY_PROFILE
  int profile_clipped_steps = 0;
#endif
#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  double basin_entry[N_SP];
#endif
  const int deterministic_reset_count =
#ifdef CARBON_OXYGEN_CHEMISTRY
      17;
#else
      3;
#endif
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
  int use_numerical_jacobian = 0;
#endif

 chem_START:
#ifdef CHEM_NEWTON_LU_REUSE
  cached_newton_lu_valid = 0;
#endif
  if(flag_reset){
    basin_iteration_start = itr;
    memset(species_clip_streak, 0, sizeof(species_clip_streak));
    convergence_hold_after_clipping = 0;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
    use_numerical_jacobian = 0;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
    if(flag_reset >= 1 && flag_reset <= deterministic_reset_count)
      memcpy(yi, y_init, sizeof(yi));
    if(flag_reset == 1) {
      yi[sHII] = 0.0;
      yi[sH2I] = 0.5 * (1.0 - yi[sHI]);
#ifdef HD
      yi[sDII] = 0.0;
      yi[sHDI] = D_number_fraction / XH - yi[sDI];
#endif
    }
    if(flag_reset == 4) {
      /* Root basin found for the saved Z=1, T=3127 K non-equilibrium
       * intermediate failure; C/O abundances scale with metallicity. */
      yi[sELEC] = 0.25179906689485865;
      yi[sHI] = 0.21111795614880233;
      yi[sHII] = 0.25179091941005283;
      yi[sH2I] = 0.268540173142787;
      yi[sHM] = 2.8871723966669333e-7;
#ifdef HD
      yi[sDI] = 1.6078502270673519e-5;
      yi[sDII] = 7.9679348332395664e-6;
      yi[sHDI] = 5.5588260539816532e-6;
#endif
      yi[sCI] = 2.39504403722442e-5 * environment->metallicity;
      yi[sCII] = 4.682672122055309e-7 * environment->metallicity;
      yi[sCO] = 1.424602206249225e-5 * environment->metallicity;
      yi[sCO2] = 2.77335270353058e-4 * environment->metallicity;
      yi[sOI] = 7.6656784138147242e-8 * environment->metallicity;
      yi[sOH] = 1.6190277220737118e-6 * environment->metallicity;
      yi[sH2O] = 1.6557922775014028e-6 * environment->metallicity;
      yi[sO2] = 1.3659802238392435e-6 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = yi[sCH] = yi[sCH2] = 0.0;
      yi[sOII] = yi[sOHII] = yi[sH2OII] = 0.0;
      yi[sH3OII] = yi[sO2II] = 0.0;
#endif
    }
    if(flag_reset == 5) {
      /* High-temperature basin found for the saved z=5, Z=0.01,
       * T=32465 K non-equilibrium intermediate failure. */
      yi[sELEC] = 0.99990537140681268;
      yi[sHI] = 1.2607213620080047e-4;
      yi[sHII] = 0.99987311647556509;
      yi[sH2I] = 1.75285745333539e-15;
      yi[sHM] = 2.0672061115436203e-12;
#ifdef HD
      yi[sDI] = 3.7373480562468117e-9;
      yi[sDII] = 2.9601525809838385e-5;
      yi[sHDI] = 1.10545931196041e-19;
#endif
      yi[sCI] = 9.8954371995761087e-5 * environment->metallicity;
      yi[sCII] = 3.0644731129628816e-14 * environment->metallicity;
      yi[sCO] = 1.2445469056970178e-4 * environment->metallicity;
      yi[sCO2] = 9.2590917520743629e-5 * environment->metallicity;
      yi[sOI] = 2.2723925332810629e-8 * environment->metallicity;
      yi[sOH] = 5.0314506352684828e-18 * environment->metallicity;
      yi[sH2O] = 6.3241931343626324e-19 * environment->metallicity;
      yi[sO2] = 6.1939284267970712e-18 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 1.255623341210446e-29 * environment->metallicity;
      yi[sCH] = 1.9883148389096883e-11 * environment->metallicity;
      yi[sCH2] = 3.7372262679145738e-19 * environment->metallicity;
      yi[sOII] = 1.8420215400720396e-4 * environment->metallicity;
      yi[sOHII] = 8.1138596456252812e-5 * environment->metallicity;
      yi[sH2OII] = 1.3106030161376541e-19 * environment->metallicity;
      yi[sH3OII] = 2.9572631077328755e-25 * environment->metallicity;
      yi[sO2II] = 1.6821257889333632e-18 * environment->metallicity;
#endif
    }
    if(flag_reset == 6) {
      /* High-temperature root found for the saved z=20, Z=0.1,
       * T=34463 K non-equilibrium intermediate failure. */
      yi[sELEC] = 1.00003264289242;
      yi[sHI] = 1.03013406908217e-5;
      yi[sHII] = 0.99998204201121632;
      yi[sH2I] = 1.2172160777026613e-17;
      yi[sHM] = 1.1701898781268036e-14;
#ifdef HD
      yi[sDI] = 3.0535703651096367e-10;
      yi[sDII] = 2.9604957800857676e-5;
      yi[sHDI] = 5.5231265654500275e-19;
#endif
      yi[sCI] = 1.2597966882885965e-4 * environment->metallicity;
      yi[sCII] = 8.3018830260684831e-14 * environment->metallicity;
      yi[sCO] = 1.5001076262865832e-5 * environment->metallicity;
      yi[sCO2] = 1.750191666936573e-4 * environment->metallicity;
      yi[sOI] = 1.3560747480040021e-9 * environment->metallicity;
      yi[sOH] = 2.1353875395991124e-13 * environment->metallicity;
      yi[sH2O] = 6.3809201639573046e-18 * environment->metallicity;
      yi[sO2] = 1.0436809651017052e-17 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 1.6945907606749465e-14 * environment->metallicity;
      yi[sCH] = 8.0572677676953679e-11 * environment->metallicity;
      yi[sCH2] = 7.5419748150528964e-12 * environment->metallicity;
      yi[sOII] = 1.3339284910343401e-4 * environment->metallicity;
      yi[sOHII] = 7.6566384940889594e-5 * environment->metallicity;
      yi[sH2OII] = 2.3006907515112232e-16 * environment->metallicity;
      yi[sH3OII] = 1.405326460540793e-19 * environment->metallicity;
      yi[sO2II] = 2.9627906658264367e-18 * environment->metallicity;
#endif
    }
    if(flag_reset == 7) {
      /* Root found for the saved z=20, Z=0.1, T=30811 K failure. */
      yi[sELEC] = 1.0000288624085911;
      yi[sHI] = 1.3988119506476293e-5;
      yi[sHII] = 0.99997788016567746;
      yi[sH2I] = 1.1813936041141837e-14;
      yi[sHM] = 2.9075900000133324e-14;
#ifdef HD
      yi[sDI] = 4.1470368754092607e-10;
      yi[sDII] = 2.9604848454206449e-5;
      yi[sHDI] = 7.4587283296514796e-19;
#endif
      yi[sCI] = 1.2979521907697177e-4 * environment->metallicity;
      yi[sCII] = 7.7535774052324431e-14 * environment->metallicity;
      yi[sCO] = 1.1185265893047512e-5 * environment->metallicity;
      yi[sCO2] = 1.7501948788490025e-4 * environment->metallicity;
      yi[sOI] = 1.8135284293301643e-9 * environment->metallicity;
      yi[sOH] = 1.3117562321818905e-17 * environment->metallicity;
      yi[sH2O] = 2.5589089726834677e-18 * environment->metallicity;
      yi[sO2] = 1.0430683407407197e-17 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 5.0009305674336627e-19 * environment->metallicity;
      yi[sCH] = 2.7067407023932438e-11 * environment->metallicity;
      yi[sCH2] = 1.37168559988535e-16 * environment->metallicity;
      yi[sOII] = 1.3245682424280528e-4 * environment->metallicity;
      yi[sOHII] = 8.1317120565874268e-5 * environment->metallicity;
      yi[sH2OII] = 5.2601891224465922e-19 * environment->metallicity;
      yi[sH3OII] = 4.4300118965131421e-24 * environment->metallicity;
      yi[sO2II] = 2.7377339204103489e-18 * environment->metallicity;
#endif
    }
    if(flag_reset == 8) {
      /* Root found for the saved z=0, Z=0.1, T=1261 K failure. */
      yi[sELEC] = 0.46299502150158972;
      yi[sHI] = 0.23374610242570271;
      yi[sHII] = 0.46296887826087041;
      yi[sH2I] = 0.15163314890128671;
      yi[sHM] = 5.9271967626337258e-7;
#ifdef HD
      yi[sDI] = 1.683184927734549e-5;
      yi[sDII] = 1.2773409639027497e-5;
      yi[sHDI] = 4.2415217524219631e-12;
#endif
      yi[sCI] = 9.8936286841851588e-5 * environment->metallicity;
      yi[sCII] = 4.966451336116868e-5 * environment->metallicity;
      yi[sCO] = 8.674984114183423e-5 * environment->metallicity;
      yi[sCO2] = 8.064859901241815e-5 * environment->metallicity;
      yi[sOI] = 3.5921356531547067e-5 * environment->metallicity;
      yi[sOH] = 4.5881503536558188e-5 * environment->metallicity;
      yi[sH2O] = 3.1222407874484122e-5 * environment->metallicity;
      yi[sO2] = 6.1983349342868853e-5 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 0.0;
      yi[sCH] = 6.6189123140422269e-10 * environment->metallicity;
      yi[sCH2] = 9.7751495971518061e-11 * environment->metallicity;
      yi[sOII] = 1.7050635464271139e-5 * environment->metallicity;
      yi[sOHII] = 7.2860046485103603e-5 * environment->metallicity;
      yi[sH2OII] = 5.029057419330885e-8 * environment->metallicity;
      yi[sH3OII] = 2.168143434270531e-11 * environment->metallicity;
      yi[sO2II] = 0.0;
#endif
    }
    if(flag_reset == 9) {
      /* Root found for the saved z=20, Z=1, T=92675 K failure. */
      yi[sELEC] = 1.0000969021535941;
      yi[sHI] = 1.9564190488324027e-6;
      yi[sHII] = 0.99991421835454997;
      yi[sH2I] = 4.4505255218988363e-16;
      yi[sHM] = 3.4963947708440967e-17;
#ifdef HD
      yi[sDI] = 5.7952040282222372e-11;
      yi[sDII] = 2.960520520585443e-5;
      yi[sHDI] = 2.7857180798267994e-20;
#endif
      yi[sCI] = 5.5686779410032298e-5 * environment->metallicity;
      yi[sCII] = 1.3672326857913847e-13 * environment->metallicity;
      yi[sCO] = 9.870514390705527e-5 * environment->metallicity;
      yi[sCO2] = 1.6160806521633251e-4 * environment->metallicity;
      yi[sOI] = 1.3195881492927567e-10 * environment->metallicity;
      yi[sOH] = 2.4340338772299575e-17 * environment->metallicity;
      yi[sH2O] = 1.0298898956331353e-18 * environment->metallicity;
      yi[sO2] = 8.6677703979452596e-18 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 1.7290999112878472e-24 * environment->metallicity;
      yi[sCH] = 1.1329856501115037e-11 * environment->metallicity;
      yi[sCH2] = 1.2358743697202434e-19 * environment->metallicity;
      yi[sOII] = 6.9253378631097865e-5 * environment->metallicity;
      yi[sOHII] = 8.3825215070314075e-5 * environment->metallicity;
      yi[sH2OII] = 3.4874029612235333e-19 * environment->metallicity;
      yi[sH3OII] = 2.8939600484293469e-26 * environment->metallicity;
      yi[sO2II] = 4.9171363067507986e-18 * environment->metallicity;
#endif
    }
    if(flag_reset == 10) {
      /* Root found for the saved z=0, Z=1, T=1501 K failure. */
      yi[sELEC] = 9.5276870514691089e-9;
      yi[sHI] = 0.011505425276457534;
      yi[sHII] = 6.8250328029227238e-10;
      yi[sH2I] = 0.4942257502637345;
      yi[sHM] = 5.2955689998561437e-15;
#ifdef HD
      yi[sDI] = 2.8704494388345964e-7;
      yi[sDII] = 1.4227412396213185e-14;
      yi[sHDI] = 2.9318218199783867e-5;
#endif
      yi[sCI] = 3.8783656922327507e-8 * environment->metallicity;
      yi[sCII] = 1.3763227740127845e-9 * environment->metallicity;
      yi[sCO] = 6.3451249345320581e-5 * environment->metallicity;
      yi[sCO2] = 2.5239240133078359e-4 * environment->metallicity;
      yi[sOI] = 4.2364060284967586e-11 * environment->metallicity;
      yi[sOH] = 3.2276692703190392e-9 * environment->metallicity;
      yi[sH2O] = 6.7532087215730819e-6 * environment->metallicity;
      yi[sO2] = 1.9285299022076825e-13 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 3.964022484995573e-16 * environment->metallicity;
      yi[sCH] = 9.134910977838556e-9 * environment->metallicity;
      yi[sCH2] = 1.0705443282522379e-7 * environment->metallicity;
      yi[sOII] = 3.581165004835212e-18 * environment->metallicity;
      yi[sOHII] = 2.7651575745559313e-17 * environment->metallicity;
      yi[sH2OII] = 7.7561751953795306e-14 * environment->metallicity;
      yi[sH3OII] = 7.4687736386544986e-9 * environment->metallicity;
      yi[sO2II] = 4.3727921331630416e-16 * environment->metallicity;
#endif
    }
    if(flag_reset == 11) {
      /* Root found for the saved z=20, Z=1, T=1515 K failure. */
      yi[sELEC] = 0.30236727032410249;
      yi[sHI] = 0.5795951115982052;
      yi[sHII] = 0.30233947754074625;
      yi[sH2I] = 0.059027888873290288;
      yi[sHM] = 9.7375063629401053e-8;
#ifdef HD
      yi[sDI] = 1.9322771279217555e-5;
      yi[sDII] = 1.0282489050219449e-5;
      yi[sHDI] = 2.8284577354283037e-12;
#endif
      yi[sCI] = 1.7626131693844425e-5 * environment->metallicity;
      yi[sCII] = 1.1724828034265817e-5 * environment->metallicity;
      yi[sCO] = 1.583620530423108e-5 * environment->metallicity;
      yi[sCO2] = 2.7081155884091612e-4 * environment->metallicity;
      yi[sOI] = 4.751595719675113e-6 * environment->metallicity;
      yi[sOH] = 2.2991851135850131e-6 * environment->metallicity;
      yi[sH2O] = 1.3739858274135889e-6 * environment->metallicity;
      yi[sO2] = 1.6165345089195457e-6 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 0.0;
      yi[sCH] = 1.2334072756372903e-9 * environment->metallicity;
      yi[sCH2] = 4.2719466865926816e-11 * environment->metallicity;
      yi[sOII] = 1.3990193853007009e-6 * environment->metallicity;
      yi[sOHII] = 4.480383830739206e-6 * environment->metallicity;
      yi[sH2OII] = 3.4372276239549154e-9 * environment->metallicity;
      yi[sH3OII] = 8.9175999873494229e-13 * environment->metallicity;
      yi[sO2II] = 0.0;
#endif
    }
    if(flag_reset == 12) {
      /* Molecular root found for the saved z=20, Z=1, T=1308 K failure. */
      yi[sELEC] = 3.3816627299155693e-8;
      yi[sHI] = 0.090029236423097111;
      yi[sHII] = 5.0973636985557156e-9;
      yi[sH2I] = 0.45495950168833826;
      yi[sHM] = 1.5124297197929668e-14;
#ifdef HD
      yi[sDI] = 2.3271032762120112e-6;
      yi[sDII] = 1.0732861297366062e-13;
      yi[sHDI] = 2.7278159774354116e-5;
#endif
      yi[sCI] = 3.3904397142241463e-8 * environment->metallicity;
      yi[sCII] = 2.2304322241103668e-9 * environment->metallicity;
      yi[sCO] = 6.9165214121594593e-5 * environment->metallicity;
      yi[sCO2] = 2.4679800839385235e-4 * environment->metallicity;
      yi[sOI] = 2.7774254938040671e-9 * environment->metallicity;
      yi[sOH] = 2.251429314960598e-8 * environment->metallicity;
      yi[sH2O] = 1.2186987979377911e-5 * environment->metallicity;
      yi[sO2] = 3.2604709313020105e-13 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 5.7265524463316557e-16 * environment->metallicity;
      yi[sCH] = 2.9756172098530313e-10 * environment->metallicity;
      yi[sCH2] = 3.4509289303746159e-10 * environment->metallicity;
      yi[sOII] = 5.2325901030760667e-17 * environment->metallicity;
      yi[sOHII] = 8.0564714123643722e-16 * environment->metallicity;
      yi[sH2OII] = 1.1360788922976149e-12 * environment->metallicity;
      yi[sH3OII] = 2.6487600249823491e-8 * environment->metallicity;
      yi[sO2II] = 1.412829757993675e-15 * environment->metallicity;
#endif
    }
    if(flag_reset == 13) {
      /* Molecular root found for the saved z=15, Z=0.1, T=1625 K
       * energy-trial failure; C/O abundances scale with metallicity. */
      yi[sELEC] = 4.5382337802913465e-6;
      yi[sHI] = 0.13212707000017435;
      yi[sHII] = 4.5367772036671654e-6;
      yi[sH2I] = 0.43392009380269708;
      yi[sHM] = 6.6903775210040103e-11;
#ifdef HD
      yi[sDI] = 3.2407443278998604e-6;
      yi[sDII] = 9.6176186589400081e-16;
      yi[sHDI] = 2.6364518829033115e-5;
#endif
      yi[sCI] = 2.9418491714483299e-8 * environment->metallicity;
      yi[sCII] = 1.4083976547487872e-9 * environment->metallicity;
      yi[sCO] = 7.269600596065384e-5 * environment->metallicity;
      yi[sCO2] = 2.4318552331419869e-4 * environment->metallicity;
      yi[sOI] = 1.1383006223979571e-9 * environment->metallicity;
      yi[sOH] = 6.3461958554081927e-8 * environment->metallicity;
      yi[sH2O] = 9.0685148668328659e-6 * environment->metallicity;
      yi[sO2] = 3.393002943944483e-6 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 3.0580833826225916e-16 * environment->metallicity;
      yi[sCH] = 6.8974084235700438e-9 * environment->metallicity;
      yi[sCH2] = 8.0746427048903346e-8 * environment->metallicity;
      yi[sOII] = 9.6314986665154312e-15 * environment->metallicity;
      yi[sOHII] = 2.8133263461932552e-12 * environment->metallicity;
      yi[sH2OII] = 4.0537034318655057e-11 * environment->metallicity;
      yi[sH3OII] = 1.3783036094629452e-8 * environment->metallicity;
      yi[sO2II] = 3.2894633155236885e-16 * environment->metallicity;
#endif
    }
    if(flag_reset == 14) {
      /* High-ionization grid root found for the saved z=0, Z=0.1,
       * T=13492 K energy trial; C/O abundances scale with metallicity. */
      yi[sELEC] = 0.94198017890404606;
      yi[sHI] = 0.05805441285023856;
      yi[sHII] = 0.94192033322456736;
      yi[sH2I] = 4.3520211046943833e-10;
      yi[sHM] = 1.3618147972046997e-8;
#ifdef HD
      yi[sDI] = 1.7239264614029546e-6;
      yi[sDII] = 2.7881336626535426e-5;
      yi[sHDI] = 6.9956361354992132e-14;
#endif
      yi[sCI] = 1.4553278290871734e-4 * environment->metallicity;
      yi[sCII] = 5.1832200636743317e-14 * environment->metallicity;
      yi[sCO] = 8.9818522821354362e-5 * environment->metallicity;
      yi[sCO2] = 8.0648662842740691e-5 * environment->metallicity;
      yi[sOI] = 4.1045415185868388e-6 * environment->metallicity;
      yi[sOH] = 2.2327059181335104e-14 * environment->metallicity;
      yi[sH2O] = 2.5110321146178189e-16 * environment->metallicity;
      yi[sO2] = 4.5117405652592979e-21 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 3.5097788223583212e-19 * environment->metallicity;
      yi[sCH] = 3.1373455903840561e-11 * environment->metallicity;
      yi[sCH2] = 1.8991432439961332e-15 * environment->metallicity;
      yi[sOII] = 6.7385275627295196e-5 * environment->metallicity;
      yi[sOHII] = 2.5239433431743299e-4 * environment->metallicity;
      yi[sH2OII] = 1.272010818433527e-15 * environment->metallicity;
      yi[sH3OII] = 5.5310536451555839e-20 * environment->metallicity;
      yi[sO2II] = 2.999343130095177e-15 * environment->metallicity;
#endif
    }
    if(flag_reset == 15) {
      /* Partly ionized root found for the saved z=0, Z=0.1,
       * T=2817 K energy trial; C/O abundances scale with metallicity. */
      yi[sELEC] = 0.12300980603699925;
      yi[sHI] = 0.85730009464603152;
      yi[sHII] = 0.12300585195318735;
      yi[sH2I] = 0.0098465571490409364;
      yi[sHM] = 3.0749809613547881e-8;
#ifdef HD
      yi[sDI] = 2.5412102424257962e-5;
      yi[sDII] = 3.5886705848171222e-6;
      yi[sHDI] = 6.0449014881965254e-7;
#endif
      yi[sCI] = 1.36278132667118e-11 * environment->metallicity;
      yi[sCII] = 4.1463444817672722e-7 * environment->metallicity;
      yi[sCO] = 6.5280680148379464e-5 * environment->metallicity;
      yi[sCO2] = 2.5030467177313828e-4 * environment->metallicity;
      yi[sOI] = 5.3695176281078972e-6 * environment->metallicity;
      yi[sOH] = 7.7137316555990897e-8 * environment->metallicity;
      yi[sH2O] = 1.0546855146720842e-7 * environment->metallicity;
      yi[sO2] = 5.4284453534823787e-9 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 2.4922831628171824e-15 * environment->metallicity;
      yi[sCH] = 0.0;
      yi[sCH2] = 0.0;
      yi[sOII] = 7.9884557638556961e-7 * environment->metallicity;
      yi[sOHII] = 2.7457484316224926e-6 * environment->metallicity;
      yi[sH2OII] = 2.4011672237160019e-9 * environment->metallicity;
      yi[sH3OII] = 7.4078187054886648e-13 * environment->metallicity;
      yi[sO2II] = 0.0;
#endif
    }
    if(flag_reset == 16) {
      /* Partly ionized root found for the saved z=0, Z=1,
       * T=3137 K energy trial. */
      yi[sELEC] = 0.24368227658044167;
      yi[sHI] = 0.74734107188133525;
      yi[sHII] = 0.24364822196622987;
      yi[sH2I] = 0.0044765617199466444;
      yi[sHM] = 2.3168149884005082e-8;
#ifdef HD
      yi[sDI] = 1.4037980521944389e-5;
      yi[sDII] = 5.8307747828085568e-6;
      yi[sHDI] = 9.7365078531417911e-6;
#endif
      yi[sCI] = 5.5117392151274453e-5 * environment->metallicity;
      yi[sCII] = 1.9187888638307844e-5 * environment->metallicity;
      yi[sCO] = 7.9794004867361318e-5 * environment->metallicity;
      yi[sCO2] = 1.6187818205856807e-4 * environment->metallicity;
      yi[sOI] = 8.0972686877595746e-5 * environment->metallicity;
      yi[sOH] = 4.7800103509553789e-5 * environment->metallicity;
      yi[sH2O] = 1.6282207832983737e-10 * environment->metallicity;
      yi[sO2] = 1.680877943286188e-5 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 0.0;
      yi[sCH] = 2.2457184144262032e-8 * environment->metallicity;
      yi[sCH2] = 7.5100344065396327e-11 * environment->metallicity;
      yi[sOII] = 9.0591189405508754e-6 * environment->metallicity;
      yi[sOHII] = 0.0;
      yi[sH2OII] = 0.0;
      yi[sH3OII] = 0.0;
      yi[sO2II] = 0.0;
#endif
    }
    if(flag_reset == 17) {
      /* Molecular root found for the saved z=0, Z=0.1,
       * T=2055 K energy trial; C/O abundances scale with metallicity. */
      yi[sELEC] = 1.7292161219077955e-4;
      yi[sHI] = 0.15260158073989594;
      yi[sHII] = 1.7291719422397354e-4;
      yi[sH2I] = 0.42359976420960388;
      yi[sHM] = 7.1611168879531515e-10;
#ifdef HD
      yi[sDI] = 3.7373343442930741e-6;
      yi[sDII] = 4.0629493336566334e-9;
      yi[sHDI] = 2.5863865864268007e-5;
#endif
      yi[sCI] = 1.6482734171834373e-8 * environment->metallicity;
      yi[sCII] = 4.7025756586749794e-11 * environment->metallicity;
      yi[sCO] = 5.7563937441663791e-5 * environment->metallicity;
      yi[sCO2] = 2.5841534867819711e-4 * environment->metallicity;
      yi[sOI] = 2.0715754337117436e-8 * environment->metallicity;
      yi[sOH] = 9.3338730701021626e-8 * environment->metallicity;
      yi[sH2O] = 4.806464506895454e-7 * environment->metallicity;
      yi[sO2] = 4.9502143698811777e-19 * environment->metallicity;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
      yi[sCOII] = 2.9752211192564415e-18 * environment->metallicity;
      yi[sCH] = 2.0528614763825189e-9 * environment->metallicity;
      yi[sCH2] = 2.1312587313702209e-9 * environment->metallicity;
      yi[sOII] = 2.159307716875683e-12 * environment->metallicity;
      yi[sOHII] = 8.2736084026850442e-11 * environment->metallicity;
      yi[sH2OII] = 2.1208936678149553e-9 * environment->metallicity;
      yi[sH3OII] = 8.4584764328809408e-9 * environment->metallicity;
      yi[sO2II] = 3.5896005208731092e-16 * environment->metallicity;
#endif
    }
#else
    if(flag_reset == 1){
      yi[sHI]    = 0.99999;
      yi[sELEC]  = 1.e-5;
      yi[sHII]   = 1.e-5;
      yi[sH2I]   = 1.e-10;
      yi[sHM]    = 1.e-15;
    }

    if(flag_reset == 2){
      yi[sHII]   = 0.99999;
      yi[sELEC]  = 0.99999;
      yi[sHI]    = 1.e-5;
      yi[sH2I]   = 1.e-10;
      yi[sHM]    = 1.e-15;
    }

    if(flag_reset == 3){
      yi[sHII]   = 1.e-10;
      yi[sELEC]  = 1.e-10;
      yi[sHI]    = 1.e-5;
      yi[sH2I]   = 0.499995;
      yi[sHM]    = 1.e-15;
    }
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
    if(flag_reset >= 2 && flag_reset <= 3){
      double carbon_total = A_C_sun * environment->metallicity;
      double oxygen_total = A_O_sun * environment->metallicity;
      yi[sCI] = yi[sCII] = yi[sCO] = yi[sCO2] = 0.0;
      yi[sOI] = yi[sOH] = yi[sH2O] = yi[sO2] = 0.0;
      if(flag_reset == 2) {
        yi[sCO] = dmin(carbon_total, oxygen_total);
        yi[sCII] = carbon_total - yi[sCO];
        yi[sO2] = 0.5 * (oxygen_total - yi[sCO]);
      } else if(flag_reset == 3) {
        yi[sCO2] = dmin(carbon_total, 0.5 * oxygen_total);
        yi[sCII] = carbon_total - yi[sCO2];
        yi[sO2] = 0.5 * (oxygen_total - 2.0 * yi[sCO2]);
      }
    }
#endif

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
  }

  const double convergence_abundance_threshold = CHEMISTRY_NEWTON_TRACE_ABUNDANCE;
  double ErrMaxCrit=1.e-6;

  if(FlagPrint)
    printf("T_gas=%g y=%g %g %g %g %g\n", T_gas, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM]);

  if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
    if(allow_restarts && flag_reset < deterministic_reset_count) {
      flag_reset++;
      flag_conv = 0;
      err_max = 100.0;
      goto chem_START;
    }
    if(allow_restarts && itr < max_iterations) {
      if(chemistry_find_random_newton_guess(nH, environment->metallicity,
                                            y_init, yi, &work->rates,
                                            &conservation, &random_state,
                                            &itr, max_iterations,
                                            &random_restart_count)) {
        flag_reset++;
        flag_conv = 0;
        err_max = 100.0;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
        use_numerical_jacobian = 0;
#endif
        goto chem_START;
      }
    }
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_set_failure(3);
#endif
    *state = state_old;
    return 0;
  }
#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  memcpy(basin_entry, yi, sizeof(basin_entry));
#endif
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    if(chemistry_species_is_algebraic(i_sp, &conservation))
      dy[i_sp] = 0.0;

  while(itr<=max_iterations && flag_conv < 2) {  // iteration
#ifdef CHEMISTRY_PROFILE
    double rhs_start = chemistry_profile_time();
#endif
    chemistry_network_rhs(CHEM_REACTION_IMPLICIT, nH, yi, work->rates.forward, work->rates.inverse, rf_tot);
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_NEWTON_RHS, chemistry_profile_time() - rhs_start);
#endif

    if(FlagPrint){
      printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g rf_tot=%g %g %g %g %g dt=%g\n"
	     , itr, T_gas, err_max, err_flag, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM], rf_tot[sELEC], rf_tot[sHI], rf_tot[sHII], rf_tot[sH2I], rf_tot[sHM], dt);
      printf("[%d] T_gas=%g y_H2=%g dy=%g resid=%g\n", itr, T_gas, yi[sH2I], yi[sH2I]-y_init[sH2I], rf_tot[sH2I]*dt);
    }
#ifdef ID_CHEM_DEBUG
    if(FlagPrint){
#ifdef HD
      printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g %g %g %g dt=%g\n", itr, T_gas, err_max, err_flag, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM], yi[sDI], yi[sDII], yi[sHDI], dt);
#else
      printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g\n", itr, T_gas, err_max, err_flag, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM]);
#endif
      printf("[%d] T_gas=%g y_H2=%g dy=%g resid=%g\n", itr, T_gas, yi[sH2I], yi[sH2I]-y_init[sH2I], rf_tot[sH2I]*dt);
    }
#endif

    int rebuild_newton_matrix = 1;
#ifdef CHEM_NEWTON_LU_REUSE
    /* Reuse the previous factorization only for the mandatory second
     * convergence-confirmation step, and only after a correction already
     * below 1e-12.  The stricter guard keeps the full non-equilibrium path
     * on the same thermochemical branch as a freshly rebuilt Jacobian. */
    rebuild_newton_matrix = !cached_newton_lu_valid || flag_conv != 1 ||
                            err_max >= 1.0e-12;
#endif
    memset(drfdy, 0, sizeof(drfdy));
    int use_analytic_jacobian = 0;
#ifdef CHEMISTRY_PROFILE
    double jacobian_start = chemistry_profile_time();
#endif
    if(rebuild_newton_matrix) {
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
    use_analytic_jacobian = !use_numerical_jacobian;
    if(use_analytic_jacobian) {
#ifdef CHEMISTRY_PROFILE
      double solver_columns_start = chemistry_profile_time();
#endif
      chemistry_network_jacobian(CHEM_REACTION_IMPLICIT, nH, yi, work->rates.forward, work->rates.inverse, drfdy);
#ifdef CHEMISTRY_PROFILE
      chemistry_profile_add(CHEM_PROFILE_JACOBIAN_SOLVER_COLUMNS, chemistry_profile_time() - solver_columns_start);
#endif
      if(!chemistry_reduce_jacobian_for_algebraic_species(nH, yi, &work->rates, solver_species, n_solver_species, &conservation, drfdy)) {
        use_numerical_jacobian = 1;
        use_analytic_jacobian = 0;
      }
    }
#endif

    if(!use_analytic_jacobian) {
      for(int solver_j=0;solver_j<n_solver_species;solver_j++) {
        int j_sp = solver_species[solver_j];
        double Dyi, epsilon_y = 1.0e-10, epsilon = 1.0e-4;
        if(fabs(yi[j_sp])<1.0e-100) Dyi = epsilon_y;
        else                       Dyi = epsilon*yi[j_sp];

        double rf_f[n_sp],rf_b[n_sp];
        int forward_valid = 0, backward_valid = 0;
        for(int attempt = 0; attempt < 20 && !forward_valid && !backward_valid; attempt++) {
          double y_tmp[n_sp];
          memcpy(y_tmp, yi, sizeof(y_tmp));
          y_tmp[j_sp] += 0.5 * Dyi;
          forward_valid = chemistry_update_implicit_algebraic_species(nH, y_tmp, &work->rates, &conservation);
          if(forward_valid)
            chemistry_network_rhs(CHEM_REACTION_IMPLICIT, nH, y_tmp, work->rates.forward, work->rates.inverse, rf_f);

          memcpy(y_tmp, yi, sizeof(y_tmp));
          y_tmp[j_sp] -= 0.5 * Dyi;
          backward_valid = y_tmp[j_sp] >= 0.0 &&
              chemistry_update_implicit_algebraic_species(nH, y_tmp, &work->rates, &conservation);
          if(backward_valid)
            chemistry_network_rhs(CHEM_REACTION_IMPLICIT, nH, y_tmp, work->rates.forward, work->rates.inverse, rf_b);

          if(!forward_valid && !backward_valid)
            Dyi *= 0.5;
        }
        if(!forward_valid && !backward_valid) {
#ifdef CHEMISTRY_PROFILE
          chemistry_profile_set_failure(5);
#endif
          *state = state_old;
          return 0;
        }

        for(int i_sp=0;i_sp<n_sp;i_sp++) {
          double drf, rf_big;
          if(forward_valid && backward_valid) {
            drfdy[j_sp][i_sp]=(rf_f[i_sp]-rf_b[i_sp])/Dyi;
            rf_big=dmax(fabs(rf_f[i_sp]),fabs(rf_b[i_sp]));
            drf=fabs(rf_f[i_sp]-rf_b[i_sp]);
          } else if(forward_valid) {
            drfdy[j_sp][i_sp]=2.0*(rf_f[i_sp]-rf_tot[i_sp])/Dyi;
            rf_big=dmax(fabs(rf_f[i_sp]),fabs(rf_tot[i_sp]));
            drf=fabs(rf_f[i_sp]-rf_tot[i_sp]);
          } else {
            drfdy[j_sp][i_sp]=2.0*(rf_tot[i_sp]-rf_b[i_sp])/Dyi;
            rf_big=dmax(fabs(rf_tot[i_sp]),fabs(rf_b[i_sp]));
            drf=fabs(rf_tot[i_sp]-rf_b[i_sp]);
          }
          if(rf_big!=0.0) {
            if(drf/rf_big<1.0e-15)
              drfdy[j_sp][i_sp]=0.0;
          }
        }
      }
    }
    }
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_NEWTON_JACOBIAN, chemistry_profile_time() - jacobian_start);
#endif

#ifdef ID_CHEM_DEBUG
    if(FlagPrint){
      for(j_sp=0;j_sp<n_sp;j_sp++)
	printf("[%d] drfdy[%d]=%g %g %g %g %g %g %g %g\n"
	       , itr, j_sp, drfdy[j_sp][0], drfdy[j_sp][1], drfdy[j_sp][2], drfdy[j_sp][3], drfdy[j_sp][4], drfdy[j_sp][5], drfdy[j_sp][6], drfdy[j_sp][7]);

    }
#endif

    if(rebuild_newton_matrix) {
      for(int solver_j=0;solver_j<n_solver_species;solver_j++) {
        int j_sp = solver_species[solver_j];
        for(int solver_i=0;solver_i<n_solver_species;solver_i++) {
          int i_sp = solver_species[solver_i];
	  if(solver_i!=solver_j) AA[solver_j * n_solver_species + solver_i] =    -dt*drfdy[j_sp][i_sp];
	  else                   AA[solver_j * n_solver_species + solver_i] = 1.0-dt*drfdy[j_sp][i_sp];
        }
      }
    }

    for(int solver_i=0;solver_i<n_solver_species;solver_i++) {
      int i_sp = solver_species[solver_i];
      solver_rhs[solver_i] = rf_tot[i_sp]*dt - dy[i_sp];
    }

#ifdef CHEMISTRY_PROFILE
    double linear_solve_start = chemistry_profile_time();
#endif
#ifdef CHEM_NEWTON_LU_REUSE
    if(rebuild_newton_matrix) {
      dgetrf_(&n_solver_species, &n_solver_species, AA, &n_solver_species,
              ipiv, &info);
      cached_newton_lu_valid = info == 0;
    }
    if(cached_newton_lu_valid) {
      char trans = 'N';
      dgetrs_(&trans, &n_solver_species, &inc, AA, &n_solver_species, ipiv,
              solver_rhs, &n_solver_species, &info);
    }
#else
    dgesv_(&n_solver_species, &inc, AA, &n_solver_species, ipiv, solver_rhs, &n_solver_species, &info);    //compile with -llapack -lblas
#endif
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_NEWTON_LINEAR_SOLVE, chemistry_profile_time() - linear_solve_start);
#endif

    memset(ddy, 0, sizeof(ddy));
    for(int solver_i=0;solver_i<n_solver_species;solver_i++) {
      ddy[solver_species[solver_i]] = solver_rhs[solver_i];
    }

    err_max=0.0;
#ifdef CHEMISTRY_PROFILE
    double newton_step_start = chemistry_profile_time();
#endif
    int invalid_species, step_was_clipped, clipped_species[N_SP];
    if(!chemistry_apply_implicit_newton_step(nH, yi, ddy, y_min, y_max,
                                             &work->rates, &conservation,
                                             !allow_restarts, &invalid_species,
                                             &step_was_clipped,
                                             clipped_species)) {
      if(allow_restarts && flag_reset < deterministic_reset_count) {
        flag_reset++;
        flag_conv = 0;
        err_max = 100.0;
        goto chem_START;
      }
#ifdef CHEMISTRY_PROFILE
      if(chemistry_profile_last_failure() == 0)
        chemistry_profile_set_failure(6);
#endif
      if(allow_restarts && itr < max_iterations) {
        if(chemistry_find_random_newton_guess(nH, environment->metallicity,
                                              y_init, yi, &work->rates,
                                              &conservation, &random_state,
                                              &itr, max_iterations,
                                              &random_restart_count)) {
          flag_reset++;
          flag_conv = 0;
          err_max = 100.0;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
          use_numerical_jacobian = 0;
#endif
          goto chem_START;
        }
      }
      *state = state_old;
      return 0;
    }
#ifdef CHEM_NEWTON_LU_REUSE
    if(step_was_clipped) {
      cached_newton_lu_valid = 0;
    }
#endif
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_add(CHEM_PROFILE_NEWTON_STEP, chemistry_profile_time() - newton_step_start);
#endif

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] += ddy[i_sp]; //dy[] is recalculated at the end of the loop

    for(int i_sp=0;i_sp<n_sp;i_sp++) {
      if(chemistry_species_is_algebraic(i_sp, &conservation)) continue;
#if defined(CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION) && !defined(CHEM_EQUILIBRIUM_SPECIES_CONVERGENCE)
      if(chemistry_is_full_implicit_initialized_intermediate(i_sp)) continue;
#endif
      double err = 0.;
      if(yi[i_sp]!=0.0) err=fabs(ddy[i_sp]/yi[i_sp]);
      if(err > err_max && fabs(yi[i_sp]) >= convergence_abundance_threshold)
        err_max = err, err_flag = i_sp;
    }
#ifdef CHEMISTRY_PROFILE
    if(err_max > 0.0)
      chemistry_profile_record_newton_limiter(err_flag, nH);
#endif

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
    for(int i_sp = 0; i_sp < n_sp; i_sp++)
      if(chemistry_species_is_algebraic(i_sp, &conservation))
        dy[i_sp] = 0.0;

    int maximum_species_clip_streak = 0;
    int convergence_relevant_step_was_clipped = 0;
    for(int i_sp = 0; i_sp < n_sp; i_sp++) {
      double abundance_before_step = yi[i_sp] - ddy[i_sp];
      int convergence_relevant_clip = clipped_species[i_sp] &&
        dmax(fabs(abundance_before_step), fabs(yi[i_sp])) >= convergence_abundance_threshold;
#if defined(CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION) && !defined(CHEM_EQUILIBRIUM_SPECIES_CONVERGENCE)
      if(chemistry_is_full_implicit_initialized_intermediate(i_sp))
        convergence_relevant_clip = 0;
#endif
      species_clip_streak[i_sp] = convergence_relevant_clip ? species_clip_streak[i_sp] + 1 : 0;
      if(convergence_relevant_clip)
        convergence_relevant_step_was_clipped = 1;
#ifdef CHEMISTRY_PROFILE
      if(convergence_relevant_clip)
        chemistry_profile_record_clipped_species(i_sp);
#endif
      if(species_clip_streak[i_sp] > maximum_species_clip_streak)
        maximum_species_clip_streak = species_clip_streak[i_sp];
    }
    if(convergence_relevant_step_was_clipped) {
      convergence_hold_after_clipping = 2;
#ifdef CHEMISTRY_PROFILE
      profile_clipped_steps++;
#endif
    }

    /* Reset 1 is a known physical basin whose normal approach to a root may
     * contain more than five limited Newton steps. */
    if(allow_restarts && maximum_species_clip_streak >= 5 &&
       flag_reset != 1
      ) {
      flag_conv = 0;
      err_max = 100.0;
      if(flag_reset < deterministic_reset_count) {
        flag_reset++;
        goto chem_START;
      }
      if(itr < max_iterations) {
        if(chemistry_find_random_newton_guess(nH, environment->metallicity,
                                              y_init, yi, &work->rates,
                                              &conservation, &random_state,
                                              &itr, max_iterations,
                                              &random_restart_count)) {
          flag_reset++;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
          use_numerical_jacobian = 0;
#endif
          goto chem_START;
        }
      }
    }

    if(itr >= 20) ErrMaxCrit = 1.e-4;
    if(convergence_relevant_step_was_clipped) {
      flag_conv = 0;
    } else if(convergence_hold_after_clipping > 0) {
      convergence_hold_after_clipping--;
      flag_conv = 0;
    } else {
      if(err_max < ErrMaxCrit) flag_conv++;
      else flag_conv = 0;
    }

    /* Do not spend the full Newton budget in a stagnant deterministic basin.
     * The initial state and deterministic resets each get 70 iterations.  A random
     * basin is not discarded on this clock; it may need a longer approach. */
    if(allow_restarts && flag_conv == 0 &&
       itr - basin_iteration_start + 1 >= 70 && flag_reset <= deterministic_reset_count) {
      flag_conv = 0;
      err_max = 100.0;
      if(flag_reset < deterministic_reset_count) {
        flag_reset++;
        goto chem_START;
      }
      if(itr < max_iterations) {
        if(chemistry_find_random_newton_guess(nH, environment->metallicity,
                                              y_init, yi, &work->rates,
                                              &conservation, &random_state,
                                              &itr, max_iterations,
                                              &random_restart_count)) {
          flag_reset++;
          goto chem_START;
        }
      }
    }

    itr++;
  }   // iteration

  if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
#ifdef CHEMISTRY_PROFILE
    chemistry_profile_set_failure(7);
#endif
    *state = state_old;
    return 0;
  }

  double convergence_error = flag_conv >= 2 ? err_max : 1.0;
#ifdef CHEMISTRY_PROFILE
  if(flag_conv >= 2)
    chemistry_profile_record_newton_solve(itr - 1, flag_reset,
                                          profile_clipped_steps);
#endif

#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  if(flag_conv >= 2) {
    printf("CHEM_RANDOM_SUCCESS iterations=%d reset=%d random_restarts=%d\n",
           itr, flag_reset, random_restart_count);
    for(int diagnostic_species = 0; diagnostic_species < N_SP;
        diagnostic_species++)
      printf("guess_%s=%.17g\n",
             chemistry_species_definition(diagnostic_species)->name,
             basin_entry[diagnostic_species]);
  }
#endif
  if(convergence_error<1.0e-6) {
    double t_chem1 = 1.0e20;
    for(int i_sp=0;i_sp<n_sp;i_sp++) {
      if(chemistry_species_is_algebraic(i_sp, &conservation))
        continue;
      if(yi[i_sp]-y_init[i_sp]!=0.0) {
	double t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
	if(t_chem2<t_chem1) t_chem1 = t_chem2;
      }

      work->evaluation.t_chem = t_chem1;
    }
  }else{
    if(convergence_error>1.e-4){
      if(!allow_restarts) {
        work->evaluation.err_chem = convergence_error;
        *state = state_old;
        return 0;
      }
      if(!report_failure) {
        work->evaluation.err_chem = convergence_error;
#ifdef CHEMISTRY_PROFILE
        chemistry_profile_set_failure(8);
#endif
        *state = state_old;
        return 0;
      }
      printf("iteration is not converged err=%le(in %d),itr=%d,"
	     " id=%d, reset=%d, random_restarts=%d, allow_restarts=%d.\n dt=%g, rho=%g, T_gas=%g, n_nuclei=%g;\n"
	     "y[%d]=%g, y[%d]=%g, y[%d]=%g, ",
	     err_max, err_flag, itr,
	     environment->ID, flag_reset, random_restart_count, allow_restarts,
	     dt, rho_gas, T_gas, n_nuclei,
	     sELEC, yi[sELEC], sHI, yi[sHI], sHII, yi[sHII]);
#ifdef HELIUM
      printf("y[%d]=%g, y[%d]=%g, y[%d]=%g, ", sHeI, yi[sHeI], sHeII, yi[sHeII], sHeIII, yi[sHeIII]);
#endif
      printf("y[%d]=%g, ", sH2I, yi[sH2I]);
#ifdef H2_PLUS
      printf("y[%d]=%g, ", sH2II, yi[sH2II]);
#endif
      printf("y[%d]=%g\n", sHM, yi[sHM]);

#ifdef HD
      printf("y[%d]=%g, y[%d]=%g y[%d]=%g", sDI, yi[sDI], sDII, yi[sDII], sHDI, yi[sHDI]);
#ifdef HD_PLUS
      printf(" y[%d]=%g", sHDII, yi[sHDII]);
#endif
#ifdef HD_MINUS
      printf(" y[%d]=%g", sDM, yi[sDM]);
#endif
      printf("\n");
#endif

      printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g, ",
	     sELEC, y_init[sELEC], sHI, y_init[sHI], sHII, y_init[sHII]);
#ifdef HELIUM
      printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g, ", sHeI, y_init[sHeI], sHeII, y_init[sHeII], sHeIII, y_init[sHeIII]);
#endif
      printf("yi[%d]=%g, ", sH2I, y_init[sH2I]);
#ifdef H2_PLUS
      printf("yi[%d]=%g, ", sH2II, y_init[sH2II]);
#endif
      printf("yi[%d]=%g\n", sHM, y_init[sHM]);

#ifdef HD
      printf("yi[%d]=%g, yi[%d]=%g yi[%d]=%g", sDI, y_init[sDI], sDII, y_init[sDII], sHDI, y_init[sHDI]);
#ifdef HD_PLUS
      printf(" yi[%d]=%g", sHDII, y_init[sHDII]);
#endif
#ifdef HD_MINUS
      printf(" yi[%d]=%g", sDM, y_init[sDM]);
#endif
      printf("\n");
#endif

      printf("CHEM_FAILURE_INPUT_BEGIN\n");
      printf("n_nuclei_cm3=%.17g\n", n_nuclei);
      printf("T_K=%.17g\n", T_gas);
      printf("T_CMB_K=%.17g\n", environment->cmb_temperature);
      printf("Z=%.17g\n", environment->metallicity);
      printf("dt_s=%.17g\n", dt);
      for(int failure_species = 0; failure_species < N_SP; failure_species++)
        printf("x0_%s=%.17g\n",
               chemistry_species_definition(failure_species)->name,
               state_old.abund[failure_species] / XH);
      printf("CHEM_FAILURE_INPUT_END\n");

      {
        double final_rhs[N_SP];
        chemistry_network_rhs(CHEM_REACTION_IMPLICIT, nH, yi,
                              work->rates.forward, work->rates.inverse,
                              final_rhs);
        printf("CHEM_FAILURE_SPECIES_BEGIN\n");
        for(int failure_species = 0; failure_species < N_SP; failure_species++) {
          double implicit_residual = final_rhs[failure_species] * dt -
              (yi[failure_species] - y_init[failure_species]);
          double relative_step = fabs(ddy[failure_species]) /
              dmax(fabs(yi[failure_species]), convergence_abundance_threshold);
          printf("CHEM_FAILURE_SPECIES name=%s abundance=%.17g initial=%.17g"
                 " newton_step=%.17g relative_step=%.17g residual=%.17g"
                 " algebraic=%d\n",
                 chemistry_species_definition(failure_species)->name,
                 yi[failure_species], y_init[failure_species],
                 ddy[failure_species], relative_step, implicit_residual,
                 chemistry_species_is_algebraic(failure_species, &conservation));
        }
        printf("CHEM_FAILURE_SPECIES_END\n");
      }

      work->evaluation.err_chem = convergence_error;
#ifdef CHEMISTRY_PROFILE
      chemistry_profile_set_failure(8);
#endif
      *state = state_old;
      return 0;
    }else{
      double t_chem1 = 1.0e20;
      for(int i_sp=0;i_sp<n_sp;i_sp++) {
	if(chemistry_species_is_algebraic(i_sp, &conservation))
	  continue;
	if(yi[i_sp]-y_init[i_sp]!=0.0) {
	  double t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
	  if(t_chem2<t_chem1) t_chem1 = t_chem2;
	}

	work->evaluation.t_chem = t_chem1;
      }
    }
  }

  work->evaluation.err_chem = convergence_error;

  for(int i_sp=0;i_sp<n_sp;i_sp++){
    if(chemistry_species_is_algebraic(i_sp, &conservation))
      continue;
    if(yi[i_sp] < y_min[i_sp]){
      yi[i_sp] = y_min[i_sp];
    }
  }

  /* Reapply floors after convergence, then reconstruct algebraic species
   * before exporting the state or evaluating chemical heating/cooling. */
  if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
    *state = state_old;
    return 0;
  }
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    if(chemistry_species_is_conservation_dependent(i_sp, &conservation) &&
       yi[i_sp] < 0.0 &&
       yi[i_sp] >= -CHEMISTRY_NEWTON_TRACE_ABUNDANCE)
      yi[i_sp] = 0.0;
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    state->abund[i_sp] = yi[i_sp] * XH;

  chemical_cooling(environment, state, work, dt, abund_init);

#ifdef CALC_CONTINUUM_OPACITY
  double ni[N_SP];
  for(int i=0;i<N_SP;i++){
    ni[i] = yi[i] * nH;
  }

#ifdef UPDATE_L_CORE_IN_SOLVER
  double l_core = get_core_radius(rho_gas, T_gas);
#else
  double l_core = step->jeans_length;
#endif
  double t_ff =  compute_freefall_time(rho_gas);
  solve_level_population_H_new(T_gas, ni, l_core, t_ff);
#endif

  return 1;
}

int chemistry_noneq_chem_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt)
{
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            NULL, MAX_ITER_CHEM, 1, 1);
}

static int chemistry_noneq_chem_implicit_trial(const ChemEnvironment *environment,
                                               const ChemStepData *step,
                                               ChemState *state,
                                               ChemWorkspace *work,
                                               double dt)
{
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            NULL, MAX_ITER_CHEM, 1, 0);
}

int chemistry_noneq_chem_implicit_with_guess(const ChemEnvironment *environment, const ChemStepData *step,
                                             ChemState *state, ChemWorkspace *work, double dt,
                                             const double *newton_guess, int max_iterations)
{
  return chemistry_noneq_chem_implicit_impl(environment, step, state, work, dt,
                                            newton_guess, max_iterations, 0, 0);
}


#pragma optimize("", on)
#undef MAX_ITER_CHEM


static void correct_initial_abundances(double *yi){
  /********** correct the total H abundances **********/
  if(yi[sHM] > 1.e-3) yi[sHM] = 1.e-10;

  double yHtot   = yi[sHI]   + yi[sHII]  + 2*yi[sH2I] + yi[sHM];
#ifdef H2_PLUS
  yHtot  += 2*yi[sH2II];
#endif

  if( fabs(1. - yHtot) > 1.e-3){
    double ratio1 = 1./yHtot;

    yi[sHI]   *= ratio1;
    yi[sHII]  *= ratio1;
    yi[sH2I]  *= ratio1;
#ifdef H2_PLUS
    yi[sH2II] *= ratio1;
#endif
    yi[sHM]   *= ratio1;
  }
  /****************************************************/

#ifdef HELIUM
  if(yi[sHeII] > 1.e10 * yi[sHeIII]){
    yi[sHeIII] = 1.e-5*yi[sHeII];
  }

  /********** correct the total He abundances **********/
  double yHetot  = yi[sHeI]  + yi[sHeII] + yi[sHeIII];

  if(yHetot > XHe/XH){
    double ratio1 = (XHe/XH) / yHetot;

    yi[sHeI]   *= ratio1;
    yi[sHeII]  *= ratio1;
    yi[sHeIII] *= ratio1;
  }
  /****************************************************/
#endif

#ifdef HD
  double yDtot = yi[sDI]   + yi[sDII]  + yi[sHDI];

#ifdef HD_PLUS
  yDtot += yi[sHDII];
#endif
#ifdef DMINUS
  yDtot += yi[sDM];
#endif

  /********** correct the total D abundances **********/
  if(yDtot > D_number_fraction /XH){
    double ratio1 = (D_number_fraction/XH) / yDtot;

    yi[sDI]   *= ratio1;
    yi[sDII]  *= ratio1;
    yi[sHDI]  *= ratio1;
#ifdef HD_PLUS
    yi[sHDII] *= ratio1;
#endif
#ifdef DMINUS
    yi[sDM]   *= ratio1;
#endif
  }
  /****************************************************/
#endif

  /*********  assure the charge neutrality  ***********/
  double y_minus = 0.0;
  double y_plus = 0.0;
  for(int species = 0; species < N_SP; species++) {
    int charge = chemistry_species_definition(species)->charge;
    if(charge < 0)
      y_minus -= charge * yi[species];
    else if(charge > 0)
      y_plus += charge * yi[species];
  }

  if(y_minus > 1.001 * y_plus){
    double minus_charge = y_minus - y_plus;

    if(minus_charge < 0.9 * yi[sELEC]){
      yi[sELEC] -= minus_charge;
    }else{
      minus_charge -= (0.9 * yi[sELEC]);
      yi[sELEC] *= 0.1;

#ifdef HELIUM
      if(yi[sHI] > 1.05 * minus_charge)
#endif
        {
          if(yi[sHI] > minus_charge){
            yi[sHI]  -= minus_charge;
            yi[sHII] += minus_charge;
          }else{
            yi[sH2I] -= 0.5*minus_charge;
            yi[sHII] += minus_charge;
          }
        }
#ifdef HELIUM
      else if(yi[sHI] + yi[sHeI] > 1.05*minus_charge){
        double diff1, diff2;

        diff1 = 0.97*yi[sHI];
        diff2 = minus_charge - diff1;

	yi[sHI]  -= diff1;
        yi[sHeI] -= diff2;

	yi[sHII]  += diff1;
        yi[sHeII] += diff2;
      }else if(yi[sHI] + 2*yi[sHeI] + yi[sHeII] > 1.05*minus_charge){
        double diff1, diff2, diff3;

        diff1 = 0.97*yi[sHI];
        diff2 = dmin(0.9*yi[sHeI], 0.5*(minus_charge - diff1)) ;
        diff3 = minus_charge - diff1 - 2*diff2;

	yi[sHI]   -= diff1;
        yi[sHeI]  -= diff2;
        yi[sHeII] -= diff3;

	yi[sHII]   += diff1;
        yi[sHeIII] += (diff3+diff2);
      }else{
        double diff1, diff2, diff3;

        diff1 = 0.97 * (yi[sHI] + 2*yi[sHeI] + yi[sHeII]);
        diff2 = minus_charge - diff1;

	yi[sELEC] -= diff2;

	diff1 = 0.9*yi[sHI];
        diff2 = 0.9*yi[sHeI];
        diff3 = 0.9*yi[sHeII];

	yi[sHI]   -= diff1;
        yi[sHeI]  -= diff2;
        yi[sHeII] -= diff3;

	yi[sHII]   += diff1;
        yi[sHeIII] += (diff3+diff2);
      }
#endif
    }
  }else if(y_plus > 1.00001 * y_minus){
    double plus_charge = y_plus - y_minus;
    yi[sELEC] += plus_charge;
  }
  /****************************************************/


  /********** Again, correct the total H, He, D abundances **********/
  if(yi[sHI] > 0.1){
    yi[sHI]     = 1. - yi[sHII] - 2*yi[sH2I] - yi[sHM];
#ifdef H2_PLUS
    yi[sHI]    -= 2*yi[sH2II];
#endif
  }else if(yi[sHII] > 0.1){
    yi[sHII]    = 1. - yi[sHI] - 2*yi[sH2I] - yi[sHM];
#ifdef H2_PLUS
    yi[sHII]    -= 2*yi[sH2II];
#endif
  }else{
    yi[sH2I]    = 0.5 * (1. - yi[sHI] - yi[sHII] - yi[sHM]);
#ifdef H2_PLUS
    yi[sHII]    -= yi[sH2II];
#endif
  }


#ifdef HELIUM
  if(yi[sHeI] > 0.1 * XHe/XH)
    yi[sHeI]   = XHe/XH - yi[sHeII] - yi[sHeIII];
  else if(yi[sHeII] > 0.1 * XHe/XH)
    yi[sHeII]  = XHe/XH - yi[sHeI] - yi[sHeIII];
  else
    yi[sHeII]  = XHe/XH - yi[sHeI] - yi[sHeIII];
#endif

#ifdef HD
  if(yi[sDI] > 0.1 * D_number_fraction/XH){
    yi[sDI]     = D_number_fraction/XH - yi[sDII] - yi[sHDI];
#ifdef HD_PLUS
    yi[sDI]    -= yi[sHDII];
#endif
#ifdef DMINUS
    yi[sDI]    -= yi[sDM];
#endif
  }else if(yi[sDII] > 0.1 * D_number_fraction/XH){
    yi[sDII]    = D_number_fraction/XH - yi[sDI] - yi[sHDI];
#ifdef HD_PLUS
    yi[sDII]   -= yi[sHDII];
#endif
#ifdef DMINUS
    yi[sDII]   -= yi[sDM];
#endif
  }else{
    yi[sHDI]    = D_number_fraction/XH - yi[sDI] - yi[sDII];
#ifdef HD_PLUS
    yi[sHDI]   -= yi[sHDII];
#endif
#ifdef DMINUS
    yi[sHDI]   -= yi[sDM];
#endif
  }
#endif
  /******************************************************************/
}

static void give_initial_guess(double T_gas, double nH, double dt, double *yi, double *y_init, double *k){
  if(T_gas > 1.e4){
    yi[sHI]   = 1.e-3;
    yi[sHII]  = 1. - yi[sHI];
    yi[sELEC] = yi[sHII];
    yi[sH2I]  = 1.e-5 * yi[sHI];
    yi[sHM]   = 1.e-10 * yi[sHI];
  }

#ifdef IONIZATION_PHOTON
  if(k[22] * dt > 0.1){
    double yHI_new = (y_init[sHI] + k[2] * dt * nH) / (1. + (k[22] + (2*k[2]+k[1]) * nH)*dt);

    if(yHI_new < yi[sHI]){
      yi[sHI]   = (y_init[sHI] + k[2] * dt * nH) / (1. + (k[22]+2*k[2]*nH+k[1]*nH)*dt);
      yi[sHII]  = 1. - yi[sHI];
      yi[sELEC] = yi[sHII];
      yi[sH2I]  = 1.e-5 * yi[sHI];
      yi[sHM]   = 1.e-10 * yi[sHI];
    }
  }
#endif
}


#endif /* CHEMISTRY_SC */
