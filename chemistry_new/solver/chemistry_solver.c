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
#include "chemistry_types.h"
#include "cooling_func.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "reaction_rates.h"
#include "chemistry_solver.h"
#include "saha.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

static void correct_initial_abundances(double *yi);
static void give_initial_guess(double temperature, double nH, double dt, double *yi, double *initial_guess, double *rates);
int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

#ifdef ENERGY_IMPLICIT
double give_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature, double *energy, double *residual);
double solve_bisection(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature_min, double residual_min, double temperature_max, double residual_max, double *temperature, double *energy, int max_iter);
static double commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double temperature, double *energy, double *residual);
#endif

#ifdef ENERGY_IMPLICIT
#define MAX_ITER 500
#define BETA_ENERGY 1.0

int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt){
  const double rho_gas = environment->rho_cgs;
  const double T_gas = work->thermodynamics.temperature;

  int iter = 0, ITER_reset = 30;
  double denom=0., err, err_chem;
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

      denom = g - g0;
      if(fabs(denom) < 1.e-10){
	if(denom > 0)
	  denom = 0.01;
	else
	  denom = -0.01;
      }

      dT = BETA_ENERGY*g*(T_old - T_new)/denom;

#ifdef ID_CHEM_DEBUG
      if(FlagPrint)
	printf("[%d] T=%g dT=%g err_chem=%g\n", j, T_new, dT, err_chem);
#endif

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
      /*else if(T_new + dT > 30000) T_new /= 1.1;*/
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

  commit_g(environment, step, state, work, chemistry_solver, dt, T_new, &e_new, &g);
  state->utherm_cgs = e_new;

 END:;
  if(isnan(err) || isnan(T_new))
    raise(11);

  if(err > 1.0e-4 || FlagPrint == 1){
    printf("not converged(in energy). err=%g id=%d T_old=%g T_new=%g Task=%d\n"
	   "All.ComovingIntegrationOn=%d, All.Time=%g, dt=%g, n_nuclei=%g, T_gas=%g, gamma=%g;\n"
	   "yi[0]=%g, yi[1]=%g, yi[2]=%g;\n"
	   "yi[%d]=%g, yi[%d]=%g;\n",
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
    *work = work_old;
    *state = state_old;
    return 1;
  }


  return 0;
}
#pragma opt


static double evaluate_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g, int commit)
{
  ChemWorkspace work_trial = *work;
  ChemState state_trial = *state;
  double e_old = state_trial.utherm_cgs;

  work_trial.thermodynamics.temperature = T_gas_new;
#ifdef SAHA
  if(chemistry_solver == CHEMISTRY_SOLVER_SAHA)
    chemistry_saha(environment, &state_trial, &work_trial, dt);
  else
#endif
    chemistry_noneq_chem_implicit(environment, step, &state_trial, &work_trial, dt);
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

  if(FlagPrint)
    printf("resid(e_new-e_old-de)=%g e_new=%g old=%g de=%g Tgas=%g rhoU_dot_chem=%g cool=%g rho_gas=%g dt=%g xmu=%g gamma=%g\n"
	   , *g, *e_new, e_old, de, T_gas_new, work_trial.evaluation.rhoU_dot_chem, work_trial.evaluation.rhoU_dot_rad, environment->rho_cgs, dt, xmu_new, gamma_new);

  if(commit){
    *state = state_trial;
    *work = work_trial;
  }

  return work_trial.evaluation.err_chem;
}

double give_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt, T_gas_new, e_new, g, 0);
}

static double commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, ChemistrySolver chemistry_solver, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(environment, step, state, work, chemistry_solver, dt, T_gas_new, e_new, g, 1);
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


static const enum SPECIES intermediate_equilibrium_species[] = {sHM};
#ifdef CHEM_INTERMEDIATE_EQUILIBRIUM
static const size_t intermediate_equilibrium_species_count = sizeof(intermediate_equilibrium_species) / sizeof(intermediate_equilibrium_species[0]);
#else
static const size_t intermediate_equilibrium_species_count = 0;
#endif

static int chemistry_species_is_intermediate_equilibrium(enum SPECIES species)
{
  for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
    if(intermediate_equilibrium_species[i] == species)
      return 1;
  return 0;
}

static int chemistry_update_intermediate_equilibrium_species(unsigned char mode, double nH, double *abundance, const ChemRates *rates)
{
  double previous[N_SP];
  double updated[sizeof(intermediate_equilibrium_species) / sizeof(intermediate_equilibrium_species[0])];

  if(intermediate_equilibrium_species_count == 0)
    return 1;

  if(intermediate_equilibrium_species_count == 1) {
    enum SPECIES species = intermediate_equilibrium_species[0];
    if(!chemistry_network_intermediate_equilibrium_abundance(species, mode, nH, abundance, rates->forward, rates->inverse, &updated[0], NULL))
      return 0;
    abundance[species] = updated[0];
    return 1;
  }

  for(int iteration = 0; iteration < 32; iteration++) {
    memcpy(previous, abundance, sizeof(previous));
    double max_change = 0.0;

    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++) {
      enum SPECIES species = intermediate_equilibrium_species[i];
      if(!chemistry_network_intermediate_equilibrium_abundance(species, mode, nH, previous, rates->forward, rates->inverse, &updated[i], NULL))
        return 0;
      double change = fabs(updated[i] - previous[species]) / (fabs(updated[i]) + 1.e-80);
      if(change > max_change)
        max_change = change;
    }

    for(size_t i = 0; i < intermediate_equilibrium_species_count; i++)
      abundance[intermediate_equilibrium_species[i]] = updated[i];

    if(max_change < 1.e-12)
      return 1;
  }
  return 0;
}

#define CHEM_MAX_CONSERVATION_CONSTRAINTS 4

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

static int chemistry_build_conservation_system(ChemConservationSystem *conservation)
{
  memset(conservation, 0, sizeof(*conservation));

  chemistry_add_conservation_constraint(conservation, sHI, 1.0);
  conservation->coefficient[0][sHI] = 1.0;
  conservation->coefficient[0][sHII] = 1.0;
  conservation->coefficient[0][sH2I] = 2.0;

#ifdef HELIUM
  int helium_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, sHeI, XHe / XH);
  conservation->coefficient[helium_constraint][sHeI] = 1.0;
  conservation->coefficient[helium_constraint][sHeII] = 1.0;
  conservation->coefficient[helium_constraint][sHeIII] = 1.0;
#endif

#ifdef HD
  int deuterium_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, sDI, D_number_fraction / XH);
  conservation->coefficient[deuterium_constraint][sDI] = 1.0;
  conservation->coefficient[deuterium_constraint][sDII] = 1.0;
  conservation->coefficient[deuterium_constraint][sHDI] = 1.0;
#endif

  int charge_constraint = conservation->count;
  chemistry_add_conservation_constraint(conservation, sELEC, 0.0);
  conservation->coefficient[charge_constraint][sELEC] = -1.0;
  conservation->coefficient[charge_constraint][sHII] = 1.0;
#ifdef HELIUM
  conservation->coefficient[charge_constraint][sHeII] = 1.0;
  conservation->coefficient[charge_constraint][sHeIII] = 2.0;
#endif

  for(int constraint = 0; constraint < conservation->count; constraint++) {
    enum SPECIES dependent = conservation->dependent_species[constraint];
    double diagonal = conservation->coefficient[constraint][dependent];
    if(diagonal == 0.0)
      return 0;
    conservation->dependent_inverse[constraint * conservation->count + constraint] = 1.0 / diagonal;
  }
  return 1;
}

static int chemistry_species_is_conservation_dependent(enum SPECIES species, const ChemConservationSystem *conservation)
{
  return conservation->is_dependent[species];
}

static int chemistry_species_is_algebraic(enum SPECIES species, const ChemConservationSystem *conservation)
{
  return chemistry_species_is_intermediate_equilibrium(species) || chemistry_species_is_conservation_dependent(species, conservation);
}

static int chemistry_update_implicit_algebraic_species(double nH, double *abundance, const ChemRates *rates, const ChemConservationSystem *conservation)
{
  abundance[sHI] = 1.0 - abundance[sHII] - 2.0 * abundance[sH2I];
#ifdef HELIUM
  abundance[sHeI] = XHe / XH - abundance[sHeII] - abundance[sHeIII];
#endif
#ifdef HD
  abundance[sDI] = D_number_fraction / XH - abundance[sDII] - abundance[sHDI];
#endif
  abundance[sELEC] = abundance[sHII];
#ifdef HELIUM
  abundance[sELEC] += abundance[sHeII] + 2.0 * abundance[sHeIII];
#endif

  if(!isfinite(abundance[sHI]) || !isfinite(abundance[sELEC])
#ifdef HELIUM
     || !isfinite(abundance[sHeI])
#endif
#ifdef HD
     || !isfinite(abundance[sDI])
#endif
    )
    return 0;

  return chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, abundance, rates);
}

static int chemistry_apply_implicit_newton_step(double nH, double *abundance, double *delta, const double *upper_bound, const ChemRates *rates, const ChemConservationSystem *conservation)
{
  double initial[N_SP];
  memcpy(initial, abundance, sizeof(initial));

  double step_scale = 1.0;
  for(int attempt = 0; attempt < 24; attempt++) {
    double trial[N_SP];
    for(int species = 0; species < N_SP; species++)
      trial[species] = initial[species] + step_scale * delta[species];

    int valid = chemistry_update_implicit_algebraic_species(nH, trial, rates, conservation);
    for(int species = 0; valid && species < N_SP; species++)
      if(!isfinite(trial[species]) || trial[species] < 0.0 || trial[species] > upper_bound[species])
        valid = 0;

    if(valid) {
      for(int species = 0; species < N_SP; species++) {
        delta[species] = trial[species] - initial[species];
        abundance[species] = trial[species];
      }
      return 1;
    }
    step_scale *= 0.5;
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

  for(int algebraic = 0; algebraic < n_algebraic; algebraic++) {
    enum SPECIES species = algebraic_species[algebraic];
    chemistry_network_jacobian_column(species, nH, abundance, rates->forward, rates->inverse, drfdy[species]);
  }

  double algebraic_response[N_SP * N_SP];
  memset(algebraic_response, 0, sizeof(algebraic_response));

  for(int solver_j = 0; solver_j < n_solver_species; solver_j++) {
    enum SPECIES column_species = solver_species[solver_j];
    for(int dependent = 0; dependent < conservation->count; dependent++)
      for(int constraint = 0; constraint < conservation->count; constraint++)
        algebraic_response[solver_j * n_algebraic + n_intermediate + dependent] +=
            conservation->dependent_inverse[constraint * conservation->count + dependent] *
            conservation->coefficient[constraint][column_species];
  }

  if(n_intermediate > 0) {
    double intermediate_jacobian[N_SP * N_SP];
    double intermediate_response[N_SP * N_SP];
    int intermediate_pivots[N_SP];

    for(int intermediate_j = 0; intermediate_j < n_intermediate; intermediate_j++) {
      enum SPECIES column_species = intermediate_equilibrium_species[intermediate_j];
      for(int intermediate_i = 0; intermediate_i < n_intermediate; intermediate_i++) {
        enum SPECIES row_species = intermediate_equilibrium_species[intermediate_i];
        intermediate_jacobian[intermediate_j * n_intermediate + intermediate_i] = drfdy[column_species][row_species];
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

    int info;
    dgesv_(&n_intermediate, &n_solver_species, intermediate_jacobian, &n_intermediate, intermediate_pivots, intermediate_response, &n_intermediate, &info);
    if(info != 0)
      return 0;

    for(int solver_j = 0; solver_j < n_solver_species; solver_j++)
      for(int intermediate = 0; intermediate < n_intermediate; intermediate++)
        algebraic_response[solver_j * n_algebraic + intermediate] = intermediate_response[solver_j * n_intermediate + intermediate];
  }

  for(int solver_j = 0; solver_j < n_solver_species; solver_j++) {
    enum SPECIES column_species = solver_species[solver_j];
    for(int solver_i = 0; solver_i < n_solver_species; solver_i++) {
      enum SPECIES row_species = solver_species[solver_i];
      for(int algebraic_k = 0; algebraic_k < n_algebraic; algebraic_k++)
        drfdy[column_species][row_species] -= drfdy[algebraic_species[algebraic_k]][row_species] * algebraic_response[solver_j * n_algebraic + algebraic_k];
    }
  }

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
};

static void chemistry_explicit_enforce_constraints(ChemState *state)
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
}

int chemistry_noneq_exp_nometal(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt)
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

    chemistry_explicit_enforce_constraints(state);
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

#define MAX_ITER_CHEM 1000
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

int chemistry_noneq_chem_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt)
{
  enum SPECIES species;
  const ChemState state_old = *state;

  int n_sp = N_SP;
  int inc  = 1, info;

  const double rho_gas = environment->rho_cgs;
  const double n_nuclei  = rho_gas/PROTONMASS;
  const double nH = XH * n_nuclei;
  const double T_gas = work->thermodynamics.temperature;

  evaluate_reaction_rates(environment, step, state, work);

  double yi[n_sp];

  yi[sELEC]  = state->abund[sELEC]/XH;
  yi[sHI]    = state->abund[sHI]/XH;
  yi[sHII]   = state->abund[sHII]/XH;
#ifdef HELIUM
  yi[sHeI]   = state->abund[sHeI]/XH;
  yi[sHeII]  = state->abund[sHeII]/XH;
  yi[sHeIII] = state->abund[sHeIII]/XH;
#endif

  yi[sH2I]   = state->abund[sH2I]/XH;
#ifdef H2_PLUS
  yi[sH2II]  = state->abund[sH2II]/XH;
#endif
  yi[sHM]    = state->abund[sHM]/XH;
#ifdef HD
  yi[sDI]    = state->abund[sDI]/XH;
  yi[sDII]   = state->abund[sDII]/XH;
  yi[sHDI]   = state->abund[sHDI]/XH;
#ifdef HD_PLUS
  yi[sHDII]  = state->abund[sHDII]/XH;
#endif
#ifdef DMINUS
  yi[sDM]    = state->abund[sDM]/XH;
#endif
#endif

  if(!chemistry_update_intermediate_equilibrium_species(CHEM_REACTION_IMPLICIT, nH, yi, &work->rates)) {
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

  y_min[sELEC] = 1.e-80,  y_mean[sELEC]  = 1.e-15   ,  y_max[sELEC]  = 2.0;

  if(yi[sELEC] > y_min[sELEC])
    yi_ion = yi[sELEC]/y_mean[sELEC];
  else
    yi_ion = 1.;

  y_min[sHI]   = 1.e-50,  y_mean[sHI]    = 0.76          ,  y_max[sHI]    = 1.001;
  y_min[sHII]  = 1.e-30,  y_mean[sHII]   = 1.e-15        ,  y_max[sHII]   = 1.001;
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
  y_min[sDII]  = 1.e-30,  y_mean[sDII]   = 6.e-09 *yi_ion,  y_max[sDII]   = D_number_fraction/XH * 1.001;
  y_min[sHDI]  = 1.e-30,  y_mean[sHDI]   = 1.e-08        ,  y_max[sHDI]   = D_number_fraction/XH * 1.001;
#ifdef HD_PLUS
  y_min[sHDII] = 1.e-30,  y_mean[sHDII]  = 2.e-18 *yi_ion,  y_max[sHDII]  = 0.5*D_number_fraction/XH;
#endif
#ifdef DMINUS
  y_min[sDM]   = 1.e-30,  y_mean[sDM]    = 5.e-16 *yi_ion,  y_max[sDM]    = D_number_fraction/XH;
#endif
#endif

  for(int i_sp=0;i_sp<n_sp;i_sp++) {
    y_mean[i_sp] /= XH;
    y_init[i_sp] = yi[i_sp];

    if(yi[i_sp] < y_min[i_sp]/* || yi[i_sp] > y_max[i_sp]*/){
      y_init[i_sp] = y_mean[i_sp];
      yi[i_sp]     = y_mean[i_sp];
    }
  }

  correct_initial_abundances(yi);
  if(!chemistry_build_conservation_system(&conservation) ||
     !chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
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

  //give_initial_guess(T_gas, nH, dt, yi, y_init, k);
  for(int i_sp=0;i_sp<n_sp;i_sp++){
    dy[i_sp] = yi[i_sp] - y_init[i_sp];
  }

  double err_max= 100.0;
  int flag_conv = 0, flag_reset = 0, itr = 1;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
  int use_numerical_jacobian = 0;
#endif

 chem_START:
  if(flag_reset){
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

    if(flag_reset >= 4){
    }

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
  }

  double ErrMaxCrit=1.e-6;

  if(FlagPrint)
    printf("T_gas=%g y=%g %g %g %g %g\n", T_gas, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM]);

  if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
    *state = state_old;
    return 0;
  }
  for(int i_sp = 0; i_sp < n_sp; i_sp++)
    if(chemistry_species_is_algebraic(i_sp, &conservation))
      dy[i_sp] = 0.0;

  while(itr<=MAX_ITER_CHEM && flag_conv < 2) {  // iteration
    chemistry_network_rhs(nH, yi, work->rates.forward, work->rates.inverse, rf_tot);

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

    memset(drfdy, 0, sizeof(drfdy));
    int use_analytic_jacobian = 0;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
    use_analytic_jacobian = !use_numerical_jacobian;
    if(use_analytic_jacobian) {
      for(int solver_j=0;solver_j<n_solver_species;solver_j++) {
        int j_sp = solver_species[solver_j];
        chemistry_network_jacobian_column(j_sp, nH, yi, work->rates.forward, work->rates.inverse, drfdy[j_sp]);
      }
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

        double y_tmp[n_sp];
        for(int i_sp=0;i_sp<n_sp;i_sp++){
          if(i_sp==j_sp) y_tmp[i_sp]=yi[i_sp]+0.5*Dyi;
          else           y_tmp[i_sp]=yi[i_sp];
        }

        double rf_f[n_sp],rf_b[n_sp];
        if(!chemistry_update_implicit_algebraic_species(nH, y_tmp, &work->rates, &conservation)) {
          *state = state_old;
          return 0;
        }
        chemistry_network_rhs(nH, y_tmp, work->rates.forward, work->rates.inverse, rf_f);

        y_tmp[j_sp]=yi[j_sp]-0.5*Dyi;
        if(!chemistry_update_implicit_algebraic_species(nH, y_tmp, &work->rates, &conservation)) {
          *state = state_old;
          return 0;
        }
        chemistry_network_rhs(nH, y_tmp, work->rates.forward, work->rates.inverse, rf_b);

        for(int i_sp=0;i_sp<n_sp;i_sp++) {
          double drf, rf_big;
          drfdy[j_sp][i_sp]=(rf_f[i_sp]-rf_b[i_sp])/Dyi;
          if(fabs(rf_f[i_sp])>fabs(rf_b[i_sp])) rf_big=fabs(rf_f[i_sp]);
          else                                rf_big=fabs(rf_b[i_sp]);
          if(rf_big!=0.0) {
            drf=fabs(rf_f[i_sp]-rf_b[i_sp]);
            if(drf/rf_big<1.0e-15)
              drfdy[j_sp][i_sp]=0.0;
          }
        }
      }
    }

#ifdef ID_CHEM_DEBUG
    if(FlagPrint){
      for(j_sp=0;j_sp<n_sp;j_sp++)
	printf("[%d] drfdy[%d]=%g %g %g %g %g %g %g %g\n"
	       , itr, j_sp, drfdy[j_sp][0], drfdy[j_sp][1], drfdy[j_sp][2], drfdy[j_sp][3], drfdy[j_sp][4], drfdy[j_sp][5], drfdy[j_sp][6], drfdy[j_sp][7]);

    }
#endif

    for(int solver_j=0;solver_j<n_solver_species;solver_j++) {
      int j_sp = solver_species[solver_j];
      for(int solver_i=0;solver_i<n_solver_species;solver_i++) {
        int i_sp = solver_species[solver_i];
	if(solver_i!=solver_j) AA[solver_j * n_solver_species + solver_i] =    -dt*drfdy[j_sp][i_sp];
	else                   AA[solver_j * n_solver_species + solver_i] = 1.0-dt*drfdy[j_sp][i_sp];
      }
    }

    for(int solver_i=0;solver_i<n_solver_species;solver_i++) {
      int i_sp = solver_species[solver_i];
      solver_rhs[solver_i] = rf_tot[i_sp]*dt - dy[i_sp];
    }

    dgesv_(&n_solver_species, &inc, AA, &n_solver_species, ipiv, solver_rhs, &n_solver_species, &info);    //compile with -llapack -lblas

    memset(ddy, 0, sizeof(ddy));
    for(int solver_i=0;solver_i<n_solver_species;solver_i++)
      ddy[solver_species[solver_i]] = solver_rhs[solver_i];

    err_max=0.0;
    if(!chemistry_apply_implicit_newton_step(nH, yi, ddy, y_max, &work->rates, &conservation)) {
      *state = state_old;
      return 0;
    }

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] += ddy[i_sp]; //dy[] is recalculated at the end of the loop

    for(int i_sp=0;i_sp<n_sp;i_sp++) {
      if(chemistry_species_is_algebraic(i_sp, &conservation)) continue;
      double err = 0.;
      if(yi[i_sp]!=0.0) err=fabs(ddy[i_sp]/yi[i_sp]);
      if(err>err_max && yi[i_sp] > 1.0e-10) err_max=err, err_flag = i_sp;
    }

    for(int i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
    for(int i_sp = 0; i_sp < n_sp; i_sp++)
      if(chemistry_species_is_algebraic(i_sp, &conservation))
        dy[i_sp] = 0.0;

    if(itr >= 20) ErrMaxCrit = 1.e-4;
#if defined(CHEM_DIFF_ANALYTIC) && !defined(CHEM_DIFF_NUMERICAL)
    if(intermediate_equilibrium_species_count > 0 && itr >= 20 && err_max > ErrMaxCrit)
      use_numerical_jacobian = 1;
#endif
    if(err_max < ErrMaxCrit) flag_conv++;
    else flag_conv = 0;

    if(flag_reset <=3 && itr - 70*(flag_reset+1) > 0 && err_max > 1){
      flag_reset++;
      goto chem_START;
    }

    itr++;
  }   // iteration

  if(!chemistry_update_implicit_algebraic_species(nH, yi, &work->rates, &conservation)) {
    *state = state_old;
    return 0;
  }

  if(err_max<1.0e-6) {
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
    if(err_max>1.e-4){
      printf("iteration is not converged err=%le(in %d),itr=%d,"
	     " id=%d.\n dt=%g, rho=%g, T_gas=%g, n_nuclei=%g;\n"
	     "y[%d]=%g, y[%d]=%g, y[%d]=%g, ",
	     err_max, err_flag, itr,
	     environment->ID, dt, rho_gas, T_gas, n_nuclei,
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

      work->evaluation.err_chem = err_max;
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

  work->evaluation.err_chem = err_max;

  state->abund[sELEC] =yi[sELEC] *XH;
  state->abund[sHI]   =yi[sHI]   *XH;
  state->abund[sHII]  =yi[sHII]  *XH;
#ifdef HELIUM
  state->abund[sHeI]  =yi[sHeI]  *XH;
  state->abund[sHeII] =yi[sHeII] *XH;
  state->abund[sHeIII]=yi[sHeIII]*XH;
#endif

  state->abund[sH2I]  =yi[sH2I]  *XH;
#ifdef H2_PLUS
  state->abund[sH2II] =yi[sH2II] *XH;
#endif
  state->abund[sHM]   =yi[sHM]   *XH;
#ifdef HD
  state->abund[sDI]   =yi[sDI]   *XH;
  state->abund[sDII]  =yi[sDII]  *XH;
  state->abund[sHDI]  =yi[sHDI]  *XH;
#ifdef HD_PLUS
  state->abund[sHDII] =yi[sHDII] *XH;
#endif
#ifdef DMINUS
  state->abund[sDM]   =yi[sDM]   *XH;
#endif
#endif

  /* The Newton solver evaluates chemical heating/cooling from the converged,
     unclamped Newton solution.  Keep that ordering while making state carry
     the new abundances consumed by chemical_cooling(). */
  chemical_cooling(environment, state, work, dt, abund_init);

  for(int i_sp=0;i_sp<n_sp;i_sp++){
    if(chemistry_species_is_algebraic(i_sp, &conservation))
      continue;
    if(yi[i_sp] < y_min[i_sp]){
      yi[i_sp] = y_min[i_sp];
      state->abund[i_sp] = yi[i_sp] * XH;
    }
  }

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
  double y_minus = yi[sELEC] + yi[sHM];
  double y_plus  = yi[sHII];

#ifdef H2_PLUS
  y_plus += yi[sH2II];
#endif

#ifdef HELIUM
  y_plus += yi[sHeII] + 2*yi[sHeIII];
#endif
#ifdef HD
  y_plus  += yi[sDII];
#ifdef DMINUS
  y_minus += yi[sDM];
#endif
#ifdef HD_PLUS
  y_plus  += yi[sHDII];
#endif
#endif

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
