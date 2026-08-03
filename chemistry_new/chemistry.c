/*>=== Includes ===<*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

#include "../allvars.h"
#include "../proto.h"
#include "../tags.h"

#include "chemistry.h"
#include "chemistry_internal.h"
#include "chemistry_types.h"

#ifdef CHEMISTRY_SC
#define DEBUG_CHEMISTRY

#ifdef CHEM_HIRANO
#ifndef O05_CHEM
#define O05_CHEM 
#endif
#endif

/*>=== Definitions ===<*/
#define N_THREEBODY 1.0e14
#define N_EQCHEM    1.0e18
#define N_NOLINE    5.0e19

#define MAX_ITERATION 1000
#define HIGH_TEMPERATURE_THRESHOLD 1.0e5

#ifdef LYMANWERNER
/* ----- Photodissociation of H2 */
#define hv_ref 12.4   /* eV */
#define hv_Lym 13.6   /* eV */
#define T_eff  1.0e4  /*  K */
#endif

/* ----- Some grid bounds. */
#define T_start    0.1
#define T_end      1.0e6

/*=== Variables ===<*/

/* ----- Grid spacings. */
static double dln_T;

double chem_e[15]; //unit in kJ/mol and mainly from KIDA data base
int FlagPrint=0;

/*=== Subroutines ===<*/
int radiative_cooling_high_T(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *work);
int chemical_cooling(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *work, double dt, double *y_old);

int chemistry_noneq_exp_nometal(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
int chemistry_noneq_chem_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, double dt);
#ifdef SAHA
static int chemistry_saha(const ChemEnvironment *environment, ChemState *state, ChemWorkspace *workspace, double dt);
#endif
int chemistry_noneq_energy_implicit(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt);

#ifdef ENERGY_IMPLICIT
double give_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double T_gas, double *e_new, double *g);
double solve_bisection(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double T_min, double g_min, double T_max, double g_max, double *T_new0, double *e_new0, int max_iter);
static double commit_g(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *workspace, ChemistrySolver chemistry_solver, double dt, double T_gas, double *e_new, double *g);
#endif

static void initialize_step_data(const ChemEnvironment *environment, const ChemState *initial_state, double initial_temperature, double duration_s, ChemStepData *step);
static int evolve_high_temperature(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, double *remaining_time_s, ChemWorkspace *workspace);
static int evaluate_chemistry_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *workspace);
static int evaluate_reaction_rates(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
static int evaluate_radiative_cooling(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
static int evaluate_timesteps(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace);
static ChemSolverPlan select_solvers(double n_nuclei, double time_cooling, double time_chemistry, double time_remaining);
int chemistry_advance(int mode, const ChemEnvironment *environment, double duration_s, ChemState *state, ChemWorkspace *workspace);

static inline int ReactionRates_nometal(int itr, double nH, double *yi,double *k, double *k_inverse, double *rf_tot);

double get_core_radius(double rho_cgs, double T_gas);

#ifdef CHEM_DIFF_ANALYTIC
int (*switch_func[256])(double, double *, double *, double *, double *);

static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *k, double *k_inverse, double *rf_tot);
int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HI(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HII(double nH,double *ni,double *k,double *k_inverse,double *rf);
#ifdef HELIUM
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *k,double *k_inverse,double *rf);
#endif
int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *k,double *k_inverse,double *rf);
#ifdef H2_PLUS
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *k,double *k_inverse,double *rf);
#endif
int ReactionRate_diff_nometal_HM(double nH,double *ni,double *k,double *k_inverse,double *rf);
#ifdef HD
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_DII(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *k,double *k_inverse,double *rf);
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *k,double *k_inverse,double *rf);
#endif
int func_other(double nH,double *ni,double *k,double *k_inverse,double *rf);
static void init_switch_func(void);
#endif

#if defined(CHEMISTRY_IMPLICIT) && !defined(CHEM_DIFF_ANALYTIC) && !defined(ALLOW_CHEM_DIFF_NUMERIC)
#error "Finite-difference chemistry Jacobian is deprecated and very slow; define CHEM_DIFF_ANALYTIC or explicitly opt in with ALLOW_CHEM_DIFF_NUMERIC."
#endif

int chem_tabulate(void);
double compute_gamma_from_abund(const double *abund, double T_gas);
double chemistry_molecular_weight(const double *abund);
double chemistry_temperature(const double *abund, double gamma, double utherm_cgs);

double (*H2_cooling)(double, double, double, double, double, double, double);

double H2_cooling_O00(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_Y06(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_GP12(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_G15(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);

#ifdef HD
double (*HD_cooling)(double, double, double, double);

double HD_cooling_O00(double nH, double T_gas, double yHI, double yHD);
double HD_cooling_L05(double nH, double T_gas, double yHI, double yHDI);
double HD_cooling_GP12(double nH, double T_gas, double yHI, double yHDI);

#ifdef H2_OPACITY
double (*HD_line_escape_fraction)(double, double);
double fesc_HD_fit(double T_K, double xNc_HD);
#endif
#endif

double continuum_cooling_rate_thin(double T_gas, double rho_cgs, const double *abund);
double continuum_cooling_rate_thin_Lambda(double T_gas);

#ifdef H2_OPACITY
double (*H2_line_escape_fraction)(double, double);

double H2line_reduction(double T);
double fesc_fit(double T_K,double xNc_H2);
double fesc_unity(double T_K, double xNc_H2);
#endif

#ifdef CONT_OPACITY
double (*continuum_escape_fraction)(double, double, double, struct opacities *);

double calc_continuum_opacities(double T_gas, double rho_cgs, double l_core, struct opacities *opac);
double calc_continuum_opacities_CIE(double T_gas, double rho_cgs, double l_core, struct opacities *opac);
#endif

void tabulate_opacities_in_LW_band(void);

double (*partition_H2)(double);

double partition_H2_Rossi1985(double T_gas);
double partition_H2_Irwin1987(double T_gas);
double partition_H2_EBM1989(double T_K);
double partition_H2_EBM1989_fitted_by_Matsukoba(double T_K);
void check_partition_functions_H2(void);

double pf_z_H(double T_K);
double partition_H(double T_gas); 

double (*partition_H2p)(double);

double (*OxygenI_cooling)(double, double, double, double);
double (*CarbonI_cooling)(double, double, double, double);
double (*CarbonII_cooling)(double, double, double, double);

double OxygenI_fivelevel_cooling(double nOI, double ne, double nH, double T_gas);
double OxygenI_twolevel_cooling(double nOI, double ne, double nH, double T_gas);
double OxygenI_twolevel_cooling_thin(double nOI, double ne, double nH, double T_gas);
double CarbonI_twolevel_cooling(double nCII, double ne, double nH, double T_gas);
double CarbonI_threelevel_cooling(double nCI, double ne, double nH, double T_gas);
double CarbonII_twolevel_cooling(double nCII, double ne, double nH, double T_gas);
double CarbonII_twolevel_cooling_thin(double nCII, double ne, double nH, double T_gas);

#ifdef METAL_LINE_COOLING_SD93
void read_continuum_cooling_rate_table_SD93(void);
double get_cooling_rate_SD93(double T_gas);
#endif

double get_PlanckOpacity(double T_gas, double rho_gas);

#ifdef SHIELD_LW
double SelfShield(double ColumnDens, double T_gas, double b);
double EachShield(double ColumnDens, double alpha, double delta, double eta);
#endif

int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

void correct_initial_abundances(double *yi);
void give_initial_guess(double T_gas, double nH, double dt, double *yi, double *y_init, double *k);
void bound_ddy(double *yi, double *ddy, double *y_max, int itr, int flag_reset, double *err_max);

#ifdef CONTINUUM_COOLING_MATSUKOBA 
double continuum_cooling_rate(double T_gas, double rho_cgs, double l_core);
double H_criticaldensity(double xnH, double T_K,double xNc_H, double y_H,double y_e);
#endif
double fesc_Lya(double xnH,double T_K,double y_H,double y_e,double rr,double vel,double *cooling_rate_Lya,double *cooling_rate_2ph);

#ifdef SAHA
double chem_saha(double xnH,double T_K,double *y_H,double *y_H2,double *y_e,double *y_Hp);
double chem_saha_5sp(double xnH, double T_K, double *y);
#endif

double dust_temperature(double rho_cgs, double T_gas, double metallicity, double l_core, double stellar_flux, struct opacities *opac, double *edot_dust);
double H2form_on_grains(double rho_cgs, double T_gas, double T_dust, double metallicity);

#ifdef CONTINUUM_COOLING_MATSUKOBA
void read_opacity_tables(void);
void tabulate_cooling_rate_cont_thin(void);
#endif
#ifdef DUST_COOLING
void read_dust_opacity_tables(void);
#endif

#ifdef LYMANWERNER
void make_rate_tables_LW_radiation(void);
#ifdef METALLICITY
double dust_PEheating(double G0, double nH, double y_e, double T_gas, double metallicity);
#endif
#endif

enum REACTIONS{
	RE1, RE2,
#ifdef HELIUM
	RE3, RE4, RE5, RE6,
#endif
	RE7, RE8, 
#ifdef H2_PLUS
	RE9, RE10, 
#endif
	RE11, 
#ifdef H2_PLUS
	RE12, 
#endif
	RE13, 
#ifdef HELIUM
	RE14, RE15,
#endif
	RE16, RE17, RE18
#ifdef H2_PLUS
	,RE19, RE20, RE21
#endif
#ifdef IONIZATION_PHOTON
	,RE22
#endif

#ifdef LYMANWERNER
	,RE25, RE29
#endif	
	,RE30
	
#ifdef HD
	,
	RE31, RE34, 
	RE32, RE33, RE35,
	RE37, RE38, RE39
#ifdef HD_PLUS
	, RE36, RE40, RE41, RE42
#endif
#ifdef DMINUS
	, RE43, RE44, RE45, RE46, RE47, RE48
#endif
#ifdef LYMANWERNER
	, RE53
#endif
#endif

#ifdef H2_FORMATION_ON_DUST
	, RE51
#endif

#ifdef THREEBODY
	,RE55, RE57, RE58
#endif
	
#ifdef INVERSE_REACTIONS
	,RE1_INV , RE2_INV , RE7_INV , RE8_INV 
	,RE13_INV, RE16_INV, RE17_INV, RE18_INV
	,RE30_INV
#endif

	,REEND
};

int N_SP, N_RE;
int count_print;
const double XH  = 0.76;
const double XHe = 0.0625;
const double XHD = 3.e-5;

/* ----- Tables ------- */
double T_table[N_T];
double k_table[k_MAX][N_T];
double k_table_inverse[k_MAX][N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T], H2_CIE[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

FILE *open_chemistry_table(const char *filename)
{
  char path[2048];

  snprintf(path, sizeof(path), "%s/%s", All.ChemistryTableDir, filename);

  FILE *fp = fopen(path, "r");
  if(fp == NULL)
    printf("Can't open chemistry table %s (All.ChemistryTableDir=%s).\n", path, All.ChemistryTableDir);

  return fp;
}

double chemistry_temperature(const double *abund, double gamma, double utherm_cgs)
{
  double molecular_weight = chemistry_molecular_weight(abund);
  double egyspec_to_temp  = (gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;

  return utherm_cgs * egyspec_to_temp;
}

double compute_gamma_from_abund(const double *abund, double T_gas)
{
  double n_H= abund[sHI] + abund[sHII] + abund[sHM] + 2*abund[sH2I];
  
  double y_HI   = abund[sHI]   /n_H;
  double y_HII  = abund[sHII]  /n_H;
  double y_He   = XHe;

  double y_HM   = abund[sHM]   /n_H;
  double y_H2   = abund[sH2I]  /n_H;
  double y_elec = abund[sELEC] /n_H;
  
  double x = 6100.0/T_gas;

  double gamma_H2_minus1_inv;
#ifdef GAMMA_NEW
  if(T_gas < 550.){
    gamma_H2_minus1_inv
      = 1./(0.133333333333 * tanh((2.08661-log10(T_gas))*4.165) + 0.533333333333);
  }else{
    gamma_H2_minus1_inv = 0.5*(5.0 + 2.0*x*x * exp(x)/(exp(x)-1.0)/(exp(x)-1.0));
  }
#else
  if(x > 100.0)
    gamma_H2_minus1_inv = 2.5;
  else
    gamma_H2_minus1_inv = 0.5*(5.0 + 2.0*x*x * exp(x)/(exp(x)-1.0)/(exp(x)-1.0));
#endif
  
  double sum_y= y_HI + y_HII + y_HM + y_He + y_elec + y_H2;  
  double denominator=1.5*(y_HI + y_HII + y_HM + y_He + y_elec) + gamma_H2_minus1_inv * y_H2;
  
  double gamma = 1.0+sum_y/denominator;

  if(gamma > 1.8 || gamma < 1.2){
    printf("Gamma odd![task%03d] gamma=%g T_gas=%g\n"
	   , ThisTask, gamma, T_gas);
    fflush(stdout);
  }
    
  return gamma;
}

double chemistry_molecular_weight(const double *abund){
  double sum_fraction =  abund[sELEC] + abund[sHI] + abund[sHII] + abund[sHM] + abund[sH2I]
#ifdef H2_PLUS
    + abund[sH2II]
#endif
#ifdef HD
    + abund[sDI] + abund[sDII] + abund[sHDI]
#ifdef HD_PLUS
    + abund[sHDII]
#endif
#ifdef DMINUS
    + abund[sDM]
#endif
#endif
#ifdef HELIUM
    + abund[sHeI] + abund[sHeII] + abund[sHeIII]
#else
    + XHe
#endif
    ;
  
  return 1.0 / sum_fraction;
}

double compute_freefall_time(double rho_gas){
  return sqrt(3.0*M_PI/32.0/GRAVITY/rho_gas);
}

double get_core_radius(double rho_cgs, double T_gas){
#ifdef JEANS_LENGTH_TMAX
  if(T_gas > JEANS_LENGTH_TMAX)
    T_gas = JEANS_LENGTH_TMAX;
#endif

  return  0.5 * M_PI * sqrt(BOLTZMANN * T_gas/ PROTONMASS / 1.2)/ sqrt(GRAVITY * rho_cgs);
}

static void update_thermodynamic_state(const ChemState *state, ChemThermodynamicState *thermodynamics)
{
  thermodynamics->molecular_weight = chemistry_molecular_weight(state->abund);
  thermodynamics->temperature = state->utherm_cgs * (state->gamma - 1.0) *
      PROTONMASS * thermodynamics->molecular_weight / BOLTZMANN;
}

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

static void evaluate_chemical_energy_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *work)
{
  const double nH = XH * environment->rho_cgs / PROTONMASS;
  const double temperature = work->thermodynamics.temperature;
  
  const double *a = state->abund;
  const double n_elec = nH * a[sELEC] / XH;
  const double n_HI = nH * a[sHI] / XH;
  const double n_HM = nH * a[sHM] / XH;
  const double n_H2I = nH * a[sH2I] / XH;
#ifdef H2_PLUS
  const double n_H2II = nH * a[sH2II] / XH;
#endif
#ifdef HELIUM
  const double n_HeI = nH * a[sHeI] / XH;
  const double n_HeII = nH * a[sHeII] / XH;
#endif

  const double n_crit = (1.e6 / sqrt(temperature)) /
      (1.6 * a[sHI] / XH * exp(-pow(400.0 / temperature, 2.0)) +
       1.4 * a[sH2I] / XH * exp(-12000.0 / (temperature + 1200.0)));

  double edot_H2form = eV_to_erg * work->rates.forward[8] * n_HM * n_HI *
                       3.73 / (1.0 + n_crit / n_HI);
#ifdef H2_PLUS
  edot_H2form += eV_to_erg * work->rates.forward[10] * n_H2II * n_HI *
                 1.83 / (1.0 + n_crit / n_HI);
#endif
#ifdef H2_FORMATION_ON_DUST
  edot_H2form += eV_to_erg * work->rates.forward[51] * n_HI *
                 (0.2 + 4.28 / (1.0 + n_crit / n_HI));
#endif
#ifdef THREEBODY
  edot_H2form += 4.48 * eV_to_erg *
      (((work->rates.forward[55] * n_HI + work->rates.forward[57] * n_H2I) * n_HI * n_HI) -
       ((work->rates.forward[11] * n_HI + work->rates.forward[58] * n_H2I) * n_H2I));
#endif

  double energy_rate = edot_H2form -
                       13.598 * eV_to_erg * work->rates.forward[1] * n_HI * n_elec;
#ifdef HELIUM
  energy_rate -= 24.587 * eV_to_erg * work->rates.forward[3] * n_HeI * n_elec;
  energy_rate -= 54.416 * eV_to_erg * work->rates.forward[5] * n_HeII * n_elec;
#endif

  work->evaluation.rhoU_dot_chem = energy_rate;
  if(fabs(energy_rate) < 1.e-40) {
    work->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  } else {
    work->evaluation.t_cool_chem =
        state->utherm_cgs * environment->rho_cgs / energy_rate;
  }
}

static int evaluate_chemistry_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *workspace)
{
  const double *a = state->abund;
  const double n_nuclei = environment->rho_cgs / PROTONMASS;
  
  double electron_dot = (workspace->rates.forward[1] * a[sHI]
#ifdef HELIUM
      + workspace->rates.forward[3] * a[sHeI] + workspace->rates.forward[5] * a[sHeII]
      - workspace->rates.forward[4] * a[sHeII] - workspace->rates.forward[6] * a[sHeIII]
#endif
      - workspace->rates.forward[2] * a[sHII]) * a[sELEC] * n_nuclei * n_nuclei
      + (workspace->rates.forward[22] * a[sHI]
#ifdef HELIUM
         + workspace->rates.forward[23] * a[sHeII] + workspace->rates.forward[24] * a[sHeI]
#endif
         + workspace->rates.forward[25] * a[sHM]) * n_nuclei
      + (workspace->rates.forward[8] * a[sHI] * a[sHM]
         + workspace->rates.forward[16] * a[sHM] * a[sELEC]
         + workspace->rates.forward[17] * a[sHI] * a[sHM]
#ifdef H2_PLUS
         + workspace->rates.forward[19] * a[sHM] * a[sHII]
#endif
#ifdef DMINUS
         + workspace->rates.forward[48] * a[sDM] * a[sHI]
#endif
        ) * n_nuclei * n_nuclei
      - (workspace->rates.forward[7] * a[sHI]
#ifdef H2_PLUS
         + workspace->rates.forward[20] * a[sH2II]
#endif
#ifdef HD_PLUS
         + workspace->rates.forward[42] * a[sHDII]
#endif
#ifdef DMINUS
         + workspace->rates.forward[43] * a[sDI]
#endif
        ) * a[sELEC] * n_nuclei * n_nuclei;

  double HI_dot = (workspace->rates.forward[2] * a[sHII] * a[sELEC]
      + 2.0 * workspace->rates.forward[11] * a[sH2I] * a[sHI]
      + workspace->rates.forward[12] * a[sH2I] * a[sHII]
      + 2.0 * workspace->rates.forward[13] * a[sH2I] * a[sELEC]
#ifdef HELIUM
      + workspace->rates.forward[15] * a[sHeI] * a[sHII]
#endif
      ) * n_nuclei * n_nuclei
#ifdef THREEBODY
      + 2.0 * workspace->rates.forward[58] * a[sH2I] * a[sH2I] * n_nuclei * n_nuclei
#endif
      + workspace->rates.forward[25] * a[sHM] * n_nuclei
#ifdef H2_PLUS
      + workspace->rates.forward[26] * a[sH2II] * n_nuclei
#endif
      + 2.0 * workspace->rates.forward[29] * a[sH2I] * n_nuclei
      - ((workspace->rates.forward[1] * a[sELEC]
#ifdef H2_PLUS
          + workspace->rates.forward[10] * a[sH2II]
#endif
#ifdef HELIUM
          + workspace->rates.forward[14] * a[sHeII]
#endif
         ) * n_nuclei
#ifdef THREEBODY
         + (2.0 * workspace->rates.forward[55] * a[sHI] + 2.0 * workspace->rates.forward[57] * a[sH2I]) *
           a[sHI] * n_nuclei * n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
         + 2.0 * workspace->rates.forward[51]
#endif
         + workspace->rates.forward[22]) * a[sHI] * n_nuclei;

  double H2_dot = (workspace->rates.forward[8] * a[sHM] * a[sHI]
#ifdef H2_PLUS
      + workspace->rates.forward[10] * a[sH2II] * a[sHI]
      + workspace->rates.forward[21] * a[sH2II] * a[sHM]
#endif
      ) * n_nuclei * n_nuclei
#ifdef H2_FORMATION_ON_DUST
      + workspace->rates.forward[51] * a[sHI] * n_nuclei
#endif
#ifdef THREEBODY
      + (workspace->rates.forward[55] * a[sHI] * a[sHI] * a[sHI] * n_nuclei
         + workspace->rates.forward[57] * a[sHI] * a[sHI] * a[sH2I] * n_nuclei) *
        n_nuclei * n_nuclei
#endif
      - (workspace->rates.forward[11] * a[sHI] * n_nuclei
         + workspace->rates.forward[12] * a[sHII] * n_nuclei
         + workspace->rates.forward[13] * a[sELEC] * n_nuclei
#ifdef THREEBODY
         + workspace->rates.forward[58] * a[sH2I] * n_nuclei
#endif
         + workspace->rates.forward[27] + workspace->rates.forward[29]) * a[sH2I] * n_nuclei;

  double timescale = dmin(fabs(a[sELEC] * n_nuclei / electron_dot),
                          fabs(a[sHI] * n_nuclei / HI_dot));
  if(n_nuclei > 1.0 || n_nuclei > N_THREEBODY)
    timescale = dmin(timescale, fabs(a[sH2I] * n_nuclei / H2_dot));

  workspace->evaluation.t_chem = timescale;
  evaluate_chemical_energy_timescale(environment, state, workspace);
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
  memcpy(state->abund, output->Abund, sizeof(state->abund));
  state->utherm_cgs = output->Utherm *
                      (All.UnitEnergy_in_cgs / All.UnitMass_in_g);
  state->gamma = output->Gamma;
}

static void export_state(const ChemState *state, chemdata_out *output)
{
  memcpy(output->Abund, state->abund, sizeof(state->abund));
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
}

int compute_abundances_d(int mode, const chemdata_in *input, chemdata_out *output, double dtime)
{
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
  int ifunc;
  int status = 0;
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

  double t_now = 0.;    /*    set current time (redshift) */
    
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
	    ifunc = chemistry_noneq_exp_nometal(environment, &step, state, workspace, dt);
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


/* ----- Compute cooling time of gas and
   ----- solve for thermal energy.  */
int evaluate_radiative_cooling(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *work)
{  
  const double T_gas = work->thermodynamics.temperature;
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas/PROTONMASS;  
#ifdef METALLICITY
  const double metallicity = environment->metallicity;
#endif
#ifdef CONT_OPACITY
  const double tau_cont = work->continuum.tau;
  const double continuum_reduction_factor = work->continuum.reduction_factor;
#endif
  
  const double ln_T0 = log(T_start);
  const double ln_T9 = log(T_end);
  
  /* --------- Compute temperature index.  */
  double ln_T = log(T_gas);
  if (ln_T < ln_T0) ln_T = ln_T0;
  if (ln_T > ln_T9) ln_T = ln_T9;
  
  volatile int i_T = (int)( (ln_T-ln_T0)/dln_T ) ;
  
  if (i_T > N_T-2) i_T = N_T - 2;     /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  const double c_T = (T_gas - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);
  
  /* --------- Interpolate temperature dependent coefficients. */
#ifdef HLINE_COOLING_THIN
  double ceHI  = ceHIa[i_T]   + c_T*(ceHIa[i_T+1] -ceHIa[i_T]);  
#endif

#ifndef CONTINUUM_COOLING_MATSUKOBA  
  double reHII = reHIIa[i_T]  + c_T*(reHIIa[i_T+1]-reHIIa[i_T]);
  double brem  = brema[i_T]   + c_T*(brema[i_T+1] -brema[i_T]);
  double cie   = H2_CIE[i_T]  + c_T*(H2_CIE[i_T+1]-H2_CIE[i_T]);
#endif

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  double ceHeI  = ceHeIa[i_T]  + c_T*(ceHeIa[i_T+1] - ceHeIa[i_T]);
  double ceHeII = ceHeIIa[i_T] + c_T*(ceHeIIa[i_T+1]- ceHeIIa[i_T]);

#ifdef HELIUM
  double ciHeI  = ciHeIa[i_T]  + c_T*(ciHeIa[i_T+1]  -ciHeIa[i_T]);
  double ciHeII = ciHeIIa[i_T] + c_T*(ciHeIIa[i_T+1] -ciHeIIa[i_T]);
  double ciHeIS = ciHeISa[i_T] + c_T*(ciHeISa[i_T+1] -ciHeISa[i_T]);

  double reHeII1= reHeII1a[i_T] + c_T*(reHeII1a[i_T+1]-reHeII1a[i_T]);
  double reHeII2= reHeII2a[i_T] + c_T*(reHeII2a[i_T+1]-reHeII2a[i_T]);
  double reHeIII= reHeIIIa[i_T] + c_T*(reHeIIIa[i_T+1]-reHeIIIa[i_T]);
#endif 
#endif 
  
  /* ========= Compton heating =========================================================*/  
  double comp1 = 5.65e-36*pow(environment->scale_factor,-4);    
  double edot_comp = -comp1*(T_gas-environment->cmb_temperature)*state->abund[sELEC]*n_nuclei;
  
  /**********************************************************************************************/ 
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/ 
  /**********************************************************************************************/ 

#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_gas, T_gas);
#else
  const double l_core = step->jeans_length;
#endif
  
  /******************************* Hydrogen lines ***********************************************/
  double edot_ceHI = 0.;
#ifdef HLINE_COOLING_THIN
  edot_ceHI   = -ceHI  *state->abund[sHI]  *state->abund[sELEC]     *n_nuclei*n_nuclei;
#else
  if(n_nuclei < 1.e8){
    double edot_Lya, edot_2ph;
    double v_thermal = sqrt(BOLTZMANN * T_gas/PROTONMASS);
    double fesc = fesc_Lya(n_nuclei, T_gas, state->abund[sHI], state->abund[sELEC], l_core, v_thermal, &edot_Lya, &edot_2ph);
    
    edot_ceHI = edot_Lya;
  }
#endif

  /******************************* Helium lines **************************************************/
#ifdef HELIUM
  double edot_ceHeI  = -ceHeI *state->abund[sHeII]*state->abund[sELEC]*state->abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
  double edot_ceHeII = -ceHeII*state->abund[sHeII]*state->abund[sELEC]     *n_nuclei*n_nuclei;

  double edot_ciHeI  = -ciHeI *state->abund[sHeI] *state->abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeII = -ciHeII*state->abund[sHeII]*state->abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeIS  = -ciHeIS *state->abund[sHeII]*state->abund[sELEC]*state->abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
#elif defined(HELIUM_EQ_COOLING)
  double denom = work->rates.forward[4]*work->rates.forward[6] + work->rates.forward[3]*work->rates.forward[6] + work->rates.forward[3]*work->rates.forward[5];
  if(denom < 1.e-100)
    denom = 1.e-100;
  
  double fHeI_eq   = XHe * work->rates.forward[4]*work->rates.forward[6] / denom;
  double fHeII_eq  = XHe * work->rates.forward[3]*work->rates.forward[6] / denom;
  double fHeIII_eq = XHe * work->rates.forward[3]*work->rates.forward[5] / denom;
  double felec_new = state->abund[sELEC] + fHeII_eq + 2 * fHeIII_eq;

  double edot_ceHeI  = -ceHeI * fHeI_eq  * felec_new * n_nuclei * n_nuclei;
  double edot_ceHeII = -ceHeII* fHeII_eq * felec_new * n_nuclei * n_nuclei;
#endif      


  /******************************* H2 & HD lines ***********************************************/  
  double nH  = n_nuclei * XH;
  double yHI = state->abund[sHI]    / XH;
  double yHII = state->abund[sHII]  / XH;
  double yH2I = state->abund[sH2I]  / XH;
  double y_e  = state->abund[sELEC] / XH;  
#ifdef HELIUM
  double yHeI = state->abund[sHeI]  / XH;  
#else
  double yHeI = XHe / XH;
#endif
#ifdef HD
  double yHDI = state->abund[sHDI]/XH;
#endif

  double edot_H2 = H2_cooling(nH, T_gas, yHI, yH2I, yHeI, yHII, y_e);
#ifdef HD  
  double edot_HD = HD_cooling(nH, T_gas, yHI, yHDI);
#endif


#ifdef H2_OPACITY  /* accounting for the local optical depth to H2 lines   */ 
  if(n_nuclei > 1.0e4){
    double N_H2 = l_core * n_nuclei * state->abund[sH2I];
    double f_reduce_H2 = H2_line_escape_fraction(T_gas, N_H2);
    work->diagnostics.H2_escape_fraction = f_reduce_H2;
    
    if(n_nuclei >= N_EQCHEM) f_reduce_H2 = 0.;
    if(f_reduce_H2 > 1.) f_reduce_H2 = 1.;
    
    if(f_reduce_H2 < 0){ /* something wrong */
      fflush(stdout);
      endrun();
    }
    
    edot_H2 *= f_reduce_H2;

    if(FlagPrint)
      printf("f_reduce=%g\n", f_reduce_H2);

#ifdef HD
    double N_HD = l_core * n_nuclei * state->abund[sHDI];
    double f_reduce_HD = HD_line_escape_fraction(T_gas, N_HD);
    work->diagnostics.HD_escape_fraction = f_reduce_HD;
    
    if(n_nuclei >= N_EQCHEM) f_reduce_HD = 0.;
    if(f_reduce_HD > 1.) f_reduce_HD = 1.;
    
    if(f_reduce_HD < 0){ /* something wrong */
      fflush(stdout);
      endrun();
    }
    
    edot_HD *= f_reduce_HD;
#endif
  } 
#endif

  /******************************* CII & OI lines ***********************************************/  
#ifdef METALLICITY
  if(FlagPrint)
    printf("f_metal=%g\n", metallicity);
  
  double nOI  = n_nuclei * XH * Y_O_sun * metallicity;  
  double nCI  = n_nuclei * XH * Y_C_sun * metallicity;  
  double nCII = n_nuclei * XH * Y_C_sun * metallicity;  
  double ne   = state->abund[sELEC] * n_nuclei;

  double cool_rate_OI  = OxygenI_cooling(nOI, ne, nH, T_gas);
  double cool_rate_CI  = CarbonI_cooling(nCI, ne, nH, T_gas);
  double cool_rate_CII = CarbonII_cooling(nCII, ne, nH, T_gas);

  double edot_OI  = -1.0 * cool_rate_OI;
  double edot_CI  = -1.0 * cool_rate_CI;
  double edot_CII = -1.0 * cool_rate_CII;
#endif  

  /******************************************************************************************/ 
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/ 
  /******************************************************************************************/ 

#ifdef CONTINUUM_COOLING_MATSUKOBA  
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, state->abund);
#else
    /* --------- Cooling by CIE  */
  /* --------- Bremsstrahlung  */
  /* --------- Recombinations  */
  double edot_reHI    = -reHII*  state->abund[sHII] * state->abund[sELEC] * n_nuclei*n_nuclei;
  double edot_CIE     = -cie  *  state->abund[sH2I] * state->abund[sH2I]  * n_nuclei*n_nuclei;
  double edot_brem    = -brem * state->abund[sHII]  *  state->abund[sELEC]*n_nuclei*n_nuclei;
  
#ifdef HELIUM
  edot_brem    -= brem*(state->abund[sHeII]+4.0*state->abund[sHeIII]) * state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII1 = -reHeII1*state->abund[sHeII] *state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII2 = -reHeII2*state->abund[sHeII] *state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeIII = -reHeIII*state->abund[sHeIII]*state->abund[sELEC]*n_nuclei*n_nuclei;
#endif  
#endif  

  /********************** opacity for the continuum process *******************************/

  double edot_dust = 0.;
#ifdef DUST_COOLING
  edot_dust = work->continuum.dust_cooling_rate;
#endif
#if defined(DUST_COOLING) && defined(LYMANWERNER)
  double edot_PE = 0.;
  if(work->continuum.dust_temperature < 1000.){
    double G0 = environment->G0;
    edot_PE = dust_PEheating(G0, nH, y_e, T_gas, metallicity);
    if(FlagPrint)
      printf("edot_dust=%g cont=%g T_dust=%g\n", edot_dust, edot_cont, work->continuum.dust_temperature);
  }
#endif
  
#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * metallicity * n_nuclei * n_nuclei;
#endif

#ifdef CONT_OPACITY  
  edot_ceHI *= exp(-tau_cont);
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  edot_ceHeI *= exp(-tau_cont);
  edot_ceHeII *= exp(-tau_cont);
#endif

  edot_H2   *= exp(-tau_cont);
#ifdef HD
  edot_HD   *= exp(-tau_cont);
#endif
#ifdef METALLICITY
  edot_CI   *= exp(-tau_cont);
  edot_CII  *= exp(-tau_cont);
  edot_OI   *= exp(-tau_cont);
#endif

#ifdef METAL_LINE_COOLING_SD93
  edot_metal_SD93 *= exp(-tau_cont);
#endif
  
  edot_cont *= continuum_reduction_factor;
#endif  

	       
  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = step->piHI  *state->abund[sHI]  *n_nuclei;
#ifdef HELIUM
  double edot_radHeI  = step->piHeI *state->abund[sHeI] *n_nuclei;
  double edot_radHeII = step->piHeII*state->abund[sHeII]*n_nuclei;
#endif
    
  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  double edot = edot_ceHI + edot_comp + edot_H2     
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    + edot_ceHeI + edot_ceHeII 
#ifdef HELIUM
    + edot_ciHeI + edot_ciHeII + edot_ciHeIS
    + edot_reHeII1 + edot_reHeII2 + edot_reHeIII
#endif
#endif
#ifdef HD
    + edot_HD
#endif
#ifdef METALLICITY
    + edot_CI + edot_CII + edot_OI
#endif
#ifdef METAL_LINE_COOLING_SD93
    + edot_metal_SD93
#endif
#ifdef DUST_COOLING
    + edot_dust
#ifdef LYMANWERNER
    + edot_PE
#endif
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
    + edot_cont
#else
    + edot_reHI + edot_brem + edot_CIE
#endif
    + edot_radHI
#ifdef HELIUM
    + edot_radHeI + edot_radHeII
#endif
    ;


  if(FlagPrint){
#ifdef METALLICITY
    printf("cool_rate:tot=%g HI=%g H2=%g HD=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g T_gas=%g radHI=%g\n",
	   edot, edot_ceHI, edot_H2, edot_HD, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont,
	   work->continuum.reduction_factor, work->continuum.dust_temperature, T_gas, edot_radHI);
#else
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g cont=%g f_reduce=%g T_gas=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_cont, f_reduce, T_gas);
#endif
  }

#ifdef ID_CHEM_DEBUG
  if(FlagPrint)
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont, f_reduce, work->continuum.dust_temperature);
#else
#if defined(THERMAL_FEEDBACK) && defined(METALLICITY)
  if(FlagPrint)
    printf("cool_rate:T_gas=%g n=%g tot=%g rad=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g f_line=%g tau=%g Tdust=%g Trad=%g\n",
	   T_gas, n_nuclei, edot, edot_radHI, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust, edot_cont,
	   f_reduce, work->diagnostics.H2_escape_fraction->EscapeProb, tau_cont,
	   work->continuum.dust_temperature, pow(environment->stellar_rad_flux / STEFAN / 4., 0.25));
#endif
#endif
  
#if defined(DEBUG_CHEMISTRY_SC)
  work->diagnostics.edot_ceHI = edot_ceHI;	       
  work->diagnostics.edot_comp = edot_comp;
  work->diagnostics.edot_H2 = edot_H2;

#ifdef HELIUM
  work->diagnostics.edot_ceHeI = edot_ceHeI;
  work->diagnostics.edot_ceHeII = edot_ceHeII;
  work->diagnostics.edot_ciHeI = edot_ciHeI;
  work->diagnostics.edot_ciHeII = edot_ciHeII;
  work->diagnostics.edot_ciHeIS = edot_ciHeIS;

  work->diagnostics.edot_reHeII1 = edot_reHeII1;
  work->diagnostics.edot_reHeII2 = edot_reHeII2;
  work->diagnostics.edot_reHeIII = edot_reHeIII;
#endif

#ifdef CONTINUUM_COOLING_MATSUKOBA
  work->diagnostics.edot_cont = edot_cont;
#else
  work->diagnostics.edot_reHI = edot_reHI;
  work->diagnostics.edot_brem = edot_brem;
  work->diagnostics.edot_CIE = edot_CIE;
#endif
  
#ifdef HD
  work->diagnostics.edot_HD = edot_HD;
#endif
#ifdef METALLICITY
  work->diagnostics.edot_CI  = edot_CI;
  work->diagnostics.edot_CII = edot_CII;
  work->diagnostics.edot_OI  = edot_OI;
#endif
#ifdef DUST_COOLING
  work->diagnostics.edot_dust = edot_dust;
#ifdef LYMANWERNER
  work->diagnostics.edot_PE = edot_PE;
#endif		
#endif		
#endif		
    
  work->evaluation.rhoU_dot_rad = edot;
  
  double rhoU = state->utherm_cgs * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.
  
  work->evaluation.t_cool_rad = (rhoU/edot); /* in seconds */  
  if(fabs(edot) < 1.0e-40)
    work->evaluation.t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}


int radiative_cooling_high_T(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *work)
{    
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei  = rho_gas/PROTONMASS;
  const double T_gas = work->thermodynamics.temperature;
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_gas, T_gas);
#else
  const double l_core = step->jeans_length;
#endif
  
  /******************************************************************************************/ 
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/ 
  /******************************************************************************************/ 

#ifdef CONTINUUM_COOLING_MATSUKOBA  
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, state->abund);
#endif  
  double edot_reHI = -8.418e-26*pow(T_gas,0.11) * n_nuclei*n_nuclei;
  double edot_brem = -1.43e-27*sqrt(T_gas) * (1. + 4.*XHe/XH) * n_nuclei * n_nuclei;

  if(FlagPrint)
    printf("T_gas=%g edot_reHI=%g brem=%g cont=%g\n", T_gas, edot_reHI, edot_brem, edot_cont);

  /********************** opacity for the continuum process *******************************/

#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * environment->metallicity * n_nuclei * n_nuclei;
#endif

  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = step->piHI  *state->abund[sHI] *n_nuclei;
    
  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  double edot = edot_reHI + edot_brem 
#ifdef METAL_LINE_COOLING_SD93
    + edot_metal_SD93
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA  
    + edot_cont
#endif
    + edot_radHI;
  
  work->evaluation.rhoU_dot_rad = edot;
  
  double rhoU = state->utherm_cgs * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.
  
  work->evaluation.t_cool_rad = (rhoU/edot); /* in seconds */

  if(FlagPrint){
#ifdef METAL_LINE_COOLING_SD93
    printf("cool_rate:tot=%g reHI=%g brem=%g cont_metal=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, edot_metal_SD93, T_gas, rhoU, work->evaluation.t_cool_rad);
#else
    printf("cool_rate:tot=%g reHI=%g brem=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, T_gas, rhoU, work->evaluation.t_cool_rad);
#endif
  }
  
  if(fabs(edot) < 1.0e-40)
    work->evaluation.t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}

int chemical_cooling(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *work, double dt, double *abund_old){
  const double T_gas = work->thermodynamics.temperature;
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas / PROTONMASS;
  
  double edot_H2form,edot_H2form_3body;
  double edot_ciHI, edot_ciHeI, edot_ciHeII;
  double edot_HMform;

  enum SPECIES species;
  double d_abundance[N_SP];
  
  /* calculate chemical heating/cooling */
  double n_elec = n_nuclei * state->abund[sELEC];
  double n_HI   = n_nuclei * state->abund[sHI];
  double n_HM   = n_nuclei * state->abund[sHM];
  double n_H2I  = n_nuclei * state->abund[sH2I];
#ifdef H2_PLUS
  double n_H2II = n_nuclei * state->abund[sH2II];
#endif
#ifdef HELIUM
  double n_HeI  = n_nuclei * state->abund[sHeI];
  double n_HeII = n_nuclei * state->abund[sHeII];
#endif
#ifdef HD
  double n_DI = n_nuclei * state->abund[sDI];
  double n_DII = n_nuclei * state->abund[sDII];
  double n_HDI = n_nuclei * state->abund[sHDI];
#endif

  for(int j=0;j<N_SP;j++)
    d_abundance[j] = state->abund[j] - abund_old[j];
  
  double n_crit = (1.e6 / sqrt(T_gas)) 
    / (1.6 * state->abund[sHI] * exp(-(400.0/T_gas)*(400.0/T_gas)) + 1.4 * state->abund[sH2I] * exp(-12000.0/(T_gas + 1200.0)));

  /********** Mode 0 is just for evaluating the timescale chemical heating / cooling ********************/
  /********** This is only used in the explicit time-stepping or hybrid time-stepping *******************/
  double edot_H2form_HM, edot_H2form_dust, edot_ciHI_EUV;
  double edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV, edot_H2form_3body_loss;

  if(dt == 0.){
    edot_H2form_HM = eV_to_erg * ( work->rates.forward[8] * n_HM * n_HI * 3.73/ (1. + n_crit/n_HI)
#ifdef H2_PLUS
				   + work->rates.forward[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/n_HI)) 
#endif
				   ); //1.38?? modified to 1.83

    
    edot_H2form_dust = 0.;
#ifdef H2_FORMATION_ON_DUST
    edot_H2form_dust = eV_to_erg * work->rates.forward[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/n_HI));
#endif
    
    edot_H2form_3body = 0.;
#ifdef THREEBODY
    edot_H2form_3body = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI
					     - ( work->rates.forward[11]*n_HI + work->rates.forward[58]*n_H2I) * n_H2I );
#endif    
    edot_H2form = edot_H2form_HM + edot_H2form_3body + edot_H2form_dust;

    edot_ciHI   = -13.598 * eV_to_erg * work->rates.forward[1] * n_HI * n_elec; 
#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * work->rates.forward[3] * n_HeI  * n_elec;
    edot_ciHeII = -54.416 * eV_to_erg * work->rates.forward[5] * n_HeII * n_elec;
#endif

    edot_HMform = 0.;
  }else{
    if(XH * n_nuclei < N_THREEBODY){
      edot_H2form_HM_loss = eV_to_erg * ( work->rates.forward[8]  * n_HM   * n_HI * 3.73 / (1. + n_nuclei/n_crit)
#ifdef H2_PLUS
					  + work->rates.forward[10] * n_H2II * n_HI * 1.83 / (1. + n_nuclei/n_crit)
#endif
					  );
 
      edot_H2form_dust_loss = 0.;
#ifdef H2_FORMATION_ON_DUST
      edot_H2form_dust_loss = eV_to_erg * work->rates.forward[51] * n_HI * 4.28 / (1. + n_nuclei/n_crit);
#endif
      edot_H2diss_FUV = eV_to_erg * work->rates.forward[29] * n_H2I * 4.48;
        
      edot_H2form_3body = 4.48 * eV_to_erg * d_abundance[sH2I] * n_nuclei / dt;
      edot_H2form_3body_loss = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI ) / (1. + n_nuclei/n_crit);


      edot_H2form = edot_H2form_3body - edot_H2form_HM_loss - edot_H2form_dust_loss + edot_H2diss_FUV;    

      edot_H2form_HM = eV_to_erg * ( work->rates.forward[8] * n_HM * n_HI * 3.73/ (1. + n_crit/n_nuclei)
#ifdef H2_PLUS
				   + work->rates.forward[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/n_nuclei)) 
#endif
				   ); //1.38?? modified to 1.83
      edot_H2form_dust = eV_to_erg * work->rates.forward[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/n_nuclei));
      edot_H2form_3body = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI
					       - ( work->rates.forward[11]*n_HI + work->rates.forward[58]*n_H2I) * n_H2I );      
      edot_HMform = 0.78 * eV_to_erg * d_abundance[sHM]  * n_nuclei / dt;
    }else{
      edot_H2form = edot_H2form_3body = 4.48 * eV_to_erg * d_abundance[sH2I] * n_nuclei / dt;
      edot_HMform = 0.78 * eV_to_erg * d_abundance[sHM]  * n_nuclei / dt;
    }
  
    edot_ciHI   = -13.598 * eV_to_erg * d_abundance[sHII]    * n_nuclei / dt;
#ifdef IONIZATION_PHOTON
    edot_ciHI_EUV = 13.598 * eV_to_erg * work->rates.forward[22] * n_HI;
    edot_ciHI  += edot_ciHI_EUV;
#endif

#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * d_abundance[sHeII]   * n_nuclei / dt;
    edot_ciHeII = -54.416 * eV_to_erg * d_abundance[sHeIII]  * n_nuclei / dt;
#endif
  }

  if(FlagPrint){
    printf("edot_chem: n_nuclei=%g H2=%g 3b=%g 3b_l=%g HM_l=%g dust_l=%g FUV=%g ciHI=%g EUV=%g HM=%g ncrit=%g dy_H2=%g\n"
	   , n_nuclei, edot_H2form, edot_H2form_3body, edot_H2form_3body_loss, edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV
	   , edot_ciHI, edot_ciHI_EUV, edot_HMform, n_crit, d_abundance[sH2I]);
  }
  
  work->evaluation.rhoU_dot_chem = edot_H2form + edot_ciHI + edot_HMform;  
  if(dt < 1.e-30){
    work->evaluation.rhoU_dot_chem = 0.;
    edot_H2form = edot_ciHI = edot_HMform = 0.;
  }

  double rhoU = state->utherm_cgs * environment->rho_cgs; 
  //modified[0814] in case for SphP[ithis].Entropy is not updated.

  work->evaluation.t_cool_chem = (rhoU/work->evaluation.rhoU_dot_chem); /* in seconds */

  return 0;
}


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


int chemistry_noneq_exp_nometal(const ChemEnvironment *environment, const ChemStepData *step, ChemState *state, ChemWorkspace *work, double dt)
{
  enum SPECIES species;

  double Ccoef,Dcoef,correction;
  double HIp,HIIp,elecp;
  double H2Ip,HMp;
	
  double t_min, elec_dot, HI_dot, H2I_dot;	

#ifdef HD
  double DIp, DIIp, HDIp, HDIIp, DMp;
#endif

  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas/PROTONMASS;
  const double T_gas = work->thermodynamics.temperature;
  
  evaluate_reaction_rates(environment, step, state, work);
  
  double abund_old[N_SP];
  for(int j=0;j<N_SP;j++)
    abund_old[j]  = state->abund[j];

  /* ----------- HI */  
  //added reaction 16,17,18,20,21[chon140112]

  Ccoef = (work->rates.forward[2]  * state->abund[sHII]  * state->abund[sELEC] 
	   + work->rates.forward[11] * state->abund[sH2I]  * state->abund[sHI]   *2.
	   + work->rates.forward[12] * state->abund[sH2I]  * state->abund[sHII]  
	   + work->rates.forward[13] * state->abund[sH2I]  * state->abund[sELEC] *2.
#ifdef HELIUM
	   + work->rates.forward[15] * state->abund[sHeI]  * state->abund[sHII]
#endif
	   + work->rates.forward[16] * state->abund[sHM] * state->abund[sELEC]
	   + work->rates.forward[17] * state->abund[sHM] * state->abund[sHI]
	   + 2 * work->rates.forward[18] * state->abund[sHM] * state->abund[sHII]
#ifdef H2_PLUS
	   + work->rates.forward[20] * state->abund[sH2II] * state->abund[sELEC] *2.
	   + work->rates.forward[21] * state->abund[sH2II] * state->abund[sHM] 
#endif
#ifdef HD
	   + work->rates.forward[32] * state->abund[sDI] * state->abund[sHII] 
#endif
	   )*n_nuclei*n_nuclei
#ifdef THREEBODY
    + work->rates.forward[58] * state->abund[sH2I]  * state->abund[sH2I]  * n_nuclei*n_nuclei *2.
#endif
    + work->rates.forward[25] * state->abund[sHM]   * n_nuclei
#ifdef H2_PLUS
    + work->rates.forward[26] * state->abund[sH2II] * n_nuclei
#endif
    + work->rates.forward[29] * state->abund[sH2I]  * n_nuclei * 2.
    ;

    
  Dcoef = (  work->rates.forward[1]  * state->abund[sELEC]
#ifdef H2_PLUS
	     + work->rates.forward[10] * state->abund[sH2II]
#endif
#ifdef HELIUM
	     + work->rates.forward[14] * state->abund[sHeII] 
#endif
#ifdef HD
	     + work->rates.forward[33] * state->abund[sDII]
	     + work->rates.forward[38] * state->abund[sHDI]
#endif
	     ) * n_nuclei
#ifdef THREEBODY
    + (work->rates.forward[55] * state->abund[sHI]  * 2.
       +  work->rates.forward[57] * state->abund[sH2I] * 2. )* state->abund[sHI] * n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
    + 2. * work->rates.forward[51]
#endif
    + work->rates.forward[22];


  HI_dot = fabs(Ccoef - Dcoef);  
  HIp   = (Ccoef*dt + state->abund[sHI]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


  /* ----------- HII */
  Ccoef = work->rates.forward[1]  * HIp * state->abund[sELEC] * n_nuclei*n_nuclei
#ifdef H2_PLUS
    + work->rates.forward[10] * HIp * state->abund[sH2II] * n_nuclei*n_nuclei 
#endif
#ifdef HELIUM
    + work->rates.forward[14] * HIp * state->abund[sHeII] * n_nuclei*n_nuclei
#endif
    + work->rates.forward[22] * HIp * n_nuclei
#ifdef H2_PLUS
    + work->rates.forward[26] * state->abund[sH2II] * n_nuclei
    + work->rates.forward[28] * state->abund[sH2II] * n_nuclei * 2
#endif
    /*added by chon[131221] */
#ifdef HD
    + work->rates.forward[33] * HIp * state->abund[sDII] * n_nuclei * n_nuclei
    + work->rates.forward[37] * state->abund[sH2I] * state->abund[sDII] * n_nuclei * n_nuclei
#ifdef HD_PLUS
    + work->rates.forward[36] * HIp * state->abund[sHDII] * n_nuclei * n_nuclei
    + work->rates.forward[52] * state->abund[sHDII] * n_nuclei
#endif
#endif
    /*[131221]*/
    ;

  Dcoef = (work->rates.forward[2]  * state->abund[sELEC]
#ifdef H2_PLUS
	   + work->rates.forward[9]  * HIp             
	   + work->rates.forward[12] * state->abund[sH2I]
#endif
#ifdef HELIUM 
	   + work->rates.forward[15] * state->abund[sHeI]
#endif
	   )* n_nuclei
    /*added by chon[131221] */
    + work->rates.forward[18] * state->abund[sHM] * n_nuclei
#ifdef H2_PLUS
    + work->rates.forward[19] * state->abund[sHM] * n_nuclei
#endif
#ifdef HD
    + (work->rates.forward[32] * state->abund[sDI] 
       + work->rates.forward[39] * state->abund[sHDI]
#ifdef HD_PLUS
       + work->rates.forward[40] * state->abund[sDI]
#endif
#ifdef DMINUS
       + work->rates.forward[45] * state->abund[sDM] 
#endif
       ) *n_nuclei /* modified state->abund[sHDI] -> state->abund[sDM] [0926]*/ 
#endif
    /*[131221]*/
    ;

    
  HIIp  = (Ccoef*dt + state->abund[sHII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
    
#ifdef HELIUM
  /* ----------- HeI */
  Ccoef = (work->rates.forward[4]  * state->abund[sHeII] * state->abund[sELEC]
	   + work->rates.forward[14] * HIp              * state->abund[sHeII])*n_nuclei*n_nuclei;
  
  Dcoef = work->rates.forward[3]  * state->abund[sELEC] * n_nuclei 
    + work->rates.forward[15] * HIIp             * n_nuclei
    + work->rates.forward[24];
  
  double HeIp  = (Ccoef*dt + state->abund[sHeI]*n_nuclei) / ( 1. + Dcoef*dt ) / n_nuclei;


  /* ----------- HeII */
  Ccoef = work->rates.forward[3]  * HeIp              * state->abund[sELEC] * n_nuclei*n_nuclei
    + work->rates.forward[6]  * state->abund[sHeIII] * state->abund[sELEC] * n_nuclei*n_nuclei
    + work->rates.forward[15] * HeIp              * HIIp             * n_nuclei*n_nuclei
    + work->rates.forward[24] * HeIp*n_nuclei;

  Dcoef = work->rates.forward[4]  * state->abund[sELEC] * n_nuclei
    + work->rates.forward[5]  * state->abund[sELEC] * n_nuclei
    + work->rates.forward[14] * HIp              * n_nuclei
    + work->rates.forward[23]; 

  double HeIIp  = (Ccoef*dt + state->abund[sHeII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

      
  /* ----------- HeIII */
  Ccoef = work->rates.forward[5]  * HeIIp * state->abund[sELEC] * n_nuclei*n_nuclei
    + work->rates.forward[23] * HeIIp * n_nuclei;
      
  Dcoef = work->rates.forward[6] * state->abund[sELEC] * n_nuclei;
    
  double HeIIIp = (Ccoef*dt + state->abund[sHeIII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


  if(HeIp   < 0) HeIp   = 0.;
  if(HeIIp  < 0) HeIIp  = 0.;
  if(HeIIIp < 0) HeIIIp = 0.;

  state->abund[sHeI]   = HeIp;
  state->abund[sHeII]  = HeIIp;
  state->abund[sHeIII] = HeIIIp;
#endif
  
  /* --------- elec */
  Ccoef = (work->rates.forward[22] *HIp  
#ifdef HELIUM	
	   + work->rates.forward[23] *HeIIp
	   + work->rates.forward[24] *HeIp 
#endif
	   )*n_nuclei;
  
  Dcoef = - (work->rates.forward[1]*HIp   - work->rates.forward[2]*HIIp
#ifdef HELIUM	
	     + work->rates.forward[3]*HeIp  - work->rates.forward[6]*HeIIIp
	     + work->rates.forward[5]*HeIIp - work->rates.forward[4]*HeIIp
#endif
	     )*n_nuclei;

  elec_dot = fabs(Ccoef - Dcoef);  
  elecp = (Ccoef*dt + state->abund[sELEC]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


  if(HIp   < 0) HIp   = 0;
  if(HIIp  < 0) HIIp  = 0;
  if(elecp < 0) elecp = 0;
    

  state->abund[sHI]   = HIp;
  state->abund[sHII]  = HIIp;    
  state->abund[sELEC] = elecp;


  /* --------- H2 */
  /* [SH131015 : crazyyyyyyy] */
  
  /* HM fast reactions */
  /*
    denom = (work->rates.forward[8] + work->rates.forward[17]) * state->abund[sHI] 
    + (work->rates.forward[18]+ work->rates.forward[19]) * state->abund[sHII] 
    +  work->rates.forward[16]         * state->abund[sELEC] 
    +  work->rates.forward[25]/n_nuclei;
  
    state->abund[sHM] = (work->rates.forward[7]*state->abund[sHI]*state->abund[sELEC]) / denom;
  */
  Ccoef = work->rates.forward[7]*state->abund[sHI]*state->abund[sELEC] * n_nuclei*n_nuclei
#ifdef HD
#ifdef DMINUS
    + work->rates.forward[47] * HIp * state->abund[sDM] *n_nuclei*n_nuclei
#endif
#endif
    ;

  Dcoef = ((work->rates.forward[8] + work->rates.forward[17]) * state->abund[sHI] 
	   + work->rates.forward[18] * state->abund[sHII] 
#ifdef H2_PLUS
	   + work->rates.forward[19] * state->abund[sHII] 
#endif
	   + work->rates.forward[16] * state->abund[sELEC]  
	   ) * n_nuclei
    +  work->rates.forward[25]
#ifdef HD	       
#ifdef DMINUS
    +  work->rates.forward[46] * state->abund[sDI] *n_nuclei
#endif
#endif
    ;
  
  HMp = (Ccoef*dt + state->abund[sHM]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


  /* H2II fast reactions */
  /*  denom = work->rates.forward[10]*state->abund[sHI] 
      + work->rates.forward[20]*state->abund[sELEC] 
      + work->rates.forward[21]*state->abund[sHM] 
      +(work->rates.forward[26]+work->rates.forward[28])/n_nuclei;
    
      state->abund[sH2II] = (  work->rates.forward[9] * state->abund[sHI]  * state->abund[sHII]
      + work->rates.forward[12] * state->abund[sH2I] * state->abund[sHII]
      + work->rates.forward[19] * state->abund[sHM]  * state->abund[sHII]
      + work->rates.forward[27] * state->abund[sH2I]/n_nuclei         )/ denom; */

  /* noneq H2II added by chon.[131228]*/

#ifdef H2_PLUS
  Ccoef = ( work->rates.forward[9] * state->abund[sHI]  * state->abund[sHII]
	    +work->rates.forward[12] * state->abund[sH2I] * state->abund[sHII] 
	    +work->rates.forward[19] * HMp  * state->abund[sHII] )*n_nuclei *n_nuclei
    + work->rates.forward[27]* state->abund[sH2I]*n_nuclei;

  Dcoef = (   work->rates.forward[10]*state->abund[sHI] 
	      + work->rates.forward[20]*state->abund[sELEC]
	      + work->rates.forward[21]*HMp ) *n_nuclei
    +work->rates.forward[26] + work->rates.forward[28];

  double H2IIp = (Ccoef*dt + state->abund[sH2II]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif


  Ccoef = ( work->rates.forward[8]  * HMp   * state->abund[sHI]
#ifdef H2_PLUS
	    + work->rates.forward[10] * H2IIp * state->abund[sHI]
	    + work->rates.forward[21] * H2IIp * HMp
#endif
#ifdef HD
	    + work->rates.forward[38] * HIp * state->abund[sHDI] 
	    + work->rates.forward[39] * HIIp * state->abund[sHDI]
#endif
	    ) *n_nuclei*n_nuclei
#ifdef THREEBODY
    +( work->rates.forward[55] * state->abund[sHI]   * state->abund[sHI] * state->abund[sHI]
       + work->rates.forward[57] * state->abund[sHI]   * state->abund[sHI] * state->abund[sH2I] )*n_nuclei*n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
    + work->rates.forward[51] * state->abund[sHI]*n_nuclei
#endif
    ;
    
  Dcoef = (work->rates.forward[11] * state->abund[sHI] 
#ifdef H2_PLUS
	   + work->rates.forward[12] * HIIp 
#endif
	   + work->rates.forward[13] * state->abund[sELEC]
#ifdef HD
	   + work->rates.forward[35] * state->abund[sDI]  
	   + work->rates.forward[37] * state->abund[sDII]  
#endif
	   )*n_nuclei
#ifdef THREEBODY
    + work->rates.forward[58] * state->abund[sH2I] * n_nuclei
#endif
    + work->rates.forward[27] 
    + work->rates.forward[29];

  H2I_dot = fabs(Ccoef - Dcoef);  
  H2Ip = (Ccoef*dt + state->abund[sH2I]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

  if(HMp < 0) HMp = 0;
#ifdef H2_PLUS
  if(H2IIp < 0) H2IIp = 0;
#endif
  if(H2Ip < 0) H2Ip = 0;

  state->abund[sHM] = HMp;
#ifdef H2_PLUS
  state->abund[sH2II] = H2IIp;
#endif
  state->abund[sH2I] = H2Ip;
  
  
#ifdef HD

  DIp=DIIp=HDIp=HDIIp=DMp=0;

  /* ----------- DI */
  Ccoef  = //work->rates.forward[31] * state->abund[sDII] * state->abund[sELEC] * n_nuclei*n_nuclei
    work->rates.forward[33] * state->abund[sDII] * HIp   * n_nuclei*n_nuclei
    + work->rates.forward[38] * state->abund[sHDI] * HIp   * n_nuclei*n_nuclei
#ifdef HD_PLUS
    + work->rates.forward[42] * state->abund[sHDII]* state->abund[sELEC] * n_nuclei*n_nuclei
    + work->rates.forward[52] * state->abund[sHDII]* n_nuclei
#endif
#ifdef DMINUS
    + work->rates.forward[44] * state->abund[sDII] * state->abund[sDM]   * n_nuclei*n_nuclei * 2.0 
    + work->rates.forward[45] * HIIp * state->abund[sDM]   * n_nuclei*n_nuclei
    + work->rates.forward[47] * HIp  * state->abund[sDM]   * n_nuclei*n_nuclei
#endif
    + work->rates.forward[53] * state->abund[sHDI] * n_nuclei;
	     
  Dcoef = ( work->rates.forward[32] * HIIp
	    //+ work->rates.forward[34] * state->abund[sHI] 
	    + work->rates.forward[35] * H2Ip
#ifdef HD_PLUS
	    + work->rates.forward[40] * HIIp
#endif
#ifdef DMINUS
	    + work->rates.forward[43] * state->abund[sELEC]
#endif
	    )* n_nuclei
    + work->rates.forward[50];

  DIp   = (Ccoef*dt + state->abund[sDI]*n_nuclei)  / (1. + Dcoef*dt) / n_nuclei;
  
  
  /* ----------- DII */
  Ccoef = work->rates.forward[32] * DIp                 * HIIp * n_nuclei*n_nuclei
    + work->rates.forward[39] * state->abund[sHDI] * HIIp * n_nuclei*n_nuclei
    + work->rates.forward[50] * DIp                 * n_nuclei
#ifdef HD_PLUS
    + work->rates.forward[51] * state->abund[sHDII]* n_nuclei
#endif
    ;

  Dcoef = ( //work->rates.forward[31] * state->abund[sELEC]
	    work->rates.forward[33] * HIp
	    + work->rates.forward[37] * H2Ip
#ifdef DMINUS
	    + work->rates.forward[41] * HIp
	    + work->rates.forward[44] * state->abund[sDM]  
#endif
	    ) * n_nuclei;
	    
  DIIp  = (Ccoef*dt + state->abund[sDII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


  /* ----------- HDI */
  Ccoef = ( //work->rates.forward[34] * DIp            * state->abund[sHI]  
	    work->rates.forward[35] * DIp            * H2Ip
	    + work->rates.forward[37] * DIIp           * H2Ip
#ifdef HD_PLUS
	    + work->rates.forward[36] * HIp * state->abund[sHDII]
#endif	      
#ifdef DMINUS
	    + work->rates.forward[48] * HIp * state->abund[sDM]  
#endif
	    )* n_nuclei*n_nuclei;

  Dcoef = work->rates.forward[38] * HIp * n_nuclei
    + work->rates.forward[39] * HIIp * n_nuclei
    + work->rates.forward[53];

  HDIp  = (Ccoef * dt + state->abund[sHDI] * n_nuclei) / (1. +  Dcoef * dt) / n_nuclei;
    
#ifdef HD_PLUS
  /* ----------- HDII */
  Ccoef = ( work->rates.forward[40] * DIp  * state->abund[sHII]
	    + work->rates.forward[41] * DIIp * state->abund[sHI] ) * n_nuclei*n_nuclei;
      
  Dcoef = work->rates.forward[36] * state->abund[sHI]   * n_nuclei
    + work->rates.forward[42] * state->abund[sELEC] * n_nuclei
    + work->rates.forward[51]
    + work->rates.forward[52];

  HDIIp = (Ccoef*dt + state->abund[sHDII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
    
#ifdef DMINUS    
  /* ----------- DM */
  Ccoef = ( work->rates.forward[43] * DIp  * state->abund[sELEC]
	    + work->rates.forward[46] * DIp  * state->abund[sHM]  ) * n_nuclei*n_nuclei;
      
  Dcoef = ( work->rates.forward[44] * DIIp           
	    + work->rates.forward[45] * state->abund[sHII]
	    + work->rates.forward[47] * state->abund[sHI] 
	    + work->rates.forward[48] * state->abund[sHI]  ) * n_nuclei;
    
  DMp   = (Ccoef*dt + state->abund[sDM]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
      
  if(DIp    < 0) DIp   =0;
  if(DIIp   < 0) DIIp  =0;
  if(HDIp   < 0) HDIp  =0;
  state->abund[sDI]    = DIp;
  state->abund[sDII]   = DIIp;
  state->abund[sHDI]   = HDIp;   

#ifdef HD_PLUS
  if(HDIIp  < 0) HDIIp =0;
  state->abund[sHDII]  = HDIIp;
#endif
#ifdef DMINUS
  if(DMp    < 0) DMp   =0;
  state->abund[sDM]    = DMp;    
#endif  
#endif       
    
  /* --------- Undertake mass conservation corrections. */
  
  //state->abund[sELEC] = 0.0;
  
  correction = H_number_fraction / (state->abund[sHI] 
				    +   state->abund[sHII] 
				    +   state->abund[sHM] 
				    + 2*state->abund[sH2I]
#ifdef H2_PLUS
				    + 2*state->abund[sH2II]
#endif
#ifdef HD
				    +   state->abund[sHDI] 
#ifdef HD_PLUS
				    +   state->abund[sHDII]
#endif
#endif
				    );
  
  if(0)
    if(correction > 2.0 || correction < 0.5){
      printf("H! ThisTask %d correction=%g\n"
	     , ThisTask, correction);
      fflush(stdout);
      endrun();
    }

  
  state->abund[sHI]  *= correction;
  state->abund[sHII] *= correction;
  //state->abund[sELEC]+= state->abund[sHII];
    
  state->abund[sHM]   *= correction;
  state->abund[sH2I]  *= correction;
#ifdef H2_PLUS
  state->abund[sH2II] *= correction;
#endif
  //state->abund[sELEC]  = state->abund[sELEC] - state->abund[sHM] + state->abund[sH2II];
  

#ifdef HD
  correction = D_number_fraction / (state->abund[sDI] 
				    + state->abund[sDII]
				    + state->abund[sHDI]
#ifdef HD_PLUS
				    + state->abund[sHDII]
#endif
#ifdef DMINUS
				    + state->abund[sDM]
#endif
				    );

  if(0)
    if(correction > 2.0 || correction < 0.5){
      printf("D! ThisTask %d correction=%g\n"
	     , ThisTask, correction);
      fflush(stdout);
      endrun();
    } 

  state->abund[sDI]   *= correction;
  state->abund[sDII]  *= correction;
  state->abund[sHDI]  *= correction;
#ifdef HD_PLUS
  state->abund[sHDII] *= correction;
#endif
#ifdef DMINUS
  state->abund[sDM]   *= correction;
#endif
  //state->abund[sELEC]  = state->abund[sELEC] + state->abund[sDII] +  state->abund[sHDII] - state->abund[sDM];
#endif 
  
  
#ifdef HELIUM
  correction = He_number_fraction / (state->abund[sHeI] 
				     + state->abund[sHeII] 
				     + state->abund[sHeIII]);

  if(0)
    if(correction > 1.1 || correction < 0.9){
      printf("He! ThisTask %d correction=%g %g %g %g\n"
	     , ThisTask, correction, state->abund[sHeI], state->abund[sHeII], state->abund[sHeIII]);
      fflush(stdout);
      endrun();
    }
  
  state->abund[sHeI]   *= correction;
  state->abund[sHeII]  *= correction;
  state->abund[sHeIII] *= correction;
  //state->abund[sELEC]  = state->abund[sELEC] + state->abund[sHeII] + 2.*state->abund[sHeIII];
#endif

  /*[SH130423: elec calculation]*/

  state->abund[sELEC] = 0.0;

  state->abund[sELEC] += state->abund[sHII];

  state->abund[sELEC]  = state->abund[sELEC] - state->abund[sHM];

#ifdef H2_PLUS
  state->abund[sELEC] += state->abund[sH2II];
#endif

#ifdef HD
  state->abund[sELEC]  += state->abund[sDII]; 
#ifdef HD_PLUS
  state->abund[sELEC]  += state->abund[sHDII];
#endif
#ifdef DMINUS
  state->abund[sELEC]  -= state->abund[sDM];
#endif
#endif

#ifdef HELIUM
  state->abund[sELEC]  = state->abund[sELEC] + state->abund[sHeII] + 2.*state->abund[sHeIII];
#endif

  if(isnan(state->abund[sHI])){
    printf("[task%03d]HI is nan (ID=%d, dt=%g, rho_gas=%g, T_gas=%g)\n", ThisTask, environment->ID, dt, rho_gas, T_gas);
    return 2;
  }
  if(isnan(T_gas)){
    printf("[task%03d]line=%d T_gas is nan (dt=%g, rho_gas=%g, T_gas=%g)\n", ThisTask, __LINE__, dt, rho_gas, T_gas);
    return 2;
  }

  /* added[0813] */
  t_min = dmin(fabs(state->abund[sELEC] * n_nuclei/elec_dot), fabs(state->abund[sHI] * n_nuclei/HI_dot));
  t_min = dmin(t_min, fabs(state->abund[sH2I] * n_nuclei/H2I_dot));

  work->evaluation.t_chem = t_min;

  chemical_cooling(environment, state, work, dt, abund_old);  

  return 1; 
}


#define MAX_ITER_CHEM 1000
#define ITER_RESET 40
#define BETA 1.0

#ifdef SAHA
static int chemistry_saha(const ChemEnvironment *environment, ChemState *state, ChemWorkspace *workspace, double dt)
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
  int n_sp;
  int inc, info;
	
  enum SPECIES species;
  const ChemState state_old = *state;
  
  n_sp = N_SP;
  inc  = 1   ;

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
		
  int i_sp,j_sp,itr,ifunc;
  int err_flag = 0, ipiv[n_sp];
  double drfdy[n_sp][n_sp],dy[n_sp];
  double y_init[n_sp], abund_init[n_sp];
  double y_mean[n_sp],y_max[n_sp],y_min[n_sp];
  double yi_ion;
  double rf_tot[n_sp];
  double AA[n_sp*n_sp],ddy[n_sp];
  double err,err_max,t_chem1,t_chem2;

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

  for(i_sp=0;i_sp<n_sp;i_sp++) {
    y_mean[i_sp] /= XH;    
    y_init[i_sp] = yi[i_sp];
    
    if(yi[i_sp] < y_min[i_sp]/* || yi[i_sp] > y_max[i_sp]*/){
      y_init[i_sp] = y_mean[i_sp];
      yi[i_sp]     = y_mean[i_sp];
    }
  }

  correct_initial_abundances(yi);
  for(i_sp=0;i_sp<n_sp;i_sp++){
    y_init[i_sp] = yi[i_sp];
    abund_init[i_sp] = y_init[i_sp] * XH;
  }

  //give_initial_guess(T_gas, nH, dt, yi, y_init, k);
  for(i_sp=0;i_sp<n_sp;i_sp++){
    dy[i_sp] = yi[i_sp] - y_init[i_sp];
  }

  int flag_conv;
  err_max= 100.0, itr=1;
  flag_conv = 0;

  int flag_reset=0;

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
    
    for(i_sp=0;i_sp<n_sp;i_sp++){
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
    }

  }

  double ErrMaxCrit=1.e-6;

  if(FlagPrint)
    printf("T_gas=%g y=%g %g %g %g %g\n", T_gas, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM]);
  
  while(itr<=MAX_ITER_CHEM && flag_conv < 2) {  // iteration
#ifdef HMINUS_EQ
    yi[sHM] =
      (work->rates.forward[7]*yi[sHI]*yi[sELEC]+k_inverse[8]*yi[sH2I]*yi[sELEC]+k_inverse[18]*yi[sHI]*yi[sHI])
      /(k_inverse[7]+work->rates.forward[8]*yi[sHI]+work->rates.forward[18]*yi[sHII]+work->rates.forward[25]);
#endif
    
    ifunc=ReactionRates_nometal(itr,nH,yi,work->rates.forward,work->rates.inverse,rf_tot);
    
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
      printf("[%d] T_gas=%g y_H2=%g dy=%g resid=%g\n", itr, T_gas, yi[sH2I], yi[sH2I]-y_init[sH2I], rf_tot[i_sp]*dt);
    }
#endif


    for(j_sp=0;j_sp<n_sp;j_sp++) {
#ifdef CHEM_DIFF_ANALYTIC
      ifunc=ReactionRates_diff_nometal_fp(itr, j_sp, nH, yi, work->rates.forward, work->rates.inverse, drfdy[j_sp]);
#else
      double Dyi, epsilon_y = 1.0e-10, epsilon = 1.0e-4;      
      if(fabs(yi[j_sp])<1.0e-100) Dyi = epsilon_y;
      else                       Dyi = epsilon*yi[j_sp];

      double y_tmp[n_sp];
      for(i_sp=0;i_sp<n_sp;i_sp++){
	if(i_sp==j_sp) y_tmp[i_sp]=yi[i_sp]+0.5*Dyi;
	else           y_tmp[i_sp]=yi[i_sp];
      }

      double rf_f[n_sp],rf_b[n_sp];
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,work->rates.forward,work->rates.inverse,rf_f);
      
      y_tmp[j_sp]=yi[j_sp]-0.5*Dyi;
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,work->rates.forward,work->rates.inverse,rf_b);
      
      for(i_sp=0;i_sp<n_sp;i_sp++) {
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
#endif
    }

#ifdef ID_CHEM_DEBUG
    if(FlagPrint){
      for(j_sp=0;j_sp<n_sp;j_sp++) 
	printf("[%d] drfdy[%d]=%g %g %g %g %g %g %g %g\n"
	       , itr, j_sp, drfdy[j_sp][0], drfdy[j_sp][1], drfdy[j_sp][2], drfdy[j_sp][3], drfdy[j_sp][4], drfdy[j_sp][5], drfdy[j_sp][6], drfdy[j_sp][7]);      
      
    }
#endif
		
    for(j_sp=0;j_sp<n_sp;j_sp++) {
      for(i_sp=0;i_sp<n_sp;i_sp++) {
	if(i_sp!=j_sp) AA[j_sp * n_sp + i_sp] =    -dt*drfdy[j_sp][i_sp];
	else           AA[j_sp * n_sp + i_sp] = 1.0-dt*drfdy[j_sp][i_sp];
      }
    }
    
    for(i_sp=0;i_sp<n_sp;i_sp++) {
      ddy[i_sp]=rf_tot[i_sp]*dt - dy[i_sp];
    }
		    
    dgesv_(&n_sp, &inc, AA, &n_sp, ipiv, ddy, &n_sp, &info);    //compile with -llapack -lblas

    /****************************************/  
    for(i_sp=0;i_sp<n_sp;i_sp++) {
      ddy[i_sp] *= BETA;
    }

    err_max=0.0;
    bound_ddy(yi, ddy, y_max, itr, flag_reset, &err_max);

    for(i_sp=0;i_sp<n_sp;i_sp++) {
      dy[i_sp] += ddy[i_sp]; //dy[] is recalculated at the end of the loop
      yi[i_sp] += ddy[i_sp];      
    }
    
    for(i_sp=0;i_sp<n_sp;i_sp++) {
#ifdef HMINUS_EQ
      if(i_sp == sHM) continue;	
#endif
      
      if(yi[i_sp]!=0.0) err=fabs(ddy[i_sp]/yi[i_sp]);
      else              err=0.0;
      if(err>err_max && yi[i_sp] > 1.0e-10) err_max=err, err_flag = i_sp;
    }
    
    for(i_sp=0;i_sp<n_sp;i_sp++)
      dy[i_sp] = yi[i_sp] - y_init[i_sp];
    
    if(itr >= 20) ErrMaxCrit = 1.e-4;
    if(err_max < ErrMaxCrit) flag_conv++;
    else flag_conv = 0;

    if(flag_reset <=3 && itr - 70*(flag_reset+1) > 0 && err_max > 1){
      flag_reset++;
      goto chem_START;
    }

    itr++;
  }   // iteration
 
  if(err_max<1.0e-6) {
    t_chem1 = 1.0e20;
    for(i_sp=0;i_sp<n_sp;i_sp++) {
      if(yi[i_sp]-y_init[i_sp]!=0.0) {
	t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
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
      t_chem1 = 1.0e20;
      for(i_sp=0;i_sp<n_sp;i_sp++) {
	if(yi[i_sp]-y_init[i_sp]!=0.0) {
	  t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
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

  /* The legacy solver evaluates chemical heating/cooling from the converged,
     unclamped Newton solution.  Keep that ordering while making state carry
     the new abundances consumed by chemical_cooling(). */
  chemical_cooling(environment, state, work, dt, abund_init);

  for(i_sp=0;i_sp<n_sp;i_sp++){
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

static inline int ReactionRates_nometal(int itr, double nH, double *yi, double *k, double *k_inverse, double *rf_tot)
{
  int i_sp, n_sp = N_SP;
  int i_re, n_re = N_RE;
  double rate;
  double rf[n_re][n_sp],ni[n_sp];
	
  enum REACTIONS reactions;
  enum SPECIES species;

  for(i_sp=0;i_sp<n_sp;i_sp++) {
    ni[i_sp]=yi[i_sp]*nH;
  }
  
  for(i_re=0;i_re<n_re;i_re++)
    for(i_sp=0;i_sp<n_sp;i_sp++) 
      rf[i_re][i_sp]=0.0;
	
  rate = k[1]*ni[sHI]*ni[sELEC];
  rf[RE1][sELEC] =  rate;
  rf[RE1][sHI] = -rate;
  rf[RE1][sHII] =  rate;

  rate = k[2]*ni[sHII]*ni[sELEC];
  rf[RE2][sELEC]  = -rate;
  rf[RE2][sHI]  =  rate;
  rf[RE2][sHII]  = -rate;
  
  rate = k[30]*ni[sHI]*ni[sHI];
  rf[RE30][sELEC] = rate;
  rf[RE30][sHI] = -rate;
  rf[RE30][sHII] =  rate;

#ifdef HELIUM
  rate = k[3]*ni[sHeI]*ni[sELEC];
  rf[RE3][sELEC]  =  rate;
  rf[RE3][sHeI]  = -rate;
  rf[RE3][sHeII]  =  rate;
  
  rate = k[4]*ni[sHeII]*ni[sELEC];
  rf[RE4][sELEC] = -rate;
  rf[RE4][sHeI] =  rate;
  rf[RE4][sHeII] = -rate;
  
  rate = k[5]*ni[sHeII]*ni[sELEC];
  rf[RE5][sELEC] =  rate;
  rf[RE5][sHeII] = -rate;
  rf[RE5][sHeIII] =  rate;
  
  rate = k[6]*ni[sHeIII]*ni[sELEC];
  rf[RE6][sELEC] = -rate;
  rf[RE6][sHeII] =  rate;
  rf[RE6][sHeIII] = -rate;
#endif  

  rate = k[7]*ni[sHI]*ni[sELEC];
  rf[RE7][sELEC] = -rate;
  rf[RE7][sHI] = -rate;
  rf[RE7][sHM] = rate;
  
  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[HI]=%g ni[sELEC]=%g ni[sHM]=%g\n", RE7, rf[RE7][sH2I], k[7], ni[sHI], ni[sELEC], ni[sHM]);     

  rate = k[8]*ni[sHM]*ni[sHI];
  rf[RE8][sHM] = -rate;
  rf[RE8][sHI] = -rate;
  rf[RE8][sH2I] =  rate;
  rf[RE8][sELEC] =  rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sHI]=%g ni[sHM]=%g\n", RE8, rf[RE8][sH2I], k[8], ni[sH2I], ni[sHI], ni[sHM]);     


#ifdef H2_PLUS
  rate = k[9]*ni[sHI]*ni[sHII];
  rf[RE9][sHI] = -rate;
  rf[RE9][sHII] = -rate;
  rf[RE9][sH2II] =  rate;
  
  rate = k[10]*ni[sHI]*ni[sH2II];
  rf[RE10][sHI] = -rate;
  rf[RE10][sH2II] = -rate;
  rf[RE10][sH2I] =  rate;
  rf[RE10][sHII] =  rate;
#endif
  
  rate = k[11]*ni[sH2I]*ni[sHI];
  rf[RE11][sH2I] = -rate;
  rf[RE11][sHI] = 2.0*rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sHI]=%g\n", RE11, rf[RE11][sH2I], k[11], ni[sH2I], ni[sHI]);     
  
#ifdef H2_PLUS  
  rate = k[12]*ni[sH2I]*ni[sHII];
  rf[RE12][sH2I] = -rate;
  rf[RE12][sHII] = -rate;
  rf[RE12][sH2II] =  rate;
  rf[RE12][sHI] =  rate;
#endif  

  rate = k[13]*ni[sH2I]*ni[sELEC];
  rf[RE13][sH2I] = -rate;
  rf[RE13][sHI] = 2.0*rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sELEC]=%g\n", RE13, rf[RE13][sH2I], k[13], ni[sH2I], ni[sELEC]);     
  
#ifdef HELIUM
  rate = k[14]*ni[sHeII]*ni[sHI];
  rf[RE14][sHeII] = -rate;
  rf[RE14][sHI] = -rate;
  rf[RE14][sHeI] =  rate;
  rf[RE14][sHII] =  rate;
	
  rate = k[15]*ni[sHeI]*ni[sHII];
  rf[RE15][sHeI] = -rate;
  rf[RE15][sHII] = -rate;
  rf[RE15][sHeII] =  rate;
  rf[RE15][sHI] =  rate;
#endif	

  rate = k[16]*ni[sHM]*ni[sELEC];
  rf[RE16][sHM] = -rate;
  rf[RE16][sELEC] =  rate;
  rf[RE16][sHI] =  rate;
	
  rate = k[17]*ni[sHM]*ni[sHI];
  rf[RE17][sHM] = -rate;
  rf[RE17][sHI] =  rate;
  rf[RE17][sELEC] =  rate;
	
  rate = k[18]*ni[sHM]*ni[sHII];
  rf[RE18][sHM] = -rate;
  rf[RE18][sHII] = -rate;
  rf[RE18][sHI] = 2.0*rate;

#ifdef H2_PLUS	
  rate = k[19]*ni[sHM]*ni[sHII];
  rf[RE19][sHM] = -rate;
  rf[RE19][sHII] = -rate;
  rf[RE19][sH2II] =  rate;
  rf[RE19][sELEC] =  rate;
	
  rate = k[20]*ni[sH2II]*ni[sELEC];
  rf[RE20][sH2II] = -rate;
  rf[RE20][sELEC] = -rate;
  rf[RE20][sHI] = 2.0*rate;
	
  rate = k[21]*ni[sH2II]*ni[sHM];
  rf[RE21][sH2II] = -rate;
  rf[RE21][sHM] = -rate;
  rf[RE21][sHI] =  rate;
  rf[RE21][sH2I] =  rate;
#endif

#ifdef IONIZATION_PHOTON
  rate = k[22] * ni[sHI];
  rf[RE22][sHI]   = -rate;
  rf[RE22][sHII]  = rate;
  rf[RE22][sELEC] = rate;
#endif

#ifdef LYMANWERNER
  rate = k[25]*ni[sHM];
  rf[RE25][sHM] = -rate;
  rf[RE25][sHI] = rate;
  rf[RE25][sELEC] = rate;
	
  rate = k[29]*ni[sH2I];
  rf[RE29][sH2I] = -rate;
  rf[RE29][sHI]  = 2*rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g\n", RE29, rf[RE29][sH2I], k[29], ni[sH2I]);     
#endif
	
#ifdef HD

  /*rate = k[31]*ni[sDII]*ni[sELEC];
  rf[RE31][sDII] = -rate;
  rf[RE31][sELEC] = -rate;
  rf[RE31][sDI] =  rate;
	
  rate = k[34]*ni[sDI]*ni[sHI];
  rf[RE34][sDI] = -rate;
  rf[RE34][sHI] = -rate;
  rf[RE34][sHDI] =  rate;
  */

  rate = k[32]*ni[sDI]*ni[sHII];
  rf[RE32][sDI] = -rate;
  rf[RE32][sHII] = -rate;
  rf[RE32][sDII] =  rate;
  rf[RE32][sHI] =  rate;
	
  rate = k[33]*ni[sDII]*ni[sHI];
  rf[RE33][sDII] = -rate;
  rf[RE33][sHI] = -rate;
  rf[RE33][sDI] =  rate;
  rf[RE33][sHII] =  rate;
	
  rate = k[35]*ni[sDI]*ni[sH2I];
  rf[RE35][sDI] = -rate;
  rf[RE35][sH2I] = -rate;
  rf[RE35][sHI] =  rate;
  rf[RE35][sHDI] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sDI]=%g\n", RE35, rf[RE35][sH2I], k[35], ni[sH2I], ni[sDI]);     
  
  rate = k[37]*ni[sDII]*ni[sH2I];
  rf[RE37][sDII] = -rate;
  rf[RE37][sH2I] = -rate;
  rf[RE37][sHII] =  rate;
  rf[RE37][sHDI] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sDII]=%g\n", RE37, rf[RE37][sH2I], k[37], ni[sH2I], ni[sDII]);     
  
  rate = k[38]*ni[sHDI]*ni[sHI];
  rf[RE38][sHDI] = -rate;
  rf[RE38][sHI] = -rate;
  rf[RE38][sH2I] =  rate;
  rf[RE38][sDI] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[HDI]=%g ni[sHI]=%g\n", RE38, rf[RE38][sH2I], k[38], ni[sHDI], ni[sHI]);     
  
  rate = k[39]*ni[sHDI]*ni[sHII];
  rf[RE39][sHDI] = -rate;
  rf[RE39][sHII] = -rate;
  rf[RE39][sH2I] =  rate;
  rf[RE39][sDII] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[HDI]=%g ni[sHII]=%g\n", RE39, rf[RE39][sH2I], k[39], ni[sHDI], ni[sHII]);     
  
#ifdef HD_PLUS
  rate = k[36]*ni[sHDII]*ni[sHI];
  rf[RE36][sHDII] = -rate;
  rf[RE36][sHI] = -rate;
  rf[RE36][sHII] =  rate;
  rf[RE36][sHDI] =  rate;
	
  rate = k[40]*ni[sDI]*ni[sHII];
  rf[RE40][sDI] = -rate;
  rf[RE40][sHII] = -rate;
  rf[RE40][sHDII] =  rate;
  
  rate = k[41]*ni[sDII]*ni[sHI];
  rf[RE41][sDII] = -rate;
  rf[RE41][sHI] = -rate;
  rf[RE41][sHDII] =  rate;
  
  rate = k[42]*ni[sHDII]*ni[sELEC];
  rf[RE42][sHDII] = -rate;
  rf[RE42][sELEC] = -rate;
  rf[RE42][sHI] =  rate;
  rf[RE42][sDI] =  rate;
#endif  

#ifdef DMINUS
  rate = k[43]*ni[sDI]*ni[sELEC];
  rf[RE43][sDI] = -rate;
  rf[RE43][sELEC] = -rate;
  rf[RE43][sDM] =  rate;
  
  rate = k[44]*ni[sDII]*ni[sDM];
  rf[RE44][sDII] = -rate;
  rf[RE44][sDM] = -rate;
  rf[RE44][sDI] =  2.0*rate;
	
  rate = k[45]*ni[sHII]*ni[sDM];
  rf[RE45][sHII] = -rate;
  rf[RE45][sDM] = -rate;
  rf[RE45][sHI] =  rate;
  rf[RE45][sDI] =  rate;
  
  rate = k[46]*ni[sHM]*ni[sDI];
  rf[RE46][sHM] = -rate;
  rf[RE46][sDI] = -rate;
  rf[RE46][sHI] =  rate;
  rf[RE46][sDM] =  rate;
  
  rate = k[47]*ni[sDM]*ni[sHI];
  rf[RE47][sDM] = -rate;
  rf[RE47][sHI] = -rate;
  rf[RE47][sHM] =  rate;
  rf[RE47][sDI] =  rate;
  
  rate = k[48]*ni[sDM]*ni[sHI];
  rf[RE48][sDM] = -rate;
  rf[RE48][sHI] = -rate;
  rf[RE48][sHDI] =  rate;
  rf[RE48][sELEC] =  rate;
#endif
	
#ifdef LYMANWERNER
  rate = k[53]*ni[sHDI];
  rf[RE53][sHDI] = -rate;
  rf[RE53][sHI] = rate;
  rf[RE53][sDI] = rate;
#endif
#endif
	
#ifdef THREEBODY
  rate = k[55]*ni[sHI]*ni[sHI]*ni[sHI];
  rf[RE55][sHI] = -2.0*rate;
  rf[RE55][sH2I] =  rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[HI]=%g\n", RE55, rf[RE55][sH2I], k[55], ni[sHI]);     
  
  rate = k[57]*ni[sHI]*ni[sHI]*ni[sH2I];
  rf[RE57][sHI] = -2.0*rate;
  rf[RE57][sH2I] =  rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sHI]=%g\n", RE57, rf[RE57][sH2I], k[57], ni[sHI], ni[sH2I]);     
  
  rate = k[58]*ni[sH2I]*ni[sH2I];
  rf[RE58][sH2I] = -rate;
  rf[RE58][sHI] = 2.0*rate;

  if(FlagPrint)
    printf("rf[%d]=%g k=%g ni[H2I]=%g\n", RE58, rf[RE58][sH2I], k[58], ni[sH2I]);     
#endif

#ifdef H2_FORMATION_ON_DUST
  rate = k[51]*ni[sHI];
  rf[RE51][sHI] = -2.0*rate;
  rf[RE51][sH2I] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[sHI]=%g\n", RE51, rf[RE51][sH2I], k[51], ni[sHI]);     
#endif

#ifdef INVERSE_REACTIONS
  rate = k_inverse[1]*ni[sHII]*ni[sELEC]*ni[sELEC];
  rf[RE1_INV][sHII] = -rate;
  rf[RE1_INV][sELEC] = -rate;
  rf[RE1_INV][sHI] = rate;

  //  if(FlagPrint)
  //printf("rf[%d]=%g k_inv=%g ni[H2I]=%g ni[sELEC]=%g\n", RE1_INV, rf[RE1_INV][sH2I], k_inverse[1], ni[sH2I], ni[sELEC]);     

  rate = k_inverse[2]*ni[sHI];
  rf[RE2_INV][sHI] = -rate;
  rf[RE2_INV][sHII] = rate;
  rf[RE2_INV][sELEC] = rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k_inv=%g ni[H2I]=%g ni[sELEC]=%g\n", RE2_INV, rf[RE2_INV][sH2I], k_inverse[2], ni[sH2I], ni[sELEC]);     

  rate = k_inverse[7]*ni[sHM];
  rf[RE7_INV][sHM] = -rate;
  rf[RE7_INV][sHI] = rate;
  rf[RE7_INV][sELEC] = rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k_inv=%g ni[H2I]=%g ni[sELEC]=%g\n", RE7_INV, rf[RE7_INV][sH2I], k_inverse[7], ni[sH2I], ni[sELEC]);     

  rate = k_inverse[8]*ni[sH2I]*ni[sELEC];
  rf[RE8_INV][sH2I] = -rate;
  rf[RE8_INV][sELEC] = -rate;
  rf[RE8_INV][sHI] = rate;
  rf[RE8_INV][sHM] = rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k_inv=%g ni[H2I]=%g ni[sELEC]=%g\n", RE8_INV, rf[RE8_INV][sH2I], k_inverse[8], ni[sH2I], ni[sELEC]);     
  
  rate = k_inverse[13]*ni[sHI]*ni[sHI]*ni[sELEC];
  rf[RE13_INV][sHI] = -2*rate;
  rf[RE13_INV][sH2I] = rate;

  //  if(FlagPrint)
  //printf("rf[%d]=%g k_inv=%g ni[HI]=%g ni[sELEC]=%g\n", RE13_INV, rf[RE13_INV][sH2I], k_inverse[13], ni[sHI], ni[sELEC]);     

  rate = k_inverse[16]*ni[sHI]*ni[sELEC]*ni[sELEC];
  rf[RE16_INV][sHI] = -rate;
  rf[RE16_INV][sELEC] = -rate;
  rf[RE16_INV][sHM] = rate;

  rate = k_inverse[17]*ni[sHI]*ni[sHI]*ni[sELEC];
  rf[RE17_INV][sHI] = -rate;
  rf[RE17_INV][sELEC] = -rate;
  rf[RE17_INV][sHM] = rate;

  rate = k_inverse[18]*ni[sHI]*ni[sHI];
  rf[RE18_INV][sHI] = -2.*rate;
  rf[RE18_INV][sHII] = rate;
  rf[RE18_INV][sHM] = rate;

  rate = k_inverse[30]*ni[sHI]*ni[sHII]*ni[sELEC];
  rf[RE30_INV][sHII] = -rate;
  rf[RE30_INV][sELEC] = -rate;
  rf[RE30_INV][sHI] = rate;
#endif
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    rf_tot[i_sp]=0.0;
  }
  
  for(i_re=0;i_re<n_re;i_re++)
    for(i_sp=0;i_sp<n_sp;i_sp++){
      rf_tot[i_sp] += rf[i_re][i_sp];

      if(FlagPrint && i_sp == sH2I && fabs(rf[i_re][i_sp]) > 1.e-100)
	printf("rf[%d]=%g tot=%g tot/nH=%g\n", i_re, rf[i_re][i_sp], rf_tot[i_sp], rf_tot[i_sp]/nH);
    }
  
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    rf_tot[i_sp] /= nH;
  }
  // rf_tot in unit of dy/dt [/sec]
  
  return 1;
}



#ifdef CHEM_DIFF_ANALYTIC
static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *k, double *k_inverse, double *rf_tot)
{
  int i_sp, n_sp = N_SP, ifunc;
  double ni[n_sp];
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    ni[i_sp]=yi[i_sp]*nH;
  }
  	
  for(i_sp=0;i_sp<n_sp;i_sp++) 
    rf_tot[i_sp]=0.0;
  
  ifunc = switch_func[i_sp_diff](nH, ni, k, k_inverse, rf_tot);

  return 1;
}

int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[1]*ni[sHI];
  rf[sELEC] +=  rate;
  rf[sHI]   += -rate;
  rf[sHII]  +=  rate;

  rate = k[2]*ni[sHII];
  rf[sELEC] += -rate;
  rf[sHI]   +=  rate;
  rf[sHII]  += -rate;

#ifdef HELIUM
  rate = k[3]*ni[sHeI];
  rf[sELEC] +=  rate;
  rf[sHeI]  += -rate;
  rf[sHeII] +=  rate;
  
  rate = k[4]*ni[sHeII];
  rf[sELEC] += -rate;
  rf[sHeI]  +=  rate;
  rf[sHeII] += -rate;
    
  rate = k[5]*ni[sHeII];
  rf[sELEC]  +=  rate;
  rf[sHeII]  += -rate;
  rf[sHeIII] +=  rate;
  
  rate = k[6]*ni[sHeIII];
  rf[sELEC]  += -rate;
  rf[sHeII]  +=  rate;
  rf[sHeIII] += -rate;
#endif  

  rate = k[7]*ni[sHI];
  rf[sELEC] += -rate;
  rf[sHI]   += -rate;
  rf[sHM]   += rate;

  rate = k[13]*ni[sH2I];
  rf[sH2I] += -rate;
  rf[sHI]  += 2.0*rate;

  rate = k[16]*ni[sHM];
  rf[sHM]   += -rate;
  rf[sELEC] +=  rate;
  rf[sHI]   +=  rate;

#ifdef H2_PLUS
  rate = k[20]*ni[sH2II];
  rf[sH2II] += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += 2.0*rate;
#endif

#ifdef HD
  /*rate = k[31]*ni[sDII];
  rf[sDII]  += -rate;
  rf[sELEC] += -rate;
  rf[sDI]   +=  rate;*/

#ifdef HD_PLUS
  rate = k[42]*ni[sHDII];
  rf[sHDII] += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   +=  rate;
  rf[sDI]   +=  rate;
#endif
  
#ifdef DMINUS
  rate = k[43]*ni[sDI];
  rf[sDI]   += -rate;
  rf[sELEC] += -rate;
  rf[sDM]   +=  rate;
#endif
#endif


#ifdef INVERSE_REACTIONS
  rate = 2*k_inverse[1]*ni[sHII]*ni[sELEC];
  rf[sHII]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;

  rate = k_inverse[8]*ni[sH2I];
  rf[sH2I]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;
  rf[sHM]   += rate;

  rate = k_inverse[13]*ni[sHI]*ni[sHI];
  rf[sHI]   += -2*rate;
  rf[sH2I]  += rate;

  rate = 2*k_inverse[16]*ni[sHI]*ni[sELEC];
  rf[sHI]   += -rate;
  rf[sELEC] += -rate;
  rf[sHM]   += rate;

  rate = k_inverse[17]*ni[sHI]*ni[sHI];
  rf[sHI]   += -rate;
  rf[sELEC] += -rate;
  rf[sHM]   += rate;

  rate = k_inverse[30]*ni[sHI]*ni[sHII];
  rf[sHII]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;
#endif
	

  return 0;
}

int ReactionRate_diff_nometal_HI(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[1]*ni[sELEC];
  rf[sELEC] +=  rate;
  rf[sHI]   += -rate;
  rf[sHII]  +=  rate;

  rate = k[7]*ni[sELEC];
  rf[sELEC] += -rate;
  rf[sHI]   += -rate;
  rf[sHM]   += rate;
  
  rate = k[8]*ni[sHM];
  rf[sHM]   += -rate;
  rf[sHI]   += -rate;
  rf[sH2I]  +=  rate;
  rf[sELEC] +=  rate;
  
  rate = 2*k[30]*ni[sHI];
  rf[sELEC] += rate;
  rf[sHI]   += -rate;
  rf[sHII]  +=  rate;


#ifdef H2_PLUS
  rate = k[9]*ni[sHII];
  rf[sHI]   += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
  
  rate = k[10]*ni[sH2II];
  rf[sHI]   += -rate;
  rf[sH2II] += -rate;
  rf[sH2I]  +=  rate;
  rf[sHII]  +=  rate;
#endif
  
  rate = k[11]*ni[sH2I];
  rf[sH2I] += -rate;
  rf[sHI]  += 2.0*rate;

#ifdef HELIUM
  rate = k[14]*ni[sHeII];
  rf[sHeII] += -rate;
  rf[sHI]   += -rate;
  rf[sHeI]  +=  rate;
  rf[sHII]  +=  rate;
#endif
    
  rate = k[17]*ni[sHM];
  rf[sHM]   += -rate;
  rf[sHI]   +=  rate;
  rf[sELEC] +=  rate;

#ifdef IONIZATION_PHOTON
  rate = k[22];
  rf[sHI]   += -rate;
  rf[sHII]  += rate;
  rf[sELEC] += rate;
#endif

#ifdef HD
  rate = k[33]*ni[sDII];
  rf[sDII] += -rate;
  rf[sHI]  += -rate;
  rf[sDI]  +=  rate;
  rf[sHII] +=  rate;
	
  /*rate = k[34]*ni[sDI];
  rf[sDI]  += -rate;
  rf[sHI]  += -rate;
  rf[sHDI] +=  rate;*/
  
  rate = k[38]*ni[sHDI];
  rf[sHDI] += -rate;
  rf[sHI]  += -rate;
  rf[sH2I] +=  rate;
  rf[sDI]  +=  rate;

#ifdef HD_PLUS
  rate = k[36]*ni[sHDII];
  rf[sHDII] += -rate;
  rf[sHI]   += -rate;
  rf[sHII]  +=  rate;
  rf[sHDI]  +=  rate;       
  
  rate = k[41]*ni[sDII];
  rf[sDII]  += -rate;
  rf[sHI]   += -rate;
  rf[sHDII] +=  rate;
#endif
  
#ifdef DMINUS
  rate = k[47]*ni[sDM];
  rf[sDM] += -rate;
  rf[sHI] += -rate;
  rf[sHM] +=  rate;
  rf[sDI] +=  rate;
  
  rate = k[48]*ni[sDM];
  rf[sDM]   += -rate;
  rf[sHI]   += -rate;
  rf[sHDI]  +=  rate;
  rf[sELEC] +=  rate;
#endif
#endif

#ifdef THREEBODY
  rate = 3.*k[55]*ni[sHI]*ni[sHI];
  rf[sHI]  += -2.0*rate;
  rf[sH2I] +=  rate;
    
  rate = 2.*k[57]*ni[sHI]*ni[sH2I];
  rf[sHI]  += -2.0*rate;
  rf[sH2I] +=  rate;
#endif

#ifdef H2_FORMATION_ON_DUST
  rate = k[51];
  rf[sHI]  += -2.0*rate;
  rf[sH2I] +=  rate;
#endif

#ifdef INVERSE_REACTIONS
  rate = k_inverse[2];
  rf[sHI]   += -rate;
  rf[sHII]  += rate;
  rf[sELEC] += rate;

  rate = 2*k_inverse[13]*ni[sHI]*ni[sELEC];
  rf[sHI]   += -2*rate;
  rf[sH2I]  += rate;

  rate = k_inverse[16]*ni[sELEC]*ni[sELEC];
  rf[sHI]   += -rate;
  rf[sELEC] += -rate;
  rf[sHM]   += rate;

  rate = 2*k_inverse[17]*ni[sHI]*ni[sELEC];
  rf[sHI]   += -rate;
  rf[sELEC] += -rate;
  rf[sHM]   += rate;

  rate = 2*k_inverse[18]*ni[sHI];
  rf[sHI]   += -2.*rate;
  rf[sHII]  += rate;
  rf[sHM]   += rate;

  rate = k_inverse[30]*ni[sHII]*ni[sELEC];
  rf[sHII]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;
#endif
	

  return 0;
}

int ReactionRate_diff_nometal_HII(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[2]*ni[sELEC];
  rf[sELEC] += -rate;
  rf[sHI]   +=  rate;
  rf[sHII]  += -rate;

#ifdef H2_PLUS
  rate = k[9]*ni[sHI];
  rf[sHI]   += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
    
  rate = k[12]*ni[sH2I];
  rf[sH2I]  += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
  rf[sHI]   +=  rate;
#endif

#ifdef HELIUM
  rate = k[15]*ni[sHeI];
  rf[sHeI]  += -rate;
  rf[sHII]  += -rate;
  rf[sHeII] +=  rate;
  rf[sHI]   +=  rate;
#endif	

  rate = k[18]*ni[sHM];
  rf[sHM]  += -rate;
  rf[sHII] += -rate;
  rf[sHI]  += 2.0*rate;
	
#ifdef H2_PLUS
  rate = k[19]*ni[sHM];
  rf[sHM]   += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
  rf[sELEC] +=  rate;
#endif

#ifdef HD
  rate = k[32]*ni[sDI];
  rf[sDI]  += -rate;
  rf[sHII] += -rate;
  rf[sDII] +=  rate;
  rf[sHI]  +=  rate;
	
  rate = k[39]*ni[sHDI];
  rf[sHDI] += -rate;
  rf[sHII] += -rate;
  rf[sH2I] +=  rate;
  rf[sDII] +=  rate;

#ifdef HD_PLUS  
  rate = k[40]*ni[sDI];
  rf[sDI]   += -rate;
  rf[sHII]  += -rate;
  rf[sHDII] +=  rate;
#endif

#ifdef DMINUS	
  rate = k[45]*ni[sDM];
  rf[sHII] += -rate;
  rf[sDM]  += -rate;
  rf[sHI]  +=  rate;
  rf[sDI]  +=  rate;
#endif
#endif


#ifdef INVERSE_REACTIONS
  rate = k_inverse[1]*ni[sELEC]*ni[sELEC];
  rf[sHII]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;

  rate = k_inverse[30]*ni[sHI]*ni[sELEC];
  rf[sHII]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;
#endif
	

  return 0;
}

#ifdef HELIUM
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[3]*ni[sELEC];
  rf[sELEC] +=  rate;
  rf[sHeI]  += -rate;
  rf[sHeII] +=  rate;
    
  rate = k[15]*ni[sHII];
  rf[sHeI]  += -rate;
  rf[sHII]  += -rate;
  rf[sHeII] +=  rate;
  rf[sHI]   +=  rate;

  return 0;
}

int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[4]*ni[sELEC];
  rf[sELEC] += -rate;
  rf[sHeI]  +=  rate;
  rf[sHeII] += -rate;
  
  rate = k[5]*ni[sELEC];
  rf[sELEC]  +=  rate;
  rf[sHeII]  += -rate;
  rf[sHeIII] +=  rate;

  rate = k[14]*ni[sHI];
  rf[sHeII] += -rate;
  rf[sHI]   += -rate;
  rf[sHeI]  +=  rate;
  rf[sHII]  +=  rate;

  return 0;
}

int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[6]*ni[sELEC];
  rf[sELEC]  += -rate;
  rf[sHeII]  +=  rate;
  rf[sHeIII] += -rate;

  return 0;
}
#endif

int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;
 
  rate = k[11]*ni[sHI];
  rf[sH2I] += -rate;
  rf[sHI]  += 2.0*rate;
  
#ifdef H2_PLUS
  rate = k[12]*ni[sHII];
  rf[sH2I]  += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
  rf[sHI]   +=  rate;
#endif
  
  rate = k[13]*ni[sELEC];
  rf[sH2I] += -rate;
  rf[sHI]  += 2.0*rate;

#ifdef LYMANWERNER
  rate = k[29];
  rf[sH2I] += -rate;
  rf[sHI]  += 2*rate;
#endif
  
#ifdef HD
  rate = k[35]*ni[sDI];
  rf[sDI]  += -rate;
  rf[sH2I] += -rate;
  rf[sHI]  +=  rate;
  rf[sHDI] +=  rate;

  rate = k[37]*ni[sDII];
  rf[sDII] += -rate;
  rf[sH2I] += -rate;
  rf[sHII] +=  rate;
  rf[sHDI] +=  rate;  
#endif
	
#ifdef THREEBODY
  rate = k[57]*ni[sHI]*ni[sHI];
  rf[sHI]  += -2.0*rate;
  rf[sH2I] +=  rate;
  
  rate = 2 * k[58]*ni[sH2I];
  rf[sH2I] += -rate;
  rf[sHI]  += 2.0*rate;
#endif


#ifdef INVERSE_REACTIONS
  rate = k_inverse[8]*ni[sELEC];
  rf[sH2I]  += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += rate;
  rf[sHM]   += rate;
#endif
	

  return 0;
}

#ifdef H2_PLUS
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[10]*ni[sHI];
  rf[sHI]   += -rate;
  rf[sH2II] += -rate;
  rf[sH2I]  +=  rate;
  rf[sHII]  +=  rate;

  rate = k[20]*ni[sELEC];
  rf[sH2II] += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   += 2.0*rate;
	
  rate = k[21]*ni[sHM];
  rf[sH2II] += -rate;
  rf[sHM]   += -rate;
  rf[sHI]   +=  rate;
  rf[sH2I]  +=  rate;

  return 0;
}
#endif

int ReactionRate_diff_nometal_HM(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[8]*ni[sHI];
  rf[sHM]   += -rate;
  rf[sHI]   += -rate;
  rf[sH2I]  +=  rate;
  rf[sELEC] +=  rate;

  rate = k[16]*ni[sELEC];
  rf[sHM]   += -rate;
  rf[sELEC] +=  rate;
  rf[sHI]   +=  rate;
	
  rate = k[17]*ni[sHI];
  rf[sHM]   += -rate;
  rf[sHI]   +=  rate;
  rf[sELEC] +=  rate;
	
  rate = k[18]*ni[sHII];
  rf[sHM]  += -rate;
  rf[sHII] += -rate;
  rf[sHI]  += 2.0*rate;
	
#ifdef H2_PLUS
  rate = k[19]*ni[sHII];
  rf[sHM]   += -rate;
  rf[sHII]  += -rate;
  rf[sH2II] +=  rate;
  rf[sELEC] +=  rate;
	
  rate = k[21]*ni[sH2II];
  rf[sH2II] += -rate;
  rf[sHM]   += -rate;
  rf[sHI]   +=  rate;
  rf[sH2I]  +=  rate;
#endif

#ifdef LYMANWERNER
  rate = k[25];
  rf[sHM]   += -rate;
  rf[sHI]   += rate;
  rf[sELEC] += rate;
#endif

#ifdef HD
#ifdef DMINUS
  rate = k[46]*ni[sDI];
  rf[sHM] += -rate;
  rf[sDI] += -rate;
  rf[sHI] +=  rate;
  rf[sDM] +=  rate;
#endif
#endif


#ifdef INVERSE_REACTIONS
  rate = k_inverse[7];
  rf[sHM]   += -rate;
  rf[sHI]   += rate;
  rf[sELEC] += rate;
#endif


  return 0;
}

#ifdef HD
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;
	
  rate = k[32]*ni[sHII];
  rf[sDI]  += -rate;
  rf[sHII] += -rate;
  rf[sDII] +=  rate;
  rf[sHI]  +=  rate;
  
  /*rate = k[34]*ni[sHI];
  rf[sDI]  += -rate;
  rf[sHI]  += -rate;
  rf[sHDI] +=  rate;*/
  
  rate = k[35]*ni[sH2I];
  rf[sDI]  += -rate;
  rf[sH2I] += -rate;
  rf[sHI]  +=  rate;
  rf[sHDI] +=  rate;

#ifdef HD_PLUS	
  rate = k[40]*ni[sHII];
  rf[sDI]   += -rate;
  rf[sHII]  += -rate;
  rf[sHDII] +=  rate;
#endif  

#ifdef DMINUS
  rate = k[43]*ni[sELEC];
  rf[sDI]   += -rate;
  rf[sELEC] += -rate;
  rf[sDM]   +=  rate;
  
  rate = k[46]*ni[sHM];
  rf[sHM] += -rate;
  rf[sDI] += -rate;
  rf[sHI] +=  rate;
  rf[sDM] +=  rate;
#endif

  return 0;
}

int ReactionRate_diff_nometal_DII(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  /*rate = k[31]*ni[sELEC];
  rf[sDII]  += -rate;
  rf[sELEC] += -rate;
  rf[sDI]   +=  rate;*/
	
  rate = k[33]*ni[sHI];
  rf[sDII] += -rate;
  rf[sHI]  += -rate;
  rf[sDI]  +=  rate;
  rf[sHII] +=  rate;
	
  rate = k[37]*ni[sH2I];
  rf[sDII] += -rate;
  rf[sH2I] += -rate;
  rf[sHII] +=  rate;
  rf[sHDI] +=  rate;
  
#ifdef HD_PLUS
  rate = k[41]*ni[sHI];
  rf[sDII]  += -rate;
  rf[sHI]   += -rate;
  rf[sHDII] +=  rate;
#endif

#ifdef DMINUS  
  rate = k[44]*ni[sDM];
  rf[sDII] += -rate;
  rf[sDM]  += -rate;
  rf[sDI]  +=  2.0*rate;
#endif

  return 0;
}

int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[38]*ni[sHI];
  rf[sHDI] += -rate;
  rf[sHI]  += -rate;
  rf[sH2I] +=  rate;
  rf[sDI]  +=  rate;
  
  rate = k[39]*ni[sHII];
  rf[sHDI] += -rate;
  rf[sHII] += -rate;
  rf[sH2I] +=  rate;
  rf[sDII] +=  rate;
  
#ifdef LYMANWERNER
  rate = k[53];
  rf[sHDI] += -rate;
  rf[sHI]  += rate;
  rf[sDI]  += rate;
#endif

  return 0;
}

#ifdef HD_PLUS
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[36]*ni[sHI];
  rf[sHDII] += -rate;
  rf[sHI]   += -rate;
  rf[sHII]  +=  rate;
  rf[sHDI]  +=  rate;
	
  rate = k[42]*ni[sELEC];
  rf[sHDII] += -rate;
  rf[sELEC] += -rate;
  rf[sHI]   +=  rate;
  rf[sDI]   +=  rate;

  return 0;
}
#endif

#ifdef DMINUS
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *k,double *k_inverse,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[44]*ni[sDII];
  rf[sDII] += -rate;
  rf[sDM]  += -rate;
  rf[sDI]  +=  2.0*rate;
	
  rate = k[45]*ni[sHII];
  rf[sHII] += -rate;
  rf[sDM]  += -rate;
  rf[sHI]  +=  rate;
  rf[sDI]  +=  rate;
  
  rate = k[47]*ni[sHI];
  rf[sDM] += -rate;
  rf[sHI] += -rate;
  rf[sHM] +=  rate;
  rf[sDI] +=  rate;
  
  rate = k[48]*ni[sHI];
  rf[sDM]   += -rate;
  rf[sHI]   += -rate;
  rf[sHDI]  +=  rate;
  rf[sELEC] +=  rate;

  return 0;
}
#endif
#endif

int func_other(double nH,double *ni,double *k,double *k_inverse,double *rf){
  return 1;
}

static void init_switch_func(void){
  for(int jj = 0;jj < 256;jj++)
    switch_func[jj] = &func_other; 
  
  switch_func[sELEC] = &ReactionRate_diff_nometal_Elec;
  switch_func[sHI]   = &ReactionRate_diff_nometal_HI;
  switch_func[sHII]  = &ReactionRate_diff_nometal_HII;
#ifdef HELIUM
  switch_func[sHeI]   = &ReactionRate_diff_nometal_HeI;
  switch_func[sHeII]  = &ReactionRate_diff_nometal_HeII;
  switch_func[sHeIII] = &ReactionRate_diff_nometal_HeIII;
#endif
  switch_func[sH2I]   = &ReactionRate_diff_nometal_H2I;
#ifdef H2_PLUS
  switch_func[sH2II]  = &ReactionRate_diff_nometal_H2II; 
#endif
  switch_func[sHM]    = &ReactionRate_diff_nometal_HM; 
#ifdef HD
  switch_func[sDI]    = &ReactionRate_diff_nometal_DI;
  switch_func[sDII]   = &ReactionRate_diff_nometal_DII; 
  switch_func[sHDI]   = &ReactionRate_diff_nometal_HDI; 
#ifdef HD_PLUS
  switch_func[sHDII]  = &ReactionRate_diff_nometal_HDII;
#endif
#ifdef DMINUS
  switch_func[sDM]    = &ReactionRate_diff_nometal_DM; 
#endif
#endif
}
#endif

static int evaluate_reaction_rates(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace){
  const double rho_cgs = environment->rho_cgs;
  const double n_nuclei = rho_cgs / PROTONMASS;
  const double nH = XH * n_nuclei;

  const double T_gas = workspace->thermodynamics.temperature;
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_cgs, T_gas);
#else
  const double l_core = step->jeans_length;
#endif
#ifdef METALLICITY
  const double metallicity = environment->metallicity;
#endif
  
  //update_radiation_state(sph_in, sph_out, sph_state, T_gas);
  
  /*********** initialize ************/
  for(int j=0; j<k_MAX; j++)
    workspace->rates.forward[j] = 0.0;
  /***********************************/

  /* ------- Compute temperature index. */
  const double ln_T0 = log(T_start);
  const double ln_T9 = log(T_end);
  
  const double log_T  = log10(T_gas);
  const double log_T2 = log_T * log_T;
  const double log_T3 = log_T * log_T * log_T;
  
  double ln_T  = log(T_gas);
  if(ln_T < ln_T0) ln_T = ln_T0;
  if(ln_T > ln_T9) ln_T = ln_T9;
  int i_T = (int) ((ln_T-ln_T0)/dln_T);
  
  if(i_T<0) i_T=0;
  if(i_T>N_T-2) i_T=N_T-2;

  if(i_T < 0){
    printf("i_T=%d (ln_T=%g,ln_T0=%g,dln_T=%g)\n",i_T,ln_T,ln_T0,dln_T);

    printf("ThisTask %d ID %d T %g n %g H2I %g HI %g HII %g\n"
	   , ThisTask, environment->ID, T_gas, n_nuclei, state->abund[sH2I], state->abund[sHI], state->abund[sHII]);
    fflush(stdout);
    
    printf("\nStop.(in chemistry_noneq_chem_implicit)\n");
    /* savepositions(All.SnapshotFileCount++); */
    
    endrun(); /* just debugging */
  }
  
  
  if (i_T > N_T-2) i_T = N_T-2;  /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  const double c_T = (T_gas - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);
  
  for(int j=1; j<=21; j++){
    workspace->rates.forward[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
  }
  
  workspace->rates.forward[30] = k_table[30][i_T] + c_T*(k_table[30][i_T+1] - k_table[30][i_T]);
  workspace->rates.forward[22] = workspace->rates.forward[23] = workspace->rates.forward[24] = 0.;

#ifdef IONIZATION_PHOTON
  workspace->rates.forward[22] = step->k22;
  workspace->rates.forward[23] = workspace->rates.forward[24] = 0.;
#endif
  
#ifdef LYMANWERNER    /* with Photo-Dissociation */
  workspace->rates.forward[25] = 0.;
  workspace->rates.forward[29] = 0.;
  workspace->rates.forward[53] = 0.;

  double J21 = 0.;
#ifdef LW_RADIATION_LOCAL
  workspace->rates.forward[25] = environment->rate_Hminus_pd;
  workspace->rates.forward[29] = environment->rate_H2_diss;
#ifdef RYTR_NHD
  workspace->rates.forward[53] = environment->rate_HD_diss;
#else
  workspace->rates.forward[53] = 1.55/ 1.39 * workspace->rates.forward[29];      /* optically-thin limit found in Wolcott-Green & Haiman 2011 */
#endif
#endif // LW_RADIATION_LOCAL

#ifdef BACKGROUND_LW_RADIATION
  workspace->rates.forward[25] += All.rate_Hminus_pd;
  
  double k29_ext = All.rate_H2_diss;
#ifdef HD
  double k53_ext = 1.55/1.39 * All.rate_H2_diss;
#endif
#ifdef SHIELD_LW
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double column_H2 = l_core * n_nuclei * state->abund[sH2I];
#ifdef HD
  const double column_HI = l_core * n_nuclei * state->abund[sHI];
  const double column_HD = l_core * n_nuclei * state->abund[sHDI];
  const double shield_HI = EachShield(column_HI, 1.49e-1, 1.62e+0, 2.85e23);
  const double shield_H2 = EachShield(column_H2, 5.20e-3, 2.38e-1, 2.34e19);
#endif
  const double shield_H2self = SelfShield(column_H2, T_gas,
                                          9.12 * sqrt(T_gas / 1.e4));
#ifdef HD
  const double shield_HDself = SelfShield(column_HD, T_gas,
                                          7.45 * sqrt(T_gas / 1.e4));
#endif
#else
  const double shield_H2self = step->shield_H2self;
#ifdef HD
  const double shield_HI = step->shield_HI;
  const double shield_H2 = step->shield_H2;
  const double shield_HDself = step->shield_HDself;
#endif
#endif
  const double f_print_H2self = shield_H2self;
#else
  const double f_print_H2self = 1.0;
#endif
#ifdef SHIELD_LW
  k29_ext *= shield_H2self /* SphP[ithis].f_byHI*/;
#ifdef HD
  k53_ext *= shield_HDself * shield_HI * shield_H2;
#endif
#endif //SHIELD_LW

  workspace->rates.forward[29] += k29_ext;
#ifdef HD
  workspace->rates.forward[53] += k53_ext;
#endif

  if(FlagPrint)
    printf("k25=%g %g k29=%g %g f=%g\n", workspace->rates.forward[25], All.rate_Hminus_pd, workspace->rates.forward[29], k29_ext, f_print_H2self);
  
#endif //BACKGROUND_LW_RADIATION
#endif // LYMANWERNER

#ifdef HD
  for(int j=31; j<=48; j++)
    workspace->rates.forward[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif  

#ifdef THREEBODY
  for(int j=55; j<=58; j++)
    workspace->rates.forward[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif

  /****************************************************************************************************************/
  /******************************** We evaluate continuum opacities here ******************************************/
  /****************************************************************************************************************/

#ifdef DUST_COOLING
  double edot_dust = 0.;
#endif

  double tau_cont=0.;
#ifdef CONT_OPACITY
  struct opacities opac;

  double dummy = calc_continuum_opacities(T_gas, rho_cgs, l_core, &opac);
  tau_cont = opac.tau;
  
#ifdef DUST_COOLING
  double stellar_flux = 0.;
#ifdef THERMAL_FEEDBACK
  stellar_flux = environment->stellar_rad_flux;
#endif
  workspace->continuum.dust_temperature = dust_temperature(rho_cgs, T_gas, metallicity, l_core, stellar_flux, &opac, &edot_dust);
  workspace->continuum.dust_cooling_rate = edot_dust;

#ifdef RYTR_TAU_DUST
  workspace->continuum.kappa_dust_Planck = opac.kappa_dust_Planck;
  if(workspace->continuum.kappa_dust_Planck > 2.06 * metallicity)
    workspace->continuum.kappa_dust_Planck = 2.06 * metallicity;
#endif

  tau_cont = opac.tau;
#endif
  
  workspace->continuum.tau = tau_cont;
  workspace->continuum.reduction_factor = 1./(1. + 3./4.*tau_cont*tau_cont);
#endif
  /****************************************************************************************************************/
  /****************************************************************************************************************/
  /****************************************************************************************************************/

#ifdef INVERSE_REACTIONS
  double n_crit = H_criticaldensity(nH, T_gas, nH*state->abund[sHI]*l_core, state->abund[sHI], state->abund[sELEC]);

  workspace->rates.inverse[1] = k_table_inverse[1][i_T] + c_T*(k_table_inverse[1][i_T+1] - k_table_inverse[1][i_T]);
  workspace->rates.inverse[2] = k_table_inverse[2][i_T] + c_T*(k_table_inverse[2][i_T+1] - k_table_inverse[2][i_T]);
  workspace->rates.inverse[7] = k_table_inverse[7][i_T] + c_T*(k_table_inverse[7][i_T+1] - k_table_inverse[7][i_T]);
  workspace->rates.inverse[8] = k_table_inverse[8][i_T] + c_T*(k_table_inverse[8][i_T+1] - k_table_inverse[8][i_T]);
  workspace->rates.inverse[13] = k_table_inverse[13][i_T] + c_T*(k_table_inverse[13][i_T+1] - k_table_inverse[13][i_T]);
  workspace->rates.inverse[16] = k_table_inverse[16][i_T] + c_T*(k_table_inverse[16][i_T+1] - k_table_inverse[16][i_T]);
  workspace->rates.inverse[17] = k_table_inverse[17][i_T] + c_T*(k_table_inverse[17][i_T+1] - k_table_inverse[17][i_T]);
  workspace->rates.inverse[18] = k_table_inverse[18][i_T] + c_T*(k_table_inverse[18][i_T+1] - k_table_inverse[18][i_T]);

  workspace->rates.inverse[30] = k_table_inverse[30][i_T] + c_T*(k_table_inverse[30][i_T+1] - k_table_inverse[30][i_T]);
  workspace->rates.inverse[55] = k_table_inverse[55][i_T] + c_T*(k_table_inverse[55][i_T+1] - k_table_inverse[55][i_T]);

  //workspace->rates.inverse[1] *= nH / (n_crit + nH); //temporally
  workspace->rates.inverse[18] *= nH / (n_crit + nH);

  //printf("tau=%g n_crit=%g\n", tau_cont, n_crit);
  
  workspace->rates.inverse[2] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  workspace->rates.inverse[7] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  //workspace->rates.inverse[2] *= (1 - exp(-tau_cont));
  //workspace->rates.inverse[7] *= (1 - exp(-tau_cont));
  workspace->rates.forward[25]        *= exp(-tau_cont);

  if(FlagPrint){
    printf("k2=%g k_inv=%g %g %g %g %g %g %g %g %g %g tau_cont=%g n=%g n_c=%g l_core=%g\n"
	   , workspace->rates.forward[2]
	   , workspace->rates.inverse[1], workspace->rates.inverse[2], workspace->rates.inverse[7], workspace->rates.inverse[8]
	   , workspace->rates.inverse[13], workspace->rates.inverse[16], workspace->rates.inverse[17], workspace->rates.inverse[18]
	   , workspace->rates.inverse[30], workspace->rates.inverse[55], tau_cont, n_nuclei, n_crit, l_core);
  }


#ifdef ID_CHEM_DEBUG
  if(FlagPrint){
    printf("k_inv=%g %g tau_cont=%g n=%g n_c=%g l_core=%g i_T=%d c_T=%g\n"
	   , workspace->rates.inverse[2], workspace->rates.inverse[7], tau_cont, n_nuclei, n_crit, l_core, i_T, c_T);
  }
#endif
#endif
  
#ifdef H2_FORMATION_ON_DUST
  workspace->rates.forward[51] = H2form_on_grains(rho_cgs, T_gas, workspace->continuum.dust_temperature, metallicity);
#endif

  if(T_gas > 10.){
#ifdef MARTIN97	
    double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21;
    double p, n_c1, n_c2, log_g_h1, log_g_l1, log_g_h2, log_g_l2, log_n_c1, log_n_c2;
    double log_gamma_cd, log_gamma_dt;

    a1  = -1.784239e2;
    a2  = -6.842243e1; 
    a3  =  4.320243e1;
    a4  = -4.633167;
    a5  =  6.970086e1;
    a6  =  4.087038e4;
    a7  = -2.370570e4;
    a8  =  1.288953e2;
    a9  = -5.391334e1;
    a10 =  5.315517;
    a11 = -1.973427e1;
    a12 =  1.678905e4;
    a13 = -2.578611e4;
    a14 =  1.482123e1;
    a15 = -4.890915;
    a16 =  4.749030e-1;
    a17 = -1.338283e2;
    a18 = -1.164408;
    a19 =  8.227443e-1; 
    a20 =  5.864073e-1;
    a21 = -2.056313;


    log_g_h1 = a1 + a2 * log_T + a3 * log_T2 + a4 * log_T3 + a5 * log10(1 + a6/T_gas);
    log_g_h2 = a7/T_gas;

    log_g_l1 = a8 + a9 * log_T + a10 * log_T2 + a11 * log10(1 + a12/T_gas);
    log_g_l2 = a13/T_gas;
  
    log_n_c1 = a14 + a15*log_T + a16*log_T2 + a17/T_gas;
    n_c1 = pow(10.0, log_n_c1);

    log_n_c2 = a18 + log_n_c1;
    n_c2 = pow(10.0, log_n_c2);

    p = a19 + a20 *exp(-T_gas/1850.0) + a21 * exp(-T_gas/440.0); 
  
    log_gamma_cd = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+pow(nH/n_c1, p)) + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+pow(nH/n_c2, p));
  
    a1  = -1.427664e2;
    a2  =  4.270741e1; 
    a3  = -2.027365;
    a4  = -2.582097e-1;
    a5  =  2.136094e1;
    a6  =  2.753531e4;
    a7  = -2.146779e4;
    a8  =  6.034928e1;
    a9  = -2.743096e1;
    a10  = 2.676150;
    a11 = -1.128215e1;
    a12 =  1.425455e4;
    a13 = -2.312520e4;
    a14 =  9.305564;
    a15 = -2.464009;
    a16 =  1.985955e-1;
    a17 =  7.430600e2;
    a18 = -1.174242;
    a19 =  7.502286e-1; 
    a20 =  2.358848e-1;
    a21 =  2.937507;

  
    log_g_h1 = a1 + a2*log_T + a3*log_T2 + a4*log_T3 + a5*log10(1+a6/T_gas);
    log_g_h2 = a7/T_gas;

    log_g_l1 = a8 + a9*log_T + a10*log_T2 + a11*log10(1+a12/T_gas);
    log_g_l2 = a13/T_gas;

    log_n_c1 = a14 + a15*log_T + a16*log_T2 + a17/T_gas;
    n_c1 = pow(10.0, log_n_c1);

    log_n_c2 = a18 + log_n_c1;
    n_c2 = pow(10.0, log_n_c2);
  
    p = a19 + a20*exp(-T_gas/1850.0) + a21*exp(-T_gas/440.0); 
  
    log_gamma_dt = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+pow(nH/n_c1, p)) + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+pow(nH/n_c2, p));  

    workspace->rates.forward[11] = pow(10.0, log_gamma_cd) + pow(10.0, log_gamma_dt);

    if(FlagPrint)
      printf("k11=%g gamma=%g %g\n", workspace->rates.forward[11], log_gamma_cd, log_gamma_dt);
   

#ifdef INVERSE_REACTIONS
    double xk_L = workspace->rates.forward[11];
    double xk_H  = workspace->rates.inverse[55];

    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    double a = 1.0/(1.0 + nH/rn_cr);

    if(a == 1.0)       workspace->rates.forward[11] = xk_L;
    else if(a == 0.0)  workspace->rates.forward[11] = xk_H;
    else               workspace->rates.forward[11] = pow(xk_H,1.0-a)*pow(xk_L,a);      
#endif
#else
    double xk_L  = 1.12e-10*exp(-7.035e+4/T_gas);
    double xk_H  = 6.5e-7/sqrt(T_gas)*exp(-5.2e+4/T_gas)*(1.0-exp(-6.e+3/T_gas));
    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    double a = 1.0/(1.0 + nH/rn_cr);
  
    if(a == 1.0)       workspace->rates.forward[11] = xk_L;
    else if(a == 0.0)  workspace->rates.forward[11] = xk_H;
    else               workspace->rates.forward[11] = pow(xk_H,1.0-a)*pow(xk_L,a);  
#endif
  }else
    workspace->rates.forward[11] = 0.;


#if defined(STROMGREN_TEST)
  for(int j=1; j<=21; j++){
#ifdef IONIZATION_PHOTON
    if(j == 2)  continue;
#endif
    workspace->rates.forward[j] = 0.;
  }

  workspace->rates.forward[30] = 0.;
#endif
  return 0;
}


void correct_initial_abundances(double *yi){
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
      else if(yi[sHI] + yi[sHeI] > 1.05*diff){
        diff1 = 0.97*yi[sHI];
        diff2 = diff - diff1;

	yi[sHI]  -= diff1;
        yi[sHeI] -= diff2;

	yi[sHII]  += diff1;
        yi[sHeII] += diff2;
      }else if(yi[sHI] + 2*yi[sHeI] + yi[sHeII] > 1.05*diff){
        diff1 = 0.97*yi[sHI];
        diff2 = dmin(0.9*yi[sHeI], 0.5*(diff - diff1)) ;
        diff3 = diff - diff1 - 2*diff2;

	yi[sHI]   -= diff1;
        yi[sHeI]  -= diff2;
        yi[sHeII] -= diff3;

	yi[sHII]   += diff1;
        yi[sHeIII] += (diff3+diff2);
      }else{
        diff1 = 0.97 * (yi[sHI] + 2*yi[sHeI] + yi[sHeII]);
        diff2 = diff - diff1;

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

void give_initial_guess(double T_gas, double nH, double dt, double *yi, double *y_init, double *k){
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


void bound_ddy(double *yi, double *ddy, double *y_max, int itr, int flag_reset, double *err_max){
  for(int i_sp=0;i_sp<N_SP;i_sp++) {
    if(yi[i_sp] + ddy[i_sp] < 0.1 * yi[i_sp]){       
      ddy[i_sp] = -0.9 * yi[i_sp];       
      *err_max += (double) (i_sp + 1);	
    }else if(yi[i_sp]+ddy[i_sp] > y_max[i_sp]){
      ddy[i_sp] = dmin(0.9*yi[i_sp], 0.9 * (y_max[i_sp] - yi[i_sp]));
      *err_max += 100. * (double) (i_sp + 1);
    }
  }    
}



#ifdef BACKGROUND_LW_RADIATION
void update_background_J21(void){
#ifdef TYPE_B
  double Teff = 1.e4;
#else
  double Teff = 1.e5;
#endif
  double J21 = All.J21 * 1.e-21;

  double hc = PLANCK * CLIGHT_REAL;
  const double E_LW = 13.6;
  double J_13eV = 2. * pow(E_LW * eV_to_erg, 3) / hc / hc / (exp(E_LW*eV_to_erg/BOLTZMANN/Teff) - 1);
  double Norm = J21 / J_13eV;

  double PD_rate, g0;
  get_Hminus_PD_rate(Teff, &PD_rate, &g0);

  All.rate_Hminus_pd = Norm * PD_rate;
#ifdef DUST_COOLING
  All.G0           = Norm * g0 / 1.6e-3;
#endif
  All.rate_H2_diss = 1.39e9 * J21;
}
#endif

int InitChem(void)
{
  enum SPECIES species;
  enum REACTIONS reactions;
	
  /* ----- cross-sections, etc ... */
  int ifunc = chem_tabulate();

#ifdef CHEM_DIFF_ANALYTIC
  init_switch_func();
#endif
 
#ifdef CONTINUUM_COOLING_MATSUKOBA
  read_opacity_tables();
  tabulate_cooling_rate_cont_thin();
#endif
#ifdef DUST_COOLING
  read_dust_opacity_tables();
#endif

#ifdef RYTR_TAU_LW
  tabulate_opacities_in_LW_band();
#endif

#ifdef LYMANWERNER
  make_rate_tables_LW_radiation();
#ifdef BACKGROUND_LW_RADIATION
  update_background_J21();
#endif
#endif
  
#ifdef METAL_LINE_COOLING_SD93
  read_continuum_cooling_rate_table_SD93();
#endif

#if defined(INIT_TEMPERATURE) || defined(FIX_TEMPERATURE)
  compute_init_gamma();
#endif

#ifdef CHEMISTRY_CORRECT_NEGATIVE_ABUNDANCE
  All.AbundMinimum[sELEC] = 1.e-10;
  All.AbundMinimum[sHI]  = 1.e-10;
  All.AbundMinimum[sHII] = 1.e-10;
  All.AbundMinimum[sH2I] = 1.e-10;
  All.AbundMinimum[sHM]  = 1.e-20;
#ifdef HD
  All.AbundMinimum[sDI]  = 1.e-20;
  All.AbundMinimum[sDII] = 1.e-20;
  All.AbundMinimum[sHDI] = 1.e-20;
#endif
#endif

  N_SP = sEND - sELEC;
  N_RE = REEND - RE1;

  return(ifunc);  
}


  /*--- [List of Chemistry] ---*/
  /*
    ------- 1:	HI    + e   -> HII   + 2e    hydrogen collisional ionization
    ------- 2:	HII   + e   -> HI    + p     hydrogen radiative recombination
    ------- 3:	HeI   + e   -> HeII  + 2e    He collisional ionization
    ------- 4:	HeII  + e   -> HeI   + p     He recombination 
    ------- 5:	HeII  + e   -> HeIII + 2e    HeII collisional ionization
    ------- 6:	HeIII + e   -> HeII  + p     HeII recombination 
    ------- 7:	HI    + e   -> HM    + p     photo-attachment
    ------- 8:	HM    + HI  -> H2I*  + e     H2 formation path (H-)
    ------- 9:  HI    + HII -> H2II  + p     H2+ channel
    ------- 10:	H2II  + HI  -> H2I*  + HII   H2+ channel
    ------- 11:	H2I   + HI  -> 3HI
    ------- 12:	H2I   + HII -> H2II  + HI    charge transfer
    ------- 13:	H2I   + e   -> 2HI   + e     collisional dissociation
    ------- 14: HeII  + HI  -> HeI   + HII + p  charge transfer   
    ------- 15: HeI   + HII -> HeII  + HI  
    ------- 16:	HM    + e   -> HI    + 2e
    ------- 17:	HM    + HI  -> 2H    + e
    ------- 18:	HM    + HII -> 2HI
    ------- 19:	HM    + HII -> H2II  + e
    ------- 20:	H2II  + e   -> 2HI
    ------- 21:	H2II  + HM  -> HI    + H2I

    ------- 22:	HI    + p   -> HII   + e     photo-ionization H
    ------- 23:	HeII  + p   -> HeIII + e     HeII ionization
    ------- 24:	HeI   + p   -> HeII  + e     photo-ionization He
    ------- 25:	HM    + p   -> HI    + e     photo-detachment
    ------- 26:	H2II  + p   -> HI    + HII   
    ------- 27:	H2I   + p   -> H2II  + e     
    ------- 28:	H2II  + p   -> 2HII  + e
    ------- 29:	H2I   + p   -> 2HI           photo-dissociation      

    ------- 60,61: H2    + e   -> HM    + HI
    ------- 62,63: HI    + HI  -> H2II  + e

    HD (Nakamura-Umemura complete model)
     
    ------- 31:  DII  + e    -> DI   + p      radiative recombination 
    ------- 32:  DI   + HII  -> DII  + HI     charge transfer
    ------- 33:  DII  + HI   -> DI   + HII    reverse reaction of 32
    ------- 34:  DI   + HI   -> HDI  + p      * neutral formation path
    ------- 35:  DI   + H2I  -> HI   + HDI    * neutral path by H2
    ------- 36:  HDII + HI   -> HII  + HDI    charge transfer
    ------- 37:  DII  + H2I  -> HII  + HDI    * D+ formation path
    ------- 38:  HDI  + HI   -> H2I  + DI     
    ------- 39:  HDI  + HII  -> H2I  + DII    
    ------- 40:  DI   + HII  -> HDII + p 
    ------- 41:  DII  + HI   -> HDII + p 
    ------- 42:  HDII + e    -> HI   + DI     collisional dissociation 
    ------- 43:  DI   + e    -> DM   + p      D- formation/photo-attachment
    ------- 44:  DII  + DM   -> 2DI           mutual neutrization
    ------- 45:  HII  + DM   -> DI   + HI     mutual neutrization
    ------- 46:  HM   + DI   -> HI   + DM 
    ------- 47:  DM   + HI   -> DI   + HM 
    ------- 48:  DM   + HI   -> HDI  + e      * D- formation path

    ------- 50:  DI   + p    -> DII  + e      photo-ionization D   
    ------- 51:  HDII + p    -> HI   + DII    photo-dissociation HD+  
    ------- 52:  HDII + p    -> HII  + DI     photo-dissociation HD+  
    ------- 53:  HDI  + p    -> HI   + DI     photo-dissociation HD  

    3 body reactions
   
    ------- 55:  HI + HI + HI -> H2I + HI  
    ------- 56:  H2I + HI     -> HI + HI + HI  same as reaction 11
    ------- 57:  HI + HI + H2I-> H2I + H2I
    ------- 58:  H2I + H2I    -> HI + HI + H2I  reverse reaction
  */ 


  /* probably this part should be parallelized in the near future */
int chem_tabulate(void)
{
 
  /* ----- Declare other values. */
  int i,j;
  double ln_T,ln_T_eV,T_eV,log10_T;
  double ln_T_eV2,ln_T_eV3,ln_T_eV4,ln_T_eV5;
  double ln_T_eV6,ln_T_eV7,ln_T_eV8,ln_T_eV9;
  double z_HI, z_HII, z_H2I, z_HM, z_elec;

  enum SPECIES species;

  chem_e[sELEC] = 0.;
  chem_e[sHI]   = 216.034;
  chem_e[sHII]  = 1528.084;
  chem_e[sH2I]  = 0.;
#ifdef H2_PLUS
  chem_e[sH2II] = 1488.365;
#endif
  chem_e[sHM]   = 143.246;
#ifdef HD
  chem_e[sDI]   = 219.804;
  chem_e[sDII]  = 1532.211;
  chem_e[sHDI]  = 0.328;
#endif

  /*************** choose functional forms of the H2 & HD cooling rates ***************/
#ifdef H2_COOLING_G15
  H2_cooling = H2_cooling_G15;
#elif defined(H2_COOLING_GP12)
  H2_cooling = H2_cooling_GP12;
#elif defined(H2_COOLING_Y06)
  H2_cooling = H2_cooling_Y06;
#elif defined(H2_COOLING_O00)
  H2_cooling = H2_cooling_O00;
#else
  H2_cooling = H2_cooling_G15;
  //printf("you should choose which H2 cooling function will be used in the simulation!\n");
  //endrun();
#endif

#ifdef HD
  HD_cooling = HD_cooling_GP12;
#endif

#ifdef METALLICITY
  OxygenI_cooling = OxygenI_twolevel_cooling;
  CarbonI_cooling = CarbonI_twolevel_cooling;
  CarbonII_cooling = CarbonII_twolevel_cooling;
#endif

  /*************** H2 opacity (fitting based on F18 or recalculate it as O00 do) ******/    
#ifdef H2_OPACITY
#ifdef H2_COOLING_O00
  H2_line_escape_fraction = fesc_unity;
#elif defined(H2_OPACITY_FIT)
  H2_line_escape_fraction = fesc_fit;  
#else
  H2_line_escape_fraction = H2line_reduction;
#endif
#endif

#if defined(HD) && defined(H2_OPACITY)
#ifdef H2_COOLING_O00
  HD_line_escape_fraction = fesc_unity;
#else
  HD_line_escape_fraction = fesc_HD_fit;
#endif
#endif

  /*************** continuum opacity (tables in M19 or Yoshida-san's formulation) *****/    
#ifdef CONT_OPACITY
  continuum_escape_fraction = calc_continuum_opacities;
#endif

  /*************** partition function of H2 ************************************** *****/    

#ifdef PARTITION_H2_EBM1989
  partition_H2 = partition_H2_EBM1989;
#elif defined(PARTITION_H2_I87)
  partition_H2 = partition_H2_Irwin1987;
#elif defined(PARTITION_H2_R85)
  partition_H2 = partition_H2_Rossi1985;
#else
  partition_H2 = partition_H2_EBM1989_fitted_by_Matsukoba;
#endif

  check_partition_functions_H2();

#ifdef H2_PLUS
  partition_H2p = pf_z_H2p;
#endif

  /* ----- Compute temperature grid. */
  dln_T= (log(T_end) - log(T_start))/N_T;
  
  for( i = 0;i< N_T;i++){
    ln_T = log(T_start) + i*dln_T;
    T_table[i] = exp(ln_T);
  }

  /* ----- Initialize arrays. */
  for(i=0;i<N_T;i++){
    
    for(j=0;j<k_MAX;j++) k_table[j][i] = 0;
    for(j=0;j<k_MAX;j++) k_table_inverse[j][i] = 0;

    //for(j=0;j<NREACTION_HI;j++) k_table_HI[j][i] = 0;
    //for(j=0;j<NREACTION_HI;j++) k_table_HI_inverse[j][i] = 0;
        
    ceHIa[i]    = 0.;
    ceHeIa[i]   = 0.;
    ceHeIIa[i]  = 0.;
    ciHIa[i]    = 0.;
    ciHeIa[i]   = 0.;
    ciHeISa[i]  = 0.;
    ciHeIIa[i]  = 0.;
#ifndef CONTINUUM_COOLING_MATSUKOBA  
    reHIIa[i]   = 0.;
#ifdef HELIUM
    reHeII1a[i] = 0.;
    reHeII2a[i] = 0.;
    reHeIIIa[i] = 0.;
#endif
    brema[i]    = 0.;
    H2_CIE[i]   = 0.;
#endif
  }
  
      
  /* ----- Fill in tables over the range T_start to T_end. */
  for(i=0 ; i<N_T ; i++){
    z_H2I = partition_H2(T_table[i]);
    z_HI  = partition_H(T_table[i]);
#ifdef H2_PLUS
    double z_H2II= partition_H2p(T_table[i]);
#endif

    z_HM  = 1.;
    z_elec = 2.;
    z_HII = 1.;

    /* ------- Compute various values of T. */
    ln_T = log(T_table[i]);
    log10_T = log10(T_table[i]);
    T_eV = T_table[i]/eV_to_K;
    ln_T_eV = log(T_eV);
    
    ln_T_eV2=ln_T_eV*ln_T_eV;
    ln_T_eV3=ln_T_eV*ln_T_eV2;
    ln_T_eV4=ln_T_eV*ln_T_eV3;
    ln_T_eV5=ln_T_eV*ln_T_eV4;
    ln_T_eV6=ln_T_eV*ln_T_eV5;
    ln_T_eV7=ln_T_eV*ln_T_eV6;
    ln_T_eV8=ln_T_eV*ln_T_eV7;
    ln_T_eV9=ln_T_eV*ln_T_eV8;
    
#ifdef CASE_B
    double lambda_HI   = 2 * 157807.0/T_table[i];
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    double lambda_HeI  = 2 * 285335.0/T_table[i];
    double lambda_HeII = 2 * 631515.0/T_table[i];
#endif
#endif

    
    /************* H collisional ionization ****************/
    k_table[1][i] = exp(-32.71396786375 
			+ 13.53655609057     *ln_T_eV
			- 5.739328757388     *ln_T_eV2 
			+ 1.563154982022     *ln_T_eV3
			- 0.2877056004391    *ln_T_eV4
			+ 0.03482559773736999*ln_T_eV5
			- 0.00263197617559   *ln_T_eV6
			+ 0.0001119543953861 *ln_T_eV7
			- 2.039149852002e-6  *ln_T_eV8);



    /************** hydrogen recombination *****************/  
#ifdef CASE_B
    k_table[2][i] = 2.753e-14 * pow(315614.0/T_table[i], 1.5)/pow(1.0+pow(115188.0/T_table[i], 0.406), 2.242);
    /* Ferland et al. 1992 fitted by Hui & Gnedin 1997 */
#else
    k_table[2][i] = exp(-28.61303380689232
			- 0.7241125657826851    *ln_T_eV
			- 0.02026044731984691   *ln_T_eV2
			- 0.002380861877349834  *ln_T_eV3
			- 0.0003212605213188796 *ln_T_eV4
			- 0.00001421502914054107*ln_T_eV5
			+ 4.989108920299513e-6  *ln_T_eV6
			+ 5.755614137575758e-7  *ln_T_eV7
			- 1.856767039775261e-8  *ln_T_eV8
			- 3.071135243196595e-9  *ln_T_eV9);
#endif

#ifdef CHEM_HIRANO
    if(T_table[i] < 5500.)
      k_table[2][i] = 0.;
#endif

    k_table[30][i] = 0.;
#ifdef O05_CHEM
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[30][i] = 1.2e-17 * pow(T_table[i], 1.2) * exp(-1.578e5/T_table[i]);
#endif

    /*************** Reactions related to Helium ****************************/
    /*************** k3, 4, 5, 6, 14, 15          ****************************/

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    /* He collisional ionization */
    k_table[3][i] = exp(-44.09864886561001
			+ 23.91596563469      *ln_T_eV
			- 10.75323019821      *ln_T_eV2
			+ 3.058038757198      *ln_T_eV3
			- 0.5685118909884001  *ln_T_eV4
			+ 0.06795391233790001 *ln_T_eV5
			- 0.005009056101857001*ln_T_eV6
			+ 0.0002067236157507  *ln_T_eV7
			- 3.649161410833e-6   *ln_T_eV8);


    /* He recombination (HeII -> HeI)  */
#ifdef CASE_B
    k_table[4][i] = 1.26e-14 * pow(lambda_HeI, 0.75) /* Case B radiative recombination */
      + 1.9e-3 * pow(T_table[i],-1.5)*exp( -470000./T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] )); /* dielectric recombination */
    /* Ferland et al. 1992 fitted by Hui & Gnedin 1997 */
#else
    //k_table[4][i] =  3.92e-13/pow(T_eV,0.6353); /* Case A radiative recombination */
    //1.9e-3 * pow(T_table[i],-1.5)*exp( -470000./T_table[i] )
    //	*(1.+0.3*exp( -94000./T_table[i] )); /* dielectric recombination */
			
    k_table[4][i] = 1.54e-9*(1.+0.3/exp(8.099328789667/T_eV))
      / (exp(40.49664394833662/T_eV)*pow(T_eV,1.5)) /* Dielectric */
      + 3.92e-13/pow(T_eV,0.6353);  /* radiative */
#endif
      
    /* HeII collisional ionization */
    k_table[5][i] = exp(-68.71040990212001
			+ 43.93347632635      *ln_T_eV
			- 18.48066993568      *ln_T_eV2
			+ 4.701626486759002   *ln_T_eV3
			- 0.7692466334492     *ln_T_eV4
			+ 0.08113042097303    *ln_T_eV5
			- 0.005324020628287001*ln_T_eV6
			+ 0.0001975705312221  *ln_T_eV7
			- 3.165581065665e-6   *ln_T_eV8);

    /* HeII recombinatiom (HeIII -> HeII) */
#ifdef CASE_B
    k_table[6][i] = 2.0 * 2.753e-14 * pow(315614.0/(T_table[i]/4.0), 1.5)/pow(1.0+pow(115188.0/(T_table[i]/4.0), 0.406), 2.242);
    /* 2 x k2(T/4) */ 
#else
    k_table[6][i] = 3.36e-10/sqrt(T_table[i])/pow(T_table[i]/1.e3,0.2)/(1+pow(T_table[i]/4.e6,0.7)); /* Fukugita-Kawasaki 94 */
#endif

    
    /* 14: HeII  + HI  -> HeI   + HII + p  charge transfer */
#ifdef ENZO_CHEM
    k_table[14][i] = 0.;
#else
    k_table[14][i] = 1.2e-15 * pow(T_table[i]/300.0, 0.25);
#endif    
    
    /* 15: HeI   + HII -> HeII  + HI */  
#ifdef ENZO_CHEM
    k_table[15][i] = 0.;
#else
    if ( T_table[i] > 10000.0)
      k_table[15][i] = 4.0e-37*pow(T_table[i], 4.74);
    else
      k_table[15][i] = 1.26e-9*pow(T_table[i], -0.75) * exp(-127500.0/T_table[i]);
#endif        
#endif


    if(T_eV < 0.01) {
      k_table[1][i] = 0;
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
      k_table[3][i] = 0;
#ifdef CASE_B
      k_table[4][i] = 1.26e-14 * pow(lambda_HeI, 0.75); // He radiative recombination only 
#else
      k_table[4][i] = 3.92e-13/pow(T_eV,0.6353);
#endif
      k_table[5][i] = 0;
#endif
    }
      
      
    /*************** Reactions related to H-      ****************************/
    /*************** k7, 8, 16, 17, 18            ****************************/
    
    /* H- formation by photo-attachment */
    /* old one k_table[7][i] = 6.77e-15*pow(T_eV,0.8779); */
#ifdef ENZO_CHEM
    k_table[7][i] = 6.775e-15*pow(T_eV,0.8779);
#elif defined(O05_CHEM)
    k_table[7][i] = 1.4e-18 * pow(T_table[i], 0.928) * exp(-T_table[i]/1.62e4);
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    if(T_table[i] < 6000.){
      k_table[7][i] = pow(10, -17.845 
			  + 0.762*log10_T 
			  + 0.1523*log10_T*log10_T 
			  - 0.03274*log10_T*log10_T *log10_T);
    }else{
      k_table[7][i] = pow(10, -16.4199 
			  + 0.1998   *log10_T*log10_T 
			  - 5.447e-3 *log10_T*log10_T*log10_T*log10_T 
			  + 4.0415e-5*log10_T*log10_T*log10_T*log10_T*log10_T*log10_T);
    } //Wishart1979
#else
    k_table[7][i] = 1.4e-18 * pow(T_table[i], 0.928) * exp(-T_table[i]/16200.0);  /* from GP98 */
#endif 
    
    /* H2 formation path H- channel */
#ifdef ENZO_CHEM
    /* old one based on Brwon \& Dalgarno */
    if (T_eV > 0.1) {
      k_table[8][i] = exp(-20.06913897587003
			  + 0.2289800603272916   *ln_T_eV
			  + 0.03599837721023835  *ln_T_eV2
			  - 0.004555120027032095 *ln_T_eV3
			  - 0.0003105115447124016*ln_T_eV4
			  + 0.0001073294010367247*ln_T_eV5
			  - 8.36671960467864e-6  *ln_T_eV6
			  + 2.238306228891639e-7 *ln_T_eV7);
    }else{
      k_table[8][i] = 1.43e-9;
    }
#elif defined(O05_CHEM)
    if(T_table[i] > 300.0) /* GP98, based on Launey et al. 1991 */
      k_table[8][i] = 4.0e-9 * pow(T_table[i], -0.17);
    else
      k_table[8][i] = 1.516e-9;    
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[8][i] = 1.35e-9* ( pow(T_table[i], 0.098493) + 0.32852 * pow(T_table[i],0.5561) + 2.771e-7 * pow(T_table[i],2.1826) )/
      (1. + 6.191e-3*pow(T_table[i],1.0461) + 8.9712e-11*pow(T_table[i], 3.0424) + 3.2576e-14*pow(T_table[i],3.7741)); //Kreckel+2010
#else
    if(T_table[i] > 300.0) /* GP98, based on Launey et al. 1991 */
      k_table[8][i] = 4.0e-9 * pow(T_table[i], -0.17);
    else
      k_table[8][i] = 1.516e-9;
#endif



    /* HM    + e   -> HI    + 2e */
    if(T_eV > 0.04){
      k_table[16][i] = exp(-18.01849334273
			   + 2.360852208681     *ln_T_eV
			   - 0.2827443061704    *ln_T_eV2
			   + 0.01623316639567   *ln_T_eV3
			   - 0.03365012031362999*ln_T_eV4
			   + 0.01178329782711   *ln_T_eV5
			   - 0.001656194699504  *ln_T_eV6
			   + 0.0001068275202678 *ln_T_eV7
			   - 2.631285809207e-6  *ln_T_eV8);
    }else
      k_table[16][i] =  0;
    

    
    /* HM    + HI  -> 2H    + e  */
    if(T_eV > 0.1){
      k_table[17][i] = exp(-20.37260896533324
			   + 1.139449335841631  *ln_T_eV
			   - 0.1421013521554148 *ln_T_eV2
			   + 0.00846445538663   *ln_T_eV3
			   - 0.0014327641212992 *ln_T_eV4
			   + 0.0002012250284791 *ln_T_eV5
			   + 0.0000866396324309 *ln_T_eV6
			   - 0.00002585009680264*ln_T_eV7
			   + 2.4555011970392e-6 *ln_T_eV8
			   - 8.06838246118e-8   *ln_T_eV9);
    }else{
      k_table[17][i] = 2.56e-9*pow(T_eV,1.78186);
    }

    
    /* HM    + HII -> 2HI */
#ifdef ENZO_CHEM
    k_table[18][i] = 6.5e-9/sqrt(T_eV);
#elif defined(CHEM_HIRANO)
    k_table[18][i] = 6.3e-8 + 5.7e-6/sqrt(T_table[i]) - 9.2e-11 * sqrt(T_table[i]) + 4.4e-13 * T_table[i];
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[16][i] = 0.;
    k_table[17][i] = 0.;
    k_table[18][i] = 2.4e-6*pow(T_table[i], -0.5) * (1. + T_table[i]/2.e4);       //Croft+1999
#else
    k_table[18][i] = 6.92e-14*pow(T_table[i], 0.5) + 8.00e-13*pow(T_table[i], 0.83) + 6.18e-7 *pow(T_table[i],-0.27);
#endif    

    /*************** Reactions related to H2II ****************************/
    /*************** k9, 10, 12, 19, 20, 21    ****************************/

#ifdef H2_PLUS
    /* 9: HI    + HII -> H2II  + p     H+ channel  */
    k_table[9][i] = 1.85e-23*pow(T_table[i],1.8);
    if (T_table[i] > 6.7e3) 
      k_table[9][i] = 5.81e-16*pow(T_table[i]/56200,-0.6657*log10(T_table[i]/56200));


    /* 10: H2II  + HI  -> H2I*  + HII   H+ channel */
    k_table[10][i] = 6.0e-10;

    /* 12: H2I   + HII -> H2II  + HI    charge transfer */
#ifdef ENZO_CHEM
    k_table[12][i] = exp(-24.249146877095
			 + 3.400824447095291     *ln_T_eV
			 - 3.898003964650152     *ln_T_eV2
			 + 2.045587822403071     *ln_T_eV3
			 - 0.5416182856220388    *ln_T_eV4
			 + 0.0841077503763412    *ln_T_eV5
			 - 0.007879026154483455  *ln_T_eV6
			 + 0.0004138398421504563 *ln_T_eV7
			 - 9.36345888928611e-6   *ln_T_eV8);
#else
    if(T_table[i] < 10000.0)
      k_table[12][i]=3.0e-10*exp(-21050.0/T_table[i]);
    else
      k_table[12][i]=1.5e-10*exp(-14000.0/T_table[i]);
#endif

    /* HM    + HII -> H2II  + e */
    if (T_table[i] > 1.0e4)
      k_table[19][i] = 4.0e-4*pow(T_table[i],-1.4)*exp(-15100.0/T_table[i]);
    else
      k_table[19][i] = 1.0e-8*pow(T_table[i],-0.4);

    /* 20: H2II  + e   -> 2HI */
#ifdef ENZO_CHEM
    if(T_table[i] > 617)
      k_table[20][i] = 1.32e-6 * pow(T_table[i],-0.76); 
    else
      k_table[20][i] = 1.e-8; 
#else
    /* old incorrect one k_table[20][i] = 5.56396e-8/pow(T_eV,0.6035); */
    k_table[20][i] = 2.0e-7/pow(T_table[i],0.5);  
#endif    
    
    /* 21: H2II  + HM  -> HI    + H2I  */
    k_table[21][i] = 4.64e-8/sqrt(T_eV);
#endif 


    /* k11 : H2I   + H   -> 3HI */
    /* k13 : H2I   + e   -> 2HI   + e     collisional dissociation */
    if(T_eV > 0.3) {
      k_table[11][i] = 1.0670825e-10*pow(T_eV,2.012)/
	(exp(4.463/T_eV)*pow(1+0.2472*T_eV,3.512));  
      /* low density limit, see k56 */
      
      k_table[13][i] = 4.38e-10*exp(-102000.0/T_table[i])*pow(T_table[i],0.35);
    }else{
      k_table[11][i] = 0; 
      k_table[13][i] = 0;
    }
#ifdef ENZO_CHEM
    k_table[11][i] = 1.0670825e-10*pow(T_eV,2.012)/
      (exp(4.463/T_eV)*pow(1+0.2472*T_eV,3.512));  
    k_table[13][i] = 5.6e-11*exp(-102124.0/T_table[i])*pow(T_table[i],0.5);
#elif defined(O05_CHEM)
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[13][i] = 1.91e-9 * pow(T_table[i], 0.136) * exp(-53407.1/T_table[i]);
#endif    


#ifdef INVERSE_REACTIONS
    double K_T;
    double Ene, dE, Z_pf, Mass;
    double m_HI, m_HII, m_HM, m_elec, m_H2I, m_H2II;

    m_HI = m_HII = m_HM = PROTONMASS;
    m_H2I = m_H2II = 2. * m_HI;
    m_elec = ELECTRONMASS;

    Ene = 2 * M_PI * BOLTZMANN * T_table[i] / PLANCK / PLANCK;


    /**** HI + e- -> Hp + 2e-****/
    if(T_table[i] > 1.e3){
      dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHII] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_elec/ z_HII / z_elec / z_elec;
      Mass = m_HI * m_elec/ m_HII / m_elec / m_elec; 
      
      K_T = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
      k_table_inverse[1][i] = k_table[1][i] * K_T;
    }



    /**** Hp + e- -> HI  ****/
    dE   = (chem_e[sHII] + chem_e[sELEC] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HII * z_elec / z_HI;
    Mass = m_HII * m_elec / m_HI;     
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);    
    k_table_inverse[2][i] = k_table[2][i] * K_T; //note this value is optically thick limit.


    /**** HI + e- -> HM  ****/
    dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHM]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_elec / z_HM;
    Mass = m_HI * m_elec / m_HM;     
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);    
    k_table_inverse[7][i] = k_table[7][i] * K_T; //note this value is optically thick limit.


    /**** HI + HM -> H2I + e-  ****/
    dE   = (chem_e[sHI] + chem_e[sHM] - chem_e[sH2I] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_HM / z_H2I / z_elec;
    Mass = m_HI * m_HM / m_H2I / m_elec;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);    
    k_table_inverse[8][i] = k_table[8][i] * K_T;

#ifdef H2_PLUS
    /**** H2I + Hp -> H2II + HI  ****/
    dE   = (chem_e[sH2I] + chem_e[sHII] - chem_e[sH2II] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_H2I * z_HII / z_HI / z_H2II;
    Mass = m_H2I * m_HII / m_HI / m_H2II;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);    
    //k_table_inverse[10][i] = k_table[10][i] * K_T;
    k_table[10][i] = k_table[12][i] * K_T;
#endif

    if(T_table[i]>200.){
      /**** H2I + e- -> 2HI + e-  ****/
      dE   = (chem_e[sH2I] + chem_e[sELEC] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_elec / z_HI / z_HI / z_elec;
      Mass = m_H2I * m_elec / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
      k_table_inverse[13][i] = k_table[13][i] * K_T;
    }


    if(T_table[i]>200.){
      /**** HM + e- -> HI + 2e-  ****/
      dE   = (chem_e[sHM] + chem_e[sELEC] - chem_e[sHI] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_elec / z_HI / z_elec / z_elec;
      Mass = m_HM * m_elec / m_HI / m_elec / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
#ifdef AIU
      k_table_inverse[16][i] = k_table[16][i] * K_T;
#endif    
    }



#ifdef AIU
    if(T_table[i]>200.){
      /**** HM + HI -> 2HI + e-  ****/
      dE   = (chem_e[sHM] + chem_e[sHI] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_HI / z_HI / z_HI / z_elec;
      Mass = m_HM * m_HI / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
      k_table_inverse[17][i] = k_table[17][i] * K_T;
    }
#endif

    /**** HM + Hp -> HI + HI  ****/
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);    
    k_table_inverse[18][i] = k_table[18][i] * K_T;

#ifdef CONTINUUM_COOLING_MATSUKOBA
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);    
    //k_table_HI_inverse[1][i] = k_table_HI[1][i] * K_T;
#endif

    if(T_table[i] > 1.e3){
      /***** HI + HI -> Hp + HI + e- *****/
      dE   = (2*chem_e[sHI] - chem_e[sHII] - chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI / z_HII / z_elec / z_HI;
      Mass = m_HI * m_HI / m_HII / m_elec / m_HI;     
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
      k_table_inverse[30][i] = k_table[30][i] * K_T;
    }

#endif


#ifdef HD
    /* deuterium network */
    //k_table[31][i]=3.6e-12*pow(T_table[i]/300.0, -0.75);   
    //k_table[34][i]=1.0e-25;
    /* new rate Savin (2002)  (old one may be 3.7e-10*pow(T_table[i],0.28)*exp(-43.0/T_table[i]) ) */
    //k_table[32][i]=2.0e-10*pow(T_table[i],0.402)*exp(-37.1/T_table[i]) - 3.31e-17*pow(T_table[i],1.48);     
    //k_table[33][i]=2.06e-10*pow(T_table[i],0.396)*exp(-33.0/T_table[i]) + 2.03e-9*pow(T_table[i],-0.332);  
    //k_table[37][i]=1.6e-9;/* old one may be 2.1e-9 */
    /* new rate Savin (2002)  (old one may be 3.7e-10*pow(T_table[i],0.28)  )                */

#ifdef CHEM_HIRANO
    k_table[31][i]=3.6e-12*pow(T_table[i]/300.0, -0.75);   
    k_table[32][i]=2.0e-10*pow(T_table[i],0.402)*exp(-37.1/T_table[i]) - 3.31e-17*pow(T_table[i],1.48);     
    k_table[33][i]=2.06e-10*pow(T_table[i],0.396)*exp(-33.0/T_table[i]) + 2.03e-9*pow(T_table[i],-0.332);  
    k_table[34][i]=1.0e-25;
    k_table[37][i]=1.6e-9;
#else
    k_table[31][i]=0.;
    k_table[32][i]=3.7e-10*pow(T_table[i],0.28)*exp(-43.0/T_table[i]);
    k_table[33][i]=3.7e-10*pow(T_table[i],0.28);
    k_table[34][i]=0.;
    k_table[37][i]=2.1e-9;
#endif

    k_table[35][i]=9.0e-11*exp(-3876.0/T_table[i]); 
    k_table[38][i]=3.2e-11*exp(-3624.0/T_table[i]);
    k_table[39][i]=1.0e-9*exp(-464.0/T_table[i]);


#ifdef HD_PLUS
    k_table[36][i]=6.4e-10;
    if(T_table[i] > 10.0){
      k_table[40][i]=pow(10.0, -19.38-1.523*log10_T+1.118*log10_T*log10_T-0.1269*log10_T*log10_T*log10_T);
      k_table[41][i]=pow(10.0, -19.38-1.523*log10_T+1.118*log10_T*log10_T-0.1269*log10_T*log10_T*log10_T);
    }else
      k_table[40][i]=k_table[41][i]=0;
    
    k_table[42][i]=7.2e-8/sqrt(T_table[i]);
#endif

#ifdef DMINUS
    /* D- part */
    k_table[43][i]=3.0e-16*pow(T_table[i]/300.0, 0.95)*exp(-T_table[i]/9320.0);   
    k_table[44][i]=5.7e-8*pow(T_table[i]/300.0, -0.5);   
    k_table[45][i]=4.6e-8*pow(T_table[i]/300.0, -0.5);   
    k_table[46][i]=6.4e-9*pow(T_table[i]/300.0, 0.41);   
    k_table[47][i]=6.4e-9*pow(T_table[i]/300.0, 0.41);   
    k_table[48][i]=1.5e-9*pow(T_table[i]/300.0, -0.1);   
#endif
#endif  
    /* 3-body reacions
       k_table[55][i]=5.5e-29/T_table[i];
       k_table[56][i]=6.5e-7/sqrt(T_table[i])*exp(-52000.0/T_table[i])*(1.0-exp(-6000.0/T_table[i]));
       k_table[57][i]=6.875e-30/T_table[i];
       k_table[58][i]=8.125e-8/sqrt(T_table[i])*exp(-52000.0/T_table[i])*(1.0-exp(-6000.0/T_table[i]));  */

    /* 3-body, new rate, consistent with equilibrium chemistry */    
    k_table[55][i]=4.5e-29/T_table[i];
    k_table[56][i]=6.5e-7/sqrt(T_table[i])*exp(-52000.0/T_table[i])*(1.0-exp(-6000.0/T_table[i]));
    k_table[57][i]=k_table[55][i]/40.0;
    k_table[58][i]=5.0e-8/sqrt(T_table[i])*exp(-52000.0/T_table[i])*(1.0-exp(-6000.0/T_table[i]));

    /* 3-body, rate coefficient uncertainty */
#ifdef TBR_FLOW07
    k_table[55][i] = 1.44e-26/pow(T_table[i], 1.54);
#elif defined(TBR_ABEL02)
    if(T_table[i] < 300){
      k_table[55][i] = 1.14e-31/pow(T_table[i], 0.38);
    }else{
      k_table[55][i] = 3.9e-30/T_table[i];
    }
#elif defined(TBR_FORR13)
    k_table[55][i] = 6.e-32/pow(T_table[i], 0.25) + 2.e-31/pow(T_table[i], 0.5);
#else
    k_table[55][i] = 6.e-32/pow(T_table[i], 0.25) + 2.e-31/pow(T_table[i], 0.5);
#endif


#ifdef INVERSE_REACTIONS
    if(T_table[i] > 500.){
      /****  HI + HI + HI -> HI + H2I  ****/
      dE   = (3*chem_e[sHI] - chem_e[sH2I] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI * z_HI / z_H2I / z_HI;
      Mass = m_HI * m_HI * m_HI / m_H2I / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);    
      k_table_inverse[55][i] = k_table[55][i] * K_T;


      /****  H2I + H2I -> H2I + HI + HI  ****/
      dE   = (2*chem_e[sH2I] - 2*chem_e[sHI] - chem_e[sH2I]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_H2I / z_HI / z_HI / z_H2I;
      Mass = m_H2I * m_H2I / m_HI / m_HI / m_H2I;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);    
      k_table[57][i] = k_table[58][i] * K_T;
    }
#endif        
  }
  
    
  /* coefficients for cooling and heating */
  for(i=0; i<N_T; i++){
#ifdef CASE_B
    double lambda_HI   = 2 * 157807.0/T_table[i];
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    double lambda_HeI  = 2 * 285335.0/T_table[i];
    double lambda_HeII = 2 * 631515.0/T_table[i];
#endif
#endif
    
    /* ------- Collisional excitations (Black 1981; Cen 1992) */
    ceHIa[i]   = 7.5e-19*exp( -118348./T_table[i] )
      /(1.+sqrt(T_table[i]/1.e5));
    
    ceHeIa[i]   = 9.1e-27*exp( -13179./T_table[i] )
      *pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5));
    
    ceHeIIa[i]  = 5.54e-17*exp( -473638./T_table[i] )
      *pow(T_table[i],-0.397)/(1.+sqrt(T_table[i]/1.e5));
    

    /* ------- Collisional ionizations (Tom's polynomial fits) */
    ciHIa[i]    = 2.18e-11*k_table[1][i];
    ciHeIa[i]   = 3.94e-11*k_table[3][i];
    ciHeIIa[i]  = 8.72e-11*k_table[5][i];

    ciHeISa[i]  = 5.01e-27*pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5))
      *exp( -55338./T_table[i] );
    
    /* ------- Recombinations */
#ifdef CASE_B
    reHIIa[i]   = 3.435e-30 * T_table[i] * pow(lambda_HI, 1.97)/pow(1.0+pow(lambda_HI/2.25, 0.376), 3.72);
    
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    reHeII1a[i] = BOLTZMANN * T_table[i] * 1.26e-14 * pow(lambda_HeI, 0.75); 
    
    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */
    
    reHeIIIa[i] = 8 * 3.435e-30 * T_table[i] * pow(lambda_HeII, 1.97)/pow(1.0+pow(lambda_HeII/2.25, 0.376), 3.72);
#endif
#else
    reHIIa[i]   = 6.50e-27*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/1.e6,0.7));

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    reHeII1a[i] = 1.55e-26*pow(T_table[i],0.3647);

    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */
    
    reHeIIIa[i] = 3.48e-26*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/4.e6,0.7));
#endif    


    /* ------- Bremsstrahlung (Black 1981)(Spitzer & Hart 1979) */
    brema[i]    = 1.43e-27*sqrt(T_table[i])
      *(1.1+0.34*exp(-pow(5.5-log10(T_table[i]),2.0)/3.0));
  
    log10_T = log10(T_table[i]);
    H2_CIE[i]   = pow(10.0, - 116.6 
v		      + 96.3   *log10_T 
		      - 47.15  *log10_T*log10_T 
		      + 10.744 *log10_T*log10_T*log10_T 
		      -  0.9155*log10_T*log10_T*log10_T*log10_T);    
#endif
  }
  
  return (i);
}


#ifdef SHIELD_LW
double SelfShield(double ColumnDens, double T_gas, double b)
{
  double x = ColumnDens/5.e14;
  double b5 = b/1.;  /* Since "b[km/sec]"*/
  
  return 0.965/pow(1.+x/b5, 1.1) + 0.035/pow(1.+x, 0.5)*exp(-8.5e-4*pow(1.+x, 0.5));
}

double EachShield(double ColumnDens, double alpha, double delta, double eta)
{
  double x = ColumnDens/eta;
  
  return exp(-alpha*x)/pow(1.+x, delta);
}
#endif
#endif /* CHEMISTRY */
