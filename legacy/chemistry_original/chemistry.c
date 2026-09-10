/*>=== Includes ===<*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

#include "../../src/standalone/allvars.h"
#include "../../src/standalone/standalone_host.h"
#include "../tags.h"
#include "chemistry.h"

#ifdef CHEMISTRY_SC
#define DEBUG_CHEMISTRY

#ifdef CHEM_HIRANO
#ifndef O05_CHEM
#define O05_CHEM 
#endif
#endif

//#define ID_CHEM_DEBUG 1921863

/*>=== Definitions ===<*/
#define N_THREEBODY 1.0e14
#define N_EQCHEM    1.0e18
#define N_NOLINE    5.0e19

#define MAX_ITERATION 1000

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

/* ----- Energy arrays  */
static double k[k_MAX];
static double k_inverse[k_MAX];

double chem_e[15]; //unit in kJ/mol and mainly from KIDA data base
int FlagPrint=0;

/*=== Subroutines ===<*/
int radiative_cooling_d_nometal(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas);
int radiative_cooling_high_T(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas);
int chemical_cooling(int mode, chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas, double dt, double *y_new, double *y_old);

int chemistry_noneq_exp_nometal(int solve_it, chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas);
int chemistry_noneq_chem_implicit(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas);
int chemistry_noneq_energy_implicit(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas);

#ifdef ENERGY_IMPLICIT
double give_g(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas, double *e_new, double *g);
double solve_bisection(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_min, double g_min, double T_max, double g_max,
		       double *T_new0, double *e_new0, int max_iter);
static double commit_g(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas, double *e_new, double *g);
#endif

void set_reaction_rate_coeff(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas, double *k);
static void update_radiation_state(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas);

static inline int ReactionRates_nometal(int itr, double nH, double *yi, double *rf_tot);

double get_core_radius(double rho_cgs, double T_gas);

#ifdef CHEM_DIFF_ANALYTIC
int (*switch_func[256])(); 

static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *rf_tot);
int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HI(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HII(double nH,double *ni,double *rf);
#ifdef HELIUM
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *rf);
#endif
int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *rf);
#ifdef H2_PLUS
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *rf);
#endif
int ReactionRate_diff_nometal_HM(double nH,double *ni,double *rf);
#ifdef HD
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_DII(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *rf);
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *rf);
#endif
int func_other(double nH,double *ni,double *rf);
static void init_switch_func(void);
#endif

#if defined(CHEMISTRY_IMPLICIT) && !defined(CHEM_DIFF_ANALYTIC) && !defined(ALLOW_CHEM_DIFF_NUMERIC)
#error "Finite-difference chemistry Jacobian is deprecated and very slow; define CHEM_DIFF_ANALYTIC or explicitly opt in with ALLOW_CHEM_DIFF_NUMERIC."
#endif

int chem_tabulate(void);
double compute_gamma_from_abund(double *abund, double T_gas);
static double compute_xmu_from_abund(double *abund);
static double compute_temp_from_abund(double *abund, double gamma, double Utherm);

double (*H2_cooling)();

double H2_cooling_O00(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_Y06(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_GP12(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_G15(double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);

#ifdef HD
double (*HD_cooling)();

double HD_cooling_O00(double nH, double T_gas, double yHI, double yHD);
double HD_cooling_L05(double nH, double T_gas, double yHI, double yHDI);
double HD_cooling_GP12(double nH, double T_gas, double yHI, double yHDI);

#ifdef H2_OPACITY
double (*HD_line_escape_fraction)();
double fesc_HD_fit(double T_K, double xNc_HD);
#endif
#endif

double continuum_cooling_rate_thin(double T_gas, double rho_cgs, double *abund);
double continuum_cooling_rate_thin_Lambda(double T_gas);

#ifdef H2_OPACITY
double (*H2_line_escape_fraction)();

double H2line_reduction(double T);
double fesc_fit(double T_K,double xNc_H2);
double fesc_unity(double T_K, double xNc_H2);
#endif

#ifdef CONT_OPACITY
double (*continuum_escape_fraction)();

double calc_continuum_opacities(double T_gas, double rho_cgs, double l_core, struct opacities *opac);
double calc_continuum_opacities_CIE(double T_gas, double rho_cgs, double l_core, struct opacities *opac);
#endif

void tabulate_opacities_in_LW_band(void);

double (*partition_H2)();

double partition_H2_Rossi1985(double T_gas);
double partition_H2_Irwin1987(double T_gas);
double partition_H2_EBM1989(double T_K);
double partition_H2_EBM1989_fitted_by_Matsukoba(double T_K);
void check_partition_functions_H2(void);

double pf_z_H(double T_K);
double partition_H(double T_gas); 

double (*partition_H2p)();

double (*OxygenI_cooling)();
double (*CarbonI_cooling)();
double (*CarbonII_cooling)();

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

static double compute_temp_from_abund(double *abund, double gamma, double Utherm)
{
  double molecular_weight = compute_xmu_from_abund(abund);
  double egyspec_to_temp  = (gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;

  return (All.UnitEnergy_in_cgs/All.UnitMass_in_g * Utherm) * egyspec_to_temp;
}

double compute_gamma_from_abund(double *abund, double T_gas)
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

static double compute_xmu_from_abund(double *abund){
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

#ifdef ALLVARS_H
double compute_xmu(struct sph_particle_data sph_p)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    abund[j] = sph_p.Abund[j];

  return compute_xmu_from_abund(abund);
}

double compute_temp(int ithis)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    abund[j] = SphP[ithis].Abund[j];

  return compute_temp_from_abund(abund, SphP[ithis].Gamma, SphP[ithis].Utherm);
}
#endif

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

static void update_radiation_state(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas)
{
  double rho_cgs = sph_in->rho_cgs;
  double n_nuclei = rho_cgs / PROTONMASS;

  sph_state->Jeans_Length = get_core_radius(rho_cgs, T_gas);

#ifdef SHIELD_LW
  double l_core = sph_state->Jeans_Length;
  double ColumnDens_H2 = l_core * n_nuclei * sph_out->Abund[sH2I];
  double ColumnDens_HI = l_core * n_nuclei * sph_out->Abund[sHI];
#ifdef HD
  double ColumnDens_HD = l_core * n_nuclei * sph_out->Abund[sHDI];
#endif

  sph_state->f_byH2self = SelfShield(ColumnDens_H2, T_gas, 9.12 * sqrt(T_gas / 1.e4));
  sph_state->f_byHI = EachShield(ColumnDens_HI, 1.49e-1, 1.62e+0, 2.85e23);
#ifdef HD
  sph_state->f_byHDself = SelfShield(ColumnDens_HD, T_gas, 7.45 * sqrt(T_gas / 1.e4));
#endif
  sph_state->f_byH2 = EachShield(ColumnDens_H2, 5.20e-3, 2.38e-1, 2.34e19);
#endif
}

int compute_abundances_d(int mode, chemdata_in *sph_in, chemdata_out *sph_out, double dtime)
{
  int ifunc;
  int status = 0;
  chemdata_state state = {0};
  chemdata_state *sph_state = &state;

  sph_state->t_cool_rad  = MAX_REAL_NUMBER;
  sph_state->t_cool_chem = MAX_REAL_NUMBER;
  FlagPrint = sph_in->flag_print;

  double dt_tot = dtime * All.UnitTime_in_s/All.HubbleParam;

  /* here, mean molecular weight must be taken into account */
  double rho_cgs      = sph_in->rho_cgs;
  double n_nuclei     = rho_cgs/PROTONMASS;
  double molecular_weight = compute_xmu_from_abund(sph_out->Abund);
  
  /* unit conversion */
  double egyspec_to_temp  = All.UnitEnergy_in_cgs/All.UnitMass_in_g * (sph_out->Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
  double T_gas            = sph_out->Utherm * egyspec_to_temp;  

#ifdef IONIZATION_PHOTON
  sph_state->k22 = sph_in->k22;
  sph_state->piHI = sph_in->piHI;
#endif

  if(FlagPrint) 
    printf("[Task%d, ID%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , ThisTask, __LINE__, sph_in->ID, n_nuclei, T_gas
	   , sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I], sph_out->Abund[sHM], dt_tot);

  if(T_gas > 1.e5){
    if(count_print == 0 && ThisTask == 0){
      count_print++;
    }

    if(FlagPrint) 
      printf("[Task%d, ID%d] n=%g T_gas=%g yi=%g %g %g %g %g t_cool=%g t_chem=%g\n"
	     , ThisTask, sph_in->ID, n_nuclei, T_gas
	     , sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I], sph_out->Abund[sHM], sph_out->t_cool, sph_out->t_chem);

    sph_out->Abund[sHI]  = 1.e-10 * XH;
    sph_out->Abund[sH2I] = 1.e-10 * XH;
    sph_out->Abund[sHM]  = 1.e-12 * XH;
    sph_out->Abund[sHII] = 1. * XH;
    sph_out->Abund[sELEC] = 1. * XH;
    sph_out->Gamma = 5./3.;

    double dt = dt_tot;
    double T_gas_init = T_gas;

    while(dt > 0){      
      ifunc = radiative_cooling_high_T(sph_in, sph_out, sph_state, T_gas);

      double dt0 = dt;
      if(0.1 * fabs(sph_state->t_cool_rad) < dt)
	dt0 = 0.1 * fabs(sph_state->t_cool_rad);

      dt -= dt0;
      
      double del_rhoU_rad  = sph_state->rhoU_dot_rad * dt0;
      double del_U_tot     = del_rhoU_rad / rho_cgs; /* delta U in erg/gram */
      del_U_tot    /= (All.UnitEnergy_in_cgs/All.UnitMass_in_g);

      if(sph_out->Utherm + del_U_tot < 0.)
        printf("Why this happen? thermal energy is negative...\n");

      sph_out->Utherm += del_U_tot;

      double molecular_weight0 = compute_xmu_from_abund(sph_out->Abund);
      double egyspec_to_temp0  = All.UnitEnergy_in_cgs/All.UnitMass_in_g * (sph_out->Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight0;
      T_gas            = sph_out->Utherm * egyspec_to_temp0;

      if(FlagPrint)
	printf("ID=%u dt=%g dt0=%g t_cool=%g rhoU_dot_rad=%g del_U_tot=%g Utherm=%g T_gas_new=%g mw=%g %g egsp=%g %g t_cool=%g t_chem=%g\n"
	       , sph_in->ID, dt, dt0, fabs(sph_state->t_cool_rad), sph_state->rhoU_dot_rad, del_U_tot, sph_out->Utherm, T_gas, molecular_weight, molecular_weight0, egyspec_to_temp, egyspec_to_temp0, sph_out->t_cool, sph_out->t_chem);

      if(T_gas < 1.e5)
        break;
    }

    if(FlagPrint)
      printf("COOL: dt_tot=%g->dt=%g t_cool=%g (T=%g->%g n=%g) t_cool=%g t_chem=%g\n"
	     , dt_tot, dt, sph_state->t_cool_rad, T_gas_init, T_gas, rho_cgs/PROTONMASS, sph_out->t_cool, sph_out->t_chem);    

	    if(dt < 1.e-100){
	      sph_out->t_cool = sph_state->t_cool_rad;
	      goto LOOP_OUT;
	    }

	    dt_tot = dt;
	    dtime = dt / All.UnitTime_in_s * All.HubbleParam;
	  }

	  for(int j=22; j<=29; j++)  k[j] = 0;

#ifdef IONIZATION_PHOTON
	  double k2 = 2.6e-13;
	  double x_neut_eq = k2 * XH * n_nuclei/ (sph_in->k22 + 1.e-300);
	  double x_neut_crit = 1.e-3;
	  double Lambda = continuum_cooling_rate_thin_Lambda(1.e4);

  if(FlagPrint && sph_in->k22 > 0.){
    printf("x_neut_eq=%g k22*dt_tot=%g k2/k22piHI=%g Lambda=%g edot=%g %g\n"
	   , x_neut_eq, sph_in->k22 * dt_tot
	   , k2/sph_in->k22*sph_in->piHI, Lambda, sph_in->piHI*sph_in->k22, (k2/sph_in->k22*sph_in->piHI-Lambda)*n_nuclei*XH*dt_tot);
  }

  if(sph_in->k22 * dt_tot > 1. && x_neut_eq < x_neut_crit){
    double fac = x_neut_eq / x_neut_crit;
    sph_state->k22 = fac * sph_in->k22;
    sph_state->piHI = fac * sph_in->piHI;
  }
  //piHI and k[22] will be given in set_reaction_rate_coeff()

  k[23] = k[24] = 0.;
#endif
  
  update_radiation_state(sph_in, sph_out, sph_state, T_gas);


  {
    /* recalculate the t_chem & t_cool for any particles */
    /* just compute the timescales and return */
#if !defined (ENERGY_IMPLICIT) || defined (HYBRID_IMPLICIT)
    double dt_in = 0.0;
    double T_gas_in = T_gas;

#if defined (CHEMISTRY_IMPLICIT) && !defined (HYBRID_IMPLICIT)
    sph_out->t_chem = MAX_REAL_NUMBER;
    set_reaction_rate_coeff(sph_in, sph_out, sph_state, T_gas_in, k);
#else
    ifunc = chemistry_noneq_exp_nometal(0, sph_in, sph_out, sph_state, dt_in, T_gas_in);
#endif
    ifunc = radiative_cooling_d_nometal(sph_in, sph_out, sph_state, T_gas_in);
    sph_out->t_cool = dmin(sph_state->t_cool_rad, sph_state->t_cool_chem);
#else
    sph_out->t_cool = MAX_REAL_NUMBER;
    sph_out->t_chem = MAX_REAL_NUMBER;      
#endif

#ifdef FIX_TEMPERATURE
    sph_out->t_cool = MAX_REAL_NUMBER;    
#endif

    if(FlagPrint)
      printf("(line%d) t_cool=%g t_chem=%g\n", __LINE__, sph_out->t_cool, sph_out->t_chem);
    
    if(mode==0 || mode ==2)
      return 0;
  }

  if(FlagPrint) 
    printf("[Task%d, ID%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , ThisTask, sph_in->ID, __LINE__, n_nuclei, T_gas
	   , sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I], sph_out->Abund[sHM], dt_tot);

  double t_now = 0.;    /*    set current time (redshift) */
    
  for(int iter=0;iter<MAX_ITERATION;iter++){
    /*---------  Set the internal time step for subcycling */
    double time1 = All.Epsilon * fabs(sph_out->t_cool);
    double time2 = All.Epsilon * fabs(sph_out->t_chem);
    double time_elapsed = t_now * All.UnitTime_in_s / All.HubbleParam;
    double time3 = dt_tot - time_elapsed;

    if(time3 <= 0.)
      goto LOOP_OUT;

    /* take the minimum */ 
    double dt = dmin(time1, dmin(time2,time3));

    if(iter >= MAX_ITERATION - 10){
      printf("[iter%d]time=%g %g %g => dt=%g\n",iter, time1, time2, time3, dt);
      printf("t_now=%g,dtime=%g\n",t_now, dtime);
      printf("n_nuclei=%g,Tgas=%g\n", n_nuclei, T_gas);
    }            

    int flag_implicit = 0, flag_chem_implicit = 0;
#if defined(ENERGY_IMPLICIT) && !defined(HYBRID_IMPLICIT)
    dt = time3;
    flag_implicit=1;
#endif

#if defined(CHEMISTRY_IMPLICIT) && !defined(HYBRID_IMPLICIT)
    dt = dmin(time1,time3); 
    flag_chem_implicit=1;
#endif

#ifdef HYBRID_IMPLICIT
    dt = time3;

    if(time3 > time1 || n_nuclei > N_EQCHEM)
      flag_implicit = 1;
    else if(time3 > time2){
      flag_implicit = 0;
      flag_chem_implicit = 1;
    }else{
      if(n_nuclei > 1.e9){
	if(time2 > time1) flag_implicit = 1;
	else flag_implicit = 0, flag_chem_implicit = 1;
      }else{
	flag_implicit = 0;
	flag_chem_implicit = 0;
      }
    }
#endif
  
    sph_out->flag_implicit = flag_implicit + 2*flag_chem_implicit;
    
    if(FlagPrint)
      printf("[Task%d, ID%d] line=%d iter=%d flag_implicit=%d flag_chem_implicit=%d\n"
	     , ThisTask,  sph_in->ID, __LINE__, iter, flag_implicit, flag_chem_implicit);

    double T_gas_in = T_gas;
    
#if defined(ENERGY_IMPLICIT)
    if(flag_implicit){
      sph_out->conv_flag = chemistry_noneq_energy_implicit(sph_in, sph_out, sph_state, dt, T_gas_in);
      if(sph_out->conv_flag)
        status = 1;
    }else
#endif
      {
        chemdata_out sph_old = *sph_out;

	int isnan_flag=0;

#if defined(CHEMISTRY_IMPLICIT)
	if(flag_chem_implicit){
	  ifunc = chemistry_noneq_chem_implicit(sph_in, sph_out, sph_state, dt, T_gas_in);
	  if(ifunc != 1)
	    status = 1;
	}else
#endif
	  {
	    ifunc = chemistry_noneq_exp_nometal(1, sph_in, sph_out, sph_state, dt, T_gas_in);
	    if(ifunc == 2){
	      isnan_flag = 1;
	      status = 1;
	      FlagPrint = 1;  //we allow to display debug info for this erroneous case
	    }
	  }

	sph_out->Gamma = compute_gamma_from_abund(sph_out->Abund, T_gas_in);
	
	/* Update temperature */ 
	int cooling_status = radiative_cooling_d_nometal(sph_in, sph_out, sph_state, T_gas_in);
	if(cooling_status)
	  status = 1;

	sph_out->t_cool  = dmin(fabs(sph_state->t_cool_rad),fabs(sph_state->t_cool_chem));
		
	double del_rhoU_rad  = sph_state->rhoU_dot_rad * dt; 		
	double del_rhoU_chem = sph_state->rhoU_dot_chem * dt; /* erg/cm3 */
	double del_rhoU_tot  = del_rhoU_rad + del_rhoU_chem;		
	double del_U_tot     = del_rhoU_tot / rho_cgs; /* delta U in erg/gram */		
	del_U_tot    /= (All.UnitEnergy_in_cgs/All.UnitMass_in_g);
	
	sph_out->Utherm += del_U_tot;		
	if(FlagPrint)
	  printf("(line%d) t_cool=%g t_chem=%g\n", __LINE__, sph_out->t_cool, sph_out->t_chem);

#if defined(HYBRID_IMPLICIT) && defined(ENERGY_IMPLICIT)
	        /************* DEBUGGING... [CHON200417] ************/
		if(0.1 * sph_old.Utherm < fabs(del_U_tot) || isnan_flag){
	          flag_implicit = 1;
	          flag_chem_implicit = 0;
	          *sph_out = sph_old;

	          sph_out->conv_flag = chemistry_noneq_energy_implicit(sph_in, sph_out, sph_state, dt, T_gas_in);
	          if(sph_out->conv_flag)
	            status = 1;
	          sph_out->flag_implicit = flag_implicit + 2*flag_chem_implicit;
	        }
#endif
	      }

    dt /= All.UnitTime_in_s; /* system units */
    dt *= All.HubbleParam;   /* h correction */      
    t_now += dt;      
            
    /* update T_gas */
    molecular_weight = compute_xmu_from_abund(sph_out->Abund);
    egyspec_to_temp  = All.UnitEnergy_in_cgs/All.UnitMass_in_g * (sph_out->Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
    T_gas            = sph_out->Utherm * egyspec_to_temp;
           
#ifdef SET_TMAX
    if(T_gas > SET_TMAX)
      sph_out->Utherm *= SET_TMAX / T_gas;
#endif

#ifdef FIX_TEMPERATURE
    sph_out->Utherm *= All.InitTemperature / T_gas;
    sph_out->t_cool = MAX_REAL_NUMBER;    
#endif

    /* stop if time = a_end */
    if(dtime - t_now < dtime/MAX_ITERATION)
      goto LOOP_OUT;      
  }  /* iter */
   
  
 LOOP_OUT:;

  if(FlagPrint) 
    printf("[Task%d, ID%d] line=%d n=%g T_gas=%g yi=%g %g %g %g %g dt=%g\n"
	   , ThisTask, sph_in->ID, __LINE__, n_nuclei, T_gas
	   , sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I], sph_out->Abund[sHM], dt_tot);
  
#ifdef CHEMISTRY_IMPLICIT
  sph_out->t_chem = MAX_REAL_NUMBER;
#else
  sph_out->t_chem = fabs(sph_out->t_chem) / All.UnitTime_in_s * All.HubbleParam;  
#endif
  
#ifdef ENERGY_IMPLICIT
  sph_out->t_cool = MAX_REAL_NUMBER;
#else
  sph_out->t_cool = fabs(sph_out->t_cool) / All.UnitTime_in_s * All.HubbleParam;  
#endif

  return status;
}



/* ----- Compute cooling time of gas and
   ----- solve for thermal energy.  */
int radiative_cooling_d_nometal(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas)
{  
  volatile int i_T;

  double ln_T0 = log(T_start);
  double ln_T9 = log(T_end);
  
  /* compute total numberdensity in units of mh
     hence n_nuclei * HI gives the number density of HI  */

  double rho_gas = sph_in->rho_cgs;
  double n_nuclei = rho_gas/PROTONMASS;  
  
  /* --------- Compute temperature index.  */
  double ln_T = log(T_gas);
  if (ln_T < ln_T0) ln_T = ln_T0;
  if (ln_T > ln_T9) ln_T = ln_T9;
  
  i_T = (int)( (ln_T-ln_T0)/dln_T ) ;
  
  if (i_T > N_T-2) i_T = N_T - 2;     /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  double c_T = (T_gas - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);
    
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
  double comp1 = 5.65e-36*pow(All.cf_atime,-4);    
  double edot_comp = -comp1*(T_gas-All.T_CMB)*sph_out->Abund[sELEC]*n_nuclei;
  
  /**********************************************************************************************/ 
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/ 
  /**********************************************************************************************/ 

  double l_core = sph_state->Jeans_Length;
  
  /******************************* Hydrogen lines ***********************************************/
  double edot_ceHI = 0.;
#ifdef HLINE_COOLING_THIN
  edot_ceHI   = -ceHI  *sph_out->Abund[sHI]  *sph_out->Abund[sELEC]     *n_nuclei*n_nuclei;
#else
  if(n_nuclei < 1.e8){
    double edot_Lya, edot_2ph;
    double v_thermal = sqrt(BOLTZMANN * T_gas/PROTONMASS);
    double fesc = fesc_Lya(n_nuclei, T_gas, sph_out->Abund[sHI], sph_out->Abund[sELEC], l_core, v_thermal, &edot_Lya, &edot_2ph);
    
    edot_ceHI = edot_Lya;
  }
#endif

  /******************************* Helium lines **************************************************/
#ifdef HELIUM
  double edot_ceHeI  = -ceHeI *sph_out->Abund[sHeII]*sph_out->Abund[sELEC]*sph_out->Abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
  double edot_ceHeII = -ceHeII*sph_out->Abund[sHeII]*sph_out->Abund[sELEC]     *n_nuclei*n_nuclei;

  double edot_ciHeI  = -ciHeI *sph_out->Abund[sHeI] *sph_out->Abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeII = -ciHeII*sph_out->Abund[sHeII]*sph_out->Abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeIS  = -ciHeIS *sph_out->Abund[sHeII]*sph_out->Abund[sELEC]*sph_out->Abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
#elif defined(HELIUM_EQ_COOLING)
  double denom = k[4]*k[6] + k[3]*k[6] + k[3]*k[5];
  if(denom < 1.e-100)
    denom = 1.e-100;
  
  double fHeI_eq   = XHe * k[4]*k[6] / denom;
  double fHeII_eq  = XHe * k[3]*k[6] / denom;
  double fHeIII_eq = XHe * k[3]*k[5] / denom;
  double felec_new = sph_out->Abund[sELEC] + fHeII_eq + 2 * fHeIII_eq;

  double edot_ceHeI  = -ceHeI * fHeI_eq  * felec_new * n_nuclei * n_nuclei;
  double edot_ceHeII = -ceHeII* fHeII_eq * felec_new * n_nuclei * n_nuclei;
#endif      


  /******************************* H2 & HD lines ***********************************************/  
  double nH  = n_nuclei * XH;
  double yHI = sph_out->Abund[sHI]    / XH;
  double yHII = sph_out->Abund[sHII]  / XH;
  double yH2I = sph_out->Abund[sH2I]  / XH;
  double y_e  = sph_out->Abund[sELEC] / XH;  
#ifdef HELIUM
  double yHeI = sph_out->Abund[sHeI]  / XH;  
#else
  double yHeI = XHe / XH;
#endif
#ifdef HD
  double yHDI = sph_out->Abund[sHDI]/XH;
#endif

  double edot_H2 = H2_cooling(nH, T_gas, yHI, yH2I, yHeI, yHII, y_e);
#ifdef HD  
  double edot_HD = HD_cooling(nH, T_gas, yHI, yHDI);
#endif


#ifdef H2_OPACITY  /* accounting for the local optical depth to H2 lines   */ 
  if(n_nuclei > 1.0e4){
    double N_H2 = sph_state->Jeans_Length * n_nuclei * sph_out->Abund[sH2I];
    double f_reduce_H2 = H2_line_escape_fraction(T_gas, N_H2);
    sph_state->EscapeProb = f_reduce_H2;
    
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
    double N_HD = sph_state->Jeans_Length * n_nuclei * sph_out->Abund[sHDI];
    double f_reduce_HD = HD_line_escape_fraction(T_gas, N_HD);
    sph_state->EscapeProb_HD = f_reduce_HD;
    
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
  double f_metal = sph_in->Metallicity;
  if(FlagPrint)
    printf("f_metal=%g\n", f_metal);
  
  double nOI  = n_nuclei * XH * Y_O_sun * f_metal;  
  double nCI  = n_nuclei * XH * Y_C_sun * f_metal;  
  double nCII = n_nuclei * XH * Y_C_sun * f_metal;  
  double ne   = sph_out->Abund[sELEC] * n_nuclei;

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
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, sph_out->Abund);
#else
    /* --------- Cooling by CIE  */
  /* --------- Bremsstrahlung  */
  /* --------- Recombinations  */
  double edot_reHI    = -reHII*  sph_out->Abund[sHII] * sph_out->Abund[sELEC] * n_nuclei*n_nuclei;
  double edot_CIE     = -cie  *  sph_out->Abund[sH2I] * sph_out->Abund[sH2I]  * n_nuclei*n_nuclei;
  double edot_brem    = -brem * sph_out->Abund[sHII]  *  sph_out->Abund[sELEC]*n_nuclei*n_nuclei;
  
#ifdef HELIUM
  edot_brem    -= brem*(sph_out->Abund[sHeII]+4.0*sph_out->Abund[sHeIII]) * sph_out->Abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII1 = -reHeII1*sph_out->Abund[sHeII] *sph_out->Abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII2 = -reHeII2*sph_out->Abund[sHeII] *sph_out->Abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeIII = -reHeIII*sph_out->Abund[sHeIII]*sph_out->Abund[sELEC]*n_nuclei*n_nuclei;
#endif  
#endif  

  /********************** opacity for the continuum process *******************************/
#ifdef CONT_OPACITY
  double tau_cont = sph_state->tau_cont;
  double f_reduce = sph_state->f_reduce;
#endif

  double edot_dust = 0.;
#ifdef DUST_COOLING
  edot_dust = sph_state->dust_cooling_rate;
#endif
#if defined(DUST_COOLING) && defined(LYMANWERNER)
  double edot_PE = 0.;
  if(sph_out->Tdust < 1000.){
#ifdef LW_RADIATION_LOCAL
    double G0 = sph_state->G0;
#else
    double G0 = All.G0;
#endif
    edot_PE = dust_PEheating(G0, nH, y_e, T_gas, sph_in->Metallicity);
    if(FlagPrint)
      printf("edot_dust=%g cont=%g T_dust=%g\n", edot_dust, edot_cont, sph_out->Tdust);
  }
#endif
  
#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * sph_in->Metallicity * n_nuclei * n_nuclei;
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

  edot_cont *= f_reduce;
#endif  

	       
  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = sph_state->piHI  *sph_out->Abund[sHI]  *n_nuclei;
#ifdef HELIUM
  double edot_radHeI  = sph_state->piHeI *sph_out->Abund[sHeI] *n_nuclei;
  double edot_radHeII = sph_state->piHeII*sph_out->Abund[sHeII]*n_nuclei;
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
	   edot, edot_ceHI, edot_H2, edot_HD, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont, f_reduce, sph_out->Tdust, T_gas, edot_radHI);
#else
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g cont=%g f_reduce=%g T_gas=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_cont, f_reduce, T_gas);
#endif
  }

#ifdef ID_CHEM_DEBUG
  if(FlagPrint)
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont, f_reduce, sph_out->Tdust);
#else
#if defined(THERMAL_FEEDBACK) && defined(METALLICITY)
  if(FlagPrint)
    printf("cool_rate:T_gas=%g n=%g ID=%d tot=%g rad=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g f_line=%g tau=%g Tdust=%g Trad=%g\n",
	   T_gas, n_nuclei, sph_in->ID, edot, edot_radHI, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust, edot_cont, f_reduce, sph_state->EscapeProb, tau_cont,
	   sph_out->Tdust, pow(sph_in->stellar_rad_flux / STEFAN / 4., 0.25));
#endif
#endif
  
#if defined(DEBUG_CHEMISTRY_SC)
  sph_state->edot_ceHI = edot_ceHI;	       
  sph_state->edot_comp = edot_comp;
  sph_state->edot_H2 = edot_H2;

#ifdef HELIUM
  sph_state->edot_ceHeI = edot_ceHeI;
  sph_state->edot_ceHeII = edot_ceHeII;
  sph_state->edot_ciHeI = edot_ciHeI;
  sph_state->edot_ciHeII = edot_ciHeII;
  sph_state->edot_ciHeIS = edot_ciHeIS;

  sph_state->edot_reHeII1 = edot_reHeII1;
  sph_state->edot_reHeII2 = edot_reHeII2;
  sph_state->edot_reHeIII = edot_reHeIII;
#endif

#ifdef CONTINUUM_COOLING_MATSUKOBA
  sph_state->edot_cont = edot_cont;
#else
  sph_state->edot_reHI = edot_reHI;
  sph_state->edot_brem = edot_brem;
  sph_state->edot_CIE = edot_CIE;
#endif
  
#ifdef HD
  sph_state->edot_HD = edot_HD;
#endif
#ifdef METALLICITY
  sph_state->edot_CI  = edot_CI;
  sph_state->edot_CII = edot_CII;
  sph_state->edot_OI  = edot_OI;
#endif
#ifdef DUST_COOLING
  sph_state->edot_dust = edot_dust;
#ifdef LYMANWERNER
  sph_state->edot_PE = edot_PE;
#endif		
#endif		
#endif		
    
  sph_state->rhoU_dot_rad = edot;
  
  double rhoU = sph_out->Utherm * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.
  rhoU *= (All.UnitEnergy_in_cgs/pow(All.UnitLength_in_cm, 3)); /* in erg/cm3 */
  
  sph_state->t_cool_rad = (rhoU/edot); /* in seconds */  
  if(fabs(edot) < 1.0e-40)
    sph_state->t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}


int radiative_cooling_high_T(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas)
{    
  double rho_gas = sph_in->rho_cgs;
  double n_nuclei  = rho_gas/PROTONMASS;  
  /**********************************************************************************************/ 
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/ 
  /**********************************************************************************************/ 

  sph_state->Jeans_Length = get_core_radius(rho_gas, T_gas);
  double l_core = sph_state->Jeans_Length;

  /******************************************************************************************/ 
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/ 
  /******************************************************************************************/ 

#ifdef CONTINUUM_COOLING_MATSUKOBA  
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, sph_out->Abund);
#endif  
  double edot_reHI = -8.418e-26*pow(T_gas,0.11) * n_nuclei*n_nuclei;
  double edot_brem = -1.43e-27*sqrt(T_gas) * (1. + 4.*XHe/XH) * n_nuclei * n_nuclei;

  if(FlagPrint)
    printf("T_gas=%g edot_reHI=%g brem=%g cont=%g\n", T_gas, edot_reHI, edot_brem, edot_cont);

  /********************** opacity for the continuum process *******************************/

#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * sph_in->Metallicity * n_nuclei * n_nuclei;
#endif

  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = sph_state->piHI  *sph_out->Abund[sHI] *n_nuclei;
    
  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  double edot = edot_reHI + edot_brem 
#ifdef METAL_LINE_COOLING_SD93
    + edot_metal_SD93
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA  
    + edot_cont
#endif
    + edot_radHI;
  
  sph_state->rhoU_dot_rad = edot;
  
  double rhoU = sph_out->Utherm * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.
  
  rhoU *= (All.UnitEnergy_in_cgs/pow(All.UnitLength_in_cm, 3)); /* in erg/cm3 */
  
  sph_state->t_cool_rad = (rhoU/edot); /* in seconds */

  if(FlagPrint){
#ifdef METAL_LINE_COOLING_SD93
    printf("cool_rate:tot=%g reHI=%g brem=%g cont_metal=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, edot_metal_SD93, T_gas, rhoU, sph_state->t_cool_rad);
#else
    printf("cool_rate:tot=%g reHI=%g brem=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, T_gas, rhoU, sph_state->t_cool_rad);
#endif
  }
  
  if(fabs(edot) < 1.0e-40)
    sph_state->t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}

/* ----- Compute cooling time of gas and
   ----- solve for thermal energy.  */

double radiative_cooling_d_nometal_test(double nH, double T_gas, double metallicity, double y_e, double yH2I, double yHDI, double *edot)
{  
  double n_nuclei = nH / XH;
  
  double yHI  = 1. - y_e - yH2I;
  double yHII = y_e;
  double yHeI = XHe / XH;
  
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  double ln_T0 = log(T_start);
  double ln_T9 = log(T_end);
  double ln_T = log(T_gas);
  if (ln_T < ln_T0) ln_T = ln_T0;
  if (ln_T > ln_T9) ln_T = ln_T9;
  
  volatile int i_T = (int)( (ln_T-ln_T0)/dln_T ) ;
  
  if (i_T > N_T-2) i_T = N_T - 2;     /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  double c_T = (dmin(T_gas, T_end) - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);  
  double ceHeI  = ceHeIa[i_T]  + c_T*(ceHeIa[i_T+1] - ceHeIa[i_T]);
  double ceHeII = ceHeIIa[i_T] + c_T*(ceHeIIa[i_T+1]- ceHeIIa[i_T]);

  double reHeII1= reHeII1a[i_T] + c_T*(reHeII1a[i_T+1]-reHeII1a[i_T]);
  double reHeII2= reHeII2a[i_T] + c_T*(reHeII2a[i_T+1]-reHeII2a[i_T]);
  double reHeIII= reHeIIIa[i_T] + c_T*(reHeIIIa[i_T+1]-reHeIIIa[i_T]);
  
  k[1] = k_table[1][i_T] + c_T*(k_table[1][i_T+1] - k_table[1][i_T]);
  k[2] = k_table[2][i_T] + c_T*(k_table[2][i_T+1] - k_table[2][i_T]);
  k[3] = k_table[3][i_T] + c_T*(k_table[3][i_T+1] - k_table[3][i_T]);
  k[4] = k_table[4][i_T] + c_T*(k_table[4][i_T+1] - k_table[4][i_T]);
  k[5] = k_table[5][i_T] + c_T*(k_table[5][i_T+1] - k_table[5][i_T]);
  k[6] = k_table[6][i_T] + c_T*(k_table[6][i_T+1] - k_table[6][i_T]);  
#endif 

  double fHI_eq = XH * k[2]/ (k[1]+k[2]);
  double fHII_eq = XH * k[1]/ (k[1]+k[2]);
  
  if(T_gas > 1.e4){
    yHI = fHI_eq / XH;
    yHII = fHII_eq / XH;
    y_e = fHII_eq / XH;
    yH2I = 1.e-10;
  }

  double felec_new = y_e * XH;
#if defined(HELIUM_EQ_COOLING)
  double denom = k[4]*k[6] + k[3]*k[6] + k[3]*k[5];
  if(denom < 1.e-100)
    denom = 1.e-100;

  double fHeI_eq   = XHe * k[4]*k[6] / denom;
  double fHeII_eq  = XHe * k[3]*k[6] / denom;
  double fHeIII_eq = XHe * k[3]*k[5] / denom;
  
  felec_new += fHeII_eq + 2 * fHeIII_eq;
  y_e = felec_new / XH;
#endif
   
  /**********************************************************************************************/ 
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/ 
  /**********************************************************************************************/ 
  double l_core = 0.;

  /******************************* Hydrogen lines ***********************************************/
  double edot_Lya;
  double v_thermal = sqrt(BOLTZMANN * T_gas/PROTONMASS);
  double edot_2ph;
  double fesc = fesc_Lya(n_nuclei, T_gas, XH * yHI, XH * y_e, l_core, v_thermal, &edot_Lya, &edot_2ph);
    
  double edot_ceHI = edot_Lya;

  /******************************* Helium lines **************************************************/
#if defined(HELIUM_EQ_COOLING)
  double edot_ceHeI  = -ceHeI * fHeI_eq * felec_new * n_nuclei * n_nuclei;
  double edot_ceHeII = -ceHeII* fHeII_eq * felec_new * n_nuclei * n_nuclei;
  double edot_reHeII1 = -reHeII1 * fHeII_eq  * felec_new * n_nuclei * n_nuclei;
  double edot_reHeII2 = -reHeII2 * fHeII_eq  * felec_new * n_nuclei * n_nuclei;
  double edot_reHeIII = -reHeIII * fHeIII_eq * felec_new * n_nuclei * n_nuclei;
#endif      

  /******************************* H2 & HD lines ***********************************************/  
  double edot_H2 = H2_cooling(nH, T_gas, yHI, yH2I, yHeI, yHII, y_e);
#ifdef HD  
  double edot_HD = HD_cooling(nH, T_gas, yHI, yHDI);
#endif

  /******************************* CII & OI lines ***********************************************/  
#ifdef METALLICITY
  double f_metal = metallicity;

  double nOI  = n_nuclei * XH * Y_O_sun * f_metal;  
  double nCI  = n_nuclei * XH * Y_C_sun * f_metal;  
  double nCII = n_nuclei * XH * Y_C_sun * f_metal;  
  double ne   = nH * y_e;

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
  double abund[N_SP];
  abund[sELEC] = felec_new;
  abund[sHI] = yHI * XH;
  abund[sHII] = yHII * XH;
  abund[sH2I] = yH2I * XH;
  abund[sHM] = 0.;
  
  double edot_cont = -continuum_cooling_rate_thin(T_gas, nH * PROTONMASS / XH, abund);

  /********************** opacity for the continuum process *******************************/
#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * metallicity * n_nuclei * n_nuclei;
#endif
  
  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  edot[0]  = -edot_ceHI;
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  edot[1]  = -edot_ceHeI - edot_ceHeII;
  edot[2]  = -edot_reHeII1 - edot_reHeII2 - edot_reHeIII;
#endif
  edot[3]  = -edot_H2;
  edot[4]  = -edot_cont;
#ifdef HD
  edot[5]  = -edot_HD;
#endif
#ifdef METALLICITY
  edot[6]  = -edot_CI;
  edot[7]  = -edot_CII;
  edot[8]  = -edot_OI;
#endif
#ifdef METAL_LINE_COOLING_SD93
  edot[9]  = -edot_metal_SD93;
#endif

  double etot = edot_ceHI + edot_H2 + edot_cont
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    + edot_ceHeI + edot_ceHeII + edot_reHeII1 + edot_reHeII2 + edot_reHeIII
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
    ;
  
  return etot; 
}


int chemical_cooling(int mode, chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas, double dt, double *y_new, double *y_old){
  double edot_H2form,edot_H2form_3body;
  double edot_ciHI, edot_ciHeI, edot_ciHeII;
  double edot_HMform;

  enum SPECIES species;
  double dyi0[N_SP];

  double rho_gas = sph_in->rho_cgs;
  double nH = XH * rho_gas / PROTONMASS;
  
  /* calculate chemical heating/cooling */
  double n_elec = nH * y_new[sELEC];
  double n_HI   = nH * y_new[sHI];
  double n_HM   = nH * y_new[sHM];
  double n_H2I  = nH * y_new[sH2I];
#ifdef H2_PLUS
  double n_H2II = nH * y_new[sH2II];
#endif
#ifdef HELIUM
  double n_HeI  = nH * y_new[sHeI];
  double n_HeII = nH * y_new[sHeII];
#endif
#ifdef HD
  double n_DI = nH * y_new[sDI];
  double n_DII = nH * y_new[sDII];
  double n_HDI = nH * y_new[sHDI];
#endif

  dyi0[sHI] = y_new[sHI] - y_old[sHI];
  dyi0[sHII] = y_new[sHII] - y_old[sHII];
#ifdef HELIUM
  dyi0[sHeII]  = y_new[sHeII]  - y_old[sHeII];
  dyi0[sHeIII] = y_new[sHeIII] - y_old[sHeIII];
#endif
  dyi0[sH2I] = y_new[sH2I] - y_old[sH2I];
  dyi0[sHM]  = y_new[sHM]  - y_old[sHM];
  
  double n_crit = (1.e6 / sqrt(T_gas)) 
    / (1.6 * y_new[sHI] * exp(-(400.0/T_gas)*(400.0/T_gas)) + 1.4 * y_new[sH2I] * exp(-12000.0/(T_gas + 1200.0)));

  /********** Mode 0 is just for evaluating the timescale chemical heating / cooling ********************/
  /********** This is only used in the explicit time-stepping or hybrid time-stepping *******************/
  double edot_H2form_HM, edot_H2form_dust, edot_ciHI_EUV;
  double edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV, edot_H2form_3body_loss;
  if(mode == 0){
    edot_H2form_HM = eV_to_erg * ( k[8] * n_HM * n_HI * 3.73/ (1. + n_crit/n_HI)
#ifdef H2_PLUS
				   + k[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/n_HI)) 
#endif
				   ); //1.38?? modified to 1.83

    
    edot_H2form_dust = 0.;
#ifdef H2_FORMATION_ON_DUST
    edot_H2form_dust = eV_to_erg * k[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/n_HI));
#endif

    edot_H2form_3body = 0.;
#ifdef THREEBODY
    edot_H2form_3body = 4.48 * eV_to_erg * ( ( k[55]*n_HI + k[57]*n_H2I ) * n_HI * n_HI
					     - ( k[11]*n_HI + k[58]*n_H2I) * n_H2I );
#endif    
    edot_H2form = edot_H2form_HM + edot_H2form_3body + edot_H2form_dust;

    edot_ciHI   = -13.598 * eV_to_erg * k[1] * n_HI * n_elec; 
#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * k[3] * n_HeI  * n_elec;
    edot_ciHeII = -54.416 * eV_to_erg * k[5] * n_HeII * n_elec;
#endif

    edot_HMform = 0.;
  }else{
    if(nH < N_THREEBODY){
      edot_H2form_HM_loss = eV_to_erg * ( k[8]  * n_HM   * n_HI * 3.73 / (1. + nH/n_crit)
#ifdef H2_PLUS
					  + k[10] * n_H2II * n_HI * 1.83 / (1. + nH/n_crit)
#endif
					  );
 
      edot_H2form_dust_loss = 0.;
#ifdef H2_FORMATION_ON_DUST
      edot_H2form_dust_loss = eV_to_erg * k[51] * n_HI * 4.28 / (1. + nH/n_crit);
#endif
      edot_H2diss_FUV = eV_to_erg * k[29] * n_H2I * 4.48;
        
      edot_H2form_3body = 4.48 * eV_to_erg * dyi0[sH2I] * nH / dt;
      edot_H2form_3body_loss = 4.48 * eV_to_erg * ( ( k[55]*n_HI + k[57]*n_H2I ) * n_HI * n_HI ) / (1. + nH/n_crit);


      edot_H2form = edot_H2form_3body - edot_H2form_HM_loss - edot_H2form_dust_loss + edot_H2diss_FUV;    

      edot_H2form_HM = eV_to_erg * ( k[8] * n_HM * n_HI * 3.73/ (1. + n_crit/nH)
#ifdef H2_PLUS
				   + k[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/nH)) 
#endif
				   ); //1.38?? modified to 1.83
      edot_H2form_dust = eV_to_erg * k[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/nH));
      edot_H2form_3body = 4.48 * eV_to_erg * ( ( k[55]*n_HI + k[57]*n_H2I ) * n_HI * n_HI
					       - ( k[11]*n_HI + k[58]*n_H2I) * n_H2I );      
      edot_HMform = 0.78 * eV_to_erg * dyi0[sHM]  * nH / dt;
    }else{
      edot_H2form = edot_H2form_3body = 4.48 * eV_to_erg * dyi0[sH2I] * nH / dt;
      edot_HMform = 0.78 * eV_to_erg * dyi0[sHM]  * nH / dt;
    }
  
    edot_ciHI   = -13.598 * eV_to_erg * dyi0[sHII]    * nH / dt;
#ifdef IONIZATION_PHOTON
    edot_ciHI_EUV = 13.598 * eV_to_erg * k[22] * n_HI;
    edot_ciHI  += edot_ciHI_EUV;
#endif

#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * dyi0[sHeII]   * nH / dt;
    edot_ciHeII = -54.416 * eV_to_erg * dyi0[sHeIII]  * nH / dt;
#endif
  }

  if(FlagPrint){
    printf("edot_chem: nH=%g H2=%g 3b=%g 3b_l=%g HM_l=%g dust_l=%g FUV=%g ciHI=%g EUV=%g HM=%g ncrit=%g dy_H2=%g\n"
	   , nH, edot_H2form, edot_H2form_3body, edot_H2form_3body_loss, edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV
	   , edot_ciHI, edot_ciHI_EUV, edot_HMform, n_crit, dyi0[sH2I]);
  }
  
  sph_state->rhoU_dot_chem = edot_H2form + edot_ciHI + edot_HMform;  
  if(dt < 1.e-30){
    sph_state->rhoU_dot_chem = 0.;
    edot_H2form = edot_ciHI = edot_HMform = 0.;
  }

  double rhoU = sph_out->Utherm * sph_in->rho_cgs; 
  //modified[0814] in case for SphP[ithis].Entropy is not updated.

  rhoU *= (All.UnitEnergy_in_cgs/pow(All.UnitLength_in_cm, 3)); /* in erg/cm3 */
  
  sph_state->t_cool_chem = (rhoU/sph_state->rhoU_dot_chem); /* in seconds */

  return 0;
}


#ifdef ENERGY_IMPLICIT
#define MAX_ITER 500
#define BETA_ENERGY 1.0

int chemistry_noneq_energy_implicit(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas){
  double rho_gas = sph_in->rho_cgs;
  
  int iter = 0, ITER_reset = 30;
  double denom=0., err, err_chem;
  double e_new;
  double T_gas_init, T_old, T_new, dT, d_T;
  const double T_max=3.e6, T_min=0.1;
  double T_large, T_small;
  double g = 0., g0, array[1000];

  chemdata_out sph_old = *sph_out;
  chemdata_state sph_state_old = *sph_state;
  T_gas_init =  T_gas;

  T_large = 5000, T_small = 1.5;
  
  T_old = T_gas_init, dT = 1.0e-3 * T_old;
  T_new = T_old;

#ifdef IONIZATION_PHOTON
  double nH = XH * rho_gas / PROTONMASS;
  double k2 = 2.6e-13;
  double x_neut_eq = MAX_REAL_NUMBER;
  if(sph_state->k22 > 0.)
    x_neut_eq = k2 * nH / sph_state->k22;
  double T_temp = sph_state->piHI * dt/ BOLTZMANN;

  if(x_neut_eq < 0.1){
    if(sph_state->k22 > 0. && sph_state->k22 * dt > 0.2)
      T_temp = dmin(T_temp, sph_state->piHI / sph_state->k22 / BOLTZMANN);
  }else{
    if(T_temp > 2.e4){
      double Lambda = continuum_cooling_rate_thin_Lambda(2.e4);
      double fneut_crit = sph_state->piHI / 2./ Lambda/ nH;
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
    printf("[2]piHI = %g, %g T_temp=%g k22=%g x_neut_eq=%g dt=%g\n", sph_state->piHI, sph_in->piHI, T_temp, sph_state->k22,x_neut_eq, dt);
  
  if(T_temp > T_gas){
    T_old = T_temp;
    dT = 0.1 * T_temp;
    T_new = T_old;
  }
#endif
  
  /***************************************************************/
  
  while(1){
    err_chem = give_g(sph_in, sph_out, sph_state, dt, T_new, &e_new, &g0);
    
    if(err_chem < 1.e-2)
      break;
    
    if(isnan(dT) || isinf(dT)){
      dT = 0.1 * T_new;
      printf("(LINE=%d) dT=nan -> %g ID=%d g=%g g0=%g denom=%g\n",
	     __LINE__, dT, sph_in->ID, g, g0, denom);
    }
    
    if(T_new + dT > T_max)
      T_new = dmax(T_large,T_min) , T_large /= 1.1;
    /*  else if(T_new + dT > 30000) T_new /= 1.1;*/
    else if(T_new + dT < 1.5)
      T_new = dmin(T_small,T_max), T_small *= 1.1;
    else T_new += dT, T_new *=1.001 ;
    
    *sph_out = sph_old;
    *sph_state = sph_state_old;

    if(isnan(T_new)){
      printf("[task%03d](iter=%d) ID=%d, T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	     ThisTask, iter, sph_in->ID, sph_out->Utherm,  sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I]);
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
  
  *sph_out = sph_old;
  *sph_state = sph_state_old;

  T_old = T_new;
  T_new = T_old + dT;

  int flag_converged = 0;
	
  int j;
  for(j=iter;j<MAX_ITER;j++){
#ifdef FLAG_BRENT
    if( flag_brent == 1){
      *sph_out = sph_old;
      *sph_state = sph_state_old;
      err = solve_bisection(sph_in, sph_out, sph_state, dt, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

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
      err_chem = give_g(sph_in, sph_out, sph_state, dt, T_new, &e_new, &g);

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
	printf("(LINE=%d) dT=nan -> %g ID=%d g=%g g0=%g denom=%g\n",
	       __LINE__, dT, sph_in->ID, g, g0, denom);
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
      
      *sph_out = sph_old;
      *sph_state = sph_state_old;
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

    *sph_out = sph_old;
    *sph_state = sph_state_old;
  }


#ifdef FLAG_BRENT			
  if(flag_converged == 0){
    T_new = T_gas_init;
		
    for(;j<MAX_ITER;j++){
      *sph_out = sph_old;
      *sph_state = sph_state_old;
	
      d_T = dmax(0.2 * T_new, 2000);
      if(sign > 0)
	T_new = dmax(T_new-d_T,T_new/10.);
      else
	T_new = dmin(T_new + d_T,T_new * 10.);

      err_chem = give_g(sph_in, sph_out, sph_state, dt, T_new, &e_new, &g);
	
      if(err_chem < 1.e-4 && sign * g < 0.){
	if(sign > 0)
	  g_minus = g, T_minus = T_new;
	else
	  g_plus = g, T_plus = T_new;

	flag_brent = 1;
	break;
      }
	
      if(isnan(T_new)){
	printf("[task%03d](iter=%d) ID=%d, T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	       ThisTask, j, sph_in->ID, sph_out->Utherm,  sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I]);
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
      *sph_out = sph_old;
      *sph_state = sph_state_old;
      err = solve_bisection(sph_in, sph_out, sph_state, dt, T_minus, g_minus, T_plus, g_plus, &T_new, &e_new, MAX_ITER + 20 - j);

      if(err < 1.e-5) flag_converged = 1;	       	
    }
  }
#endif

  commit_g(sph_in, sph_out, sph_state, dt, T_new, &e_new, &g);
  sph_out->Utherm = e_new;

 END:;
  if(isnan(err) || isnan(T_new))
    raise(11);

  if(err > 1.0e-4 || FlagPrint == 1){
    printf("not converged(in energy). err=%g id=%d T_old=%g T_new=%g Task=%d\n"
	   "All.ComovingIntegrationOn=%d, All.Time=%g, dt=%g, n_nuclei=%g, T_gas=%g, gamma=%g;\n"
	   "yi[0]=%g, yi[1]=%g, yi[2]=%g;\n"
	   "yi[%d]=%g, yi[%d]=%g;\n",	   
	   err, sph_in->ID, T_gas_init, T_new, ThisTask,
	   All.ComovingIntegrationOn, All.Time, dt, rho_gas/PROTONMASS, T_gas_init, sph_old.Gamma,
	   sph_old.Abund[sELEC], sph_old.Abund[sHI], sph_old.Abund[sHII], 
	   sH2I, sph_old.Abund[sH2I], sHM, sph_old.Abund[sHM]);
#ifdef H2_PLUS
    printf("yi[%d]=%g;\n", sH2II, sph_old.Abund[sH2II]);
#endif
#ifdef HELIUM
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n"
	   , sHeI, sph_old.Abund[sHeI], sHeII, sph_old.Abund[sHeII], sHeIII, sph_old.Abund[sHeIII]);
#endif
#ifdef HD
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n",
	   sDI,   sph_old.Abund[sDI], 
	   sDII,  sph_old.Abund[sDII], 
	   sHDI,  sph_old.Abund[sHDI]);
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
    printf("sph_in->rate_Hminus_pd=%g, sph_in->rate_H2_diss=%g;\n"
	   ,sph_in->rate_Hminus_pd, sph_in->rate_H2_diss);
#ifdef DUST_COOLING
    printf("sph_in->G0=%g;\n", sph_in->G0);
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
    printf("sph_in->k22=%g, sph_in->piHI=%g;\n"
           , sph_in->k22, sph_in->piHI);
#endif

#ifdef METALLICITY
    printf("sph_in->Metallicity=%g;\n", sph_in->Metallicity);
#endif

    printf("yf[0]=%g, yf[1]=%g, yf[2]=%g, yf[3]=%g, yf[4]=%g;\n",
	   sph_out->Abund[sELEC], sph_out->Abund[sHI], sph_out->Abund[sHII], sph_out->Abund[sH2I], sph_out->Abund[sHM]);
  }

  if(err > 1.0e-4){
    *sph_out = sph_old;
    return 1;
  }
  

  return 0;
}
#pragma opt


static double evaluate_g(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state,
			 double dt, double T_gas_new, double *e_new, double *g, int commit)
{
  chemdata_out sph_trial = *sph_out;
  chemdata_state sph_state_trial = *sph_state;
  double e_old = sph_out->Utherm;
  
  chemistry_noneq_chem_implicit(sph_in, &sph_trial, &sph_state_trial, dt, T_gas_new);
  radiative_cooling_d_nometal(sph_in, &sph_trial, &sph_state_trial, T_gas_new);
  double T_gamma = compute_temp_from_abund(sph_trial.Abund, sph_trial.Gamma, sph_trial.Utherm);
  sph_trial.Gamma = compute_gamma_from_abund(sph_trial.Abund, T_gamma);
  
  double gamma_new = sph_trial.Gamma;
  double xmu_new = compute_xmu_from_abund(sph_trial.Abund);

  *e_new = BOLTZMANN * T_gas_new/( (gamma_new - 1.0)*xmu_new*PROTONMASS )/(All.UnitEnergy_in_cgs/All.UnitMass_in_g);
  double de  = (sph_state_trial.rhoU_dot_chem + sph_state_trial.rhoU_dot_rad)/sph_in->rho_cgs *dt / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) ;

#ifdef IONIZATION_PHOTON
  double e_max = BOLTZMANN * 1.e6 /( (gamma_new - 1.0)*xmu_new*PROTONMASS )/(All.UnitEnergy_in_cgs/All.UnitMass_in_g);
  if(e_old + de > e_max && de > 0.){
    de =  dmax(e_max - e_old, 0.);
  }
#endif

  *g  = *e_new - e_old - de;

  if(FlagPrint)
    printf("resid(e_new-e_old-de)=%g e_new=%g old=%g de=%g Tgas=%g rhoU_dot_chem=%g cool=%g rho_gas=%g dt=%g xmu=%g gamma=%g\n"
	   , *g, *e_new, e_old, de, T_gas_new, sph_state_trial.rhoU_dot_chem, sph_state_trial.rhoU_dot_rad, sph_in->rho_cgs, dt, xmu_new, gamma_new);
  
  if(commit){
    *sph_out = sph_trial;
    *sph_state = sph_state_trial;
  }
  
  return sph_state_trial.err_chem;
}

double give_g(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(sph_in, sph_out, sph_state, dt, T_gas_new, e_new, g, 0);
}

static double commit_g(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas_new, double *e_new, double *g){
  return evaluate_g(sph_in, sph_out, sph_state, dt, T_gas_new, e_new, g, 1);
}


double solve_bisection(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state,
		       double dt, double T_min, double g_min, double T_max, double g_max,
		       double *T_new0, double *e_new0, int max_iter){
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
    
    err = give_g(sph_in, sph_out, sph_state, dt, T_cand, &e_new, &val_new);

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


int chemistry_noneq_exp_nometal(int solve_it, chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas)
{
  int mode;
  enum SPECIES species;

  double Ccoef,Dcoef,correction;
  double HIp,HIIp,elecp;
  double H2Ip,HMp;
	
  double t_min, elec_dot, HI_dot, H2I_dot;	

#ifdef HD
  double DIp, DIIp, HDIp, HDIIp, DMp;
#endif

  double rho_gas = sph_in->rho_cgs;
  double n_nuclei = rho_gas/PROTONMASS;

  set_reaction_rate_coeff(sph_in, sph_out, sph_state, T_gas, k);
  
  double y_old[N_SP], y_new[N_SP];
  
  y_old[sELEC]  = sph_out->Abund[sELEC]/XH;
  y_old[sHI]    = sph_out->Abund[sHI]/XH;
  y_old[sHII]   = sph_out->Abund[sHII]/XH;
#ifdef HELIUM
  y_old[sHeI]   = sph_out->Abund[sHeI]/XH;
  y_old[sHeII]  = sph_out->Abund[sHeII]/XH;
  y_old[sHeIII] = sph_out->Abund[sHeIII]/XH;
#endif
  
  y_old[sH2I]   = sph_out->Abund[sH2I]/XH;
#ifdef H2_PLUS
  y_old[sH2II]  = sph_out->Abund[sH2II]/XH;
#endif
  y_old[sHM]    = sph_out->Abund[sHM]/XH;
#ifdef HD
  y_old[sDI]    = sph_out->Abund[sDI]/XH;
  y_old[sDII]   = sph_out->Abund[sDII]/XH;
  y_old[sHDI]   = sph_out->Abund[sHDI]/XH;
#ifdef HD_PLUS
  y_old[sHDII]  = sph_out->Abund[sHDII]/XH;
#endif
#ifdef DMINUS
  y_old[sDM]    = sph_out->Abund[sDM]/XH;
#endif
#endif


  if (!solve_it){
    elec_dot = (k[1]  *sph_out->Abund[sHI]   
#ifdef HELIUM
		+ k[3]  *sph_out->Abund[sHeI]  
		+ k[5]  *sph_out->Abund[sHeII] 
		- k[4]  *sph_out->Abund[sHeII] 
		- k[6]  *sph_out->Abund[sHeIII]
#endif
		- k[2]  *sph_out->Abund[sHII]  
		)*sph_out->Abund[sELEC]*n_nuclei*n_nuclei
			
      + (k[22] *sph_out->Abund[sHI] 
#ifdef HELIUM  
	 + k[23] *sph_out->Abund[sHeII]
	 + k[24] *sph_out->Abund[sHeI]
#endif
	 + k[25] *sph_out->Abund[sHM]  )*n_nuclei
      + (k[8]  *sph_out->Abund[sHI]  * sph_out->Abund[sHM]  
	 + k[16] *sph_out->Abund[sHM]  * sph_out->Abund[sELEC]
	 + k[17] *sph_out->Abund[sHI]  * sph_out->Abund[sHM]  
#ifdef H2_PLUS
	 + k[19] *sph_out->Abund[sHM]  * sph_out->Abund[sHII]
#endif
#ifdef DMINUS
	 + k[48] *sph_out->Abund[sDM]  * sph_out->Abund[sHI] //k[58] originally
#endif
	 )*n_nuclei*n_nuclei
		
      - (k[7]  *sph_out->Abund[sHI]  
#ifdef H2_PLUS
	 + k[20] *sph_out->Abund[sH2II]
#endif
#ifdef HD
	 //+ k[31] *sph_out->Abund[sDII] 
#ifdef HD_PLUS
	 + k[42] *sph_out->Abund[sHDII]
#endif
#ifdef DMINUS
	 + k[43] *sph_out->Abund[sDI] 
#endif
#endif
	 ) *n_nuclei*sph_out->Abund[sELEC]*n_nuclei
      ;


    HI_dot =(k[2]  * sph_out->Abund[sHII]  * sph_out->Abund[sELEC] 
	     + k[11] * sph_out->Abund[sH2I]  * sph_out->Abund[sHI]   *2.
	     + k[12] * sph_out->Abund[sH2I]  * sph_out->Abund[sHII]  
	     + k[13] * sph_out->Abund[sH2I]  * sph_out->Abund[sELEC] *2.
#ifdef HELIUM
	     + k[15] * sph_out->Abund[sHeI]  * sph_out->Abund[sHII]    
#endif
	     )*n_nuclei*n_nuclei

#ifdef THREEBODY
      + k[58] * sph_out->Abund[sH2I]  * sph_out->Abund[sH2I]  * n_nuclei*n_nuclei *2
#endif
      + k[25] * sph_out->Abund[sHM]   * n_nuclei
#ifdef H2_PLUS
      + k[26] * sph_out->Abund[sH2II] * n_nuclei
#endif
      + k[29] * sph_out->Abund[sH2I]  * n_nuclei * 2
      -(
	(k[1]  * sph_out->Abund[sELEC]
#ifdef H2_PLUS
	 + k[10] * sph_out->Abund[sH2II]
#endif
#ifdef HELIUM
	 + k[14] * sph_out->Abund[sHeII]
#endif
	 ) * n_nuclei
#ifdef THREEBODY
	+ (k[55] * sph_out->Abund[sHI]  * 2.
	   +  k[57] * sph_out->Abund[sH2I] * 2. )* sph_out->Abund[sHI] * n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
	+ 2.0 * k[51]
#endif
	+ k[22]
	) * sph_out->Abund[sHI] * n_nuclei;
		
		
    H2I_dot = (k[8]  *sph_out->Abund[sHM]  *sph_out->Abund[sHI]  
#ifdef H2_PLUS
	       + k[10] *sph_out->Abund[sH2II]*sph_out->Abund[sHI]  
	       + k[21] *sph_out->Abund[sH2II]*sph_out->Abund[sHM] 
#endif
	       )*n_nuclei*n_nuclei
#ifdef H2_FORMATION_ON_DUST
      + k[51] *sph_out->Abund[sHI] * n_nuclei
#endif
      
#ifdef THREEBODY
      + (k[55] *sph_out->Abund[sHI]  *sph_out->Abund[sHI]  *sph_out->Abund[sHI]  *n_nuclei
	 + k[57] *sph_out->Abund[sHI]  *sph_out->Abund[sHI]  *sph_out->Abund[sH2I] *n_nuclei)*n_nuclei*n_nuclei
#endif 
      - (k[11] *sph_out->Abund[sHI]   * n_nuclei
	 + k[12] *sph_out->Abund[sHII]  * n_nuclei
	 + k[13] *sph_out->Abund[sELEC] * n_nuclei
#ifdef THREEBODY
	 + k[58] *sph_out->Abund[sH2I] * n_nuclei
#endif
	 + k[27] 
	 + k[29]
	 ) * sph_out->Abund[sH2I] * n_nuclei;
		
    t_min = dmin(fabs(sph_out->Abund[sELEC] * n_nuclei/elec_dot), fabs(sph_out->Abund[sHI] * n_nuclei/HI_dot));
    if(n_nuclei > 1. || n_nuclei>N_THREEBODY)  //removed[0815], modified[1009] 
      t_min = dmin(t_min, fabs(sph_out->Abund[sH2I] * n_nuclei/H2I_dot));
    
    sph_out->t_chem = t_min;

    chemical_cooling(mode=0, sph_in, sph_out, sph_state, T_gas, dt, y_old, y_old);  
  }else{  
    /* ----------- HI */  
    //added reaction 16,17,18,20,21[chon140112]

    Ccoef = (k[2]  * sph_out->Abund[sHII]  * sph_out->Abund[sELEC] 
	     + k[11] * sph_out->Abund[sH2I]  * sph_out->Abund[sHI]   *2.
	     + k[12] * sph_out->Abund[sH2I]  * sph_out->Abund[sHII]  
	     + k[13] * sph_out->Abund[sH2I]  * sph_out->Abund[sELEC] *2.
#ifdef HELIUM
	     + k[15] * sph_out->Abund[sHeI]  * sph_out->Abund[sHII]
#endif
	     + k[16] * sph_out->Abund[sHM] * sph_out->Abund[sELEC]
	     + k[17] * sph_out->Abund[sHM] * sph_out->Abund[sHI]
	     + 2 * k[18] * sph_out->Abund[sHM] * sph_out->Abund[sHII]
#ifdef H2_PLUS
	     + k[20] * sph_out->Abund[sH2II] * sph_out->Abund[sELEC] *2.
	     + k[21] * sph_out->Abund[sH2II] * sph_out->Abund[sHM] 
#endif
#ifdef HD
	     + k[32] * sph_out->Abund[sDI] * sph_out->Abund[sHII] 
#endif
	     )*n_nuclei*n_nuclei
#ifdef THREEBODY
      + k[58] * sph_out->Abund[sH2I]  * sph_out->Abund[sH2I]  * n_nuclei*n_nuclei *2.
#endif
      + k[25] * sph_out->Abund[sHM]   * n_nuclei
#ifdef H2_PLUS
      + k[26] * sph_out->Abund[sH2II] * n_nuclei
#endif
      + k[29] * sph_out->Abund[sH2I]  * n_nuclei * 2.
      ;

    
    Dcoef = (  k[1]  * sph_out->Abund[sELEC]
#ifdef H2_PLUS
	       + k[10] * sph_out->Abund[sH2II]
#endif
#ifdef HELIUM
	       + k[14] * sph_out->Abund[sHeII] 
#endif
#ifdef HD
	       + k[33] * sph_out->Abund[sDII]
	       + k[38] * sph_out->Abund[sHDI]
#endif
	       ) * n_nuclei
#ifdef THREEBODY
      + (k[55] * sph_out->Abund[sHI]  * 2.
	 +  k[57] * sph_out->Abund[sH2I] * 2. )* sph_out->Abund[sHI] * n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
      + 2. * k[51]
#endif
      + k[22];


    HI_dot = fabs(Ccoef - Dcoef);  
    HIp   = (Ccoef*dt + sph_out->Abund[sHI]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* ----------- HII */
    Ccoef = k[1]  * HIp * sph_out->Abund[sELEC] * n_nuclei*n_nuclei
#ifdef H2_PLUS
      + k[10] * HIp * sph_out->Abund[sH2II] * n_nuclei*n_nuclei 
#endif
#ifdef HELIUM
      + k[14] * HIp * sph_out->Abund[sHeII] * n_nuclei*n_nuclei
#endif
      + k[22] * HIp * n_nuclei
#ifdef H2_PLUS
      + k[26] * sph_out->Abund[sH2II] * n_nuclei
      + k[28] * sph_out->Abund[sH2II] * n_nuclei * 2
#endif
      /*added by chon[131221] */
#ifdef HD
      + k[33] * HIp * sph_out->Abund[sDII] * n_nuclei * n_nuclei
      + k[37] * sph_out->Abund[sH2I] * sph_out->Abund[sDII] * n_nuclei * n_nuclei
#ifdef HD_PLUS
      + k[36] * HIp * sph_out->Abund[sHDII] * n_nuclei * n_nuclei
      + k[52] * sph_out->Abund[sHDII] * n_nuclei
#endif
#endif
      /*[131221]*/
      ;

    Dcoef = (k[2]  * sph_out->Abund[sELEC]
#ifdef H2_PLUS
	     + k[9]  * HIp             
	     + k[12] * sph_out->Abund[sH2I]
#endif
#ifdef HELIUM 
	     + k[15] * sph_out->Abund[sHeI]
#endif
	     )* n_nuclei
      /*added by chon[131221] */
      + k[18] * sph_out->Abund[sHM] * n_nuclei
#ifdef H2_PLUS
      + k[19] * sph_out->Abund[sHM] * n_nuclei
#endif
#ifdef HD
      + (k[32] * sph_out->Abund[sDI] 
	 + k[39] * sph_out->Abund[sHDI]
#ifdef HD_PLUS
	 + k[40] * sph_out->Abund[sDI]
#endif
#ifdef DMINUS
	 + k[45] * sph_out->Abund[sDM] 
#endif
	 ) *n_nuclei /* modified sph_out->Abund[sHDI] -> sph_out->Abund[sDM] [0926]*/ 
#endif
      /*[131221]*/
      ;

    
    HIIp  = (Ccoef*dt + sph_out->Abund[sHII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
    
#ifdef HELIUM
    /* ----------- HeI */
    Ccoef = (k[4]  * sph_out->Abund[sHeII] * sph_out->Abund[sELEC]
	     + k[14] * HIp              * sph_out->Abund[sHeII])*n_nuclei*n_nuclei;
  
    Dcoef = k[3]  * sph_out->Abund[sELEC] * n_nuclei 
      + k[15] * HIIp             * n_nuclei
      + k[24];
  
    double HeIp  = (Ccoef*dt + sph_out->Abund[sHeI]*n_nuclei) / ( 1. + Dcoef*dt ) / n_nuclei;


    /* ----------- HeII */
    Ccoef = k[3]  * HeIp              * sph_out->Abund[sELEC] * n_nuclei*n_nuclei
      + k[6]  * sph_out->Abund[sHeIII] * sph_out->Abund[sELEC] * n_nuclei*n_nuclei
      + k[15] * HeIp              * HIIp             * n_nuclei*n_nuclei
      + k[24] * HeIp*n_nuclei;

    Dcoef = k[4]  * sph_out->Abund[sELEC] * n_nuclei
      + k[5]  * sph_out->Abund[sELEC] * n_nuclei
      + k[14] * HIp              * n_nuclei
      + k[23]; 

    double HeIIp  = (Ccoef*dt + sph_out->Abund[sHeII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

      
    /* ----------- HeIII */
    Ccoef = k[5]  * HeIIp * sph_out->Abund[sELEC] * n_nuclei*n_nuclei
      + k[23] * HeIIp * n_nuclei;
      
    Dcoef = k[6] * sph_out->Abund[sELEC] * n_nuclei;
    
    double HeIIIp = (Ccoef*dt + sph_out->Abund[sHeIII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    if(HeIp   < 0) HeIp   = 0.;
    if(HeIIp  < 0) HeIIp  = 0.;
    if(HeIIIp < 0) HeIIIp = 0.;

    sph_out->Abund[sHeI]   = HeIp;
    sph_out->Abund[sHeII]  = HeIIp;
    sph_out->Abund[sHeIII] = HeIIIp;
#endif
  
    /* --------- elec */
    Ccoef = (k[22] *HIp  
#ifdef HELIUM	
	     + k[23] *HeIIp
	     + k[24] *HeIp 
#endif
	     )*n_nuclei;
  
    Dcoef = - (k[1]*HIp   - k[2]*HIIp
#ifdef HELIUM	
	       + k[3]*HeIp  - k[6]*HeIIIp
	       + k[5]*HeIIp - k[4]*HeIIp
#endif
	       )*n_nuclei;

    elec_dot = fabs(Ccoef - Dcoef);  
    elecp = (Ccoef*dt + sph_out->Abund[sELEC]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    if(HIp   < 0) HIp   = 0;
    if(HIIp  < 0) HIIp  = 0;
    if(elecp < 0) elecp = 0;
    

    sph_out->Abund[sHI]   = HIp;
    sph_out->Abund[sHII]  = HIIp;    
    sph_out->Abund[sELEC] = elecp;


    /* --------- H2 */
    /* [SH131015 : crazyyyyyyy] */
  
    /* HM fast reactions */
    /*
      denom = (k[8] + k[17]) * sph_out->Abund[sHI] 
      + (k[18]+ k[19]) * sph_out->Abund[sHII] 
      +  k[16]         * sph_out->Abund[sELEC] 
      +  k[25]/n_nuclei;
  
      sph_out->Abund[sHM] = (k[7]*sph_out->Abund[sHI]*sph_out->Abund[sELEC]) / denom;
    */
    Ccoef = k[7]*sph_out->Abund[sHI]*sph_out->Abund[sELEC] * n_nuclei*n_nuclei
#ifdef HD
#ifdef DMINUS
	       + k[47] * HIp * sph_out->Abund[sDM] *n_nuclei*n_nuclei
#endif
#endif
      ;

    Dcoef = ((k[8] + k[17]) * sph_out->Abund[sHI] 
	     + k[18] * sph_out->Abund[sHII] 
#ifdef H2_PLUS
	     + k[19] * sph_out->Abund[sHII] 
#endif
	     + k[16] * sph_out->Abund[sELEC]  
	     ) * n_nuclei
	       +  k[25]
#ifdef HD	       
#ifdef DMINUS
      +  k[46] * sph_out->Abund[sDI] *n_nuclei
#endif
#endif
      ;
  
    HMp = (Ccoef*dt + sph_out->Abund[sHM]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* H2II fast reactions */
    /*  denom = k[10]*sph_out->Abund[sHI] 
	+ k[20]*sph_out->Abund[sELEC] 
	+ k[21]*sph_out->Abund[sHM] 
	+(k[26]+k[28])/n_nuclei;
    
	sph_out->Abund[sH2II] = (  k[9] * sph_out->Abund[sHI]  * sph_out->Abund[sHII]
	+ k[12] * sph_out->Abund[sH2I] * sph_out->Abund[sHII]
	+ k[19] * sph_out->Abund[sHM]  * sph_out->Abund[sHII]
	+ k[27] * sph_out->Abund[sH2I]/n_nuclei         )/ denom; */

    /* noneq H2II added by chon.[131228]*/

#ifdef H2_PLUS
    Ccoef = ( k[9] * sph_out->Abund[sHI]  * sph_out->Abund[sHII]
	      +k[12] * sph_out->Abund[sH2I] * sph_out->Abund[sHII] 
	      +k[19] * HMp  * sph_out->Abund[sHII] )*n_nuclei *n_nuclei
      + k[27]* sph_out->Abund[sH2I]*n_nuclei;

    Dcoef = (   k[10]*sph_out->Abund[sHI] 
		+ k[20]*sph_out->Abund[sELEC]
		+ k[21]*HMp ) *n_nuclei
      +k[26] + k[28];

    double H2IIp = (Ccoef*dt + sph_out->Abund[sH2II]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif


    Ccoef = ( k[8]  * HMp   * sph_out->Abund[sHI]
#ifdef H2_PLUS
	      + k[10] * H2IIp * sph_out->Abund[sHI]
	      + k[21] * H2IIp * HMp
#endif
#ifdef HD
	       + k[38] * HIp * sph_out->Abund[sHDI] 
	       + k[39] * HIIp * sph_out->Abund[sHDI]
#endif
	      ) *n_nuclei*n_nuclei
#ifdef THREEBODY
      +( k[55] * sph_out->Abund[sHI]   * sph_out->Abund[sHI] * sph_out->Abund[sHI]
	 + k[57] * sph_out->Abund[sHI]   * sph_out->Abund[sHI] * sph_out->Abund[sH2I] )*n_nuclei*n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
      + k[51] * sph_out->Abund[sHI]*n_nuclei
#endif
      ;
    
    Dcoef = (k[11] * sph_out->Abund[sHI] 
#ifdef H2_PLUS
	     + k[12] * HIIp 
#endif
	     + k[13] * sph_out->Abund[sELEC]
#ifdef HD
	     + k[35] * sph_out->Abund[sDI]  
	     + k[37] * sph_out->Abund[sDII]  
#endif
	     )*n_nuclei
#ifdef THREEBODY
	     + k[58] * sph_out->Abund[sH2I] * n_nuclei
#endif
	     + k[27] 
	     + k[29];

    H2I_dot = fabs(Ccoef - Dcoef);  
    H2Ip = (Ccoef*dt + sph_out->Abund[sH2I]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

    if(HMp < 0) HMp = 0;
#ifdef H2_PLUS
    if(H2IIp < 0) H2IIp = 0;
#endif
    if(H2Ip < 0) H2Ip = 0;

    sph_out->Abund[sHM] = HMp;
#ifdef H2_PLUS
    sph_out->Abund[sH2II] = H2IIp;
#endif
    sph_out->Abund[sH2I] = H2Ip;
  
  
#ifdef HD

    DIp=DIIp=HDIp=HDIIp=DMp=0;

    /* ----------- DI */
    Ccoef  = //k[31] * sph_out->Abund[sDII] * sph_out->Abund[sELEC] * n_nuclei*n_nuclei
	     k[33] * sph_out->Abund[sDII] * HIp   * n_nuclei*n_nuclei
	     + k[38] * sph_out->Abund[sHDI] * HIp   * n_nuclei*n_nuclei
#ifdef HD_PLUS
	     + k[42] * sph_out->Abund[sHDII]* sph_out->Abund[sELEC] * n_nuclei*n_nuclei
	     + k[52] * sph_out->Abund[sHDII]* n_nuclei
#endif
#ifdef DMINUS
	     + k[44] * sph_out->Abund[sDII] * sph_out->Abund[sDM]   * n_nuclei*n_nuclei * 2.0 
	     + k[45] * HIIp * sph_out->Abund[sDM]   * n_nuclei*n_nuclei
	     + k[47] * HIp  * sph_out->Abund[sDM]   * n_nuclei*n_nuclei
#endif
	     + k[53] * sph_out->Abund[sHDI] * n_nuclei;
	     
    Dcoef = ( k[32] * HIIp
	      //+ k[34] * sph_out->Abund[sHI] 
	      + k[35] * H2Ip
#ifdef HD_PLUS
	      + k[40] * HIIp
#endif
#ifdef DMINUS
	      + k[43] * sph_out->Abund[sELEC]
#endif
	      )* n_nuclei
      + k[50];

    DIp   = (Ccoef*dt + sph_out->Abund[sDI]*n_nuclei)  / (1. + Dcoef*dt) / n_nuclei;
  
  
    /* ----------- DII */
    Ccoef = k[32] * DIp                 * HIIp * n_nuclei*n_nuclei
      + k[39] * sph_out->Abund[sHDI] * HIIp * n_nuclei*n_nuclei
      + k[50] * DIp                 * n_nuclei
#ifdef HD_PLUS
      + k[51] * sph_out->Abund[sHDII]* n_nuclei
#endif
      ;

    Dcoef = ( //k[31] * sph_out->Abund[sELEC]
	      k[33] * HIp
	      + k[37] * H2Ip
#ifdef DMINUS
	      + k[41] * HIp
	      + k[44] * sph_out->Abund[sDM]  
#endif
	      ) * n_nuclei;
	    
    DIIp  = (Ccoef*dt + sph_out->Abund[sDII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* ----------- HDI */
    Ccoef = ( //k[34] * DIp            * sph_out->Abund[sHI]  
	      k[35] * DIp            * H2Ip
	      + k[37] * DIIp           * H2Ip
#ifdef HD_PLUS
	      + k[36] * HIp * sph_out->Abund[sHDII]
#endif	      
#ifdef DMINUS
	      + k[48] * HIp * sph_out->Abund[sDM]  
#endif
	      )* n_nuclei*n_nuclei;

    Dcoef = k[38] * HIp * n_nuclei
      + k[39] * HIIp * n_nuclei
      + k[53];

    HDIp  = (Ccoef * dt + sph_out->Abund[sHDI] * n_nuclei) / (1. +  Dcoef * dt) / n_nuclei;
    
#ifdef HD_PLUS
    /* ----------- HDII */
    Ccoef = ( k[40] * DIp  * sph_out->Abund[sHII]
	      + k[41] * DIIp * sph_out->Abund[sHI] ) * n_nuclei*n_nuclei;
      
    Dcoef = k[36] * sph_out->Abund[sHI]   * n_nuclei
      + k[42] * sph_out->Abund[sELEC] * n_nuclei
      + k[51]
      + k[52];

    HDIIp = (Ccoef*dt + sph_out->Abund[sHDII]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
    
#ifdef DMINUS    
    /* ----------- DM */
    Ccoef = ( k[43] * DIp  * sph_out->Abund[sELEC]
	      + k[46] * DIp  * sph_out->Abund[sHM]  ) * n_nuclei*n_nuclei;
      
    Dcoef = ( k[44] * DIIp           
	      + k[45] * sph_out->Abund[sHII]
	      + k[47] * sph_out->Abund[sHI] 
	      + k[48] * sph_out->Abund[sHI]  ) * n_nuclei;
    
    DMp   = (Ccoef*dt + sph_out->Abund[sDM]*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
      
    if(DIp    < 0) DIp   =0;
    if(DIIp   < 0) DIIp  =0;
    if(HDIp   < 0) HDIp  =0;
    sph_out->Abund[sDI]    = DIp;
    sph_out->Abund[sDII]   = DIIp;
    sph_out->Abund[sHDI]   = HDIp;   

#ifdef HD_PLUS
    if(HDIIp  < 0) HDIIp =0;
    sph_out->Abund[sHDII]  = HDIIp;
#endif
#ifdef DMINUS
    if(DMp    < 0) DMp   =0;
    sph_out->Abund[sDM]    = DMp;    
#endif  
#endif       
    
    /* --------- Undertake mass conservation corrections. */
  
    //sph_out->Abund[sELEC] = 0.0;
  
    correction = H_number_fraction / (sph_out->Abund[sHI] 
				      +   sph_out->Abund[sHII] 
				      +   sph_out->Abund[sHM] 
				      + 2*sph_out->Abund[sH2I]
#ifdef H2_PLUS
				      + 2*sph_out->Abund[sH2II]
#endif
#ifdef HD
				      +   sph_out->Abund[sHDI] 
#ifdef HD_PLUS
				      +   sph_out->Abund[sHDII]
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

  
    sph_out->Abund[sHI]  *= correction;
    sph_out->Abund[sHII] *= correction;
    //sph_out->Abund[sELEC]+= sph_out->Abund[sHII];
    
    sph_out->Abund[sHM]   *= correction;
    sph_out->Abund[sH2I]  *= correction;
#ifdef H2_PLUS
    sph_out->Abund[sH2II] *= correction;
#endif
    //sph_out->Abund[sELEC]  = sph_out->Abund[sELEC] - sph_out->Abund[sHM] + sph_out->Abund[sH2II];
  

#ifdef HD
    correction = D_number_fraction / (sph_out->Abund[sDI] 
				      + sph_out->Abund[sDII]
				      + sph_out->Abund[sHDI]
#ifdef HD_PLUS
				      + sph_out->Abund[sHDII]
#endif
#ifdef DMINUS
				      + sph_out->Abund[sDM]
#endif
				      );

    if(0)
      if(correction > 2.0 || correction < 0.5){
	printf("D! ThisTask %d correction=%g\n"
	       , ThisTask, correction);
	fflush(stdout);
	endrun();
      } 

    sph_out->Abund[sDI]   *= correction;
    sph_out->Abund[sDII]  *= correction;
    sph_out->Abund[sHDI]  *= correction;
#ifdef HD_PLUS
    sph_out->Abund[sHDII] *= correction;
#endif
#ifdef DMINUS
    sph_out->Abund[sDM]   *= correction;
#endif
    //sph_out->Abund[sELEC]  = sph_out->Abund[sELEC] + sph_out->Abund[sDII] +  sph_out->Abund[sHDII] - sph_out->Abund[sDM];
#endif 
  
  
#ifdef HELIUM
    correction = He_number_fraction / (sph_out->Abund[sHeI] 
				       + sph_out->Abund[sHeII] 
				       + sph_out->Abund[sHeIII]);

    if(0)
      if(correction > 1.1 || correction < 0.9){
	printf("He! ThisTask %d correction=%g %g %g %g\n"
	       , ThisTask, correction, sph_out->Abund[sHeI], sph_out->Abund[sHeII], sph_out->Abund[sHeIII]);
	fflush(stdout);
	endrun();
      }
  
    sph_out->Abund[sHeI]   *= correction;
    sph_out->Abund[sHeII]  *= correction;
    sph_out->Abund[sHeIII] *= correction;
    //sph_out->Abund[sELEC]  = sph_out->Abund[sELEC] + sph_out->Abund[sHeII] + 2.*sph_out->Abund[sHeIII];
#endif

    /*[SH130423: elec calculation]*/

    sph_out->Abund[sELEC] = 0.0;

    sph_out->Abund[sELEC] += sph_out->Abund[sHII];

    sph_out->Abund[sELEC]  = sph_out->Abund[sELEC] - sph_out->Abund[sHM];

#ifdef H2_PLUS
    sph_out->Abund[sELEC] += sph_out->Abund[sH2II];
#endif

#ifdef HD
    sph_out->Abund[sELEC]  += sph_out->Abund[sDII]; 
#ifdef HD_PLUS
    sph_out->Abund[sELEC]  += sph_out->Abund[sHDII];
#endif
#ifdef DMINUS
    sph_out->Abund[sELEC]  -= sph_out->Abund[sDM];
#endif
#endif

#ifdef HELIUM
    sph_out->Abund[sELEC]  = sph_out->Abund[sELEC] + sph_out->Abund[sHeII] + 2.*sph_out->Abund[sHeIII];
#endif

    if(isnan(sph_out->Abund[sHI])){
      printf("[task%03d]HI is nan (ID=%d, dt=%g, rho_gas=%g, T_gas=%g, mode=%d)\n", ThisTask, sph_in->ID, dt, rho_gas, T_gas, solve_it);
      return 2;
    }
    if(isnan(T_gas)){
      printf("[task%03d]line=%d T_gas is nan (dt=%g, rho_gas=%g, T_gas=%g, mode=%d)\n", ThisTask, __LINE__, dt, rho_gas, T_gas, solve_it);
      return 2;
    }

    /* added[0813] */
    t_min = dmin(fabs(sph_out->Abund[sELEC] * n_nuclei/elec_dot), fabs(sph_out->Abund[sHI] * n_nuclei/HI_dot));
    t_min = dmin(t_min, fabs(sph_out->Abund[sH2I] * n_nuclei/H2I_dot));

    sph_out->t_chem = t_min;


    y_new[sELEC]  = sph_out->Abund[sELEC]/XH;
    y_new[sHI]    = sph_out->Abund[sHI]/XH;
    y_new[sHII]   = sph_out->Abund[sHII]/XH;
#ifdef HELIUM
    y_new[sHeI]   = sph_out->Abund[sHeI]/XH;
    y_new[sHeII]  = sph_out->Abund[sHeII]/XH;
    y_new[sHeIII] = sph_out->Abund[sHeIII]/XH;
#endif
    y_new[sH2I]   = sph_out->Abund[sH2I]/XH;
#ifdef H2_PLUS
    y_new[sH2II]  = sph_out->Abund[sH2II]/XH;
#endif
    y_new[sHM]    = sph_out->Abund[sHM]/XH;
#ifdef HD
    y_new[sDI]    = sph_out->Abund[sDI]/XH;
    y_new[sDII]   = sph_out->Abund[sDII]/XH;
    y_new[sHDI]   = sph_out->Abund[sHDI]/XH;
#ifdef HD_PLUS
    y_new[sHDII]  = sph_out->Abund[sHDII]/XH;
#endif
#ifdef DMINUS
    y_new[sDM]    = sph_out->Abund[sDM]/XH;
#endif
#endif
  
    chemical_cooling(mode=1, sph_in, sph_out, sph_state, T_gas, dt, y_new, y_old);  
  }

  return 1; 
}


#define MAX_ITER_CHEM 1000
#define ITER_RESET 40
#define BETA 1.0

int chemistry_noneq_chem_implicit(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double dt, double T_gas)
{
  int n_sp;
  int inc, info;
	
  enum SPECIES species;
  chemdata_out sph_old;

  sph_old = *sph_out;
  
  n_sp = N_SP;
  inc  = 1   ;

  double rho_gas = sph_in->rho_cgs;
  double n_nuclei  = rho_gas/PROTONMASS;
  double nH = XH * n_nuclei;
    
#ifdef SAHA
  if(n_nuclei> N_EQCHEM){
    double y_H2_init = sph_out->Abund[sH2I] / XH;
    double y_Hp_init = sph_out->Abund[sHII] / XH;
    double y_Hm_init = sph_out->Abund[sHM]  / XH;

    double y_input[N_SP];
    y_input[sELEC] = sph_out->Abund[sELEC] / XH;
    y_input[sHI]   = sph_out->Abund[sHI]   / XH;
    y_input[sHII]  = sph_out->Abund[sHII]  / XH;
    y_input[sH2I]  = sph_out->Abund[sH2I]  / XH;
    y_input[sHM]   = sph_out->Abund[sHM]   / XH;

    double err = chem_saha_5sp(nH, T_gas, y_input);

    sph_out->Abund[sHI]  = y_input[sHI]  * XH;
    sph_out->Abund[sHII] = y_input[sHII] * XH;
    sph_out->Abund[sH2I] = y_input[sH2I] * XH;
    sph_out->Abund[sHM]  = y_input[sHM]  * XH;
    sph_out->Abund[sELEC] = y_input[sELEC] * XH;    
    
    double edot_H2form = 4.48    * eV_to_erg * (y_input[sH2I] - y_H2_init) * nH / dt;    
    double edot_ionize = -13.598 * eV_to_erg * (y_input[sHII] - y_Hp_init) * nH / dt;
    double edot_HMform = 0.78    * eV_to_erg * (y_input[sHM]  - y_Hm_init) * nH / dt;

    if(dt < 1.e-30)
      edot_H2form = edot_ionize = edot_HMform = 0.;

    sph_state->rhoU_dot_chem = edot_H2form + edot_ionize + edot_HMform;

#ifdef DEBUG_CHEMISTRY_SC
    sph_out->edot_ciHI = 0.;
    sph_out->edot_chem = edot_H2form;
#endif

    sph_state->err_chem = err;

    return 1;
  }
#endif


  set_reaction_rate_coeff(sph_in, sph_out, sph_state, T_gas, k);
  
  double yi[n_sp];
       
  yi[sELEC]  = sph_out->Abund[sELEC]/XH;
  yi[sHI]    = sph_out->Abund[sHI]/XH;
  yi[sHII]   = sph_out->Abund[sHII]/XH;
#ifdef HELIUM
  yi[sHeI]   = sph_out->Abund[sHeI]/XH;
  yi[sHeII]  = sph_out->Abund[sHeII]/XH;
  yi[sHeIII] = sph_out->Abund[sHeIII]/XH;
#endif
  
  yi[sH2I]   = sph_out->Abund[sH2I]/XH;
#ifdef H2_PLUS
  yi[sH2II]  = sph_out->Abund[sH2II]/XH;
#endif
  yi[sHM]    = sph_out->Abund[sHM]/XH;
#ifdef HD
  yi[sDI]    = sph_out->Abund[sDI]/XH;
  yi[sDII]   = sph_out->Abund[sDII]/XH;
  yi[sHDI]   = sph_out->Abund[sHDI]/XH;
#ifdef HD_PLUS
  yi[sHDII]  = sph_out->Abund[sHDII]/XH;
#endif
#ifdef DMINUS
  yi[sDM]    = sph_out->Abund[sDM]/XH;
#endif
#endif
		
  int i_sp,j_sp,itr,ifunc;
  int err_flag = 0, ipiv[n_sp];
  double drfdy[n_sp][n_sp],dy[n_sp];
  double y_init[n_sp], yi_init[n_sp];
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
  for(i_sp=0;i_sp<n_sp;i_sp++) 
    y_init[i_sp] = yi[i_sp];

  //give_initial_guess(T_gas, nH, dt, yi, y_init, k);
  for(i_sp=0;i_sp<n_sp;i_sp++){
    dy[i_sp] = yi[i_sp] - y_init[i_sp];
    yi_init[i_sp] = yi[i_sp]; 
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
      (k[7]*yi[sHI]*yi[sELEC]+k_inverse[8]*yi[sH2I]*yi[sELEC]+k_inverse[18]*yi[sHI]*yi[sHI])
      /(k_inverse[7]+k[8]*yi[sHI]+k[18]*yi[sHII]+k[25]);
#endif
    
    ifunc=ReactionRates_nometal(itr,nH,yi,rf_tot);
    
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
      ifunc=ReactionRates_diff_nometal_fp(itr, j_sp, nH, yi, drfdy[j_sp]);
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
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,rf_f);
      
      y_tmp[j_sp]=yi[j_sp]-0.5*Dyi;
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,rf_b);
      
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
			
      sph_out->t_chem = t_chem1;
    }
  }else{
    if(err_max>1.e-4){
      printf("iteration is not converged err=%le(in %d),itr=%d,"
	     " id=%d.\n dt=%g, rho=%g, T_gas=%g, n_nuclei=%g;\n"
	     "y[%d]=%g, y[%d]=%g, y[%d]=%g, ",
	     err_max, err_flag, itr, 
	     sph_in->ID, dt, rho_gas, T_gas, n_nuclei,
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

      sph_state->err_chem = err_max;
      *sph_out = sph_old;
      return 0;
    }else{
      t_chem1 = 1.0e20;
      for(i_sp=0;i_sp<n_sp;i_sp++) {
	if(yi[i_sp]-y_init[i_sp]!=0.0) {
	  t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
	  if(t_chem2<t_chem1) t_chem1 = t_chem2;
	}
			
	sph_out->t_chem = t_chem1;
      }
    }
  }
  
  sph_state->err_chem = err_max;

  int mode;
  chemical_cooling(mode=1, sph_in, sph_out, sph_state, T_gas, dt, yi, y_init);
  
  for(i_sp=0;i_sp<n_sp;i_sp++){
    if(yi[i_sp] < y_min[i_sp])
      yi[i_sp] = y_min[i_sp];
  }

  sph_out->Abund[sELEC] =yi[sELEC] *XH;
  sph_out->Abund[sHI]   =yi[sHI]   *XH;
  sph_out->Abund[sHII]  =yi[sHII]  *XH;
#ifdef HELIUM
  sph_out->Abund[sHeI]  =yi[sHeI]  *XH;
  sph_out->Abund[sHeII] =yi[sHeII] *XH;
  sph_out->Abund[sHeIII]=yi[sHeIII]*XH;
#endif
  
  sph_out->Abund[sH2I]  =yi[sH2I]  *XH;
#ifdef H2_PLUS
  sph_out->Abund[sH2II] =yi[sH2II] *XH;
#endif
  sph_out->Abund[sHM]   =yi[sHM]   *XH;
#ifdef HD    
  sph_out->Abund[sDI]   =yi[sDI]   *XH;
  sph_out->Abund[sDII]  =yi[sDII]  *XH;
  sph_out->Abund[sHDI]  =yi[sHDI]  *XH;
#ifdef HD_PLUS
  sph_out->Abund[sHDII] =yi[sHDII] *XH;
#endif
#ifdef DMINUS
  sph_out->Abund[sDM]   =yi[sDM]   *XH;
#endif
#endif

#ifdef CALC_CONTINUUM_OPACITY
  double ni[N_SP];
  for(int i=0;i<N_SP;i++){
    ni[i] = yi[i] * nH;
  }

  double l_core = sph_state->Jeans_Length;
  double t_ff =  compute_freefall_time(rho_gas); 
  solve_level_population_H_new(T_gas, ni, l_core, t_ff);
#endif

#ifdef DEBUG_CHEMISTRY_SC
  for(int i=0;i<70;i++){
    sph_state->k[i] = k[i];
    sph_state->k_inv[i] = k_inverse[i];
  }
#endif

  return 1;
}


#pragma optimize("", on)
#undef MAX_ITER_CHEM

static inline int ReactionRates_nometal(int itr, double nH, double *yi, double *rf_tot)
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
static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *rf_tot)
{
  int i_sp, n_sp = N_SP, ifunc;
  double ni[n_sp];
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    ni[i_sp]=yi[i_sp]*nH;
  }
  	
  for(i_sp=0;i_sp<n_sp;i_sp++) 
    rf_tot[i_sp]=0.0;
  
  ifunc = switch_func[i_sp_diff](nH, ni, rf_tot);

  return 1;
}

int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HI(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HII(double nH,double *ni,double *rf){
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
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *rf){
  double rate;
  enum SPECIES species;

  rate = k[6]*ni[sELEC];
  rf[sELEC]  += -rate;
  rf[sHeII]  +=  rate;
  rf[sHeIII] += -rate;

  return 0;
}
#endif

int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *rf){
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
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HM(double nH,double *ni,double *rf){
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
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_DII(double nH,double *ni,double *rf){
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

int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *rf){
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
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *rf){
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
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *rf){
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

int func_other(double nH,double *ni,double *rf){
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


void set_reaction_rate_coeff(chemdata_in *sph_in, chemdata_out *sph_out, chemdata_state *sph_state, double T_gas, double *k){
  update_radiation_state(sph_in, sph_out, sph_state, T_gas);

  double rho_cgs = sph_in->rho_cgs;
  double n_nuclei = rho_cgs / PROTONMASS;
  double nH = XH * n_nuclei;
  
  /*********** initialize ************/
  for(int j=0; j<k_MAX; j++)
    k[j] = 0.0;
  /***********************************/

#if defined(LW_RADIATION_LOCAL) && defined(DUST_COOLING)
  sph_state->G0 = sph_in->G0;
#endif

  /* ------- Compute temperature index. */
  double ln_T0 = log(T_start);
  double ln_T9 = log(T_end);
  
  double log_T  = log10(T_gas);
  double log_T2 = log_T * log_T;
  double log_T3 = log_T * log_T * log_T;
  
  double ln_T  = log(T_gas);
  if(ln_T < ln_T0) ln_T = ln_T0;
  if(ln_T > ln_T9) ln_T = ln_T9;
  int i_T = (int) ((ln_T-ln_T0)/dln_T);
  
  if(i_T<0) i_T=0;
  if(i_T>N_T-2) i_T=N_T-2;

  if(i_T < 0){
    printf("i_T=%d (ln_T=%g,ln_T0=%g,dln_T=%g)\n",i_T,ln_T,ln_T0,dln_T);

    printf("ThisTask %d ID %d T %g n %g H2I %g HI %g HII %g\n"
	   , ThisTask, sph_in->ID, T_gas, n_nuclei, sph_out->Abund[sH2I], sph_out->Abund[sHI], sph_out->Abund[sHII]);
    fflush(stdout);
    
    printf("\nStop.(in chemistry_noneq_chem_implicit)\n");
    /* savepositions(All.SnapshotFileCount++); */
    
    endrun(); /* just debugging */
  }
  
  
  if (i_T > N_T-2) i_T = N_T-2;  /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  double c_T = (T_gas - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);
  
  for(int j=1; j<=21; j++){
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
  }
  
  k[30] = k_table[30][i_T] + c_T*(k_table[30][i_T+1] - k_table[30][i_T]);
  k[22] = k[23] = k[24] = 0.;

#ifdef IONIZATION_PHOTON
  k[22] = sph_state->k22;
  k[23] = k[24] = 0.;
#endif
  
#ifdef LYMANWERNER    /* with Photo-Dissociation */
  k[25] = 0.;
  k[29] = 0.;
  k[53] = 0.;

  double J21 = 0.;
#ifdef LW_RADIATION_LOCAL
  k[25] = sph_in->rate_Hminus_pd;
  k[29] = sph_in->rate_H2_diss;
#ifdef RYTR_NHD
  k[53] = sph_in->rate_HD_diss;
#else
  k[53] = 1.55/ 1.39 * k[29];      /* optically-thin limit found in Wolcott-Green & Haiman 2011 */
#endif
#endif // LW_RADIATION_LOCAL

#ifdef BACKGROUND_LW_RADIATION
	  k[25] += All.rate_Hminus_pd;

	  double k29_ext = All.rate_H2_diss;
#ifdef HD
	  double k53_ext = 1.55/1.39 * All.rate_H2_diss;
#endif
#ifdef SHIELD_LW
	  double f_print_H2self = sph_state->f_byH2self;
#else
	  double f_print_H2self = 1.0;
#endif
#ifdef SHIELD_LW
	  k29_ext *= sph_state->f_byH2self /* SphP[ithis].f_byHI*/;
#ifdef HD
	  k53_ext *= sph_state->f_byHDself * sph_state->f_byHI * sph_state->f_byH2;
#endif
#endif //SHIELD_LW

  k[29] += k29_ext;
#ifdef HD
	  k[53] += k53_ext;
#endif

	  if(FlagPrint)
	    printf("k25=%g %g k29=%g %g f=%g\n", k[25], All.rate_Hminus_pd, k[29], k29_ext, f_print_H2self);

#if defined(METALLICITY) && defined(DUST_COOLING)
#ifdef LW_RADIATION_LOCAL
	  if(FlagPrint)
	    printf("G0_local=%g global=%g\n", sph_in->G0, All.G0);

	  sph_state->G0 += All.G0;
#endif
#endif
  
#endif //BACKGROUND_LW_RADIATION
#endif // LYMANWERNER

#ifdef HD
  for(int j=31; j<=48; j++)
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif  

#ifdef THREEBODY
  for(int j=55; j<=58; j++)
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif

  /****************************************************************************************************************/
  /******************************** We evaluate continuum opacities here ******************************************/
  /****************************************************************************************************************/

  double l_core = sph_state->Jeans_Length;
  double T_dust = All.T_CMB;
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
  stellar_flux = sph_in->stellar_rad_flux;
#endif
  T_dust = dust_temperature(rho_cgs, T_gas, sph_in->Metallicity, l_core, stellar_flux, &opac, &edot_dust);
  sph_out->Tdust = T_dust;
  sph_state->dust_cooling_rate = edot_dust;

#ifdef RYTR_TAU_DUST
  sph_out->kappa_dust = opac.kappa_dust_Planck;
  if(sph_out->kappa_dust > 2.06 * sph_in->Metallicity)
    sph_out->kappa_dust = 2.06 * sph_in->Metallicity;
#endif

  tau_cont = opac.tau;
#endif
  
  sph_state->tau_cont = tau_cont;
  sph_state->f_reduce = 1./(1. + 3./4.*tau_cont*tau_cont);
#endif
  /****************************************************************************************************************/
  /****************************************************************************************************************/
  /****************************************************************************************************************/

#ifdef INVERSE_REACTIONS
  double n_crit = H_criticaldensity(nH, T_gas, nH*sph_out->Abund[sHI]*l_core, sph_out->Abund[sHI], sph_out->Abund[sELEC]);

  k_inverse[1] = k_table_inverse[1][i_T] + c_T*(k_table_inverse[1][i_T+1] - k_table_inverse[1][i_T]);
  k_inverse[2] = k_table_inverse[2][i_T] + c_T*(k_table_inverse[2][i_T+1] - k_table_inverse[2][i_T]);
  k_inverse[7] = k_table_inverse[7][i_T] + c_T*(k_table_inverse[7][i_T+1] - k_table_inverse[7][i_T]);
  k_inverse[8] = k_table_inverse[8][i_T] + c_T*(k_table_inverse[8][i_T+1] - k_table_inverse[8][i_T]);
  k_inverse[13] = k_table_inverse[13][i_T] + c_T*(k_table_inverse[13][i_T+1] - k_table_inverse[13][i_T]);
  k_inverse[16] = k_table_inverse[16][i_T] + c_T*(k_table_inverse[16][i_T+1] - k_table_inverse[16][i_T]);
  k_inverse[17] = k_table_inverse[17][i_T] + c_T*(k_table_inverse[17][i_T+1] - k_table_inverse[17][i_T]);
  k_inverse[18] = k_table_inverse[18][i_T] + c_T*(k_table_inverse[18][i_T+1] - k_table_inverse[18][i_T]);

  k_inverse[30] = k_table_inverse[30][i_T] + c_T*(k_table_inverse[30][i_T+1] - k_table_inverse[30][i_T]);
  k_inverse[55] = k_table_inverse[55][i_T] + c_T*(k_table_inverse[55][i_T+1] - k_table_inverse[55][i_T]);

  //k_inverse[1] *= nH / (n_crit + nH); //temporally
  k_inverse[18] *= nH / (n_crit + nH);

  //printf("tau=%g n_crit=%g\n", tau_cont, n_crit);
  
  k_inverse[2] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  k_inverse[7] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  //k_inverse[2] *= (1 - exp(-tau_cont));
  //k_inverse[7] *= (1 - exp(-tau_cont));
  k[25]        *= exp(-tau_cont);

  if(FlagPrint){
    printf("k2=%g k_inv=%g %g %g %g %g %g %g %g %g %g tau_cont=%g n=%g n_c=%g l_core=%g\n"
	   , k[2]
	   , k_inverse[1], k_inverse[2], k_inverse[7], k_inverse[8]
	   , k_inverse[13], k_inverse[16], k_inverse[17], k_inverse[18], k_inverse[30], k_inverse[55], tau_cont, n_nuclei, n_crit, l_core);
  }


#ifdef ID_CHEM_DEBUG
  if(FlagPrint){
    printf("k_inv=%g %g tau_cont=%g n=%g n_c=%g l_core=%g i_T=%d c_T=%g\n"
	   , k_inverse[2], k_inverse[7], tau_cont, n_nuclei, n_crit, l_core, i_T, c_T);
  }
#endif

#endif
  
#ifdef H2_FORMATION_ON_DUST
  k[51] = H2form_on_grains(rho_cgs, T_gas, T_dust, sph_in->Metallicity);
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

    k[11] = pow(10.0, log_gamma_cd) + pow(10.0, log_gamma_dt);

    if(FlagPrint)
      printf("k11=%g gamma=%g %g\n", k[11], log_gamma_cd, log_gamma_dt);
   

#ifdef INVERSE_REACTIONS
    double xk_L = k[11];
    double xk_H  = k_inverse[55];

    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    double a = 1.0/(1.0 + nH/rn_cr);

    if(a == 1.0)       k[11] = xk_L;
    else if(a == 0.0)  k[11] = xk_H;
    else               k[11] = pow(xk_H,1.0-a)*pow(xk_L,a);      
#endif
#else
    double xk_L  = 1.12e-10*exp(-7.035e+4/T_gas);
    double xk_H  = 6.5e-7/sqrt(T_gas)*exp(-5.2e+4/T_gas)*(1.0-exp(-6.e+3/T_gas));
    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    double a = 1.0/(1.0 + nH/rn_cr);
  
    if(a == 1.0)       k[11] = xk_L;
    else if(a == 0.0)  k[11] = xk_H;
    else               k[11] = pow(xk_H,1.0-a)*pow(xk_L,a);  
#endif
  }else
    k[11] = 0.;


#if defined(STROMGREN_TEST)
  for(int j=1; j<=21; j++){
#ifdef IONIZATION_PHOTON
    if(j == 2)  continue;
#endif
    k[j] = 0.;
  }

  k[30] = 0.;
#endif
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


void test_cooling_table(void){
  if(ThisTask != 0) return;

  double nH = 1.;

  const int nh = 3;
  const int ne = 3;
  const int nM = 4;
  const int nT = 800;
  const double Tmin = 10.;
  const double Tmax = 1.e8;
  double dln_T = log(Tmax/Tmin)/nT;
  
  for(int ih=0;ih<nh;ih++){
    double yH2I = pow(10.,-3-ih);
    for(int ie=0;ie<ne;ie++){
      double y_e = pow(10.,-3-ie);
      for(int j=0;j<nM;j++){
	double metallicity = pow(10., -j);
	char fname[100];
	sprintf(fname, "./cooling_table_z%d_H2%d_e%d.txt", j, 3+ih, 3+ie);

	FILE *fp = fopen(fname, "w");
	
	for(int i=0;i<nT;i++){
	  double T_gas = Tmin * exp(dln_T * i);
	  double edot[10];

	  double etot = radiative_cooling_d_nometal_test(nH, T_gas, metallicity, y_e, yH2I, 0.0, edot);

	  fprintf(fp, "%g %g %g %g %g %g %g %g %g %g %g %g\n", T_gas, fabs(etot), edot[0], edot[1], edot[2], edot[3], edot[4], edot[5], edot[6], edot[7], edot[8], edot[9]);
	}
    
	fclose(fp);
      }
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
		      + 96.3   *log10_T 
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
