#ifdef CHEMISTRY
#define DEBUG_CHEMISTRY

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
#ifdef GREY_SPECTRUM
#define nu_start   13.6
#define nu_end     15.0
double nu_max, nu_min;
#else
#define nu_start   0.74
#define nu_end     0.74e3
#endif


#ifdef TRACE_PART
#define N_TRACE 0
/*int TRACE_ID[N_TRACE] = {88536120}*/;
#endif

/*=== Variables ===<*/

/* ----- Grid spacings. */
static double dln_T;

static double piHI,piHeI,piHeII;

/* ----- Energy arrays  */
static double edot_ceHI,edot_ceHeI,edot_ceHeII;
static double edot_ciHI,edot_ciHeI,edot_ciHeII;
static double edot_ciHeIS,edot_reHI,edot_reHeII1;
static double edot_reHeII2,edot_reHeIII,edot_comp;
static double edot_brem,edot_H2,edot_radHI;
static double edot_radHeI,edot_radHeII,edot_CIE;
static double edot_H2form_gas,edot_compression;
static double edot_2ph;
static double edot_dust;

#ifdef HD
static double edot_HD;
#endif

static double k[k_MAX];
static double k_inverse[k_MAX];

struct sph_data_local{
#ifdef SHIELD_LW
  FLOAT f_byH2self;
#ifdef HD
  FLOAT f_byHDself;
#endif
  FLOAT f_byHI;
  FLOAT f_byH2;
#endif

#ifdef H2_OPACITY
  double EscapeProb;
#ifdef HD
  double EscapeProb_HD;
#endif
#endif
  double f_reduce;
};

static struct sph_data_local sph_temp;

double chem_e[15]; //unit in kJ/mol and mainly from KIDA data base
int FlagPrint=0;

/*=== Subroutines ===<*/
int radiative_cooling_d_nometal(int ithis, double rho_gas, double T_gas, double *rhoU_dot, double *t_cool);
int chemical_cooling(int mode,int ithis, double nH, double T_gas, double dt, double *y_new, double *y_old, double *rhoU_dot, double *t_cool);

int chemistry_noneq_exp_nometal(int solve_it, int ithis, double dt, double rho_gas, double T_gas, double *t_chem, double *rhoU, double *t_cool);
int chemistry_noneq_chem_implicit(int solve_it, int ithis, double dt, double rho_gas, double T_gas, double *rhoU_dot, double *t_cool, double *err);
int chemistry_noneq_energy_implicit(int ithis, double dt, double rho_gas, double T_gas, double *rhoU_chem, double *rhoU_rad);

static inline int ReactionRates_nometal(int itr, double nH, double *yi, double *rf_tot, int ithis);
static inline double compute_xmu(struct sph_particle_data sph_p);

double set_core_radius(int ithis, double T_gas);
double get_core_radius(int ithis);

#ifdef CHEM_DIFF_ANALYTIC
int (*switch_func[256])(); 

static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *rf_tot, int ithis);
int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HI(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HII(double nH,double *ni,double *rf,int ithis);
#ifdef HELIUM
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *rf,int ithis);
#endif
int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *rf,int ithis);
#ifdef H2_PLUS
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *rf,int ithis);
#endif
int ReactionRate_diff_nometal_HM(double nH,double *ni,double *rf,int ithis);
#ifdef HD
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_DII(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *rf,int ithis);
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *rf,int ithis);
#endif
int func_other(double nH,double *ni,double *rf,int ithis);
#endif

int compute_abundances_d(int mode, int ithis, double a_start, double a_end, int *st_flag, int flag_print);
int chem_tabulate(void);
int heap_sort(unsigned long n, double ra[]);
int compute_gamma(int ipart);

double (*H2_cooling)();

double H2_cooling_O00(int i , double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_Y06(int i , double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_GP12(int i, double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);
double H2_cooling_G15(int i , double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e);

#ifdef HD
double (*HD_cooling)();

double HD_cooling_O00(int ithis, double nH, double T_gas, double yHI, double yHD);
double HD_cooling_L05(int i, double nH, double T_gas, double yHI, double yHDI);
double HD_cooling_GP12(int i, double nH, double T_gas, double yHI, double yHDI);

#ifdef H2_OPACITY
double (*HD_line_escape_fraction)();
double fesc_HD_fit(double T_K, double xNc_HD);
#endif
#endif

double continuum_cooling_rate_thin(int ithis, double T_gas, double rho_cgs);

int H2line_tabulate();
#ifdef H2_OPACITY
double (*H2_line_escape_fraction)();

double H2line_reduction(int ithis, double T);
double fesc_fit(int ithis,double T_K,double xNc_H2);
double fesc_unity(int ithis, double T_K, double xNc_H2);
#endif

#ifdef CONT_OPACITY
double (*continuum_escape_fraction)();

double calc_continuum_opacities(int ithis, double T_gas, double rho_cgs, double l_core, struct opacities *opac);
double calc_continuum_opacities_CIE(int ithis, double T_gas, double rho_cgs, double l_core, struct opacities *opac);
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
double pf_z_H2p(double T_K);


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

double get_PlanckOpacity(double T_gas, double rho_gas);

#ifdef IONIZATION_PHOTON
void calc_ionization_rate(int ithis, double Teff, double *k22, double *gamma_ion, double dr, double ri);
#endif

#ifdef SHIELD_LW
double SelfShield(double ColumnDens, double T_gas, double b);
double EachShield(double ColumnDens, double alpha, double delta, double eta);
#endif

#ifdef ENERGY_IMPLICIT
double give_g(int ithis, double dt, double rho_gas, double T_gas, double *e_new, double *g, double *rhoU_chem, double *rhoU_rad);
double solve_brent(int ithis, double dt, double rho_gas,double T_min, double g_min, double T_max, double g_max,
								 double *e_new0, double *rhoU_chem, double *rhoU_rad, int max_iter);
int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
#endif

void set_reaction_rate_coeff(int ithis, double T_gas, double rho_gas, double *k);
void correct_initial_abundances(double *yi);
void give_initial_guess(double T_gas, double nH, double dt, double *yi, double *y_init, double *k);
void bound_ddy(double *yi, double *ddy, double *y_max, int itr, int flag_reset, double *err_max);


#ifdef CONTINUUM_COOLING_MATSUKOBA 
double continuum_cooling_rate(int ithis, double T_gas, double rho_cgs, double l_core);
double H_criticaldensity(double xnH, double T_K,double xNc_H, double y_H,double y_e);
double fesc_Lya(int ithis,double xnH,double T_K,double y_H,double y_e,double rr,double vel,double *cooling_rate_Lya,double *cooling_rate_2ph);
#endif

#ifdef SAHA
double chem_saha(double xnH,double T_K,double *y_H,double *y_H2,double *y_e,double *y_Hp);
double chem_saha_5sp(double xnH, double T_K, double *y);
#endif

double dust_cross_section(void);
double dust_temperature(double rho_cgs, double T_gas, double metallicity, double l_core, double stellar_flux, struct opacities *opac, double *edot_dust);
double H2form_on_grains(double rho_cgs, double T_gas, double T_dust, double metallicity);
double dust_PEheating(int ithis, double nH, double y_e, double T_gas, double metallicity);
double calc_G0(void);

#ifdef CONTINUUM_COOLING_MATSUKOBA
void read_opacity_tables(void);
void tabulate_cooling_rate_cont_thin(void);
#endif
#ifdef DUST_COOLING
void read_dust_opacity_tables(void);
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
	//RE31, RE34, 
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


int N_SP;
int N_RE;
const double XH  = 0.76;
const double XHe = 0.0625;
const double XHD = 3.e-5;

/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k_table[k_MAX][N_T];
double k_table_inverse[k_MAX][N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#if defined(IONIZATION_PHOTON)
double sigma[s_MAX][N_nu];
double SIGMA_e22, SIGMA_e23, SIGMA_e24;
#endif


double compute_temp(int ithis)
{
  double sum_fraction;
  double molecular_weight;
  double egyspec_to_temp;
  double a3inv, ThermalEnergy, T_gas;

  a3inv            = All.cf_a3inv;
  molecular_weight = compute_xmu(SphP[ithis]);
  egyspec_to_temp  = (SphP[ithis].Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;

  ThermalEnergy    = SphP[ithis].Entropy / (SphP[ithis].Gamma-1.0)  * pow(SphP[ithis].Density * a3inv, SphP[ithis].Gamma-1.0);
  T_gas            = (All.UnitEnergy_in_cgs/All.UnitMass_in_g * ThermalEnergy) * egyspec_to_temp;
  
  //  SphP[ithis].Temp = T_gas;

  return T_gas;
}

int compute_Nnuclei(int ithis){
  double a3inv;
  double rho_cgs;
  double n_nuclei;

  a3inv   = All.cf_a3inv; 
  rho_cgs = All.UnitDensity_in_cgs * SphP[ithis].Density * All.HubbleParam * All.HubbleParam * a3inv; /* physical gram / cm3 */
        
  /* here, mean molecular weight must be taken into account */
    
  n_nuclei = rho_cgs / PROTONMASS;   
  SphP[ithis].n_nuclei = n_nuclei;
	
  return 1;
}

int compute_gamma(int ithis)
{
  /* variable gamma : Equations (5)-(7) in Omukai-Nishi */
  int i=1;
  double y_HI,y_HII,y_HM,y_HeI,y_HeII,y_HeIII,y_elec,y_H2;
  double n_H; /* number density of hydrogen nuclei */
  double x,denominator,sum_y;
  double gamma_H2_minus1_inv;
  double gamma0;
  double a3inv, ThermalEnergy;
  double T_gas, egyspec_to_temp;
  double sum_fraction, molecular_weight;

  gamma0 = SphP[ithis].Gamma;

  a3inv            = All.cf_a3inv;
  molecular_weight = compute_xmu(SphP[ithis]);
  egyspec_to_temp  = (SphP[ithis].Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
  
  ThermalEnergy    = SphP[ithis].Entropy / (SphP[ithis].Gamma-1.0)  * pow(SphP[ithis].Density * a3inv, SphP[ithis].Gamma-1.0);
  T_gas            = (All.UnitEnergy_in_cgs/All.UnitMass_in_g * ThermalEnergy) * egyspec_to_temp;

  n_H=SphP[ithis].HI + SphP[ithis].HII + SphP[ithis].HM + 2*SphP[ithis].H2I;
  /* number fraction of hydrogen atoms  */
  
  y_HI   =SphP[ithis].HI   /n_H;
  y_HII  =SphP[ithis].HII  /n_H;
#ifdef HELIUM
  y_HeI  =SphP[ithis].HeI  /n_H;
  y_HeII =SphP[ithis].HeII /n_H;
  y_HeIII=SphP[ithis].HeIII/n_H;
#else
  y_HeI   = XHe;
  y_HeII  = 0.;
  y_HeIII = 0.;
#endif
  y_HM   =SphP[ithis].HM   /n_H;
  y_H2   =SphP[ithis].H2I  /n_H;
  y_elec =SphP[ithis].elec /n_H;
  
  x=6100.0/T_gas;
  
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
  
  sum_y= y_HI + y_HII + y_HM + y_HeI + y_HeII + y_HeIII + y_elec + y_H2;
  
  denominator=1.5*(y_HI + y_HII + y_HM + y_HeI + y_HeII + y_HeIII + y_elec) + gamma_H2_minus1_inv * y_H2;
  
  SphP[ithis].Gamma = 1.0+sum_y/denominator;

  if(SphP[ithis].Gamma > 1.8 || SphP[ithis].Gamma < 1.2){
    printf("Gamma odd![task%03d]ithis=%d(ID=%d) gamma_old=%g gamma=%g T_gas=%g n_nucle=%g rhoU=%g Ent=%g dtEnt=%g at (x=%g y=%g z=%g)\n"
	   , ThisTask,ithis,P[ithis].ID,gamma0,SphP[ithis].Gamma, T_gas, SphP[ithis].n_nuclei
	   , SphP[ithis].ThermalEnergy, SphP[ithis].Entropy, SphP[ithis].DtEntropy	
	   , P[ithis].Pos[0], P[ithis].Pos[1], P[ithis].Pos[2]);
#ifdef IONIZATION_PHOTON
    printf("SphP[ithis].tau_ion=%g, SphP[ithis].k22=%g, SphP[ithis].piHI=%g\n"
	   , SphP[ithis].tau_ion, SphP[ithis].k22, SphP[ithis].piHI);
#endif
    fflush(stdout);
    /* endrun(332); */
  }
    
  return(i);
}


static inline double compute_xmu(struct sph_particle_data sph_p){
  double sum_fraction;
  
  sum_fraction =  sph_p.elec 
    + sph_p.HI + sph_p.HII + sph_p.HM + sph_p.H2I
#ifdef H2_PLUS
    + sph_p.H2II
#endif
#ifdef HD
    + sph_p.DI + sph_p.DII + sph_p.HDI
#ifdef HD_PLUS
    + sph_p.HDII  
#endif
#ifdef D_MINUS
    + sph_p.DM
#endif
#endif
#ifdef HELIUM
    + sph_p.HeI + sph_p.HeII + sph_p.HeIII
#else
    + XHe
#endif
    ;
  
  return 1.0 / sum_fraction;
}

double compute_freefall_time(double rho_gas){
  return sqrt(3.0*M_PI/32.0/GRAVITY/rho_gas);
}

double set_core_radius(int ithis, double T_gas){
#ifdef DENSITY_GRAD
  return SphP[ithis].l_core;
#else  
#ifdef JEANS_LENGTH_TMAX
  if(T_gas > JEANS_LENGTH_TMAX)
    T_gas = JEANS_LENGTH_TMAX;
#endif

  SphP[ithis].Jeans_Length
    = 0.5 * M_PI *
    sqrt(BOLTZMANN * T_gas/ PROTONMASS / 1.2)/
    sqrt(GRAVITY * SphP[ithis].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
  
  return SphP[ithis].Jeans_Length;
#endif
}

double get_core_radius(int ithis){
#ifdef DENSITY_GRAD
  return SphP[ithis].l_core;
#else
  return SphP[ithis].Jeans_Length;
#endif
}



int compute_abundances_d(int mode, int ithis, double a_start, double a_end, int *flag_implicit0, int flag_print)
{
  /* *******************************************************************
     Variables:
     ********************************************************************  */
  
  /*    time variable */
  double a_now;
  
  /*    temperature, density in g/cm3, pressure in K/cm3  */
  double T_gas, rho_cgs, T_gas_in;
  
  double dt=0, dt_in=0, t_cool_rad=0, t_cool_chem=0;
  double dt_min=0, time1=0, time2=0, time3=0, t_min=0;
  
  int j, ifunc, advance_flag, iter;
  int conv_flag = 0;
  int flag_implicit=0;
  int flag_chem_implicit=0;

  double n_nuclei;
  double rhoU_dot_rad, rhoU_dot_chem;
  double del_rhoU_tot, del_rhoU_rad, del_rhoU_chem;
  double del_U_rad, del_U_tot;
      
  double a3inv, da;
  double egyspec_to_temp;
  double sum_fraction, molecular_weight;
  double t_chem;
  double err;

#ifdef SHIELD_LW
  double ColumnDens_H2, ColumnDens_HI, ColumnDens_HD;
#endif

  FlagPrint = flag_print;

  double dt_tot = (a_end - a_start) * All.UnitTime_in_s/All.HubbleParam;
  if(All.ComovingIntegrationOn)
    dt_tot /= (a_start * All.cf_hubble_a);
  
  /* *******************************************************************
     Initialization
     ******************************************************************* */

  a3inv   = All.cf_a3inv;
  rho_cgs = All.UnitDensity_in_cgs * SphP[ithis].Density * All.HubbleParam * All.HubbleParam * a3inv; /* physical gram / cm3 */
  
  /* here, mean molecular weight must be taken into account */
  
  n_nuclei     = rho_cgs/PROTONMASS;
  SphP[ithis].n_nuclei = n_nuclei;

  molecular_weight = compute_xmu(SphP[ithis]);
  
  /* unit conversion */
  egyspec_to_temp           = All.UnitEnergy_in_cgs/All.UnitMass_in_g * (SphP[ithis].Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
  SphP[ithis].ThermalEnergy = SphP[ithis].Entropy / (SphP[ithis].Gamma-1.0)  * pow(SphP[ithis].Density * a3inv, SphP[ithis].Gamma-1.0);
  T_gas                     = SphP[ithis].ThermalEnergy * egyspec_to_temp;  

  if(P[ithis].ID == 1415878)
    FlagPrint = 1;

  if(FlagPrint) 
    printf("[Task%d, ID%d] n=%g T_gas=%g yi=%g %g %g %g %g\n"
	   , ThisTask, P[ithis].ID, n_nuclei, T_gas
	   , SphP[ithis].elec, SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I, SphP[ithis].HM);

#ifdef ID_CHEM_DEBUG
  if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint) 
    printf("n=%g, T_gas=%g\n", n_nuclei, T_gas);
#else
  if(FlagPrint){
#ifdef DENSITY_GRAD
    double l_J = 0.5 * M_PI *
      sqrt(BOLTZMANN * T_gas/ PROTONMASS / 1.2)/
      sqrt(GRAVITY * SphP[ithis].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);

    printf("[Task%d] n=%g T_gas=%g yi=%g %g %g %g %g %g %g %g l_c=%g l_J=%g\n"
	   , ThisTask, n_nuclei, T_gas
	   , SphP[ithis].elec, SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I, SphP[ithis].HM
	   , SphP[ithis].DI, SphP[ithis].DII, SphP[ithis].HDI, SphP[ithis].l_core, l_J);
#endif
  }
#endif

  for(j=22; j<=29; j++)  k[j] = 0;
  
  piHI = piHeI = piHeII = 0;

#ifdef IONIZATION_PHOTON
  //double k22, gamma_ion;

  // commented out just for test14
  double alpha = SphP[ithis].k22 / n_nuclei / 2.53e-13;
  double beta = SphP[ithis].k22 * n_nuclei * dt_tot;
    
  if(beta > 1.e15 && alpha > 1.e7){
    //double fac =  1.e7/alpha;    
    //SphP[ithis].k22  *= fac;
    //SphP[ithis].piHI *= fac;
  }
  
  double x_neut_crit = 1.e-15;
  double T_temp = SphP[ithis].piHI * n_nuclei * x_neut_crit * dt_tot /BOLTZMANN;
  if(T_temp > 1.e6){
    double fac = 1.e6 / T_temp;
    SphP[ithis].k22 *= fac;
    SphP[ithis].piHI *= fac;
  }
  
  k[22] = SphP[ithis].k22;
  piHI  = SphP[ithis].piHI * n_nuclei * SphP[ithis].HI;

  //piHI here is only used when flag_implicit=0. otherwise redefined at function "give_g" 
  
  k[23] = k[24] = piHeI = piHeII = 0.;
#endif
  
#ifdef SHIELD_LW
#ifdef SHIELD_6RAY
  int i;
  ColumnDens_H2 = ColumnDens_HI = ColumnDens_HD = 0.;
  for(i=0; i<6; i++){
    ColumnDens_H2 += 1./SphP[ithis].ColumnDens_H2[i];
    ColumnDens_HI += 1./SphP[ithis].ColumnDens_HI[i];
#ifdef HD		
    ColumnDens_HD += 1./SphP[ithis].ColumnDens_HD[i];
#endif		
  }
  ColumnDens_H2 = 6./ColumnDens_H2;
  ColumnDens_HI = 6./ColumnDens_HI;
#ifdef HD	
  ColumnDens_HD = 6./ColumnDens_HD;
#endif
#else

  double l_core = set_core_radius(ithis, T_gas);

  ColumnDens_H2 = l_core * n_nuclei * SphP[ithis].H2I;
  ColumnDens_HI = l_core * n_nuclei * SphP[ithis].HI;
#ifdef HD	
  ColumnDens_HD = l_core * n_nuclei * SphP[ithis].HDI;
#endif
#endif
  //>--- [H2]
  // eq-12 in Wolcott-Green et al. (2011)
  //   - Dopper broadening is defined in Ahn&Shapiro(2007)
  //        b = 1.29e4*(T/atomic-weight)**0.5 [cm/s]  (Spitzer 1978)
  //          = 1.29e4*(T/2)**0.5 [cm/s] (for H2)
  //          = 9.12e5*(T/1e4)**0.5 [cm/s]
  sph_temp.f_byH2self = SelfShield(ColumnDens_H2, T_gas, 9.12*pow(T_gas/1.e4, 0.5));

  // eq-12 & Table-1 in Wolcott-Green&Haiman(2011)
  sph_temp.f_byHI     = EachShield(ColumnDens_HI, 1.49e-1, 1.62e+0, 2.85e23);
  
#ifdef HD	
  //>--- [HD]
  // eq-12 in Wolcott-Green et al. (2011)
  //   - Dopper broadening is defined in Ahn&Shapiro(2007)
  //        b = 1.29e4*(T/atomic-weight)**0.5 [cm/s]  (Spitzer 1978)
  //          = 1.29e4*(T/3)**0.5 [cm/s] for HD
  //          = 7.45e5*(T/1e4)**0.5 [cm/s]
  sph_temp.f_byHDself = SelfShield(ColumnDens_HD, T_gas, 7.45*pow(T_gas/1.e4, 0.5));    
#endif

  // eq-12 & Table-1 in Wolcott-Green&Haiman(2011)
  sph_temp.f_byH2     = EachShield(ColumnDens_H2, 5.20e-3, 2.38e-1, 2.34e19);
#endif


  {
    /* recalculate the t_chem & t_cool for any particles */
    /* just compute the timescales and return */
#if !defined (ENERGY_IMPLICIT) || defined (HYBRID_IMPLICIT)
    dt_in = 0.0;
    T_gas_in = T_gas;

#if defined (CHEMISTRY_IMPLICIT) && !defined (HYBRID_IMPLICIT)
    t_chem = MAX_REAL_NUMBER;
#else
    ifunc = chemistry_noneq_exp_nometal(advance_flag=0, ithis, dt_in, rho_cgs, T_gas_in, &t_chem, &rhoU_dot_chem, &t_cool_chem);
#endif
    ifunc = radiative_cooling_d_nometal(ithis, rho_cgs, T_gas_in,  &rhoU_dot_rad, &t_cool_rad);    
    if(n_nuclei < N_EQCHEM)
      SphP[ithis].t_cool = dmin(fabs(t_cool_rad),fabs(t_cool_chem));

    SphP[ithis].t_chem = t_chem;				     
#else
    SphP[ithis].t_cool = MAX_REAL_NUMBER;
    SphP[ithis].t_chem = MAX_REAL_NUMBER;      
#endif
    
    if(mode==0 || mode ==2)
      return (ithis);
  }
    

  a_now = a_start;    /*    set current time (redshift) */
    
  for(iter=0;iter<MAX_ITERATION;iter++){

    if(iter >= MAX_ITERATION - 10){
      printf("[iter%d]time=%g %g %g => dt=%g\n",iter, time1, time2, time3, dt);
      printf("a_start=%g,a_now=%g,a_end=%g\n", a_start, a_now, a_end);
      printf("n_nuclei=%g,Tgas=%g,flag=%d\n", n_nuclei, T_gas, flag_implicit);
    }

    
    /*---------  Set the internal time step for subcycling */
    time1 = All.Epsilon * fabs(SphP[ithis].t_cool);
    time2 = All.Epsilon * fabs(SphP[ithis].t_chem);
    time3 = dt_tot;
        
    /* take the minimum */ 
    dt = dmin(time1, dmin(time2,time3));

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
    
    T_gas_in = T_gas;
    
#if defined(ENERGY_IMPLICIT) || defined(HYBRID_IMPLICIT)
    if(flag_implicit){
      conv_flag = chemistry_noneq_energy_implicit(ithis, dt, rho_cgs, T_gas_in, &rhoU_dot_chem, &rhoU_dot_rad);		
      SphP[ithis].t_cool = MAX_REAL_NUMBER, SphP[ithis].t_chem = MAX_REAL_NUMBER;
    }else
#endif
      {
        struct sph_particle_data sph_old;
        sph_old = SphP[ithis];

#if defined(CHEMISTRY_IMPLICIT)|| defined(HYBRID_IMPLICIT)
	if(flag_chem_implicit){
	  ifunc = chemistry_noneq_chem_implicit(advance_flag=1, ithis, dt, rho_cgs, T_gas_in, &rhoU_dot_chem, &t_cool_chem, &err);    		
	  SphP[ithis].t_chem = MAX_REAL_NUMBER;       	
	}else
#endif
	  {
	    ifunc = chemistry_noneq_exp_nometal(advance_flag=1, ithis, dt, rho_cgs, T_gas_in, &t_chem,&rhoU_dot_chem, &t_cool_chem);
	    SphP[ithis].t_chem = t_chem;
	  }

#ifdef VARIABLEGAMMA
	ifunc = compute_gamma(ithis);   /* Update gamma */
#endif
	
	/* Update temperature */ 
	ifunc = radiative_cooling_d_nometal(ithis, rho_cgs, T_gas_in, &rhoU_dot_rad, &t_cool_rad);     

	SphP[ithis].t_cool  = dmin(fabs(t_cool_rad),fabs(t_cool_chem));
		
	del_rhoU_rad  = rhoU_dot_rad * dt; 		
	del_rhoU_chem = rhoU_dot_chem * dt; /* erg/cm3 */
		
	del_rhoU_tot  = del_rhoU_rad + del_rhoU_chem;		
	del_U_tot     = del_rhoU_tot / rho_cgs; /* delta U in erg/gram */		
	del_U_tot    /= (All.UnitEnergy_in_cgs/All.UnitMass_in_g);
	
	SphP[ithis].ThermalEnergy += del_U_tot;		

#ifdef HYBRID_IMPLICIT
        /************* DEBUGGING... [CHON200417] ************/
	if(0.1 * sph_old.ThermalEnergy < fabs(del_U_tot)){
          flag_implicit = 1;
          SphP[ithis] = sph_old;

          conv_flag = chemistry_noneq_energy_implicit(ithis, dt, rho_cgs, T_gas_in, &rhoU_dot_chem, &rhoU_dot_rad);
          SphP[ithis].t_cool = MAX_REAL_NUMBER, SphP[ithis].t_chem = MAX_REAL_NUMBER;
        }
#endif
      }


    if(All.ComovingIntegrationOn){
      /* get dt to be in the system units. Up to here, dt in seconds */
      dt /= All.UnitTime_in_s; /* system units */
      dt *= All.HubbleParam; /* h correction */
        
      a_now = a_now + a_now * All.cf_hubble_a * dt;      
    }else{
      dt /= All.UnitTime_in_s; /* system units */
      dt *= All.HubbleParam;   /* h correction */      
        
      a_now = a_now + dt;  
    }
            
    /* update T_gas */
    molecular_weight = compute_xmu(SphP[ithis]);
    egyspec_to_temp  = All.UnitEnergy_in_cgs/All.UnitMass_in_g * (SphP[ithis].Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
    T_gas            = SphP[ithis].ThermalEnergy * egyspec_to_temp;
           
#ifdef SET_TMAX
    if(T_gas > SET_TMAX)
      SphP[ithis].ThermalEnergy *= SET_TMAX / T_gas;
#endif

    /* stop if time = a_end */
    if(a_end - a_now < ((a_end-a_start)/MAX_ITERATION))
      goto LOOP_OUT;      
  }  /* iter */
   
  
 LOOP_OUT:;
    
  SphP[ithis].Entropy  = SphP[ithis].ThermalEnergy * (SphP[ithis].Gamma-1.0) / pow(SphP[ithis].Density * a3inv, SphP[ithis].Gamma-1.0);
  *flag_implicit0 = flag_implicit + 2*flag_chem_implicit;
    
  if(flag_implicit)
    return conv_flag;
  else
    return (iter);
}



/* ----- Compute cooling time of gas and
   ----- solve for thermal energy.  */

int radiative_cooling_d_nometal(int ithis, double rho_gas, double T_gas, double *rhoU_dot, double *t_cool)
{  
  /* ----- Collisional excitation coefs. */
  double  ceHI,ceHeI,ceHeII;
  
  /* ----- Collisional ionization coefs. */
  double  ciHI,ciHeI,ciHeIS,ciHeII;
  
  /* ----- Recombination coefs. */
  double  reHII,reHeII1,reHeII2,reHeIII;
  
  /* ----- Compton cooling/heating and bremstrahlung. */
  double  comp1=0,brem,cie;
  
  volatile int i_T;
  int i_energy;
  double  edot=0.;
  double  ln_T,c_T,ln_T0,ln_T9;
  
  /* molecular cooling rate */
  
  double  T_CMB; 
  double  rhoU, n_nuclei, a3inv;
  
  double  f_reduce;
  
  double edot_cont=0.;
		
#if defined(METAL_COOLING) || defined(METAL_COOLING_LW) 
  double edot_CI, edot_CII, edot_OI;
#endif
	
	
  ln_T0 = log(T_start);
  ln_T9 = log(T_end);
  
  /* compute total numberdensity in units of mh
     hence n_nuclei * HI gives the number density of HI  */
  
  n_nuclei  = rho_gas/PROTONMASS;  
  a3inv     = All.cf_a3inv;

  
  /* --------- Compute temperature index.  */
  ln_T = log(T_gas);
  if (ln_T < ln_T0) ln_T = ln_T0;
  if (ln_T > ln_T9) ln_T = ln_T9;
  
  i_T = (int)( (ln_T-ln_T0)/dln_T ) ;
  
  if (i_T > N_T-2) i_T = N_T - 2;     /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  c_T = (T_gas - T[i_T])/(T[i_T+1]-T[i_T]);
  
  
  /* --------- Interpolate temperature dependent coefficients. */
#ifdef HLINE_COOLING_THIN
  ceHI  = ceHIa[i_T]   + c_T*(ceHIa[i_T+1] -ceHIa[i_T]);  
  //ciHI  = ciHIa[i_T]   + c_T*(ciHIa[i_T+1] -ciHIa[i_T]);
#endif

#ifndef CONTINUUM_COOLING_MATSUKOBA  
  reHII = reHIIa[i_T]  + c_T*(reHIIa[i_T+1]-reHIIa[i_T]);
  brem  = brema[i_T]   + c_T*(brema[i_T+1] -brema[i_T]);
  cie   = H2_CIE[i_T]  + c_T*(H2_CIE[i_T+1]-H2_CIE[i_T]);
#endif

#ifdef HELIUM
  ceHeI  = ceHeIa[i_T]  + c_T*(ceHeIa[i_T+1] - ceHeIa[i_T]);
  ceHeII = ceHeIIa[i_T] + c_T*(ceHeIIa[i_T+1]- ceHeIIa[i_T]);

  ciHeI  = ciHeIa[i_T]  + c_T*(ciHeIa[i_T+1]  -ciHeIa[i_T]);
  ciHeII = ciHeIIa[i_T] + c_T*(ciHeIIa[i_T+1] -ciHeIIa[i_T]);
  ciHeIS = ciHeISa[i_T] + c_T*(ciHeISa[i_T+1] -ciHeISa[i_T]);

  reHeII1= reHeII1a[i_T] + c_T*(reHeII1a[i_T+1]-reHeII1a[i_T]);
  reHeII2= reHeII2a[i_T] + c_T*(reHeII2a[i_T+1]-reHeII2a[i_T]);
  reHeIII= reHeIIIa[i_T] + c_T*(reHeIIIa[i_T+1]-reHeIIIa[i_T]);
#endif 

  
  /* ========= Compton heating =========================================================*/  
  if(All.ComovingIntegrationOn){
    comp1 = 5.65e-36*pow(All.Time,-4);
  }else{
    comp1 = 5.65e-36*pow(All.ScaleFactor,-4);
  }
  
  edot_comp = -comp1*(T_gas-All.T_CMB)*SphP[ithis].elec*n_nuclei;
  
  /**********************************************************************************************/ 
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/ 
  /**********************************************************************************************/ 

  double l_core = set_core_radius(ithis, T_gas);

  /******************************* Hydrogen lines ***********************************************/
#ifdef HLINE_COOLING_THIN
  edot_ceHI   = -ceHI  *SphP[ithis].HI  *SphP[ithis].elec     *n_nuclei*n_nuclei;
#else
  if(n_nuclei < 1.e8){
    double v_thermal, fesc;
    double edot_Lya;

    v_thermal = sqrt(BOLTZMANN * T_gas/PROTONMASS);
    fesc = fesc_Lya(ithis, n_nuclei, T_gas, SphP[ithis].HI, SphP[ithis].elec, l_core, v_thermal, &edot_Lya, &edot_2ph);
    
    edot_ceHI = edot_Lya;
  }else
    edot_ceHI = 0.;
#endif

  /******************************* Helium lines **************************************************/
#ifdef HELIUM
  edot_ceHeI  = -ceHeI *SphP[ithis].HeII*SphP[ithis].elec*SphP[ithis].elec*n_nuclei*n_nuclei*n_nuclei;
  edot_ceHeII = -ceHeII*SphP[ithis].HeII*SphP[ithis].elec     *n_nuclei*n_nuclei;

  edot_ciHeI  = -ciHeI *SphP[ithis].HeI *SphP[ithis].elec     *n_nuclei*n_nuclei;
  edot_ciHeII = -ciHeII*SphP[ithis].HeII*SphP[ithis].elec     *n_nuclei*n_nuclei;
  edot_ciHeIS  = -ciHeIS *SphP[ithis].HeII*SphP[ithis].elec*SphP[ithis].elec*n_nuclei*n_nuclei*n_nuclei;
#endif      


  /******************************* H2 & HD lines ***********************************************/  
  double nH;
  double yHI, yHII, yH2I, y_e, yHeI, yHDI;

  nH  = n_nuclei * XH;
  yHI = SphP[ithis].HI    / XH;
  yHII = SphP[ithis].HII  / XH;
  yH2I = SphP[ithis].H2I  / XH;
  y_e  = SphP[ithis].elec / XH;  
#ifdef HELIUM
  yHeI = SphP[ithis].HeI  / XH;  
#else
  yHeI = XHe / XH;
#endif
#ifdef HD
  yHDI = SphP[ithis].HDI/XH;
#endif

#ifdef NO_H2COOLING
  edot_H2 = 0.;
#else
  edot_H2 = H2_cooling(ithis, nH, T_gas, yHI, yH2I, yHeI, yHII, y_e);
#endif

#ifdef HD  
  edot_HD = HD_cooling(ithis, nH, T_gas, yHI, yHDI);
#endif


#ifdef H2_OPACITY  /* accounting for the local optical depth to H2 lines   */ 
  if(n_nuclei > 1.0e4){
    double N_H2;
    N_H2 = l_core * n_nuclei * SphP[ithis].H2I;

    f_reduce = H2_line_escape_fraction(ithis, T_gas, N_H2);
    sph_temp.EscapeProb = f_reduce;
    
    if(n_nuclei >= N_EQCHEM) f_reduce = 0.;
    if(f_reduce > 1.) f_reduce = 1.;
    
    if(f_reduce < 0){ /* something wrong */
      fflush(stdout);
      endrun(367523);
    }
    
    edot_H2 *= f_reduce;

    if(FlagPrint)
      printf("f_reduce=%g\n", f_reduce);

#ifdef HD
    double N_HD;
    N_HD = l_core * n_nuclei * SphP[ithis].HDI;

    f_reduce = HD_line_escape_fraction(T_gas, N_HD);
    sph_temp.EscapeProb_HD = f_reduce;
    
    if(n_nuclei >= N_EQCHEM) f_reduce = 0.;
    if(f_reduce > 1.) f_reduce = 1.;
    
    if(f_reduce < 0){ /* something wrong */
      fflush(stdout);
      endrun(367523);
    }
    
    edot_HD *= f_reduce;
#endif

  } 
#endif

  /******************************* CII & OI lines ***********************************************/  
#ifdef METAL_COOLING_LW
  double f_metal;
  f_metal = SphP[ithis].Metallicity;

  double ne;
  double nOI, nCI, nCII;
  double cool_rate_OI, cool_rate_CI, cool_rate_CII;
  
  nOI  = n_nuclei * XH * Y_O_sun * f_metal;  
  nCI  = n_nuclei * XH * Y_C_sun * f_metal;  
  nCII = n_nuclei * XH * Y_C_sun * f_metal;  
  ne   = SphP[ithis].elec * n_nuclei;
  //nH   = SphP[ithis].HI   * n_nuclei;
  
  //cool_rate_OI = OxygenI_cooling(nOI, ne, nH, T_gas);
  //cool_rate_CII = CarbonII_cooling(nCII, ne, nH, T_gas);

#ifdef DUST_COOLING
  if(SphP[ithis].Tdust < 2000.){
    //nOI  *= 0.54;
    //nCI  *= 0.28;
    //nCII *= 0.28;
  }
#endif

  cool_rate_OI  = OxygenI_cooling(nOI, ne, nH, T_gas);
  cool_rate_CI  = CarbonI_cooling(nCI, ne, nH, T_gas);
  cool_rate_CII = CarbonII_cooling(nCII, ne, nH, T_gas);

  //edot_OI = -1.0 * n_nuclei * cool_rate_OI * A_O_sun * f_metal;
  edot_OI  = -1.0 * cool_rate_OI;
  edot_CI  = -1.0 * cool_rate_CI;
  edot_CII = -1.0 * cool_rate_CII;
#endif  

  /******************************************************************************************/ 
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/ 
  /******************************************************************************************/ 

  edot_ciHI = 0.;
#ifdef CONTINUUM_COOLING_MATSUKOBA  
  edot_reHI = edot_brem = edot_CIE = 0;
  edot_cont = -continuum_cooling_rate_thin(ithis, T_gas, rho_gas);
#else
  /* --------- Cooling by CIE  */
  /* --------- Bremsstrahlung  */
  /* --------- Recombinations  */
  edot_reHI    = -reHII*  SphP[ithis].HII * SphP[ithis].elec * n_nuclei*n_nuclei;
  edot_CIE     = -cie  *  SphP[ithis].H2I * SphP[ithis].H2I  * n_nuclei*n_nuclei;

#ifndef HELIUM
  edot_brem    = -brem * SphP[ithis].HII  *  SphP[ithis].elec*n_nuclei*n_nuclei;
#else
  edot_brem    = -brem*(SphP[ithis].HII+SphP[ithis].HeII+4.0*SphP[ithis].HeIII)  *  SphP[ithis].elec*n_nuclei*n_nuclei;
  edot_reHeII1 = -reHeII1*SphP[ithis].HeII *SphP[ithis].elec*n_nuclei*n_nuclei;
  edot_reHeII2 = -reHeII2*SphP[ithis].HeII *SphP[ithis].elec*n_nuclei*n_nuclei;
  edot_reHeIII = -reHeIII*SphP[ithis].HeIII*SphP[ithis].elec*n_nuclei*n_nuclei;
#endif  
#endif  

  /********************** opacity for the continuum process *******************************/
  double tau_cont;
#ifdef CONT_OPACITY
  tau_cont = SphP[ithis].tau_cont;
  f_reduce = sph_temp.f_reduce;

#ifdef DEBUG_1ZONE
  SphP[ithis].f_reduce = f_reduce;
#endif
#endif

#ifdef DUST_COOLING
  if(SphP[ithis].Tdust < 2000.){
    //edot_CII *= (1. - DUST_TO_GAS_MASS_RATIO);
    //edot_OI  *= (1. - DUST_TO_GAS_MASS_RATIO);
    //edot_OI  *= 0.54;
    //edot_CII *= 0.28;
  }

#ifdef LYMANWERNER
  double edot_PE;
  edot_PE = dust_PEheating(ithis, nH, y_e, T_gas, SphP[ithis].Metallicity);
#endif
  if(FlagPrint)
    printf("edot_dust=%g cont=%g PE=%g T_dust=%g\n", edot_dust, edot_cont, edot_PE, SphP[ithis].Tdust);
#endif


#ifdef CONT_OPACITY  
  edot_H2   *= exp(-tau_cont);
#ifdef HD
  edot_HD   *= exp(-tau_cont);
#endif
#ifdef METAL_COOLING_LW
  edot_CI   *= exp(-tau_cont);
  edot_CII  *= exp(-tau_cont);
  edot_OI   *= exp(-tau_cont);
#endif

  edot_cont *= f_reduce;
  edot_CIE  *= f_reduce;
#endif  

	       
  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  edot_radHI   = piHI  ;//*SphP[ithis].HI  *n_nuclei;
#ifdef HELIUM
  edot_radHeI  = piHeI *SphP[ithis].HeI *n_nuclei;
  edot_radHeII = piHeII*SphP[ithis].HeII*n_nuclei;
#endif

    
  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  edot = edot_ceHI + edot_ciHI + edot_reHI 
    + edot_comp + edot_brem + edot_H2 
#ifdef HELIUM
    + edot_ceHeI + edot_ceHeII 
    + edot_ciHeI + edot_ciHeII + edot_ciHeIS
    + edot_reHeII1 + edot_reHeII2 + edot_reHeIII
#endif
#ifdef HD
    + edot_HD
#endif
#ifdef METAL_COOLING_LW
    + edot_CI + edot_CII + edot_OI
#endif
#ifdef CALC_CONTINUUM_OPACITY
    + edot_CIE +  edot_2ph
#endif
#ifdef DUST_COOLING
    + edot_dust
#ifdef LYMANWERNER
    + edot_PE
#endif
#endif
    + edot_cont
    + edot_radHI + edot_radHeI + edot_radHeII;

#ifdef ID_CHEM_DEBUG
  if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint)
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont, f_reduce, SphP[ithis].Tdust);
#else
#ifdef THERMAL_FEEDBACK
  if(FlagPrint)
    printf("cool_rate:T_gas=%g n=%g ID=%d tot=%g rad=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g PE=%g cont=%g f_reduce=%g f_line=%g tau=%g Tdust=%g Trad=%g\n",
	   T_gas, n_nuclei, P[ithis].ID, edot, edot_radHI, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust, edot_PE, edot_cont, f_reduce, sph_temp.EscapeProb, tau_cont,
	   SphP[ithis].Tdust, pow(SphP[ithis].stellar_rad_flux / STEFAN / 4., 0.25));
#else
  /*if(FlagPrint)
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g\n",
    edot, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust, edot_cont, f_reduce, SphP[ithis].Tdust);*/
#endif
#endif

#if defined(DEBUG_1ZONE) || defined(FLAG_DEBUG_CHEMISTRY)
#ifdef HELIUM
  SphP[ithis].edot_ceHeI = edot_ceHeI;
  SphP[ithis].edot_ceHeII = edot_ceHeII;
  SphP[ithis].edot_ciHeI = edot_ciHeI;
  SphP[ithis].edot_ciHeII = edot_ciHeII;
  SphP[ithis].edot_ciHeIS = edot_ciHeIS;

  SphP[ithis].edot_reHeII1 = edot_reHeII1;
  SphP[ithis].edot_reHeII2 = edot_reHeII2;
  SphP[ithis].edot_reHeIII = edot_reHeIII;
#endif

  SphP[ithis].edot_ceHI = edot_ceHI;
  //SphP[ithis].edot_ciHI = edot_ciHI;		
  //SphP[ithis].edot_reHII = edot_reHI;
		
  SphP[ithis].edot_comp = edot_comp;
  //SphP[ithis].edot_brem = edot_brem;
  SphP[ithis].edot_cont = edot_cont;
		
  SphP[ithis].edot_H2 = edot_H2;
#ifdef HD
  SphP[ithis].edot_HD = edot_HD;
#endif
#ifdef CALC_CONTINUUM_OPACITY 
  SphP[ithis].edot_CIE  = edot_CIE;
  SphP[ithis].edot_cont = edot_cont;
#endif
  SphP[ithis].edot_2ph  = edot_2ph;
#ifdef METAL_COOLING_LW
  SphP[ithis].edot_CI  = edot_CI;
  SphP[ithis].edot_CII = edot_CII;
  SphP[ithis].edot_OI  = edot_OI;
#endif
#ifdef DUST_COOLING
  SphP[ithis].edot_dust = edot_dust;
#ifdef LYMANWERNER
  SphP[ithis].edot_PE = edot_PE;
#endif		
#endif		
  
  SphP[ithis].edot_rad = edot;
#endif		

  if(n_nuclei > N_NOLINE) /* clouds optically thick to atomic lines */
    edot = edot_CIE;
    
  *rhoU_dot = edot;
  
  rhoU = SphP[ithis].ThermalEnergy * SphP[ithis].Density * a3inv * All.HubbleParam * All.HubbleParam; //modified[0814] SphP[ithis].Entropy is not updated.
  
  rhoU *= (All.UnitEnergy_in_cgs/pow(All.UnitLength_in_cm, 3)); /* in erg/cm3 */
  
  *t_cool = (rhoU/edot); /* in seconds */
  
  if(fabs(edot) < 1.0e-40)
    *t_cool = MAX_REAL_NUMBER;
  
  /* for checking the gravitational heating */    
  edot_compression = SphP[ithis].DtEntropy / (SphP[ithis].Gamma-1.0) * pow(SphP[ithis].Density, SphP[ithis].Gamma-1.0);
  /* this is the rate of change of thermal energy per unit mass */
  
  edot_compression *= (All.UnitEnergy_in_cgs/All.UnitMass_in_g/All.UnitTime_in_s); /* erg/sec/gram */
  edot_compression *= rho_gas; /* erg/sec/cm^3 */

#if defined(DEBUG_1ZONE) || defined(FLAG_DEBUG_CHEMISTRY)
  SphP[ithis].tcool_rad = *t_cool;
#endif  
  
  i_energy = ithis;

  return i_energy;  
}


int chemical_cooling(int mode, int ithis, double nH, double T_gas, double dt, double *y_new, double *y_old, double *rhoU_dot, double *t_cool){
  double n_crit,n_elec,n_HI,n_HM,n_H2I,n_H2II;
  double n_HeI, n_HeII;
  double edot_H2form,edot_H2form_3body;
  double edot_ciHI, edot_ciHeI, edot_ciHeII;
  double edot_HMform;
  double edot_chem;

  enum SPECIES species;
  double dyi0[N_SP];

  /* calculate chemical heating/cooling */

  double a3inv;
  a3inv  = All.cf_a3inv;

  n_elec = nH * y_new[sELEC];
  n_HI   = nH * y_new[sHI];
  n_HM   = nH * y_new[sHM];
  n_H2I  = nH * y_new[sH2I];
  double n_HII  = nH * y_new[sHII];
#ifdef H2_PLUS
  n_H2II = nH * y_new[sH2II];
#endif
#ifdef HELIUM
  n_HeI  = nH * y_new[sHeI];
  n_HeII = nH * y_new[sHeII];
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
  
  n_crit = (1.e6 / sqrt(T_gas)) 
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

      double edot_H2form1;
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
    //commented out for test13
    edot_ciHI_EUV = 13.598 * eV_to_erg * SphP[ithis].k22 * n_HI;
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
  
  edot_chem = edot_H2form + edot_ciHI + edot_HMform;

  *rhoU_dot = edot_chem;
  if(dt < 1.e-30){
    *rhoU_dot = 0.;
    edot_H2form = edot_ciHI = edot_HMform = 0.;
  }

  double rhoU;
  rhoU = SphP[ithis].ThermalEnergy*SphP[ithis].Density * a3inv * All.HubbleParam * All.HubbleParam; 
  //modified[0814] in case for SphP[ithis].Entropy is not updated.

  rhoU *= (All.UnitEnergy_in_cgs/pow(All.UnitLength_in_cm, 3)); /* in erg/cm3 */
  
  *t_cool = (rhoU/edot_chem); /* in seconds */

#if defined(DEBUG_1ZONE) || defined(FLAG_DEBUG_CHEMISTRY)
  SphP[ithis].n_crit = n_crit;
  
  SphP[ithis].edot_chem    = edot_chem;
  SphP[ithis].edot_H2form  = edot_H2form;
  SphP[ithis].edot_HMform  = edot_HMform;
  SphP[ithis].edot_ciHI   = edot_ciHI;		  
  SphP[ithis].edot_H2form_3body  = edot_H2form_3body;
  SphP[ithis].edot_H2form_HM_loss  = edot_H2form_HM_loss;
  SphP[ithis].edot_H2form_dust_loss  = edot_H2form_dust_loss;

  SphP[ithis].tcool_chem  = *t_cool;

#ifdef HELIUM
  SphP[ithis].edot_ciHeI  = edot_ciHeI;
  SphP[ithis].edot_ciHeII = edot_ciHeII;
#endif
#endif

  return ithis;
}


#ifdef ENERGY_IMPLICIT
#define MAX_ITER 500
#define BETA_ENERGY 1.0


int chemistry_noneq_energy_implicit(int ithis, double dt, double rho_gas, double T_gas, double *rhoU_chem, double *rhoU_rad){
  int j;
  int iter, ITER_reset = 30;
  double denom=0., err, err_chem;
  double e_new;
  double T_gas_init, T_old, T_new, dT, d_T;
  const double T_max=3.e6, T_min=0.1;
  double T_large, T_small;
  double g, g0, array[1000];
  double a3inv;
  struct sph_particle_data sph_old;
  
  a3inv  = All.cf_a3inv;

#ifdef IONIZATION_PHOTON
  // commented out just for test12
  double x_neut_crit = 1.e-15;
  double T_temp = SphP[ithis].piHI * x_neut_crit * rho_gas / PROTONMASS * dt /BOLTZMANN;
  /*if(T_temp > 1.e6)
    SphP[ithis].piHI *= (1.e6 / T_temp);*/

  if(FlagPrint)
    printf("[2]piHI = %g, %g T_temp=%g\n", SphP[ithis].piHI, piHI, T_temp);
#endif


  iter = 0;
  sph_old    = SphP[ithis];
  T_gas_init =  T_gas;

  T_large = 5000, T_small = 1.5;
  
  T_old = T_gas_init, dT = 1.0e-3 * T_old;
  T_new = T_old;

  /******************* some special case?? *************************/
  //T_old = 3000.;
  //dT    = 0.001 * T_old;                                                                                                                           
  //T_new = T_old; 

  /*
  if(T_gas < 100.){
    T_old = 2000.;
    dT    = 0.001 * T_old;
    T_new = T_old;
#ifdef IONIZATION_PHOTON
  }else if(SphP[ithis].k22 > 1.e-3){
    T_old = 20000.;
    dT    = 0.001 * T_old;
    T_new = T_old;
#endif
  }else if(T_gas > 3.e4){
    T_old = 8000.;
    dT    = 0.001 * T_old;
    T_new = T_old;
  }else{
    T_old = T_gas_init, dT = 1.0e-3 * T_old;
    T_new =T_old;
  }
  
#ifdef IONIZATION_PHOTON
  T_large = 5000, T_small = 1e7;
#else
  T_large = 5000, T_small = 1e7;
#endif
  */

#ifdef IONIZATION_PHOTON
  // commented out just for test11
  if(T_temp > 1.e6){
    T_old = 5.e5;
    dT    = 0.5 * T_old;
    T_new = T_old;
  }else if(T_temp > 1.e5 ){
    T_old = 5.e5;
    dT    = 0.5 * T_old;
    T_new = T_old;
  }
#endif
  
  /***************************************************************/
  


  while(1){
    err_chem = give_g(ithis, dt, rho_gas, T_new, &e_new, &g0, rhoU_chem, rhoU_rad);
    
    if(err_chem < 1.e-2)
      break;
    
    if(isnan(dT) || isinf(dT)){
      dT = 0.1 * T_new;
      printf("(LINE=%d) dT=nan -> %g ID=%d g=%g g0=%g denom=%g\n",
	     __LINE__, dT, P[ithis].ID, g, g0, denom);
    }
    
    if(T_new + dT > T_max)
      T_new = dmax(T_large,T_min) , T_large /= 1.1;
    /*  else if(T_new + dT > 30000) T_new /= 1.1;*/
    else if(T_new + dT < 1.5)
      T_new = dmin(T_small,T_max), T_small *= 1.1;
    else T_new += dT, T_new *=1.001 ;
    
    SphP[ithis] = sph_old;

    if(isnan(T_new)){
      printf("[task%03d](iter=%d) ID=%d, T_gas is nan! Ent=%g, elec=%g, HI=%g, HII=%g, H2I=%g\n",
	     ThisTask, j, P[ithis].ID, SphP[ithis].ThermalEnergy,  SphP[ithis].elec, SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I);
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
  
  SphP[ithis] = sph_old;

  T_old = T_new;
  T_new = T_old + dT;

  int flag_converged;
  flag_converged = 0;
	
  for(j=iter;j<MAX_ITER;j++){    
#ifdef FLAG_BRENT
    if( flag_brent == 1){
      SphP[ithis] = sph_old;
      err = solve_brent(ithis, dt, rho_gas, T_minus, g_minus, T_plus, g_plus, &e_new, rhoU_chem, rhoU_rad, MAX_ITER + 20 - j);

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
      err_chem = give_g(ithis, dt, rho_gas, T_new, &e_new, &g, rhoU_chem, rhoU_rad);

      denom = g - g0;
      if(fabs(denom) < 1.e-10){
	if(denom > 0)
	  denom = 0.01;
	else
	  denom = -0.01;
      }
      
      dT = BETA_ENERGY*g*(T_old - T_new)/denom;
      
#ifdef ID_CHEM_DEBUG
      if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint)
	printf("[%d] T=%g dT=%g err_chem=%g\n", j, T_new, dT, err_chem);
#endif

      if(err_chem < 1e-4)
	break;      
      
      if(isnan(dT) || isinf(dT)){
	dT = 0.1 * T_new;
	printf("(LINE=%d) dT=nan -> %g ID=%d g=%g g0=%g denom=%g\n",
	       __LINE__, dT, P[ithis].ID, g, g0, denom);
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
      }else T_new *=1.001 ;
      
      if(j++ > MAX_ITER){
	printf("didn't converged T_new=%g dT=%g.\n"
	       , T_new, dT), err = 1.0;
	goto END;
      }
      
      SphP[ithis] = sph_old;
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

    SphP[ithis] = sph_old;
  }


#ifdef FLAG_BRENT			
  if(flag_converged == 0){
    T_new = T_gas_init;
		
    for(;j<MAX_ITER;j++){
      SphP[ithis] = sph_old;
	
      d_T = dmax(0.2 * T_new, 2000);
      if(sign > 0)
	T_new = dmax(T_new-d_T,T_new/10.);
      else
	T_new = dmin(T_new + d_T,T_new * 10.);

      err_chem = give_g(ithis, dt, rho_gas, T_new, &e_new, &g, rhoU_chem, rhoU_rad);
	
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
	       ThisTask, j, P[ithis].ID, SphP[ithis].ThermalEnergy,  SphP[ithis].elec, SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I);
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
      SphP[ithis] = sph_old;
      err = solve_brent(ithis, dt, rho_gas, T_minus, g_minus, T_plus, g_plus, &e_new, rhoU_chem, rhoU_rad, MAX_ITER + 20 - j);

      if(err < 1.e-5) flag_converged = 1;	       	
    }
  }
#endif

  SphP[ithis].ThermalEnergy = e_new;
  
 END:;
  if(isnan(err) || isnan(T_new))
    raise(11);

  if(err > 1.0e-4 || FlagPrint == 1){
    printf("not converged(in energy). err=%g id=%d T_old=%g T_new=%g (%g %g %g) Task=%d index=%d\n"
	   "ent=%g min_ent=%g dt_ent=%g , All.ComovingIntegrationOn=%d, All.Time=a_start0=%g, dt=%g, n_nuclei=%g, rho=%g, T_gas=%g, gamma=%g;\n"
	   "yi[0]=%g, yi[1]=%g, yi[2]=%g;\n"
	   "yi[%d]=%g, yi[%d]=%g;\n",	   
	   err, P[ithis].ID, T_gas_init, T_new,P[ithis].Pos[0],P[ithis].Pos[1],P[ithis].Pos[2],ThisTask, ithis,
	   sph_old.Entropy, All.MinEgySpec * 0.666666 / pow(SphP[ithis].Density * a3inv, 0.6666666),
	   SphP[ithis].DtEntropy, All.ComovingIntegrationOn, All.Time, dt, sph_old.n_nuclei, sph_old.Density, T_gas_init, sph_old.Gamma,
	   sph_old.elec, sph_old.HI, sph_old.HII, 
	   sH2I, sph_old.H2I, sHM, sph_old.HM);
#ifdef H2_PLUS
    printf("yi[%d]=%g;\n", sH2II, sph_old.H2II);
#endif
#ifdef HELIUM
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n"
	   , sHeI, sph_old.HeI, sHeII, sph_old.HeII, sHeIII, sph_old.HeIII);
#endif
#ifdef HD
    printf("yi[%d]=%g, yi[%d]=%g, yi[%d]=%g;\n",
	   sDI,   sph_old.DI, 
	   sDII,  sph_old.DII, 
	   sHDI,  sph_old.HDI);
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
    printf("SphP[ithis].rate_Hminus_pd=%g, SphP[ithis].rate_H2_diss=%g, SphP[ithis].J21=%g;\n"
	   ,SphP[ithis].rate_Hminus_pd,SphP[ithis].rate_H2_diss, SphP[ithis].J21);
#else
    printf("All.rate_Hminus_pd=%g, All.rate_H2_diss=%g, All.J21=%g;\n"
	   ,All.rate_Hminus_pd, All.rate_H2_diss, All.J21);
#endif
    
#ifdef SHIELD_6RAY
    //printf("J21=%g,fac=%g;\n",All.J21,SphP[ithis].J21_factor);
    printf("J21=%g;\n",All.J21);
    printf("SphP[ithis].ColumnDens_H2[0]=%g,SphP[ithis].ColumnDens_H2[1]=%g,SphP[ithis].ColumnDens_H2[2]=%g;\n"
	   "SphP[ithis].ColumnDens_H2[3]=%g,SphP[ithis].ColumnDens_H2[1]=%g,SphP[ithis].ColumnDens_H2[5]=%g;\n"
	   "SphP[ithis].ColumnDens_HI[0]=%g,SphP[ithis].ColumnDens_HI[1]=%g,SphP[ithis].ColumnDens_HI[2]=%g;\n"
	   "SphP[ithis].ColumnDens_HI[3]=%g,SphP[ithis].ColumnDens_HI[4]=%g,SphP[ithis].ColumnDens_HI[5]=%g;\n"
	   ,SphP[ithis].ColumnDens_H2[0],SphP[ithis].ColumnDens_H2[1],SphP[ithis].ColumnDens_H2[2]
	   ,SphP[ithis].ColumnDens_H2[3],SphP[ithis].ColumnDens_H2[4],SphP[ithis].ColumnDens_H2[5]
	   ,SphP[ithis].ColumnDens_HI[0],SphP[ithis].ColumnDens_HI[1],SphP[ithis].ColumnDens_HI[2]
	   ,SphP[ithis].ColumnDens_HI[3],SphP[ithis].ColumnDens_HI[4],SphP[ithis].ColumnDens_HI[5]);
#ifdef HD
    printf("SphP[ithis].ColumnDens_HD[0]=%g,SphP[ithis].ColumnDens_HD[1]=%g,SphP[ithis].ColumnDens_HD[2]=%g;\n"
	   "SphP[ithis].ColumnDens_HD[3]=%g,SphP[ithis].ColumnDens_HD[4]=%g,SphP[ithis].ColumnDens_HD[5]=%g;\n"
	   ,SphP[ithis].ColumnDens_HD[0],SphP[ithis].ColumnDens_HD[1],SphP[ithis].ColumnDens_HD[2]
	   ,SphP[ithis].ColumnDens_HD[3],SphP[ithis].ColumnDens_HD[4],SphP[ithis].ColumnDens_HD[5]);
#endif
#endif
#endif        

#ifdef DENSITY_GRAD
    printf("SphP[ithis].l_core=%g;\n", SphP[ithis].l_core);
#endif
    
#ifdef IONIZATION_PHOTON
    printf("SphP[ithis].tau_ion=%g, SphP[ithis].k22=%g, SphP[ithis].piHI=%g;\n"
	   , SphP[ithis].tau_ion, sph_old.k22, sph_old.piHI );
#endif
#ifdef SOBOLEV_APP
    printf("SphP[ithis].DelV[0]=%g, SphP[ithis].DelV[1]=%g, SphP[ithis].DelV[2]=%g;\n"
	   , SphP[ithis].DelV[0], SphP[ithis].DelV[1], SphP[ithis].DelV[2] );
#endif

    printf("yf[0]=%g, yf[1]=%g, yf[2]=%g, yf[3]=%g, yf[4]=%g;\n",
	   SphP[ithis].elec, SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I, SphP[ithis].HM);
  }

  if(err > 1.0e-4){
    SphP[ithis] = sph_old;
    return 1;
  }
  

  return 0;
}
#pragma opt


double give_g(int ithis, double dt, double rho_gas, double T_gas_new, double *e_new, double *g, double *rhoU_dot_chem, double *rhoU_dot_rad){
  int ifunc;
  int solve_it = 1;
  double e_old, de;
  double gmm_cmp_l;
  double gamma_new, xmu_new;
  double t_cool_chem, t_cool_rad;
  double err;
  
  e_old = SphP[ithis].ThermalEnergy;
  
  gmm_cmp_l = -SphP[ithis].DtEntropy / (SphP[ithis].Gamma-1.0) * pow(SphP[ithis].Density, SphP[ithis].Gamma-1.0);
  gmm_cmp_l *= (All.UnitEnergy_in_cgs/All.UnitMass_in_g/All.UnitTime_in_s); /* erg/sec/gram */
  gmm_cmp_l *= rho_gas; /* erg/sec/cm^3 */
  
  // < chemistry >
  chemistry_noneq_chem_implicit(solve_it, ithis, dt, rho_gas, T_gas_new, rhoU_dot_chem, &t_cool_chem, &err);

  ifunc = compute_gamma(ithis);
  radiative_cooling_d_nometal(ithis, rho_gas, T_gas_new, rhoU_dot_rad, &t_cool_rad);
  
  SphP[ithis].t_cool = dmin(fabs(t_cool_rad),fabs(t_cool_chem));
  
  gamma_new = SphP[ithis].Gamma;
  xmu_new = compute_xmu(SphP[ithis]);

  *e_new = BOLTZMANN * T_gas_new/( (gamma_new - 1.0)*xmu_new*PROTONMASS )/(All.UnitEnergy_in_cgs/All.UnitMass_in_g);
  de  = (*rhoU_dot_chem + *rhoU_dot_rad)/rho_gas *dt / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) ;

#ifdef IONIZATION_PHOTON
  double e_max = BOLTZMANN * 3.e6 /( (gamma_new - 1.0)*xmu_new*PROTONMASS )/(All.UnitEnergy_in_cgs/All.UnitMass_in_g); 
  if(e_old + de > e_max && de > 0.){
    de =  dmax(e_max - e_old, 0.);
  }
#endif
  
  *g  = *e_new - e_old - de;

  if(FlagPrint)
    printf("resid(e_new-e_old-de)=%g e_new=%g old=%g de=%g Tgas=%g rhoU_dot_chem=%g cool=%g\n"
	   , *g, *e_new, e_old, de, T_gas_new, *rhoU_dot_chem, *rhoU_dot_rad);
  
#ifdef ID_CHEM_DEBUG
  if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint)
    printf("resid=%g e_new=%g old=%g de=%g rhoU_dot_chem=%g cool=%g\n"
	   , *g, *e_new, e_old, de, *rhoU_dot_chem, *rhoU_dot_rad);
#endif


#ifdef DEBUG_CHEMISTRY
  if(err > 1.){
    if(*g > *e_new * 100){
      *g = *e_new * 100 - e_old;
      printf("[ID=%d] g=>%g \n", P[ithis].ID, *g);
    }
  }
#endif
 
  return err;
}


double solve_brent(int ithis, double dt, double rho_gas,double T_min, double g_min, double T_max, double g_max,
		   double *e_new0, double *rhoU_chem, double *rhoU_rad, int max_iter){
  int iter, mflag;
  double err = 1.0;
  double T_a, T_b, T_c, T_d;
  double f1, f2, f3;
  double T_cand, T_crit;
  double val_a, val_b, val_c, val_d;
  double val_new;
  double xm;
  double dummy1;
  double e_new;
  struct sph_particle_data sph_old;
  
  const double Delta = 1.0e-30;
  sph_old = SphP[ithis];
  
  if(g_min * g_max > 0){
    printf("something wrong in solve_brent."
	   " g_min and g_max must have opposite sign.\n");
    return 1e30;
  }
  
  if( fabs(g_min) < fabs(g_max) ){
    T_b = T_min, val_b = g_min;
    T_a = T_max, val_a = g_max;
  }else{
    T_a = T_min, val_a = g_min;
    T_b = T_max, val_b = g_max;
  }
  
  T_c = T_a, val_c = val_a;
  
  iter = 0, mflag = 1;
  while(iter < max_iter && err > 1.0e-6){
    
    f1 = val_a - val_b;
    
    if(val_a != val_c && val_b != val_c){
      f2 = val_b - val_c, f3= val_c - val_a;
      T_cand = -( T_a *val_b*val_c/f3/f2 + T_b *val_c*val_a/f1/f3
		  + T_c *val_a*val_b/f2/f1 );
    }else{
      T_cand = T_b - val_b*(T_a - T_b)/f1;
    }
    
    T_crit = (3*T_a + T_b) * 0.25;
    if(mflag == 1) xm = fabs(T_b - T_c);
    else xm = fabs(T_c - T_d);
    
    mflag = 0;
    if(  ((T_cand - T_crit) * (T_cand - T_b) > 0) && 
	 ( (fabs(T_cand - T_b) > 0.5 * xm) || xm < Delta )  ){
      
      T_cand = 0.5 * (T_a + T_b);
      mflag = 1;
    }
    
    T_cand = 0.5 * (T_a + T_b);
    
    err = give_g(ithis,dt,rho_gas,T_cand,&e_new,&val_new,rhoU_chem,rhoU_rad);
    //  val_new = func(T_cand);

    T_d = T_c, T_c = T_b;
    val_d = val_c, val_c = val_b;
    if(val_a * val_new < 0) T_b = T_cand, val_b = val_new;
    else T_a = T_cand, val_a = val_new;
    
    //    printf("%g %g %d\n",T_cand,val_new,mflag);
    
    if( fabs(val_a) < fabs(val_b) ){
      dummy1 = val_a, val_a = val_b, val_b = dummy1;
      dummy1 = T_a, T_a = T_b, T_b = dummy1;
    }
      
    err = fabs( (T_a - T_b)/T_b );
    
    if( err < 1.0e-5) break;
    SphP[ithis] = sph_old;    
    

    iter++;
  }
	
  *e_new0 = e_new;
	
  return err;
}


#endif




int chemistry_noneq_exp_nometal(int solve_it, int ithis, double dt, double rho_gas, double T_gas,
			    double *t_chem, double *rhoU_dot_chem, double *t_cool_chem)
{
  int j, mode;
  enum SPECIES species;

  double Ccoef,Dcoef,correction;
  double HIp,HIIp,HeIp,HeIIp,HeIIIp,elecp;
  double H2Ip,H2IIp,HMp;
	
  double t_min, elec_dot, HI_dot, H2I_dot;	
  double n_nuclei;

#ifdef HD
  double DIp, DIIp, HDIp, HDIIp, DMp;
#endif

  n_nuclei = rho_gas/PROTONMASS;

  set_reaction_rate_coeff(ithis, T_gas, rho_gas, k);
  k[22] /= (SphP[ithis].HI * n_nuclei);  
  //we here redefine k[22] just for the convinience
  
  double y_old[N_SP], y_new[N_SP];
  
  y_old[sELEC]  = SphP[ithis].elec/XH;
  y_old[sHI]    = SphP[ithis].HI/XH;
  y_old[sHII]   = SphP[ithis].HII/XH;
#ifdef HELIUM
  y_old[sHeI]   = SphP[ithis].HeI/XH;
  y_old[sHeII]  = SphP[ithis].HeII/XH;
  y_old[sHeIII] = SphP[ithis].HeIII/XH;
#endif
  
  y_old[sH2I]   = SphP[ithis].H2I/XH;
#ifdef H2_PLUS
  y_old[sH2II]  = SphP[ithis].H2II/XH;
#endif
  y_old[sHM]    = SphP[ithis].HM/XH;
#ifdef HD
  y_old[sDI]    = SphP[ithis].DI/XH;
  y_old[sDII]   = SphP[ithis].DII/XH;
  y_old[sHDI]   = SphP[ithis].HDI/XH;
#ifdef HD_PLUS
  y_old[sHDII]  = SphP[ithis].HDII/XH;
#endif
#ifdef DMINUS
  y_old[sDM]    = SphP[ithis].DM/XH;
#endif
#endif


  if (!solve_it){
    elec_dot = (k[1]  *SphP[ithis].HI   
#ifdef HELIUM
		+ k[3]  *SphP[ithis].HeI  
		+ k[5]  *SphP[ithis].HeII 
		- k[4]  *SphP[ithis].HeII 
		- k[6]  *SphP[ithis].HeIII
#endif
		- k[2]  *SphP[ithis].HII  
		)*SphP[ithis].elec*n_nuclei*n_nuclei
			
      + (k[22] *SphP[ithis].HI 
#ifdef HELIUM  
	 + k[23] *SphP[ithis].HeII
	 + k[24] *SphP[ithis].HeI
#endif
	 + k[25] *SphP[ithis].HM  )*n_nuclei
      + (k[8]  *SphP[ithis].HI  * SphP[ithis].HM  
	 + k[16] *SphP[ithis].HM  * SphP[ithis].elec
	 + k[17] *SphP[ithis].HI  * SphP[ithis].HM  
#ifdef H2_PLUS
	 + k[19] *SphP[ithis].HM  * SphP[ithis].HII
#endif
#ifdef DMINUS
	 + k[48] *SphP[ithis].DM  * SphP[ithis].HI //k[58] originally
#endif
	 )*n_nuclei*n_nuclei
		
      - (k[7]  *SphP[ithis].HI  
#ifdef H2_PLUS
	 + k[20] *SphP[ithis].H2II
#endif
#ifdef HD
	 //+ k[31] *SphP[ithis].DII 
#ifdef HD_PLUS
	 + k[42] *SphP[ithis].HDII
#endif
#ifdef DMINUS
	 + k[43] *SphP[ithis].DI 
#endif
#endif
	 ) *n_nuclei*SphP[ithis].elec*n_nuclei
      ;


    HI_dot =(k[2]  * SphP[ithis].HII  * SphP[ithis].elec 
	     + k[11] * SphP[ithis].H2I  * SphP[ithis].HI   *2.
	     + k[12] * SphP[ithis].H2I  * SphP[ithis].HII  
	     + k[13] * SphP[ithis].H2I  * SphP[ithis].elec *2.
#ifdef HELIUM
	     + k[15] * SphP[ithis].HeI  * SphP[ithis].HII    
#endif
	     )*n_nuclei*n_nuclei

#ifdef THREEBODY
      + k[58] * SphP[ithis].H2I  * SphP[ithis].H2I  * n_nuclei*n_nuclei *2
#endif
      + k[25] * SphP[ithis].HM   * n_nuclei
#ifdef H2_PLUS
      + k[26] * SphP[ithis].H2II * n_nuclei
#endif
      + k[29] * SphP[ithis].H2I  * n_nuclei * 2
      -(
	(k[1]  * SphP[ithis].elec
#ifdef H2_PLUS
	 + k[10] * SphP[ithis].H2II
#endif
#ifdef HELIUM
	 + k[14] * SphP[ithis].HeII
#endif
	 ) * n_nuclei
#ifdef THREEBODY
	+ (k[55] * SphP[ithis].HI  * 2.
	   +  k[57] * SphP[ithis].H2I * 2. )* SphP[ithis].HI * n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
	+ 2.0 * k[51]
#endif
	+ k[22]
	) * SphP[ithis].HI * n_nuclei;
		
		
    H2I_dot = (k[8]  *SphP[ithis].HM  *SphP[ithis].HI  
#ifdef H2_PLUS
	       + k[10] *SphP[ithis].H2II*SphP[ithis].HI  
	       + k[21] *SphP[ithis].H2II*SphP[ithis].HM 
#endif
	       )*n_nuclei*n_nuclei
#ifdef H2_FORMATION_ON_DUST
      + k[51] *SphP[ithis].HI * n_nuclei
#endif
      
#ifdef THREEBODY
      + (k[55] *SphP[ithis].HI  *SphP[ithis].HI  *SphP[ithis].HI  *n_nuclei
	 + k[57] *SphP[ithis].HI  *SphP[ithis].HI  *SphP[ithis].H2I *n_nuclei)*n_nuclei*n_nuclei
#endif 
      - (k[11] *SphP[ithis].HI   * n_nuclei
	 + k[12] *SphP[ithis].HII  * n_nuclei
	 + k[13] *SphP[ithis].elec * n_nuclei
#ifdef THREEBODY
	 + k[58] *SphP[ithis].H2I * n_nuclei
#endif
	 + k[27] 
	 + k[29]
	 ) * SphP[ithis].H2I * n_nuclei;
		
    t_min = dmin(fabs(SphP[ithis].elec * n_nuclei/elec_dot), fabs(SphP[ithis].HI * n_nuclei/HI_dot));
    if(n_nuclei > 1. || n_nuclei>N_THREEBODY)  //removed[0815], modified[1009] 
      t_min = dmin(t_min, fabs(SphP[ithis].H2I * n_nuclei/H2I_dot));
				
    *t_chem = t_min;

    chemical_cooling(mode=0, ithis, n_nuclei*XH, T_gas, dt, y_old, y_old, rhoU_dot_chem, t_cool_chem);  
  }else{  
    /* ----------- HI */  
    //added reaction 16,17,18,20,21[chon140112]

    Ccoef = (k[2]  * SphP[ithis].HII  * SphP[ithis].elec 
	     + k[11] * SphP[ithis].H2I  * SphP[ithis].HI   *2.
	     + k[12] * SphP[ithis].H2I  * SphP[ithis].HII  
	     + k[13] * SphP[ithis].H2I  * SphP[ithis].elec *2.
#ifdef HELIUM
	     + k[15] * SphP[ithis].HeI  * SphP[ithis].HII
#endif
	     + k[16] * SphP[ithis].HM * SphP[ithis].elec
	     + k[17] * SphP[ithis].HM * SphP[ithis].HI
	     + 2 * k[18] * SphP[ithis].HM * SphP[ithis].HII
#ifdef H2_PLUS
	     + k[20] * SphP[ithis].H2II * SphP[ithis].elec *2.
	     + k[21] * SphP[ithis].H2II * SphP[ithis].HM 
#endif
	     )*n_nuclei*n_nuclei
#ifdef THREEBODY
      + k[58] * SphP[ithis].H2I  * SphP[ithis].H2I  * n_nuclei*n_nuclei *2.
#endif
      + k[25] * SphP[ithis].HM   * n_nuclei
#ifdef H2_PLUS
      + k[26] * SphP[ithis].H2II * n_nuclei
#endif
      + k[29] * SphP[ithis].H2I  * n_nuclei * 2.;

    
    Dcoef = (  k[1]  * SphP[ithis].elec
#ifdef H2_PLUS
	       + k[10] * SphP[ithis].H2II
#endif
#ifdef HELIUM
	       + k[14] * SphP[ithis].HeII 
#endif
	       ) * n_nuclei
#ifdef THREEBODY
      + (k[55] * SphP[ithis].HI  * 2.
	 +  k[57] * SphP[ithis].H2I * 2. )* SphP[ithis].HI * n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
      + 2. * k[51]
#endif
      + k[22];


    HI_dot = fabs(Ccoef - Dcoef);  
    HIp   = (Ccoef*dt + SphP[ithis].HI*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* ----------- HII */
    Ccoef = k[1]  * HIp * SphP[ithis].elec * n_nuclei*n_nuclei
#ifdef H2_PLUS
      + k[10] * HIp * SphP[ithis].H2II * n_nuclei*n_nuclei 
#endif
#ifdef HELIUM
      + k[14] * HIp * SphP[ithis].HeII * n_nuclei*n_nuclei
#endif
      + k[22] * HIp * n_nuclei
#ifdef H2_PLUS
      + k[26] * SphP[ithis].H2II * n_nuclei
      + k[28] * SphP[ithis].H2II * n_nuclei * 2
#endif
      /*added by chon[131221] */
#ifdef HD
      + k[33] * HIp * SphP[ithis].DII * n_nuclei * n_nuclei
      + k[37] * SphP[ithis].H2I * SphP[ithis].DII * n_nuclei * n_nuclei
#ifdef HD_PLUS
      + k[36] * HIp * SphP[ithis].HDII * n_nuclei * n_nuclei
      + k[52] * SphP[ithis].HDII * n_nuclei
#endif
#endif
      /*[131221]*/
      ;

    Dcoef = (k[2]  * SphP[ithis].elec
#ifdef H2_PLUS
	     + k[9]  * HIp             
	     + k[12] * SphP[ithis].H2I
#endif
#ifdef HELIUM 
	     + k[15] * SphP[ithis].HeI
#endif
	     )* n_nuclei
      /*added by chon[131221] */
      + k[18] * SphP[ithis].HM * n_nuclei
#ifdef H2_PLUS
      + k[19] * SphP[ithis].HM * n_nuclei
#endif
#ifdef HD
      + (k[32] * SphP[ithis].DI 
	 + k[39] * SphP[ithis].HDI
#ifdef HD_PLUS
	 + k[40] * SphP[ithis].DI
#endif
#ifdef DMINUS
	 + k[45] * SphP[ithis].DM 
#endif
	 ) *n_nuclei /* modified SphP[ithis].HDI -> SphP[ithis].DM [0926]*/ 
#endif
      /*[131221]*/
      ;

    
    HIIp  = (Ccoef*dt + SphP[ithis].HII*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
    
#ifdef HELIUM
    /* ----------- HeI */
    Ccoef = (k[4]  * SphP[ithis].HeII * SphP[ithis].elec
	     + k[14] * HIp              * SphP[ithis].HeII)*n_nuclei*n_nuclei;
  
    Dcoef = k[3]  * SphP[ithis].elec * n_nuclei 
      + k[15] * HIIp             * n_nuclei
      + k[24];
  
    HeIp  = (Ccoef*dt + SphP[ithis].HeI*n_nuclei) / ( 1. + Dcoef*dt ) / n_nuclei;


    /* ----------- HeII */
    Ccoef = k[3]  * HeIp              * SphP[ithis].elec * n_nuclei*n_nuclei
      + k[6]  * SphP[ithis].HeIII * SphP[ithis].elec * n_nuclei*n_nuclei
      + k[15] * HeIp              * HIIp             * n_nuclei*n_nuclei
      + k[24] * HeIp*n_nuclei;

    Dcoef = k[4]  * SphP[ithis].elec * n_nuclei
      + k[5]  * SphP[ithis].elec * n_nuclei
      + k[14] * HIp              * n_nuclei
      + k[23]; 

    HeIIp  = (Ccoef*dt + SphP[ithis].HeII*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

      
    /* ----------- HeIII */
    Ccoef = k[5]  * HeIIp * SphP[ithis].elec * n_nuclei*n_nuclei
      + k[23] * HeIIp * n_nuclei;
      
    Dcoef = k[6] * SphP[ithis].elec * n_nuclei;
    
    HeIIIp = (Ccoef*dt + SphP[ithis].HeIII*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    if(HeIp   < 0) HeIp   = 0.;
    if(HeIIp  < 0) HeIIp  = 0.;
    if(HeIIIp < 0) HeIIIp = 0.;

    SphP[ithis].HeI   = HeIp;
    SphP[ithis].HeII  = HeIIp;
    SphP[ithis].HeIII = HeIIIp;
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
    elecp = (Ccoef*dt + SphP[ithis].elec*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    if(HIp   < 0) HIp   = 0;
    if(HIIp  < 0) HIIp  = 0;
    if(elecp < 0) elecp = 0;
    

    SphP[ithis].HI   = HIp;
    SphP[ithis].HII  = HIIp;    
    SphP[ithis].elec = elecp;


    /* count the number of photo-ionization */
    /* SphP[ithis].N_ion_cum += k[24] * HIp  * dt;  */
    

    /* --------- H2 */
    /* [SH131015 : crazyyyyyyy] */
  
    /* HM fast reactions */
    /*
      denom = (k[8] + k[17]) * SphP[ithis].HI 
      + (k[18]+ k[19]) * SphP[ithis].HII 
      +  k[16]         * SphP[ithis].elec 
      +  k[25]/n_nuclei;
  
      SphP[ithis].HM = (k[7]*SphP[ithis].HI*SphP[ithis].elec) / denom;
    */
    Ccoef = k[7]*SphP[ithis].HI*SphP[ithis].elec * n_nuclei*n_nuclei
#ifdef HD
#ifdef DMINUS
      +k[47] * HIp * SphP[ithis].DM *n_nuclei*n_nuclei
#endif
#endif
      ;

    Dcoef = ((k[8] + k[17]) * SphP[ithis].HI 
	     + k[18] * SphP[ithis].HII 
#ifdef H2_PLUS
	     + k[19] * SphP[ithis].HII 
#endif
	     + k[16] * SphP[ithis].elec  ) * n_nuclei
      +  k[25]
#ifdef HD
#ifdef DMINUS
      +  k[46] * SphP[ithis].DI *n_nuclei
#endif
#endif
      ;
  
    HMp = (Ccoef*dt + SphP[ithis].HM*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* H2II fast reactions */
    /*  denom = k[10]*SphP[ithis].HI 
	+ k[20]*SphP[ithis].elec 
	+ k[21]*SphP[ithis].HM 
	+(k[26]+k[28])/n_nuclei;
    
	SphP[ithis].H2II = (  k[9] * SphP[ithis].HI  * SphP[ithis].HII
	+ k[12] * SphP[ithis].H2I * SphP[ithis].HII
	+ k[19] * SphP[ithis].HM  * SphP[ithis].HII
	+ k[27] * SphP[ithis].H2I/n_nuclei         )/ denom; */

    /* noneq H2II added by chon.[131228]*/

#ifdef H2_PLUS
    Ccoef = ( k[9] * SphP[ithis].HI  * SphP[ithis].HII
	      +k[12] * SphP[ithis].H2I * SphP[ithis].HII 
	      +k[19] * HMp  * SphP[ithis].HII )*n_nuclei *n_nuclei
      + k[27]* SphP[ithis].H2I*n_nuclei;

    Dcoef = (   k[10]*SphP[ithis].HI 
		+ k[20]*SphP[ithis].elec
		+ k[21]*HMp ) *n_nuclei
      +k[26] + k[28];

    H2IIp = (Ccoef*dt + SphP[ithis].H2II*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif


    Ccoef = ( k[8]  * HMp   * SphP[ithis].HI
#ifdef H2_PLUS
	      + k[10] * H2IIp * SphP[ithis].HI
	      + k[21] * H2IIp * HMp
#endif
	      ) *n_nuclei*n_nuclei
#ifdef HD
      + (k[38]*HIp*SphP[ithis].HDI
	 + k[39]*HIIp*SphP[ithis].HDI) *n_nuclei*n_nuclei
#endif
#ifdef THREEBODY
      +( k[55] * SphP[ithis].HI   * SphP[ithis].HI * SphP[ithis].HI
	 + k[57] * SphP[ithis].HI   * SphP[ithis].HI * SphP[ithis].H2I )*n_nuclei*n_nuclei*n_nuclei
#endif
#ifdef H2_FORMATION_ON_DUST
      + k[51] * SphP[ithis].HI*n_nuclei
#endif
      ;
    
    Dcoef = (k[11] * SphP[ithis].HI 
#ifdef H2_PLUS
	     + k[12] * HIIp 
#endif
	     + k[13] * SphP[ithis].elec)*n_nuclei
#ifdef HD
      + (k[35]*SphP[ithis].DI
	 +k[37]*SphP[ithis].DII) *n_nuclei
#endif
#ifdef THREEBODY
      + k[58] * SphP[ithis].H2I * n_nuclei
#endif
      + k[27] 
      + k[29];

    H2I_dot = fabs(Ccoef - Dcoef);  
    H2Ip = (Ccoef*dt + SphP[ithis].H2I*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;

    if(HMp < 0) HMp = 0;
#ifdef H2_PLUS
    if(H2IIp < 0) H2IIp = 0;
#endif
    if(H2Ip < 0) H2Ip = 0;

    SphP[ithis].HM = HMp;
#ifdef H2_PLUS
    SphP[ithis].H2II = H2IIp;
#endif
    SphP[ithis].H2I = H2Ip;
  
  
#ifdef HD

    DIp=DIIp=HDIp=HDIIp=DMp=0;

    /* ----------- DI */
    Ccoef  = //k[31] * SphP[ithis].DII * SphP[ithis].elec * n_nuclei*n_nuclei
      k[33] * SphP[ithis].DII * SphP[ithis].HI   * n_nuclei*n_nuclei
      + k[38] * SphP[ithis].HDI * SphP[ithis].HI   * n_nuclei*n_nuclei
#ifdef HD_PLUS
      + k[42] * SphP[ithis].HDII* SphP[ithis].elec * n_nuclei*n_nuclei
      + k[52] * SphP[ithis].HDII* n_nuclei
#endif
#ifdef DMINUS
      + k[44] * SphP[ithis].DII * SphP[ithis].DM   * n_nuclei*n_nuclei * 2.0 
      + k[45] * SphP[ithis].HII * SphP[ithis].DM   * n_nuclei*n_nuclei
      + k[47] * SphP[ithis].HI  * SphP[ithis].DM   * n_nuclei*n_nuclei
#endif
      + k[53] * SphP[ithis].HDI * n_nuclei;
	     
    Dcoef = ( k[32] * SphP[ithis].HII
	      //+ k[34] * SphP[ithis].HI 
	      + k[35] * SphP[ithis].H2I
#ifdef HD_PLUS
	      + k[40] * SphP[ithis].HII
#endif
#ifdef DMINUS
	      + k[43] * SphP[ithis].elec
#endif
	      )* n_nuclei
      + k[50];

    DIp   = (Ccoef*dt + SphP[ithis].DI*n_nuclei)  / (1. + Dcoef*dt) / n_nuclei;
  
  
    /* ----------- DII */
    Ccoef = k[32] * DIp             * SphP[ithis].HII * n_nuclei*n_nuclei
      + k[39] * SphP[ithis].HDI * SphP[ithis].HII * n_nuclei*n_nuclei
      + k[50] * DIp             * n_nuclei
#ifdef HD_PLUS
      + k[51] * SphP[ithis].HDII* n_nuclei
#endif
      ;

    Dcoef = ( //k[31] * SphP[ithis].elec
	      k[33] * SphP[ithis].HI  
	      + k[37] * SphP[ithis].H2I 
#ifdef DMINUS
	      + k[41] * SphP[ithis].HI  
	      + k[44] * SphP[ithis].DM  
#endif
	      ) * n_nuclei;
	    
    DIIp  = (Ccoef*dt + SphP[ithis].DII*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;


    /* ----------- HDI */
    Ccoef = ( //k[34] * DIp            * SphP[ithis].HI  
	      k[35] * DIp            * SphP[ithis].H2I 
	      + k[37] * DIIp           * SphP[ithis].H2I 
#ifdef HD_PLUS
	      + k[36] * SphP[ithis].HI * SphP[ithis].HDII
#endif	      
#ifdef DMINUS
	      + k[48] * SphP[ithis].HI * SphP[ithis].DM  
#endif
	      )* n_nuclei*n_nuclei;

    Dcoef = k[38] * SphP[ithis].HI  * n_nuclei
      + k[39] * SphP[ithis].HII * n_nuclei
      + k[53];

    HDIp  = (Ccoef * dt + SphP[ithis].HDI * n_nuclei) / (1. +  Dcoef * dt) / n_nuclei;
    
#ifdef HD_PLUS
    /* ----------- HDII */
    Ccoef = ( k[40] * DIp  * SphP[ithis].HII
	      + k[41] * DIIp * SphP[ithis].HI ) * n_nuclei*n_nuclei;
      
    Dcoef = k[36] * SphP[ithis].HI   * n_nuclei
      + k[42] * SphP[ithis].elec * n_nuclei
      + k[51]
      + k[52];

    HDIIp = (Ccoef*dt + SphP[ithis].HDII*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
    
#ifdef DMINUS    
    /* ----------- DM */
    Ccoef = ( k[43] * DIp  * SphP[ithis].elec
	      + k[46] * DIp  * SphP[ithis].HM  ) * n_nuclei*n_nuclei;
      
    Dcoef = ( k[44] * DIIp           
	      + k[45] * SphP[ithis].HII
	      + k[47] * SphP[ithis].HI 
	      + k[48] * SphP[ithis].HI  ) * n_nuclei;
    
    DMp   = (Ccoef*dt + SphP[ithis].DM*n_nuclei) / (1. + Dcoef*dt) / n_nuclei;
#endif
      
    if(DIp    < 0) DIp   =0;
    if(DIIp   < 0) DIIp  =0;
    if(HDIp   < 0) HDIp  =0;
    SphP[ithis].DI    = DIp;
    SphP[ithis].DII   = DIIp;
    SphP[ithis].HDI   = HDIp;   

#ifdef HD_PLUS
    if(HDIIp  < 0) HDIIp =0;
    SphP[ithis].HDII  = HDIIp;
#endif
#ifdef DMINUS
    if(DMp    < 0) DMp   =0;
    SphP[ithis].DM    = DMp;    
#endif
  
#endif 
    
    
   

    
    /* --------- Undertake mass conservation corrections. */
  
    //SphP[ithis].elec = 0.0;
  
    correction = H_number_fraction / (SphP[ithis].HI 
				      +   SphP[ithis].HII 
				      +   SphP[ithis].HM 
				      + 2*SphP[ithis].H2I
#ifdef H2_PLUS
				      + 2*SphP[ithis].H2II
#endif
#ifdef HD
				      +   SphP[ithis].HDI 
#ifdef HD_PLUS
				      +   SphP[ithis].HDII
#endif
#endif
				      );
  
    if(0)
      if(correction > 2.0 || correction < 0.5){
	printf("H! ThisTask %d ithis %d correction=%g\n"
	       , ThisTask, ithis, correction);
	fflush(stdout);
	endrun(1100);
      }

  
    SphP[ithis].HI  *= correction;
    SphP[ithis].HII *= correction;
    //SphP[ithis].elec+= SphP[ithis].HII;
    
    SphP[ithis].HM   *= correction;
    SphP[ithis].H2I  *= correction;
#ifdef H2_PLUS
    SphP[ithis].H2II *= correction;
#endif
    //SphP[ithis].elec  = SphP[ithis].elec - SphP[ithis].HM + SphP[ithis].H2II;
  

#ifdef HD
    correction = D_number_fraction / (SphP[ithis].DI 
				      + SphP[ithis].DII
				      + SphP[ithis].HDI
#ifdef HD_PLUS
				      + SphP[ithis].HDII
#endif
#ifdef DMINUS
				      + SphP[ithis].DM
#endif
				      );

    if(0)
      if(correction > 2.0 || correction < 0.5){
	printf("D! ThisTask %d ithis %d correction=%g\n"
	       , ThisTask, ithis, correction);
	fflush(stdout);
	endrun(1101);
      } 

    SphP[ithis].DI   *= correction;
    SphP[ithis].DII  *= correction;
    SphP[ithis].HDI  *= correction;
#ifdef HD_PLUS
    SphP[ithis].HDII *= correction;
#endif
#ifdef DMINUS
    SphP[ithis].DM   *= correction;
#endif
    //SphP[ithis].elec  = SphP[ithis].elec + SphP[ithis].DII +  SphP[ithis].HDII - SphP[ithis].DM;
#endif 
  
  
#ifdef HELIUM
    correction = He_number_fraction / (SphP[ithis].HeI 
				       + SphP[ithis].HeII 
				       + SphP[ithis].HeIII);

    if(0)
      if(correction > 1.1 || correction < 0.9){
	printf("He! ThisTask %d ithis %d correction=%g %g %g %g\n"
	       , ThisTask, ithis, correction, SphP[ithis].HeI, SphP[ithis].HeII, SphP[ithis].HeIII);
	fflush(stdout);
	endrun(1102);
      }
  
    SphP[ithis].HeI   *= correction;
    SphP[ithis].HeII  *= correction;
    SphP[ithis].HeIII *= correction;
    //SphP[ithis].elec  = SphP[ithis].elec + SphP[ithis].HeII + 2.*SphP[ithis].HeIII;
#endif

    /*[SH130423: elec calculation]*/

    SphP[ithis].elec = 0.0;

    SphP[ithis].elec += SphP[ithis].HII;

    SphP[ithis].elec  = SphP[ithis].elec - SphP[ithis].HM;

#ifdef H2_PLUS
    SphP[ithis].elec += SphP[ithis].H2II;
#endif

#ifdef HD
    SphP[ithis].elec  += SphP[ithis].DII; 
#ifdef HD_plus
    SphP[ithis].elec  += SphP[ithis].HDII;
#endif
#ifdef DMINUS
    SphP[ithis].elec  -= SphP[ithis].DM;
#endif
#endif

#ifdef HELIUM
    SphP[ithis].elec  = SphP[ithis].elec + SphP[ithis].HeII + 2.*SphP[ithis].HeIII;
#endif

    if(isnan(SphP[ithis].HI))
      printf("[task%03d]HI is nan (ithis=%d, dt=%g, rho_gas=%g, T_gas=%g, mode=%d)\n", ThisTask, ithis, dt, rho_gas, T_gas, solve_it);
    if(isnan(T_gas))
      printf("[task%03d]line=%d T_gas is nan (ithis=%d, dt=%g, rho_gas=%g, T_gas=%g, mode=%d)\n", ThisTask, __LINE__, ithis, dt, rho_gas, T_gas, solve_it);

    /* added[0813] */
    t_min = dmin(fabs(SphP[ithis].elec * n_nuclei/elec_dot), fabs(SphP[ithis].HI * n_nuclei/HI_dot));
    t_min = dmin(t_min, fabs(SphP[ithis].H2I * n_nuclei/H2I_dot));

    *t_chem = t_min;


    y_new[sELEC]  = SphP[ithis].elec/XH;
    y_new[sHI]    = SphP[ithis].HI/XH;
    y_new[sHII]   = SphP[ithis].HII/XH;
#ifdef HELIUM
    y_new[sHeI]   = SphP[ithis].HeI/XH;
    y_new[sHeII]  = SphP[ithis].HeII/XH;
    y_new[sHeIII] = SphP[ithis].HeIII/XH;
#endif
    y_new[sH2I]   = SphP[ithis].H2I/XH;
#ifdef H2_PLUS
    y_new[sH2II]  = SphP[ithis].H2II/XH;
#endif
    y_new[sHM]    = SphP[ithis].HM/XH;
#ifdef HD
    y_new[sDI]    = SphP[ithis].DI/XH;
    y_new[sDII]   = SphP[ithis].DII/XH;
    y_new[sHDI]   = SphP[ithis].HDI/XH;
#ifdef HD_PLUS
    y_new[sHDII]  = SphP[ithis].HDII/XH;
#endif
#ifdef DMINUS
    y_new[sDM]    = SphP[ithis].DM/XH;
#endif
#endif
  
    chemical_cooling(mode=1, ithis, n_nuclei*XH, T_gas, dt, y_new, y_old, rhoU_dot_chem, t_cool_chem);  
  }

  return 1; 
}







#define MAX_ITER_CHEM 1000
#define BETA 1.0
 
int chemistry_noneq_chem_implicit(int solve_it, int ithis, double dt, double rho_gas, double T_gas, double *rhoU_dot, double *t_cool, double *err0)
{
  int j,n_sp;
  int inc, info;
  double n_nuclei,nH;
  double a3inv;
	
  enum SPECIES species;
  struct sph_particle_data sph_old;	

  sph_old = SphP[ithis];
  a3inv   = All.cf_a3inv;
  
  n_sp = N_SP;
  inc  = 1   ;

  n_nuclei  = rho_gas/PROTONMASS;
  nH = XH * n_nuclei;
    

#ifdef SAHA
  if(n_nuclei> N_EQCHEM){
    double y_H, y_H2, y_e, y_Hp;
    double y_H_init, y_H2_init, y_e_init, y_Hp_init, y_Hm_init;
    double y_input[N_SP];

    //printf("init saha...\n");

    y_H_init = SphP[ithis].HI   / XH;
    y_e_init = SphP[ithis].elec / XH;
    y_H2_init = SphP[ithis].H2I / XH;
    y_Hp_init = SphP[ithis].HII / XH;
    y_Hm_init = SphP[ithis].HM  / XH;

    /*
    double err;
       err = chem_saha(nH, T_gas, &y_H, &y_H2, &y_e, &y_Hp);
    

    SphP[ithis].HI   = y_H  * XH;
    SphP[ithis].HII  = y_Hp * XH;
    SphP[ithis].H2I  = y_H2 * XH;
    SphP[ithis].elec = y_e  * XH;    
    
    y_input[sHI]  = y_H;
    y_input[sHII] = y_Hp;
    y_input[sH2I] = y_H2;
    y_input[sHM]  = 0.;
    y_input[sELEC] = y_e;
    */
    
    double err;
    y_input[sELEC] = SphP[ithis].elec / XH;
    y_input[sHI]   = SphP[ithis].HI   / XH;
    y_input[sHII]  = SphP[ithis].HII  / XH;
    y_input[sH2I]  = SphP[ithis].H2I  / XH;
    y_input[sHM]   = SphP[ithis].HM   / XH;

    err = chem_saha_5sp(nH, T_gas, y_input);

    SphP[ithis].HI  = y_input[sHI]  * XH;
    SphP[ithis].HII = y_input[sHII] * XH;
    SphP[ithis].H2I = y_input[sH2I] * XH;
    SphP[ithis].HM  = y_input[sHM]  * XH;
    SphP[ithis].elec = y_input[sELEC] * XH;    
    
    double edot_H2form,edot_ionize,edot_HMform;
    
    edot_H2form = 4.48    * eV_to_erg * (y_input[sH2I] - y_H2_init) * nH / dt;    
    edot_ionize = -13.598 * eV_to_erg * (y_input[sHII] - y_Hp_init) * nH / dt;
    edot_HMform = 0.78    * eV_to_erg * (y_input[sHM]  - y_Hm_init) * nH / dt;

    if(dt < 1.e-30)
      edot_H2form = edot_ionize = edot_HMform = 0.;

    *rhoU_dot = edot_H2form + edot_ionize + edot_HMform;

#ifdef DEBUG_1ZONE
    SphP[ithis].edot_ciHI = 0.;
    SphP[ithis].edot_chem = edot_H2form;
#endif

    *err0 = err;

    return 1;
  }
#endif


#ifdef CHEM_DIFF_ANALYTIC
  int jj;

  for(jj = 0;jj < 256;jj++)
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
#endif  
  
  set_reaction_rate_coeff(ithis, T_gas, rho_gas, k);
  
  double yi[n_sp];
       
  yi[sELEC]  = SphP[ithis].elec/XH;
  yi[sHI]    = SphP[ithis].HI/XH;
  yi[sHII]   = SphP[ithis].HII/XH;
#ifdef HELIUM
  yi[sHeI]   = SphP[ithis].HeI/XH;
  yi[sHeII]  = SphP[ithis].HeII/XH;
  yi[sHeIII] = SphP[ithis].HeIII/XH;
#endif
  
  yi[sH2I]   = SphP[ithis].H2I/XH;
#ifdef H2_PLUS
  yi[sH2II]  = SphP[ithis].H2II/XH;
#endif
  yi[sHM]    = SphP[ithis].HM/XH;
#ifdef HD
  yi[sDI]    = SphP[ithis].DI/XH;
  yi[sDII]   = SphP[ithis].DII/XH;
  yi[sHDI]   = SphP[ithis].HDI/XH;
#ifdef HD_PLUS
  yi[sHDII]  = SphP[ithis].HDII/XH;
#endif
#ifdef DMINUS
  yi[sDM]    = SphP[ithis].DM/XH;
#endif
#endif
	
	
  int i_sp,j_sp,itr,ifunc;
  int err_flag = 0, ipiv[n_sp];
  double epsilon_y = 1.0e-10;
  double epsilon = 1.0e-4;
  double Dyi,drfdy[n_sp][n_sp],drf,dy[n_sp];
  double y_tmp[n_sp],y_init[n_sp], yi_init[n_sp];
  double y_mean[n_sp],y_max[n_sp],y_min[n_sp];
  double yi_ion;
  double rf_tot[n_sp],rf_f[n_sp],rf_b[n_sp],rf_big;
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

  int flag_min = 0;
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    y_mean[i_sp] /= XH;    
    y_init[i_sp] = yi[i_sp];
    
    if(yi[i_sp] < y_min[i_sp]/* || yi[i_sp] > y_max[i_sp]*/){
      flag_min = 1;
      
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
    
    ifunc=ReactionRates_nometal(itr,nH,yi,rf_tot,ithis);
    
    if(FlagPrint)
      printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g\n"
	     , itr, T_gas, err_max, err_flag, yi[sELEC], yi[sHI], yi[sHII], yi[sH2I], yi[sHM]);    
#ifdef ID_CHEM_DEBUG
    if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint){
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
      ifunc=ReactionRates_diff_nometal_fp(itr, j_sp, nH, yi, drfdy[j_sp], ithis);
#else
      if(fabs(yi[j_sp])<1.0e-100) Dyi = epsilon_y;
      else                       Dyi = epsilon*yi[j_sp];
      
      for(i_sp=0;i_sp<n_sp;i_sp++) {
	if(i_sp==j_sp) y_tmp[i_sp]=yi[i_sp]+0.5*Dyi;
	else           y_tmp[i_sp]=yi[i_sp];
      }
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,rf_f,ithis);
      
      y_tmp[j_sp]=yi[j_sp]-0.5*Dyi;
      ifunc=ReactionRates_nometal(itr,nH,y_tmp,rf_b,ithis);
      
      for(i_sp=0;i_sp<n_sp;i_sp++) {
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
    if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint){
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
			
      SphP[ithis].t_chem = t_chem1;
    }
  }else{
    if(err_max>1.e-4){
      printf("iteration is not converged err=%le(in %d),itr=%d,Ent=%g,min_Ent=%g "
	     "x=%g y=%g z=%g in id=%d.\n dt=%g, rho=%g, T_gas=%g, n_nuclei=%g;\n"
	     "y[%d]=%g, y[%d]=%g, y[%d]=%g, ",
	     err_max, err_flag, itr, SphP[ithis].Entropy, All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[ithis].Density * a3inv, GAMMA_MINUS1),
	     P[ithis].Pos[0], P[ithis].Pos[1], P[ithis].Pos[2], P[ithis].ID, dt, rho_gas, T_gas, SphP[ithis].n_nuclei,
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

    }else{
      t_chem1 = 1.0e20;
      for(i_sp=0;i_sp<n_sp;i_sp++) {
	if(yi[i_sp]-y_init[i_sp]!=0.0) {
	  t_chem2 = fabs((yi[i_sp]+y_init[i_sp])/2.0 / (yi[i_sp]-y_init[i_sp])) * dt;
	  if(t_chem2<t_chem1) t_chem1 = t_chem2;
	}
			
	SphP[ithis].t_chem = t_chem1;
      }
    }
  }
  
  *err0 = err_max;

  int mode;
  chemical_cooling(mode=1, ithis, nH, T_gas, dt, yi, y_init, rhoU_dot, t_cool);
  
  int isnan_flag = 0;

  for(i_sp=0;i_sp<n_sp;i_sp++){
    if(yi[i_sp] < y_min[i_sp])
      yi[i_sp] = y_min[i_sp];
  }

  SphP[ithis].elec =yi[sELEC] *XH;
  SphP[ithis].HI   =yi[sHI]   *XH;
  SphP[ithis].HII  =yi[sHII]  *XH;
#ifdef HELIUM
  SphP[ithis].HeI  =yi[sHeI]  *XH;
  SphP[ithis].HeII =yi[sHeII] *XH;
  SphP[ithis].HeIII=yi[sHeIII]*XH;
#endif
  
  SphP[ithis].H2I  =yi[sH2I]  *XH;
#ifdef H2_PLUS
  SphP[ithis].H2II =yi[sH2II] *XH;
#endif
  SphP[ithis].HM   =yi[sHM]   *XH;
#ifdef HD    
  SphP[ithis].DI   =yi[sDI]   *XH;
  SphP[ithis].DII  =yi[sDII]  *XH;
  SphP[ithis].HDI  =yi[sHDI]  *XH;
#ifdef HD_PLUS
  SphP[ithis].HDII =yi[sHDII] *XH;
#endif
#ifdef DMINUS
  SphP[ithis].DM   =yi[sDM]   *XH;
#endif
#endif

#ifdef CALC_CONTINUUM_OPACITY
  double t_ff;
  double ni[N_SP];

  int i;
  for(i=0;i<N_SP;i++){
    ni[i] = yi[i] * nH;
  }

  double l_core = get_core_radius(ithis);
  
  t_ff =  compute_freefall_time(SphP[ithis].Density*All.UnitDensity_in_cgs); 
  solve_level_population_H_new(T_gas, ni, l_core, t_ff);
#endif

#ifdef DEBUG_1ZONE
  int i;
  for(i=0;i<70;i++){
    SphP[ithis].k[i] = k[i];
    SphP[ithis].k_inv[i] = k_inverse[i];
  }
#endif

  return 1;
}



#pragma optimize("", on)

#undef MAX_ITER_CHEM


static inline int ReactionRates_nometal(int itr, double nH, double *yi, double *rf_tot, int ithis)
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
  
  rate = k[8]*ni[sHM]*ni[sHI];
  rf[RE8][sHM] = -rate;
  rf[RE8][sHI] = -rate;
  rf[RE8][sH2I] =  rate;
  rf[RE8][sELEC] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[HM]=%g ni[sHI]=%g\n", RE8, rf[RE8][sH2I], k[8], ni[sHM], ni[sHI]);     
  
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

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sHI]=%g\n", RE11, rf[RE11][sH2I], k[11], ni[sH2I], ni[sHI]);     
  
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

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sELEC]=%g\n", RE13, rf[RE13][sH2I], k[13], ni[sH2I], ni[sELEC]);     
  
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
  // commented out just for test2a
  double k22  = SphP[ithis].k22  * ni[sHI];
  piHI = SphP[ithis].piHI * ni[sHI];

  if(FlagPrint)
    printf("[3]piHI = %g, %g ire=%d rate=%g\n", SphP[ithis].piHI, piHI, RE22, k22);  

  rate = k22;// *ni[sHI];
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
  /*
  rate = k[31]*ni[sDII]*ni[sELEC];
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

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[HI]=%g\n", RE55, rf[RE55][sH2I], k[55], ni[sHI]);     
  
  rate = k[57]*ni[sHI]*ni[sHI]*ni[sH2I];
  rf[RE57][sHI] = -2.0*rate;
  rf[RE57][sH2I] =  rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g ni[sHI]=%g\n", RE57, rf[RE57][sH2I], k[57], ni[sHI], ni[sH2I]);     
  
  rate = k[58]*ni[sH2I]*ni[sH2I];
  rf[RE58][sH2I] = -rate;
  rf[RE58][sHI] = 2.0*rate;

  //if(FlagPrint)
  //printf("rf[%d]=%g k=%g ni[H2I]=%g\n", RE58, rf[RE58][sH2I], k[58], ni[sH2I]);     
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

  rate = k_inverse[2]*ni[sHI];
  rf[RE2_INV][sHI] = -rate;
  rf[RE2_INV][sHII] = rate;
  rf[RE2_INV][sELEC] = rate;

  rate = k_inverse[7]*ni[sHM];
  rf[RE7_INV][sHM] = -rate;
  rf[RE7_INV][sHI] = rate;
  rf[RE7_INV][sELEC] = rate;

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

  //if(FlagPrint)
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

      if(FlagPrint && i_sp == sELEC && fabs(rf[i_re][i_sp]) > 1.e-100)
	printf("rf[%d]=%g tot=%g\n", i_re, rf[i_re][i_sp], rf_tot[i_sp]);     
    }
  
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    rf_tot[i_sp] /= nH;
  }
  // rf_tot in unit of dy/dt [/sec]
  
  return 1;
}



#ifdef CHEM_DIFF_ANALYTIC
static inline int ReactionRates_diff_nometal_fp(int itr, int i_sp_diff, double nH, double *yi, double *rf_tot, int ithis)
{
  int i_sp, n_sp = N_SP, ifunc;
  double ni[n_sp];
	
  for(i_sp=0;i_sp<n_sp;i_sp++) {
    ni[i_sp]=yi[i_sp]*nH;
  }
  	
  for(i_sp=0;i_sp<n_sp;i_sp++) 
    rf_tot[i_sp]=0.0;
  
  ifunc = switch_func[i_sp_diff](nH, ni, rf_tot, ithis);

  return 1;
}

int ReactionRate_diff_nometal_Elec(double nH,double *ni,double *rf,int ithis){
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
  rate = k[31]*ni[sDII];
  rf[sDII]  += -rate;
  rf[sELEC] += -rate;
  rf[sDI]   +=  rate;

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

int ReactionRate_diff_nometal_HI(double nH,double *ni,double *rf,int ithis){
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
  // commented out just for test2a
  rate = SphP[ithis].k22;
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

int ReactionRate_diff_nometal_HII(double nH,double *ni,double *rf,int ithis){
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
int ReactionRate_diff_nometal_HeI(double nH,double *ni,double *rf,int ithis){
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

int ReactionRate_diff_nometal_HeII(double nH,double *ni,double *rf,int ithis){
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

int ReactionRate_diff_nometal_HeIII(double nH,double *ni,double *rf,int ithis){
  double rate;
  enum SPECIES species;

  rate = k[6]*ni[sELEC];
  rf[sELEC]  += -rate;
  rf[sHeII]  +=  rate;
  rf[sHeIII] += -rate;

  return 0;
}
#endif

int ReactionRate_diff_nometal_H2I(double nH,double *ni,double *rf,int ithis){
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
int ReactionRate_diff_nometal_H2II(double nH,double *ni,double *rf,int ithis){
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

int ReactionRate_diff_nometal_HM(double nH,double *ni,double *rf,int ithis){
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
int ReactionRate_diff_nometal_DI(double nH,double *ni,double *rf,int ithis){
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

int ReactionRate_diff_nometal_DII(double nH,double *ni,double *rf,int ithis){
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

int ReactionRate_diff_nometal_HDI(double nH,double *ni,double *rf,int ithis){
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
int ReactionRate_diff_nometal_HDII(double nH,double *ni,double *rf,int ithis){
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
int ReactionRate_diff_nometal_DM(double nH,double *ni,double *rf,int ithis){
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

int func_other(double nH,double *ni,double *rf,int ithis){
  return 1;
}
#endif



void set_reaction_rate_coeff(int ithis, double T_gas, double rho_cgs, double *k){
  //volatile int i_T;
  int i_T;
  int j;
  double ln_T,c_T,ln_T0,ln_T9;
  double log_T, log_T2, log_T3;
    
  double n_nuclei;  

  n_nuclei = rho_cgs / PROTONMASS;

  /*********** initialize ************/
  for(j=0; j<k_MAX; j++)
    k[j] = 0.0;
  /***********************************/

  /* ------- Compute temperature index. */
  ln_T0 = log(T_start);
  ln_T9 = log(T_end);
  
  log_T  = log10(T_gas);
  log_T2 = log_T * log_T;
  log_T3 = log_T * log_T * log_T;
  
  ln_T  = log(T_gas);
  if(ln_T < ln_T0) ln_T = ln_T0;
  if(ln_T > ln_T9) ln_T = ln_T9;
  i_T = (int) ((ln_T-ln_T0)/dln_T);

  /*i_T = 0;
    while(i_T*dln_T + ln_T0 < ln_T){
    i_T++;
    }
    if(i_T>0) i_T--;*/
  
  
  if(i_T<0) i_T=0;
  if(i_T>N_T-2) i_T=N_T-2;

  if(i_T < 0){
    printf("i_T=%d (ln_T=%g,ln_T0=%g,dln_T=%g)\n",i_T,ln_T,ln_T0,dln_T);

    printf("ThisTask %d ithis %d ID %d T %g n %g H2I %g HI %g HII %g\n"
	   , ThisTask, ithis, P[ithis].ID, T_gas, n_nuclei, SphP[ithis].H2I, SphP[ithis].HI, SphP[ithis].HII);
    printf("         POS=[%1.10e,%1.10e,%1.10e]\n"
	   , P[ithis].Pos[0], P[ithis].Pos[1], P[ithis].Pos[2]);  /*--- [SH120616] ---*/
    fflush(stdout);
    
    printf("\nStop.(in chemistry_noneq_chem_implicit)\n");
    /* savepositions(All.SnapshotFileCount++); */
    
    endrun(8765); /* just debugging */
  }
  
  
  if (i_T > N_T-2) i_T = N_T-2;  /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  c_T = (T_gas - T[i_T])/(T[i_T+1]-T[i_T]);
  
  for(j=1; j<=21; j++){
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
  }
  
  j=30;
  k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);

  k[22] = 0.;
#ifdef IONIZATION_PHOTON
  
  //commented out just for test2a
  k[22] = SphP[ithis].k22  * SphP[ithis].HI * n_nuclei;
  //piHI  = SphP[ithis].piHI * SphP[ithis].HI * n_nuclei;
  
  k[23] = k[24] = piHeI = piHeII = 0.;
#endif

  
#ifdef LYMANWERNER    /* with Photo-Dissociation */
  k[25] = 0.;
  k[29] = 0.;
  k[53] = 0.;

#ifdef LW_RADIATION_LOCAL
  k[25] = SphP[ithis].rate_Hminus_pd;
  k[29] = SphP[ithis].rate_H2_diss;
#ifdef RYTR_NHD
  k[53] = SphP[ithis].rate_HD_diss;
#else
  k[53] = 1.55/ 1.39 * k[29];      /* optically-thin limit found in Wolcott-Green & Haiman 2011 */
#endif
#endif // LW_RADIATION_LOCAL

#ifdef BACKGROUND_LW_RADIATION
  k[25] += All.rate_Hminus_pd;

  double k29_ext = All.rate_H2_diss;  
#ifdef SHIELD_LW
  k29_ext *= sph_temp.f_byH2self /* SphP[ithis].f_byHI*/;
#ifdef HD
  double k53_ext = 1.55/1.39 * All.rate_H2_diss * sph_temp.f_byHDself * sph_temp.f_byHI * sph_temp.f_byH2;
#endif
#endif //SHIELD_LW

  k[29] += k29_ext;
#ifdef HD
  k[53] += k53_ext; 
#endif
  if(FlagPrint)
    printf("k25=%g %g k29=%g %g f=%g\n", k[25], All.rate_Hminus_pd, k[29], k29_ext, sph_temp.f_byH2self);
#endif //BACKGROUND_LW_RADIATION
#endif // LYMANWERNER

#ifdef HD
  for(j=31; j<=48; j++)
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif  

#ifdef THREEBODY
  for(j=55; j<=58; j++)
    k[j] = k_table[j][i_T] + c_T*(k_table[j][i_T+1] - k_table[j][i_T]);
#endif

  /****************************************************************************************************************/
  /******************************** We evaluate continuum opacities here ******************************************/
  /****************************************************************************************************************/

  double l_core = get_core_radius(ithis);

  double tau_cont;
#ifdef CONT_OPACITY
#ifdef CONTINUUM_COOLING_MATSUKOBA
  double dummy;
  struct opacities opac;

  dummy = calc_continuum_opacities(ithis, T_gas, rho_cgs, l_core, &opac);
  tau_cont = opac.tau;
#else
  tau_cont = 0.;
  opac.tau = 0.;
#endif
  
#ifdef DUST_COOLING
  double T_dust;

  double stellar_flux = 0.;
#ifdef THERMAL_FEEDBACK
  stellar_flux = SphP[ithis].stellar_rad_flux;
#endif

  T_dust = dust_temperature(rho_cgs, T_gas, SphP[ithis].Metallicity, l_core, stellar_flux, &opac, &edot_dust);
  SphP[ithis].Tdust = T_dust;
#ifdef RYTR_TAU_DUST
  SphP[ithis].kappa_dust = opac.kappa_dust_Planck;
  if(SphP[ithis].kappa_dust > 2.06 * SphP[ithis].Metallicity)
    SphP[ithis].kappa_dust = 2.06 * SphP[ithis].Metallicity;
#endif

  tau_cont = opac.tau;
#endif
  
  SphP[ithis].tau_cont = tau_cont;
  sph_temp.f_reduce = 1./(1. + 3./4.*tau_cont*tau_cont);
#else
  tau_cont = 0.;
#endif
  /****************************************************************************************************************/
  /****************************************************************************************************************/
  /****************************************************************************************************************/
	
  double nH;
#ifdef INVERSE_REACTIONS
  double n_crit;

  nH = XH * n_nuclei;
  //n_crit = H_criticaldensity(nH, T_gas, rho_cgs*SphP[ithis].HI*l_core, SphP[ithis].HI, SphP[ithis].elec);
  n_crit = H_criticaldensity(nH, T_gas, nH*SphP[ithis].HI*l_core, SphP[ithis].HI, SphP[ithis].elec);

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
    printf("k_inv=%g %g tau_cont=%g n=%g n_c=%g l_core=%g i_T=%d c_T=%g\n"
	   , k_inverse[2], k_inverse[7], tau_cont, n_nuclei, n_crit, l_core, i_T, c_T);
  }


#ifdef ID_CHEM_DEBUG
  if(P[ithis].ID == ID_CHEM_DEBUG || FlagPrint){
    printf("k_inv=%g %g tau_cont=%g n=%g n_c=%g l_core=%g i_T=%d c_T=%g\n"
	   , k_inverse[2], k_inverse[7], tau_cont, n_nuclei, n_crit, l_core, i_T, c_T);
  }
#endif

#endif
  
#ifdef H2_FORMATION_ON_DUST
  k[51] = H2form_on_grains(rho_cgs, T_gas, T_dust, SphP[ithis].Metallicity);
#endif

  
  //reaction 11 needs special treatment
  //nH = SphP[ithis].HI * n_nuclei;

  if(T_gas > 10.){
    //double nH;

    nH = XH * n_nuclei;	
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

    double xk_L, xk_H, rlgT4, rn_cr;
    double kcid_sk, a;

    xk_L = pow(10.0, log_gamma_cd) + pow(10.0, log_gamma_dt);
    xk_H  = k_inverse[55];

    rlgT4 = log10(T_gas*1.e-4);
    rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    a = 1.0/(1.0 + nH/rn_cr);

    if(a == 1.0)       k[11] = xk_L;
    else if(a == 0.0)  k[11] = xk_H;
    else               k[11] = pow(xk_H,1.0-a)*pow(xk_L,a);      
#else
    double xk_L, xk_H, rlgT4, rn_cr;
    double kcid_sk, a;
  
    xk_L  = 1.12e-10*exp(-7.035e+4/T_gas);
    xk_H  = 6.5e-7/sqrt(T_gas)*exp(-5.2e+4/T_gas)*(1.0-exp(-6.e+3/T_gas));
    rlgT4 = log10(T_gas*1.e-4);
    rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);
  
    a = 1.0/(1.0 + nH/rn_cr);
  
    if(a == 1.0)       k[11] = xk_L;
    else if(a == 0.0)  k[11] = xk_H;
    else               k[11] = pow(xk_H,1.0-a)*pow(xk_L,a);  
#endif
  }else
    k[11] = 0.;
}


void correct_initial_abundances(double *yi){
  int i_sp;
  double ratio1;

  /********** correct the total H abundances **********/
  if(yi[sHM] > 1.e-3) yi[sHM] = 1.e-10;
  
  double yHtot;  
  yHtot   = yi[sHI]   + yi[sHII]  + 2*yi[sH2I] + yi[sHM];
#ifdef H2_PLUS
  yHtot  += 2*yi[sH2II];
#endif

  if( fabs(1. - yHtot) > 1.e-3){
    ratio1 = 1./yHtot;

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
  double yHetot;
  yHetot  = yi[sHeI]  + yi[sHeII] + yi[sHeIII];

  if(yHetot > XHe/XH){
    ratio1 = (XHe/XH) / yHetot;
  
    yi[sHeI]   *= ratio1;
    yi[sHeII]  *= ratio1;
    yi[sHeIII] *= ratio1;
  }
  /****************************************************/
#endif

#ifdef HD
  double yDtot;
  yDtot   = yi[sDI]   + yi[sDII]  + yi[sHDI];

#ifdef HD_PLUS
  yDtot += yi[sHDII];
#endif
#ifdef DMINUS
  yDtot += yi[sDM];
#endif

  /********** correct the total D abundances **********/
  if(yDtot > D_number_fraction /XH){
    ratio1 = (D_number_fraction/XH) / yDtot;

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
  double y_minus,y_plus;
  y_minus = yi[sELEC] + yi[sHM];
  y_plus  = yi[sHII];

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
    double plus_charge;
    plus_charge = y_plus - y_minus;
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
  }else{
    yi[sDII]    = D_number_fraction/XH - yi[sDI] - yi[sHDI];
#ifdef HD_PLUS
    yi[sDII]   -= yi[sHDII];
#endif
#ifdef DMINUS
    yi[sDII]   -= yi[sDM];      
#endif
  }
  #endif
  /******************************************************************/
}

void give_initial_guess(double T_gas, double nH, double dt, double *yi, double *y_init, double *k){
  /************ Initial guess for the electron abundances ***********/  
  /*
  double n_e_crit, dn_e;
  n_e_crit = 2*(dt*k[1]*nH*XH-1)/(2*(k[2]+k[1])*dt);

  if(nH*XH*yi[sELEC] < n_e_crit){
    double d_yie, d_yie1, d_yie2, d_yie3;
    
    d_yie = dmin(0.8 - yi[sELEC], n_e_crit/nH/XH - yi[sELEC]);
    yi[sELEC] += d_yie;
#ifdef HELIUM
    if(yi[sHI] > 1.01 * d_yie)
#endif
      {
	yi[sHII]  += d_yie;
	yi[sHI]   -= d_yie;
      }
#ifdef HELIUM
    else if(yi[sHI]+yi[sHeI] > 1.01*d_yie){
      d_yie1 = 1./1.01 * yi[sHI];
      d_yie2 = d_yie - d_yie1;

      yi[sHII]   += d_yie1;
      yi[sHeII]  += d_yie2;

      yi[sHI]    -= d_yie1;
      yi[sHeI]   -= d_yie2;
    }else if(yi[sHI]+2*yi[sHeI]+yi[sHeII] > 1.01*d_yie){
      d_yie1 = 1./1.01 * yi[sHI];
      d_yie2 = dmin(0.9 * yi[sHeI], 0.5 * (d_yie - d_yie1));
      d_yie3 = d_yie - d_yie1 - 2*d_yie2;

      yi[sHII]   += d_yie1;
      yi[sHeIII] += (d_yie2 + d_yie3);
      yi[sHI]    -= d_yie1;
      yi[sHeI]   -= d_yie2;
      yi[sHeII]  -= d_yie3;
    }else{
      d_yie1 = 0.97 * (yi[sHI] + 2*yi[sHeI] + yi[sHeII]);
      d_yie2 = d_yie - d_yie1;

      yi[sELEC] -= d_yie2;
      
      d_yie1 = 0.9*yi[sHI];
      d_yie2 = 0.9*yi[sHeI];
      d_yie3 = 0.9*yi[sHeII];
      
      yi[sHII]   += d_yie1;
      yi[sHeIII] += (d_yie3+d_yie2);
      yi[sHI]    -= d_yie1;
      yi[sHeI]   -= d_yie2;
      yi[sHeII]  -= d_yie3;
    }
#endif

    if(yi[sHI] > 0.1){
      yi[sHI]  = 1. - yi[sHII] - yi[sHM] - 2*yi[sH2I];
#ifdef H2_PLUS
      yi[sHI] -= 2*yi[sH2II];
#endif
    }else{
      yi[sHII]  = 1. - yi[sHI]  - yi[sHM] - 2*yi[sH2I];
#ifdef H2_PLUS
      yi[sHII] -= 2*yi[sH2II];
#endif
    }
#ifdef HELIUM
    if(yi[sHeI] > 0.1 * XHe/XH)
      yi[sHeI]   = XHe/XH - yi[sHeII] - yi[sHeIII];
    else if(yi[sHeII] > 0.1 * XHe/XH)
      yi[sHeII]  = XHe/XH - yi[sHeI]  - yi[sHeIII];
    else
      yi[sHeIII] = XHe/XH - yi[sHeI]  - yi[sHeII];
#endif    
  }

#ifdef HELIUM
  yi[sHeII]  = 1./(1+dt*nH*XH*yi[sELEC]*(k[3]+k[4]+k[5])) * (dt*nH*XH*yi[sELEC]*k[3]*XHe   +y_init[sHeII]);
  yi[sHeIII] = 1./(1+dt*nH*XH*yi[sELEC]*k[6])             * (dt*nH*XH*yi[sELEC]*k[5]*yi[sHeII]+y_init[sHeIII]);
#endif
  
  //yi[sHM]    = 1./(1+dt*nH*XH*(k[8]+k[17]+yi[sELEC]*(k[16]+k[18]+k[19]))) * (dt*nH*XH*yi[sELEC]*k[7] + y_init[sHM]);
  */
#ifdef IONIZATION_PHOTON
  /*double hi;
  double k22_crit, x_neut_crit=1.e-3;  
  k22_crit = ((-k[1] + k[2]/x_neut_crit)*nH + (y_init[sHI]/x_neut_crit - 1.)/dt) * nH;

  if(T_gas > 5000 && k[22] > k22_crit){
    hi = 0.999*yi[sHI];
    yi[sHII]  += hi;
    yi[sELEC] += hi;
    yi[sHI]    = 0.001*yi[sHI];
    }

  // commented out just for test2c
  
  if(k[22] * nH * dt > 0.1 || T_gas > 3.e4){
    yi[sHII] = 0.99999;
    yi[sELEC] = 0.99999;
    yi[sH2I] = 1.e-9;
    yi[sHI] = 1.e-5;
  } 
  */
#endif
  
}


void bound_ddy(double *yi, double *ddy, double *y_max, int itr, int flag_reset, double *err_max){
  int i_sp, n_sp;
  n_sp = N_SP;

  if(itr%70 == 0){
    if( (i_sp == 0 || i_sp == 2) && fabs(ddy[i_sp]) > 1.e10*fabs(yi[i_sp]) ){
      double y_old[3];
      double d_yie;
      double yHtot0,yHctot;
	  
      y_old[0] = yi[0];
      y_old[1] = yi[1];
      y_old[2] = yi[2];

      yHtot0 = yi[sHI] + yi[sHII] + yi[sHM] + 2*yi[sH2I];
      yHctot = yi[sHII] + yi[sHM];
	  	  
      d_yie = dmax(0.99 - yi[sELEC], 0.);
	  
      double d_yie1, d_yie2, d_yie3;
	  
      yi[sELEC] += d_yie;
      {
	yi[sHII]  += d_yie;
	yi[sHI]   -= d_yie;
      }

      if(yi[sHI] > 0.1){
	yi[sHI]   = 1.  - yi[sHII] - yi[sHM] - 2*yi[sH2I];
      }else if(yi[sHII] > 0.1){
	yi[sHII]  = 1.  - yi[sHI] - yi[sHM] - 2*yi[sH2I];
      }else{
	yi[sH2I]  = 0.5 * (1.  - yi[sHI] - yi[sHII] - yi[sHM]);
      }


      ddy[sELEC] = ddy[sHI] = ddy[sHII] = 0.;
    }
#ifdef HD
    else if(i_sp == sDII  && (fabs(ddy[i_sp]) > 1.e30*fabs(yi[i_sp]) || itr%50==0) ){
      double f, sum_D;
      f = 1.0 - 0.1 * ( 1 + itr % 9);
      sum_D = yi[sDI] + yi[sDII];

      yi[sDI]  = f        * sum_D;
      yi[sDII] = (1. - f) * sum_D;

      ddy[sDI] = ddy[sDII] = 0.;
    }
#endif
  }
  
  if(itr%70 != 0){
    for(i_sp=0;i_sp<n_sp;i_sp++) {
      if(yi[i_sp] + ddy[i_sp] < 0.1 * yi[i_sp]){       
	ddy[i_sp] = -0.9 * yi[i_sp];
	
	if(yi[i_sp] > 1.e-10 * y_max[i_sp])
	  *err_max += (double)i_sp;	
      }else if(yi[i_sp]+ddy[i_sp] > y_max[i_sp]){
	ddy[i_sp] = dmin(0.9*yi[i_sp], 0.9 * (y_max[i_sp] - yi[i_sp]));
	*err_max += 100. * (double) i_sp;
      }
    }
  }
  
}


int InitChem(void)
{
  int ifunc;
  enum SPECIES species;
  enum REACTIONS reactions;

  /* ----- Set various coefs. */
	
  ifunc = H2line_tabulate();
	
  /* ----- cross-sections, etc ... */
  ifunc = chem_tabulate();
 
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
#ifdef TYPE_B
  double Teff = 1.e4;
#else
  double Teff = 1.e5;
#endif
  
  double hc = PLANCK * C_LIGHT;
  const double E_LW = 13.6;
  double J_13eV = 2. * pow(E_LW * eV_to_erg, 3) / hc / hc / (exp(E_LW*eV_to_erg/BOLTZMANN/Teff) - 1);
  double Norm = All.J21 / J_13eV;

  double PD_rate, g0;
  get_Hminus_PD_rate(Teff, &PD_rate, &g0);

  All.rate_Hminus_pd = Norm * PD_rate;
#ifdef DUST_COOLING
  All.G0           = Norm * g0 / 1.6e-3;
#endif
  All.rate_H2_diss = 1.39e9 * All.J21;
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
  double z_HI, z_HII, z_H2I, z_H2II, z_HM, z_elec;
  FILE *fd;

#ifdef CASE_B
  double lambda_HI, lambda_HeI, lambda_HeII;
#endif
  
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
  printf("you should choose which H2 cooling function will be used in the simulation!\n");
  endrun(200210);
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
#ifdef CONTINUUM_COOLING_MATSUKOBA
  continuum_escape_fraction = calc_continuum_opacities;
#else
  continuum_escape_fraction = calc_continuum_opacities_CIE;
#endif
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
    T[i] = exp(ln_T);
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
    z_H2I = partition_H2(T[i]);
    z_HI  = partition_H(T[i]);
#ifdef H2_PLUS
    z_H2II= partition_H2p(T[i]);
#endif

    z_HM  = 1.;
    z_elec = 2.;
    z_HII = 1.;

    /* ------- Compute various values of T. */
    ln_T = log(T[i]);
    log10_T = log10(T[i]);
    T_eV = T[i]/eV_to_K;
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
    lambda_HI   = 2 * 157807.0/T[i];
    lambda_HeI  = 2 * 285335.0/T[i];
    lambda_HeII = 2 * 631515.0/T[i];
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
    k_table[2][i] = 2.753e-14 * pow(315614.0/T[i], 1.5)/pow(1.0+pow(115188.0/T[i], 0.406), 2.242);
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

    
#ifdef CONTINUUM_COOLING_MATSUKOBA
    k_table[30][i] = 1.2e-17 * pow(T[i], 1.2) * exp(-1.578e5/T[i]);
#endif
      


    /*************** Reactions related to Helium ****************************/
    /*************** k3, 4, 5, 6, 14, 15          ****************************/

#ifdef HELIUM
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
      + 1.9e-3 * pow(T[i],-1.5)*exp( -470000./T[i] )
      *(1.+0.3*exp( -94000./T[i] )); /* dielectric recombination */
    /* Ferland et al. 1992 fitted by Hui & Gnedin 1997 */
#else
    //k_table[4][i] =  3.92e-13/pow(T_eV,0.6353); /* Case A radiative recombination */
    //1.9e-3 * pow(T[i],-1.5)*exp( -470000./T[i] )
    //	*(1.+0.3*exp( -94000./T[i] )); /* dielectric recombination */
			
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
    k_table[6][i] = 2.0 * 2.753e-14 * pow(315614.0/(T[i]/4.0), 1.5)/pow(1.0+pow(115188.0/(T[i]/4.0), 0.406), 2.242);
    /* 2 x k2(T/4) */ 
#else
    k_table[6][i] = 3.36e-10/sqrt(T[i])/pow(T[i]/1.e3,0.2)/(1+pow(T[i]/4.e6,0.7)); /* Fukugita-Kawasaki 94 */
#endif

    
    /* 14: HeII  + HI  -> HeI   + HII + p  charge transfer */
#ifdef ENZO_CHEM
    k_table[14][i] = 0.;
#else
    k_table[14][i] = 1.2e-15 * pow(T[i]/300.0, 0.25);
#endif    
    
    /* 15: HeI   + HII -> HeII  + HI */  
#ifdef ENZO_CHEM
    k_table[15][i] = 0.;
#else
    if ( T[i] > 10000.0)
      k_table[15][i] = 4.0e-37*pow(T[i], 4.74);
    else
      k_table[15][i] = 1.26e-9*pow(T[i], -0.75) * exp(-127500.0/T[i]);
#endif        
#endif


    if(T_eV < 0.01) {
      k_table[1][i] = 0;
#ifdef HELIUM
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
    k_table[7][i] = 1.4e-18 * pow(T[i], 0.928) * exp(-T[i]/1.62e4);
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    if(T[i] < 6000.){
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
    k_table[7][i] = 1.4e-18 * pow(T[i], 0.928) * exp(-T[i]/16200.0);  /* from GP98 */
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
    if(T[i] > 300.0) /* GP98, based on Launey et al. 1991 */
      k_table[8][i] = 4.0e-9 * pow(T[i], -0.17);
    else
      k_table[8][i] = 1.516e-9;    
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[8][i] = 1.35e-9* ( pow(T[i], 0.098493) + 0.32852 * pow(T[i],0.5561) + 2.771e-7 * pow(T[i],2.1826) )/
      (1. + 6.191e-3*pow(T[i],1.0461) + 8.9712e-11*pow(T[i], 3.0424) + 3.2576e-14*pow(T[i],3.7741)); //Kreckel+2010
#else
    if(T[i] > 300.0) /* GP98, based on Launey et al. 1991 */
      k_table[8][i] = 4.0e-9 * pow(T[i], -0.17);
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
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[18][i] = 2.4e-6*pow(T[i], -0.5) * (1. + T[i]/2.e4);       //Croft+1999

    //k_table_HI[1][i] = 6.92e-14* pow(T[i], 0.5);
    //k_table_HI[2][i] = 8.0e-13 * pow(T[i], 0.83);
    //k_table_HI[3][i] = 6.18e-7 * pow(T[i], -0.27);
#else
    k_table[18][i] = 6.92e-14*pow(T[i], 0.5) + 8.00e-13*pow(T[i], 0.83) + 6.18e-7 *pow(T[i],-0.27);
    //k_table[18][i] = 6.3e-8 + 5.7e-6/sqrt(T[i]) - 9.2e-11 * sqrt(T[i]) + 4.4e-13 * T[i];
#endif    


#ifdef CONTINUUM_COOLING_MATSUKOBA
    k_table[16][i] = 0.;
    k_table[17][i] = 0.;
    //k_table[18][i] = 0.;
#endif    



    /*************** Reactions related to H2II ****************************/
    /*************** k9, 10, 12, 19, 20, 21    ****************************/

#ifdef H2_PLUS
    /* 9: HI    + HII -> H2II  + p     H+ channel  */
    k_table[9][i] = 1.85e-23*pow(T[i],1.8);
    if (T[i] > 6.7e3) 
      k_table[9][i] = 5.81e-16*pow(T[i]/56200,-0.6657*log10(T[i]/56200));


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
    if(T[i] < 10000.0)
      k_table[12][i]=3.0e-10*exp(-21050.0/T[i]);
    else
      k_table[12][i]=1.5e-10*exp(-14000.0/T[i]);
#endif

    /* HM    + HII -> H2II  + e */
    if (T[i] > 1.0e4)
      k_table[19][i] = 4.0e-4*pow(T[i],-1.4)*exp(-15100.0/T[i]);
    else
      k_table[19][i] = 1.0e-8*pow(T[i],-0.4);

    /* 20: H2II  + e   -> 2HI */
#ifdef ENZO_CHEM
    if(T[i] > 617)
      k_table[20][i] = 1.32e-6 * pow(T[i],-0.76); 
    else
      k_table[20][i] = 1.e-8; 
#else
    /* old incorrect one k_table[20][i] = 5.56396e-8/pow(T_eV,0.6035); */
    k_table[20][i] = 2.0e-7/pow(T[i],0.5);  
#endif    
    
    /* 21: H2II  + HM  -> HI    + H2I  */
    k_table[21][i] = 4.64e-8/sqrt(T_eV);
#endif 


    /* k11 : H2I   + H   -> 3HI */
    /* k13 : H2I   + e   -> 2HI   + e     collisional dissociation */
#ifdef ENZO_CHEM
    k_table[11][i] = 1.0670825e-10*pow(T_eV,2.012)/
      (exp(4.463/T_eV)*pow(1+0.2472*T_eV,3.512));  
    k_table[13][i] = 5.6e-11*exp(-102124.0/T[i])*pow(T[i],0.5);
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
    k_table[13][i] = 1.91e-9 * pow(T[i], 0.136) * exp(-53407.1/T[i]);
#else
    if(T_eV > 0.3) {
      k_table[11][i] = 1.0670825e-10*pow(T_eV,2.012)/
	(exp(4.463/T_eV)*pow(1+0.2472*T_eV,3.512));  
      /* low density limit, see k56 */
      
      k_table[13][i] = 4.38e-10*exp(-102000.0/T[i])*pow(T[i],0.35);
    }else{
      k_table[11][i] = 0; 
      k_table[13][i] = 0;
    }
#endif    


#ifdef INVERSE_REACTIONS
    double K_T, f_elec;
    double Ene, dE, Z_pf, Mass;
    double m_HI, m_HII, m_HM, m_elec, m_H2I, m_H2II;

    m_HI = m_HII = m_HM = PROTONMASS;
    m_H2I = m_H2II = 2. * m_HI;
    m_elec = ELECTRONMASS;

    Ene = 2 * M_PI * BOLTZMANN * T[i] / PLANCK / PLANCK;


    /**** HI + e- -> Hp + 2e-****/
    if(T[i] > 1.e3){
      dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHII] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_elec/ z_HII / z_elec / z_elec;
      Mass = m_HI * m_elec/ m_HII / m_elec / m_elec; 
      
      K_T = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
      k_table_inverse[1][i] = k_table[1][i] * K_T;
    }



    /**** Hp + e- -> HI  ****/
    dE   = (chem_e[sHII] + chem_e[sELEC] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HII * z_elec / z_HI;
    Mass = m_HII * m_elec / m_HI;     
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T[i]);    
    k_table_inverse[2][i] = k_table[2][i] * K_T; //note this value is optically thick limit.


    /**** HI + e- -> HM  ****/
    dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHM]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_elec / z_HM;
    Mass = m_HI * m_elec / m_HM;     
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T[i]);    
    k_table_inverse[7][i] = k_table[7][i] * K_T; //note this value is optically thick limit.


    /**** HI + HM -> H2I + e-  ****/
    dE   = (chem_e[sHI] + chem_e[sHM] - chem_e[sH2I] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_HM / z_H2I / z_elec;
    Mass = m_HI * m_HM / m_H2I / m_elec;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T[i]);    
    k_table_inverse[8][i] = k_table[8][i] * K_T;

#ifdef H2_PLUS
    /**** H2I + Hp -> H2II + HI  ****/
    dE   = (chem_e[sH2I] + chem_e[sHII] - chem_e[sH2II] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_H2I * z_HII / z_HI / z_H2II;
    Mass = m_H2I * m_HII / m_HI / m_H2II;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T[i]);    
    //k_table_inverse[10][i] = k_table[10][i] * K_T;
    k_table[10][i] = k_table[12][i] * K_T;
#endif

    if(T[i]>200.){
      /**** H2I + e- -> 2HI + e-  ****/
      dE   = (chem_e[sH2I] + chem_e[sELEC] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_elec / z_HI / z_HI / z_elec;
      Mass = m_H2I * m_elec / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
      k_table_inverse[13][i] = k_table[13][i] * K_T;
    }


    if(T[i]>200.){
      /**** HM + e- -> HI + 2e-  ****/
      dE   = (chem_e[sHM] + chem_e[sELEC] - chem_e[sHI] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_elec / z_HI / z_elec / z_elec;
      Mass = m_HM * m_elec / m_HI / m_elec / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
#ifdef AIUEO
      k_table_inverse[16][i] = k_table[16][i] * K_T;
#endif    
    }



#ifdef AIUEO
    if(T[i]>200.){
      /**** HM + HI -> 2HI + e-  ****/
      dE   = (chem_e[sHM] + chem_e[sHI] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_HI / z_HI / z_HI / z_elec;
      Mass = m_HM * m_HI / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
      k_table_inverse[17][i] = k_table[17][i] * K_T;
    }
#endif

    //#ifdef AIUEO
    /**** HM + Hp -> HI + HI  ****/
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T[i]);    
    k_table_inverse[18][i] = k_table[18][i] * K_T;

#ifdef CONTINUUM_COOLING_MATSUKOBA
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T[i]);    
    //k_table_HI_inverse[1][i] = k_table_HI[1][i] * K_T;
#endif

    if(T[i] > 1.e3){
      /***** HI + HI -> Hp + HI + e- *****/
      dE   = (2*chem_e[sHI] - chem_e[sHII] - chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI / z_HII / z_elec / z_HI;
      Mass = m_HI * m_HI / m_HII / m_elec / m_HI;     
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
      k_table_inverse[30][i] = k_table[30][i] * K_T;
    }

#endif


#ifdef HD
    /* deuterium network */
    //k_table[31][i]=3.6e-12*pow(T[i]/300.0, -0.75);   
    //k_table[34][i]=1.0e-25;
    /* new rate Savin (2002)  (old one may be 3.7e-10*pow(T[i],0.28)*exp(-43.0/T[i]) ) */
    //k_table[32][i]=2.0e-10*pow(T[i],0.402)*exp(-37.1/T[i]) - 3.31e-17*pow(T[i],1.48);     
    //k_table[33][i]=2.06e-10*pow(T[i],0.396)*exp(-33.0/T[i]) + 2.03e-9*pow(T[i],-0.332);  
    //k_table[37][i]=1.6e-9;/* old one may be 2.1e-9 */
    /* new rate Savin (2002)  (old one may be 3.7e-10*pow(T[i],0.28)  )                */

    k_table[32][i]=3.7e-10*pow(T[i],0.28)*exp(-43.0/T[i]);
    k_table[33][i]=3.7e-10*pow(T[i],0.28);
    k_table[35][i]=9.0e-11*exp(-3876.0/T[i]); 
    k_table[37][i]=2.1e-9;
    k_table[38][i]=3.2e-11*exp(-3624.0/T[i]);
    k_table[39][i]=1.0e-9*exp(-464.0/T[i]);

#ifdef HD_PLUS
    k_table[36][i]=6.4e-10;
    if(T[i] > 10.0){
      k_table[40][i]=pow(10.0, -19.38-1.523*log10_T+1.118*log10_T*log10_T-0.1269*log10_T*log10_T*log10_T);
      k_table[41][i]=pow(10.0, -19.38-1.523*log10_T+1.118*log10_T*log10_T-0.1269*log10_T*log10_T*log10_T);
    }else
      k_table[40][i]=k_table[41][i]=0;
    
    k_table[42][i]=7.2e-8/sqrt(T[i]);
#endif

#ifdef DMINUS
    /* D- part */
    k_table[43][i]=3.0e-16*pow(T[i]/300.0, 0.95)*exp(-T[i]/9320.0);   
    k_table[44][i]=5.7e-8*pow(T[i]/300.0, -0.5);   
    k_table[45][i]=4.6e-8*pow(T[i]/300.0, -0.5);   
    k_table[46][i]=6.4e-9*pow(T[i]/300.0, 0.41);   
    k_table[47][i]=6.4e-9*pow(T[i]/300.0, 0.41);   
    k_table[48][i]=1.5e-9*pow(T[i]/300.0, -0.1);   
#endif
#endif  
    /* 3-body reacions
       k_table[55][i]=5.5e-29/T[i];
       k_table[56][i]=6.5e-7/sqrt(T[i])*exp(-52000.0/T[i])*(1.0-exp(-6000.0/T[i]));
       k_table[57][i]=6.875e-30/T[i];
       k_table[58][i]=8.125e-8/sqrt(T[i])*exp(-52000.0/T[i])*(1.0-exp(-6000.0/T[i]));  */

    /* 3-body, new rate, consistent with equilibrium chemistry */    
    k_table[55][i]=4.5e-29/T[i];
    k_table[56][i]=6.5e-7/sqrt(T[i])*exp(-52000.0/T[i])*(1.0-exp(-6000.0/T[i]));
    k_table[57][i]=k_table[55][i]/40.0;
    k_table[58][i]=5.0e-8/sqrt(T[i])*exp(-52000.0/T[i])*(1.0-exp(-6000.0/T[i]));

    /* 3-body, rate coefficient uncertainty */
#ifdef TBR_FLOW07
    k_table[55][i] = 1.44e-26/pow(T[i], 1.54);
#elif defined(TBR_ABEL02)
    if(T[i] < 300){
      k_table[55][i] = 1.14e-31/pow(T[i], 0.38);
    }else{
      k_table[55][i] = 3.9e-30/T[i];
    }
#elif defined(TBR_FORR13)
    k_table[55][i] = 6.e-32/pow(T[i], 0.25) + 2.e-31/pow(T[i], 0.5);
#else
    k_table[55][i] = 6.e-32/pow(T[i], 0.25) + 2.e-31/pow(T[i], 0.5);
#endif


#ifdef INVERSE_REACTIONS
    if(T[i] > 500.){
      /****  HI + HI + HI -> HI + H2I  ****/
      dE   = (3*chem_e[sHI] - chem_e[sH2I] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI * z_HI / z_H2I / z_HI;
      Mass = m_HI * m_HI * m_HI / m_H2I / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T[i]);    
      k_table_inverse[55][i] = k_table[55][i] * K_T;


      /****  H2I + H2I -> H2I + HI + HI  ****/
      dE   = (2*chem_e[sH2I] - 2*chem_e[sHI] - chem_e[sH2I]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_H2I / z_HI / z_HI / z_H2I;
      Mass = m_H2I * m_H2I / m_HI / m_HI / m_H2I;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T[i]);    
      k_table[57][i] = k_table[58][i] * K_T;
    }
#endif        
  }
  
  
  
  /* coefficients for cooling and heating */
  for(i=0; i<N_T; i++){
    
    /* ------- Collisional excitations (Black 1981; Cen 1992) */
    ceHIa[i]   = 7.5e-19*exp( -118348./T[i] )
      /(1.+sqrt(T[i]/1.e5));
    
    ceHeIa[i]   = 9.1e-27*exp( -13179./T[i] )
      *pow(T[i],-0.1687)/(1.+sqrt(T[i]/1.e5));
    
    ceHeIIa[i]  = 5.54e-17*exp( -473638./T[i] )
      *pow(T[i],-0.397)/(1.+sqrt(T[i]/1.e5));
    

    /* ------- Collisional ionizations (Tom's polynomial fits) */
    ciHIa[i]    = 2.18e-11*k_table[1][i];
    ciHeIa[i]   = 3.94e-11*k_table[3][i];
    ciHeIIa[i]  = 8.72e-11*k_table[5][i];

    ciHeISa[i]  = 5.01e-27*pow(T[i],-0.1687)/(1.+sqrt(T[i]/1.e5))
      *exp( -55338./T[i] );
    
    /* ------- Recombinations */
#ifndef CONTINUUM_COOLING_MATSUKOBA  
#ifdef CASE_B
    reHIIa[i]   = 3.435e-30 * T[i] * pow(lambda_HI, 1.97)/pow(1.0+pow(lambda_HI/2.25, 0.376), 3.72);
    
#ifdef HELIUM
    reHeII1a[i] = BOLTZMANN * T[i] * 1.26e-14 * pow(lambda_HeI, 0.75); 
    
    reHeII2a[i] = 1.24e-13*pow(T[i],-1.5)*exp( -470000/T[i] )
      *(1.+0.3*exp( -94000./T[i] ));     /* He dielectric recombination */
    
    reHeIIIa[i] = 8 * 3.435e-30 * T[i] * pow(lambda_HeII, 1.97)/pow(1.0+pow(lambda_HeII/2.25, 0.376), 3.72);
#endif
#else
    reHIIa[i]   = 6.50e-27*sqrt(T[i])*pow(T[i]/1000.0,-0.2)
      / (1.0 + pow(T[i]/1.e6,0.7));

#ifdef HELIUM    
    reHeII1a[i] = 1.55e-26*pow(T[i],0.3647);

    reHeII2a[i] = 1.24e-13*pow(T[i],-1.5)*exp( -470000/T[i] )
      *(1.+0.3*exp( -94000./T[i] ));     /* He dielectric recombination */
    
    reHeIIIa[i] = 3.48e-26*sqrt(T[i])*pow(T[i]/1000.0,-0.2)
      / (1.0 + pow(T[i]/4.e6,0.7));
#endif    
#endif    


    /* ------- Bremsstrahlung (Black 1981)(Spitzer & Hart 1979) */
    brema[i]    = 1.43e-27*sqrt(T[i])
      *(1.1+0.34*exp(-pow(5.5-log10(T[i]),2.0)/3.0));
  
    log10_T = log10(T[i]);
    H2_CIE[i]   = pow(10.0, - 116.6 
		      + 96.3   *log10_T 
		      - 47.15  *log10_T*log10_T 
		      + 10.744 *log10_T*log10_T*log10_T 
		      -  0.9155*log10_T*log10_T*log10_T*log10_T);    
#endif
  }

  

#if defined(BG_RADIATION) || defined(IONIZATION_PHOTON)
  for(j=0;j<s_MAX;j++)
    for(i=0;i<N_nu;i++)
      sigma[j][i]=0;

  SIGMA_e22 = 6.30e-18;  //[cm^2] evaluated at nu = e22
  SIGMA_e23 = 1.58e-18;  //[cm^2] evaluated at nu = e23
  SIGMA_e24 = 7.42e-18;  //[cm^2] evaluated at nu = e24

  /* ----- Compute cross-sections. */
  for(i = 0;i< N_nu;i++){

    /* photoionisation HI */
    if( nu[i] > e22 ){
      double epsilon = sqrt(nu[i]/e22-1);
      sigma[22][i] = 6.3e-18*pow(e22/nu[i], 4)*exp(4.0-4.0*atan(epsilon)/epsilon)/ (1-exp(-2.0*M_PI/epsilon));
    }else{
      sigma[22][i] = 0.0;
    }

    /* photoionisation HeII */
    if( nu[i] > e23 ){
      double epsilon = sqrt(nu[i]/e23-1);
      sigma[23][i] = 1.58e-18*pow(e23/nu[i], 4)*exp(4.0-4.0*atan(epsilon)/epsilon)/ (1-exp(-2.0*M_PI/epsilon));
    }else{
      sigma[23][i] = 0.0;
    }

    /* photoionisation HeI */
    if( nu[i] > e24 )
      sigma[24][i] = 7.42e-18*(1.66*pow(nu[i]/e24, -2.05)- 0.66*pow(nu[i]/e24,-3.05));
    else
      sigma[24][i] = 0.0;
	

    /* photodissociation HM */
    if( nu[i] > e25 )
      sigma[25][i] = 2.11e-16*pow(nu[i]-e25,1.5)/pow(nu[i],3);  /* unit normalization enrgy in eV */
    else
      sigma[25][i] = 0.0;
      
    
    /* photodissociation H2II  Shapiro & Kang Table 3 */
    if( nu[i] > e26a && nu[i] <= e26b )
      sigma[26][i] = pow(10,-40.97+6.03*nu[i]-0.504*nu[i]*nu[i]+1.387e-2*nu[i]*nu[i]*nu[i]);
    else 
      if( nu[i] > e26b && nu[i] < e26c )
	sigma[26][i] = pow(10,-30.26+2.79*nu[i]-0.184*nu[i]*nu[i]+3.535e-3*nu[i]*nu[i]*nu[i]); 
      else
	sigma[26][i] = 0.0;
    
    
    /* photodissociation H2I  Shapiro & Kang Table 3 */
    if( nu[i] > e27a && nu[i] <= e27b )
      sigma[27][i] = 6.2e-18*nu[i] - 9.4e-17;
    else 
      if( nu[i] > e27b && nu[i] <= e27c )
	sigma[27][i] = 1.4e-18*nu[i] - 1.48e-17;
      else 
	if( nu[i] > e27c )
	  sigma[27][i] = 2.5e-14*pow(nu[i],-2.71);
	else
	  sigma[27][i] = 0.0;
          

    /* photodissociation H2II Shapiro & Kang Table 3 */
    if( nu[i] >= e28a && nu[i] < e28b )
      sigma[28][i] = pow(10.0,-16.926-4.528e-2*nu[i]+2.238e-4*nu[i]*nu[i] +4.245e-7*nu[i]*nu[i]*nu[i]);
    else
      sigma[28][i] = 0.0;
        

    /* photodissociation H2I by LW band photons */
    sigma[29][i] = 0.0;  /* given directly as k29 */
    
  }
#endif

#ifdef CALC_CONT_OPACITY
  for(i = 0;i< N_nu;i++){
    /****** H bound-free cross section ******/
    sigma[1] = 7.909e-18 * pow(nu[i]/13.598,-3);
 
    /****** He bound-free cross section ******/
    sigma[2] = 7.83e-18 * (1.66 * pow(nu[i]/24.586, -2.05) - 0.66*pow(nu[i]/24.586),-3.05);

    /****** H- bound-free cross section ******/
    sigma[3] = 1.e-18*1;
  }
#endif


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


#if defined(LYMANWERNER) && !defined(LW_RADIATION_LOCAL)
double calc_G0(void){
  int i;
  double Norm, hc;
  double nu, nu_old, dnu, sum, J0;
  const double nu_min=5.12, nu_max= 13.6;
  const int n_bin=500;

  hc    = CLIGHT * PLANCK;

  Norm = 2.0 * pow(e22*eV_to_erg, 3) / (hc*hc) /(exp(e22*eV_to_erg/(BOLTZMANN*T_eff))-1.0);
  Norm = 1.0 / Norm;
  
  double dln_nu;
  dln_nu = log(nu_max / nu_min)/n_bin;

  for(i=0,sum=0.; i<n_bin; i++){
    nu = nu_min * exp(dln_nu * i);
    if(nu < nu_min || nu > nu_max)
      continue;
    
    J0 = 2.0 * pow(nu*eV_to_erg, 3) / (hc*hc) / (exp(nu*eV_to_erg/(BOLTZMANN*T_eff))-1.0);
    if(i>0){
      dnu = (nu - nu_old) * eV_to_Hz;
      sum += J0 * dnu;
    }

    nu_old = nu;
  }

  sum *= Norm;
  sum *= (4 * M_PI / 1.6e-3);

  /****** the result is expressed in [(intensity/nu) * nu / (intensity/nu)] = [nu] ******/
  /****** Multiplying [intensity/nu], i.e. J_21[nu=13.6eV],                        ******/                            
  /****** we can get the total intensity integrated over FUV (5.12-13.6 eV) band.  ******/

  return sum;
}
#endif

#ifdef IONIZATION_PHOTON
#define DELTA_TAU_CRIT 0.001
#define N_tau  750
#define N_Teff 300

//const int N_tau  = 250;
//const int N_Teff = 300;
const double Teff_min = 1.e3;
const double Teff_max = 1.e6;
const double tau_min  = 1.e-3;
const double tau_max  = 1.e14;
double dln_tau, dln_Teff;
double Teffs[N_Teff], Tau[N_tau];
double PHI22a_thin[N_Teff][N_tau], PHI22b_thin[N_Teff][N_tau];
double PHI22a[N_Teff][N_tau], PHI22b[N_Teff][N_tau];
double k_ion1[N_Teff], GammaIon1[N_Teff], J_nu_integ[N_Teff];

void calc_phis(double Teff, double tau, double *phi1, double *phi2, double *kion1, double *gamma1);

int init_ionization_rad(void){
  int i,j,ifunc;
  double dlog_nu,log_nu,eps,nu_bound[N_nu_bound];
  static int flag_rate=0;

  /* ----- Compute base nu grid (in eV) and spectral shape J0_nu. */
  dlog_nu = (log(nu_end) - log(nu_start))/N_nu_base;
  for(i = 0; i< N_nu_base;i++){
    log_nu = log(nu_start) + i*dlog_nu;
    nu[i] = exp(log_nu);
  }

  /* ----- Set transition points. */
  nu_bound[0]  = e22;
  nu_bound[1]  = e23;
  nu_bound[2]  = e24;
  nu_bound[3]  = e25;
  nu_bound[4]  = e26a;
  nu_bound[5]  = e26b;
  nu_bound[6]  = e26c;
  nu_bound[7]  = e27a;
  nu_bound[8]  = e27b;
  nu_bound[9]  = e27c;
  nu_bound[10] = e28a;
  nu_bound[11] = e28b;
  nu_bound[12] = e30;

  /* ----- Append boundary points to nu. */
  eps = 1.e-5;
  for(i=0;i<N_nu_bound;i++){
    j = i*N_pts_bound + N_nu_base-1;
    nu[j+1] = nu_bound[i]*(1. - 4.*eps);
    nu[j+2] = nu_bound[i]*(1. - 2.*eps);
    nu[j+3] = nu_bound[i]*(1. - 1.*eps);
    nu[j+4] = nu_bound[i];
    nu[j+5] = nu_bound[i]*(1. + 1.*eps);
    nu[j+6] = nu_bound[i]*(1. + 2.*eps);
    nu[j+7] = nu_bound[i]*(1. + 4.*eps);
  }

  /* ----- Sort nu. */
  ifunc= heap_sort(N_nu,nu-1);
  
  return i;
}


int make_ionization_rate_table(void)
{
  /*----- Declare other values. */
  int i,j;
  int i_tau, i_Teff;
  double Teff, tau;
  double phi1a, phi1b, phi2a, phi2b, Jtot;
  double J_nu_integ0;
  double delta_nu, c0, c1;
  double Norm, hc, J0_nu[N_nu];

  FILE *fp, *fp2;

  dln_tau  = log(tau_max/tau_min)   / N_tau;
  dln_Teff = log(Teff_max/Teff_min) / N_Teff;

  hc   = CLIGHT * PLANCK;

  if(ThisTask == 0){
    if( !(fp = fopen("ionization_rate_thick.txt","w")) ){
      printf("Can't open the file, ./ionization_rate_thick.txt\n");
      exit(1);
    }

    if( !(fp2 = fopen("ionization_rate_thin.txt","w")) ){
      printf("Can't open the file, ./ionization_rate_thin.txt\n");
      exit(1);
    }
  }

  for(i_Teff=0;i_Teff<N_Teff;i_Teff++){
    Teff = Teff_min * exp(i_Teff * dln_Teff);
    Teffs[i_Teff] = Teff;

    /* normalization at Lyman limit */
    /* J21 normalization (amplitude) is done seperately */
    /* unit of J0_nu: J_nu[hv=e22] erg/sec/cm^2/Hz/str (cgs.) */
    
    Norm = 2.0 * pow(e22*eV_to_erg, 3) / (hc*hc) /(exp(e22*eV_to_erg/(BOLTZMANN*Teff))-1.0);
    Norm = 1.0 / Norm;
    
    for(i=0; i<N_nu; i++){
      if(nu[i] > 0)
	J0_nu[i] = Norm * 2.0 * pow(nu[i]*eV_to_erg, 3) / (hc*hc) / (exp(nu[i]*eV_to_erg/(BOLTZMANN*Teff))-1.0);
      else
	J0_nu[i] = 0.;
    }

#ifdef GREY_SPECTRUM
    nu_max = nu_end;
    nu_min = nu_start;
    
    Norm = (nu_max - nu_min) * eV_to_Hz;

    for(i=0; i<N_nu; i++){
      if(nu[i] > nu_min && nu[i] < nu_max)
	J0_nu[i] = 1.;
      else
	J0_nu[i] = 0.;
    }
#endif

    phi1a = phi2a = Jtot = 0.;
    /*----- Compute k values integrating over frequency. */
    for(j=1; j<N_nu; j++){      
      if(nu[j - 1] > e22){
	delta_nu = (nu[j] - nu[j - 1]) * eV_to_Hz;
	c0 = J0_nu[j] / (nu[j] * eV_to_erg);
	c1 = J0_nu[j - 1] / (nu[j - 1] * eV_to_erg);
	
	phi1a += (0.5 * 0.25 * M_1_PI) * (c0 * sigma[22][j] + c1 * sigma[22][j - 1] ) * delta_nu;
	
	c0 = J0_nu[j] * (nu[j] - e22) / nu[j];
	c1 = J0_nu[j - 1] * (nu[j - 1] - e22) / nu[j - 1];
	 	
	phi2a += (0.5 * 0.25 * M_1_PI) * ( c0*sigma[22][j] + c1*sigma[22][j - 1] ) * delta_nu;
      }    

      delta_nu = (nu[j] - nu[j - 1]) * eV_to_Hz;
      c0 = J0_nu[j];
      c1 = J0_nu[j - 1];

      Jtot += 0.5 * (c0 + c1) * delta_nu; 
    }

    k_ion1[i_Teff]     = phi1a;
    GammaIon1[i_Teff]  = phi2a;

    J_nu_integ[i_Teff] = Jtot;
    J_nu_integ0 = 2*CLIGHT*pow(M_PI*BOLTZMANN*Teff, 4)/hc/hc/hc/15. * Norm;

    if(ThisTask == 0){
      fprintf(fp2, "%g %g %g\n", Teff, phi1a, phi2a);
      //printf("Teff=%g: J_nu=%g %g\n", Teff, J_nu_integ[i_Teff], J_nu_integ0);
    }

    for(i_tau=0;i_tau<N_tau;i_tau++){
      tau = tau_min * exp(i_tau * dln_tau);
      Tau[i_tau] = tau;
      
      k[22] = 0.;
    
      /*----- Compute k values integrating over frequency. */
      phi1a = phi2a = 0.;
      phi1b = phi2b = 0.;
      for(j=1; j<N_nu; j++){      

	if(nu[j - 1] > e22){
	  delta_nu = (nu[j] - nu[j - 1]) * eV_to_Hz;
	  c0 = J0_nu[j] / (nu[j] * eV_to_erg);
	  c1 = J0_nu[j - 1] / (nu[j - 1] * eV_to_erg);
	
	  phi1a += (0.5 * 0.25 * M_1_PI) * (c0 * exp( -tau * sigma[22][j]/SIGMA_e22) * sigma[22][j]
					    + c1 * exp( -tau * sigma[22][j - 1]/SIGMA_e22) * sigma[22][j-1] ) * delta_nu;
	  phi1b += (0.5 * 0.25 * M_1_PI) * (c0 * exp( -tau * sigma[22][j]/SIGMA_e22) + c1 * exp( -tau * sigma[22][j - 1]/SIGMA_e22) ) * delta_nu;
	  
	  c0 = J0_nu[j] * (nu[j] - e22) / nu[j];
	  c1 = J0_nu[j - 1] * (nu[j - 1] - e22) / nu[j - 1];
	 	
	  phi2a += (0.5 * 0.25 * M_1_PI) * (c0 * exp( -tau * sigma[22][j]/SIGMA_e22) * sigma[22][j]
					    + c1 * exp( -tau * sigma[22][j - 1]/SIGMA_e22) * sigma[22][j-1] ) * delta_nu;
	  phi2b += (0.5 * 0.25 * M_1_PI) * (c0 * exp( -tau * sigma[22][j]/SIGMA_e22) + c1 * exp( -tau * sigma[22][j - 1]/SIGMA_e22) ) * delta_nu;

	  if(tau > 10. && ThisTask == 0 && (c0>0 || c1 > 0) && (Teff > 10000 && Teff < 12000) ){
	    //printf("tau=%g phi1=%g %g phi2=%g %g, exp=%g c0 =%g c1=%g delta=%g\n",
	    //       tau, phi1a, phi1b, phi2a, phi2b, exp( -tau * sigma[22][j]/SIGMA_e22), c0, c1, delta_nu), fflush(stdout);
	    //printf("tau=%g exp=%g %g c0 =%g c1=%g delta=%g\n",
	    //      tau,  exp( -tau * sigma[22][j]/SIGMA_e22), exp( -tau * sigma[22][j - 1]/SIGMA_e22), J0_nu[j] / (nu[j] * eV_to_erg), J0_nu[j - 1] / (nu[j - 1] * eV_to_erg),
	    //	    delta_nu);
	  }
	}
    
      }

      PHI22a_thin[i_Teff][i_tau] = phi1a;
      PHI22b_thin[i_Teff][i_tau] = phi2a;

      PHI22a[i_Teff][i_tau] = phi1b;
      PHI22b[i_Teff][i_tau] = phi2b;

      if(ThisTask == 0)
	fprintf(fp, "%g %g %g %g\n", Teff, tau, phi1b, phi2b);
    }

  }

  if(ThisTask == 0)
    fclose(fp), fclose(fp2);

  return j;
}


void calc_ionization_rate(int ithis, double Teff, double *k22, double *gamma_ion, double dr, double ri){
  double phi1a, phi2a, phi1b, phi2b;
  double kion1, gamma_ion1, kion2, gamma_ion2;
  double denom;
  double tau, d_tau, n_H;
  double w1, w2;
  double trans_len_cm;

  if(All.ComovingIntegrationOn)
    trans_len_cm =  All.Time / All.HubbleParam * All.UnitLength_in_cm ;  /* [comoving kpc/h --> physical cm] */
  else
    trans_len_cm =  1. / All.HubbleParam * All.UnitLength_in_cm;         /* [physical kpc/h --> physical cm] */
  
  n_H   = SphP[ithis].n_nuclei * SphP[ithis].HI;
  tau   = SphP[ithis].tau_ion;
  d_tau = SphP[ithis].d_tau_0;

  if(ri < 0.001 * SphP[ithis].Hsml)
    ri = 0.001 * SphP[ithis].Hsml;

  if(dr < 1.e-30 * SphP[ithis].Hsml)
    dr = 1.e-30 * SphP[ithis].Hsml;
  
#ifndef TEST_MONOCHROMATIC_SPECTRUM
  calc_phis(Teff, tau        , &phi1a, &phi2a, &kion1, &gamma_ion1);
  calc_phis(Teff, tau - d_tau, &phi1b, &phi2b, &kion1, &gamma_ion1);

  denom = 4 * M_PI * ri * ri * pow(trans_len_cm, 2);

  kion1      *= (n_H/denom);
  gamma_ion1 *= (n_H/denom);

  denom = dr * (ri*ri - ri*dr + 1./3 * dr*dr) * pow(trans_len_cm, 3);

  kion2      = (phi1b - phi1a)/denom ;
  gamma_ion2 = (phi2b - phi2a)/denom ;

  if(d_tau/tau < 1.e-13)
    kion2 = gamma_ion2 = 0.;
#else
  kion1      = exp(-tau) * SIGMA_e22/(e22 * eV_to_erg);
  gamma_ion1 = 0;
  denom      = 4 * M_PI * ri * ri * pow(trans_len_cm, 2);

  kion1      *= (n_H/denom);
  gamma_ion1 *= (n_H/denom);

  phi1a      = exp(-tau)       / (e22 * eV_to_erg) / (4 * M_PI);
  phi1b      = exp(-tau+d_tau) / (e22 * eV_to_erg) / (4 * M_PI);
  denom = dr * (ri*ri - ri*dr + 1./3 * dr*dr) * pow(trans_len_cm, 3);
  
  /* didn't work!!
    phi1a      = exp(-tau-d_tau/2) / (e22 * eV_to_erg) / (4 * M_PI);
    phi1b      = exp(-tau+d_tau/2) / (e22 * eV_to_erg) / (4 * M_PI);
  
    denom = dr * (ri*ri - 1./12 * dr*dr) * pow(trans_len_cm, 3);*/

  kion2      = (phi1b - phi1a)/denom ;
  gamma_ion2 = 0;
#endif

  w1 = 1/ (1. + d_tau/DELTA_TAU_CRIT); 
  w2 = 1/ (1. + DELTA_TAU_CRIT/d_tau); 

  if(n_H > 1.e-20){
    *k22       = (w1 * kion1      + w2 * kion2)/n_H;
    *gamma_ion = (w1 * gamma_ion1 + w2 * gamma_ion2)/n_H;
  }else{
    *k22 = 0.;
    *gamma_ion = 0.;
  }

  if(P[ithis].ID == 15269986){
    printf("[ID%d] n_HI=%g k22=%g %g ", P[ithis].ID, SphP[ithis].n_nuclei * SphP[ithis].HI, *k22, SphP[ithis].k22);
    printf("dr=%g ri=%g denom=%g tau=%g d_tau=%g\n", dr, ri, denom, tau, d_tau);
  }

  if(*k22 < 0){
    printf("minus! (n_HI=%g) k22=%g %g\n", SphP[ithis].n_nuclei * SphP[ithis].HI, *k22, SphP[ithis].k22);
    printf("dr=%g ri=%g denom=%g tau=%g d_tau=%g\n", dr, ri, denom, tau, d_tau);

    printf("kion=%g %g gamma=%g %g w=%g %g\n", kion1, kion2, gamma_ion1, gamma_ion2,w1,w2);
    raise(11);
  }
}



void calc_phis(double Teff, double tau, double *phi1, double *phi2, double *kion1, double *gamma1){
  volatile int i_T, i_tau;
  double ln_Teff0, ln_tau0;
  double ln_Teff, ln_tau;
  double c_Teff, c_tau;

  ln_Teff0 = log(Teff_min);

  if (Teff < Teff_min) Teff = Teff_min;
  if (Teff > Teffs[N_Teff-1]) Teff = Teffs[N_Teff-1];
  ln_Teff = log(Teff);
  
  i_T = (int)( (ln_Teff-ln_Teff0)/dln_Teff ) ;
  
  /*i_T = 0;
  while(i_T*dln_Teff + ln_Teff0 < ln_Teff){
    i_T++;
  }
  i_T--;*/

  if (i_T > N_Teff-2) i_T = N_Teff - 2;     /* [i_T+1] appears below */  
  if (i_T < 0    ) i_T = 0;
  
  c_Teff = (Teff - Teffs[i_T])/(Teffs[i_T+1]-Teffs[i_T]);
 
  ln_tau0 = log(tau_min);

  if (tau < tau_min) tau = tau_min;
  if (tau > Tau[N_tau-1]) tau = Tau[N_tau-1];
  ln_tau = log(tau);
  
  i_tau = (int)( (ln_tau-ln_tau0)/dln_tau ) ;
  
  /*i_tau = 0;
  while(i_tau*dln_tau + ln_tau0 < ln_tau){
    i_tau++;
  }
  i_tau--;*/

  if (i_tau > N_tau-2) i_tau = N_tau - 2;     /* [i_T+1] appears below */  
  
  c_tau = (tau - Tau[i_tau])/(Tau[i_tau+1]-Tau[i_tau]);

  //printf("c_tau=%g (%g %g) tau=%g c_Teff=%g (%g %g)\n",
  //       c_tau,  Tau[i_tau],  Tau[i_tau+1], tau, c_Teff, Teffs[i_T], Teffs[i_T+1]);

  *phi1 = c_Teff * c_tau * PHI22a[i_T+1][i_tau+1] + c_Teff * (1 - c_tau) * PHI22a[i_T+1][i_tau]
    + (1 - c_Teff) * c_tau * PHI22a[i_T][i_tau+1] + (1 - c_Teff) * (1 - c_tau) * PHI22a[i_T][i_tau];

  *phi2 = c_Teff * c_tau * PHI22b[i_T+1][i_tau+1] + c_Teff * (1 - c_tau) * PHI22b[i_T+1][i_tau]
    + (1 - c_Teff) * c_tau * PHI22b[i_T][i_tau+1] + (1 - c_Teff) * (1 - c_tau) * PHI22b[i_T][i_tau];

  *kion1  = c_Teff * c_tau * PHI22a_thin[i_T+1][i_tau+1] + c_Teff * (1 - c_tau) * PHI22a_thin[i_T+1][i_tau]
    + (1 - c_Teff) * c_tau * PHI22a_thin[i_T][i_tau+1]   + (1 - c_Teff) * (1 - c_tau) * PHI22a_thin[i_T][i_tau];

  *gamma1 = c_Teff * c_tau * PHI22b_thin[i_T+1][i_tau+1] + c_Teff * (1 - c_tau) * PHI22b_thin[i_T+1][i_tau]
    + (1 - c_Teff) * c_tau * PHI22b_thin[i_T][i_tau+1] + (1 - c_Teff) * (1 - c_tau) * PHI22b_thin[i_T][i_tau];

  //*kion1  = c_Teff * k_ion1[i_T+1]    + (1 - c_Teff) * k_ion1[i_T];
  //*gamma1 = c_Teff * GammaIon1[i_T+1] + (1 - c_Teff) * GammaIon1[i_T];

  if(c_Teff > 1.001 || c_Teff < -0.001 || c_tau > 1.001 || c_tau < -0.001){
    printf("minus!? c_T=%g c_tau=%g i_T=%d i_tau=%d\n", c_Teff, c_tau, i_T, i_tau);
    printf("tau=%g Tau[i]=%g Tau[i+1]=%g\n", tau, Tau[i_tau], Tau[i_tau+1]);
    printf("Teff=%g Teff[i]=%g Teff[i+1]=%g\n", Teff, Teffs[i_T], Teffs[i_T+1]);
    raise(11);
  }
  
  if(*kion1 < 0){
    printf("minus! c_T=%g c_tau=%g i_T=%d i_tau=%d\n", c_Teff, c_tau, i_T, i_tau);
    printf("PHI22a_thin=%g %g %g %g\n"
           , PHI22a_thin[i_T+1][i_tau+1], PHI22a_thin[i_T+1][i_tau]
           , PHI22a_thin[i_T][i_tau+1]  , PHI22a_thin[i_T][i_tau]);
  }
}

#ifdef RYTR_NH2
#define DELTA_NH2_CRIT 1.e14

void calc_shielding_factor(int ithis, double *k_diss_fac, double dr, double ri){
  double trans_len_cm;
  if(All.ComovingIntegrationOn)
    trans_len_cm =  All.Time / All.HubbleParam * All.UnitLength_in_cm ;  /* [comoving kpc/h --> physical cm] */
  else
    trans_len_cm =  1. / All.HubbleParam * All.UnitLength_in_cm;         /* [physical kpc/h --> physical cm] */
  
  double n_H2  = SphP[ithis].n_nuclei * SphP[ithis].H2I;
  double N_H2  = SphP[ithis].N_H2;
  double d_NH2 = SphP[ithis].d_NH2_0;

  if(ri < 0.001 * SphP[ithis].Hsml)
    ri = 0.001 * SphP[ithis].Hsml;

  if(dr < 1.e-3 * SphP[ithis].Hsml)
    dr = 1.e-3 * SphP[ithis].Hsml;

  double x_max = N_H2 / 1.e14;
  double x_min = (N_H2 - d_NH2) / 1.e14;
  double dx    = d_NH2 / 1.e14;

  double fac;
  if(x_max < 1.){
    fac = x_max - x_min;  
  }else if(x_min > 1.){
    fac = 4. * (pow(x_max, 0.25) - pow(x_min, 0.25));
  }else{
    fac = 4. * pow(x_max, 0.25) - x_min - 3.;
  }
  
  double f_shield = fac / dx;
  double denom = (ri*ri - ri*dr + 1./3 * dr*dr) * pow(trans_len_cm, 2);  

  double k_diss_fac1 = 1./ri/ri/trans_len_cm/trans_len_cm;
  if(x_max > 1.)
    k_diss_fac1 *= pow(x_max, -0.75);
  
  double k_diss_fac2 = f_shield / denom;

  double w1 = 1/ (1. + d_NH2/DELTA_NH2_CRIT); 
  double w2 = 1/ (1. + DELTA_NH2_CRIT/d_NH2); 

  if(n_H2 > 1.e-20){
    *k_diss_fac = w1 * k_diss_fac1 + w2 * k_diss_fac2;
  }else{
    *k_diss_fac = 0.;
  }

  if(ithis%10000 == 0 || ri < 0.00001){
    //printf("[%d] k_diss_fac=%g %g w=%g %g f_shield=%g x_max=%g x_min=%g f=%g N_H2/ri=%g\n", ithis, k_diss_fac1, k_diss_fac2, w1, w2, f_shield, x_max, x_min, pow(x_max, -0.75), N_H2/ri/trans_len_cm);
  }


  if(*k_diss_fac < 0){
    printf("k_diss_fac=%g %g w=%g %g\n", k_diss_fac1, k_diss_fac2, w1, w2);
    raise(11);
  }
}
#endif
#endif
 
 
