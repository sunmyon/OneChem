#ifndef __CHEMISTRY
#define __CHEMISTRY

#include <math.h>
#include <stdio.h>

#ifdef CHEMISTRY_SC
/* ----- chemistry part ------- */
#define H_number_fraction  0.76
#define He_number_fraction 0.06
#define D_number_fraction  2.25e-5

#define N_T            10000
#define k_MAX      60
#ifndef T_CMB0
#define T_CMB0     2.728
#endif

#ifdef METALLICITY
/* number fraction of heavy elements relative to H nuclei for the SOLAR photosphere */
/* Caffau+ 2011 Nat */
/* added by G. Chiaki Apr 20 2013 */
#define Y_C_sun  0.927e-4
#define Y_O_sun  3.568e-4

#define A_C_sun  3.16e-4
#define A_O_sun  5.75e-4
#define A_Mg_sun 3.47e-5
#define A_Si_sun 3.31e-5
#define A_Fe_sun 3.31e-5
#define Zmet_sun 0.0153         
#endif


/* ----- Tables ------- */
extern double  T_table[N_T];
extern double  k_table[k_MAX][N_T];
extern double  k_table_inverse[k_MAX][N_T];
extern double  ciHIa[N_T],ciHeIa[N_T],ciHeIIa[N_T],ciHeISa[N_T],reHIIa[N_T],brema[N_T],H2_CIE[N_T];
extern double  ceHIa[N_T],ceHeIa[N_T],ceHeIIa[N_T],reHeII1a[N_T],reHeII2a[N_T],reHeIIIa[N_T];
#endif


/* ----- Declare some physical constants (cgs). */
#define eV_to_K      11606.0
#define eV_to_erg    1.60184e-12
#define eV_to_Hz     2.41838e14
#define eV_to_um_inv 0.80669

#define kJ_per_mol_to_eV 1.04e-2

#define Z_SUN 0.02
#define DUST_TO_GAS_MASS_RATIO 0.00934
#define STEFAN 5.6692e-5  /* erg cm-2 K-4 s-1 */
#ifndef CLIGHT_REAL
#define CLIGHT_REAL 2.99792458e10
#endif
#ifndef CLIGHT
#define CLIGHT CLIGHT_REAL
#endif

#ifndef dmax
#define dmax(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef dmin
#define dmin(a,b) ((a)<(b)?(a):(b))
#endif

extern const double XH;
extern const double XHe;
extern const double XHD;

extern double chem_e[15];
FILE *open_chemistry_table(const char *filename);

enum SPECIES{
	sELEC,
	sHI,
	sHII,
#ifdef HELIUM
	sHeI,
	sHeII,
	sHeIII,
#endif
	sH2I,
#ifdef H2_PLUS
	sH2II,
#endif
	sHM,
#ifdef HD
	sDI,
	sDII,
	sHDI,
#ifdef HD_PLUS
	sHDII,
#endif
#ifdef DMINUS
	sDM,
#endif
#endif
	sEND
};

struct opacities{
  double tau;
  double tau_Rossland;
  double tau_Planck;
  double kappa_Rossland;
  double kappa_Planck;
#ifdef DUST_COOLING
  double kappa_dust_Rossland;
  double kappa_dust_Planck;
#endif
  double x;
  double f_reduce;
};

typedef struct chemdata_in
{
  double Mass;
  MyIDType ID;
  double rho_cgs;
  int flag_print;
  double Density;
#ifdef METALLICITY
  double Metallicity;
#endif

#ifdef THERMAL_FEEDBACK
  double stellar_rad_flux;
#endif

#ifdef LW_RADIATION_LOCAL
  double rate_Hminus_pd;
  double rate_H2_diss;
  double J21;
#ifdef HD
  double rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double G0;
#endif
#endif

#ifdef IONIZATION_PHOTON
  double k22;
  double piHI;
#endif
} chemdata_in;

typedef struct chemdata_out
{
  /* primitive variables */
  double Utherm;
  
  double t_cool;
  double t_chem;
  double Tdust;
  
  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;

  char flag_implicit;
  char conv_flag;
  double kappa_dust;
} chemdata_out;


typedef struct chemdata_state{
  double t_cool_rad, t_cool_chem;
  double rhoU_dot_rad, rhoU_dot_chem;
  double err_chem;
  
#ifdef SHIELD_LW
  double f_byH2self;
#ifdef HD
  double f_byHDself;
#endif
  double f_byHI;
  double f_byH2;
#endif

#ifdef H2_OPACITY
  double EscapeProb;
#ifdef HD
  double EscapeProb_HD;
#endif
#endif
  double tau_cont;
  double f_reduce;
  double Jeans_Length;

#if defined(LW_RADIATION_LOCAL) && defined(DUST_COOLING)
  double G0;
#endif
#ifdef DUST_COOLING
  double dust_cooling_rate;
#endif
  double k22;
  double piHI;
  
#if defined(DEBUG_CHEMISTRY_SC)
  double edot_H2, edot_comp, edot_ceHI, edot_radHI;
  double edot_ciHI;
  double edot_H2form_gas;

#ifdef HELIUM
  double edot_ceHeI, edot_ceHeII, edot_ciHeI;
  double edot_ciHeI, edot_ciHeII, edot_ciHeIS;
  double edot_reHeII1, edot_reHeII2, edot_reHeIII;
#endif
  
#ifdef HD
  double edot_HD;
#endif
#ifdef METALLICITY
  double edot_CI, edot_CII, edot_OI;
#endif
#ifdef DUST_COOLING
  double edot_dust;
#ifdef LYMANWERNER
  double edot_PE;
#endif
#endif

#ifdef CONTINUUM_COOLING_MATSUKOBA
  double edot_cont;
#else
  double edot_reHI;
  double edot_brem;
  double edot_CIE;
#endif
  
  double k[70], k_inv[70];
#endif
} chemdata_state;



extern double (*H2_cooling)();
#ifdef HD
extern double (*HD_cooling)();
#endif
extern double (*H2_line_escape_fraction)();
#ifdef CONT_OPACITY
extern double (*continuum_escape_fraction)();
#endif
extern double (*partition_H2)();

extern int FlagPrint;
extern int count_print;
#endif
