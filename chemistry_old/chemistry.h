#ifndef __CHEMISTRY
#define __CHEMISTRY

#include <math.h>

#define N_pts_bound    7
#define N_nu_base      200
#ifdef X_RAY
#define N_nu_bound     18 /*12*/
#else
#define N_nu_bound     13
#endif
#define N_T            10000
#define N_nu           N_nu_base + N_nu_bound*N_pts_bound /*  284 */

#ifdef CHEMISTRY
/* ----- chemistry part ------- */
#define H_number_fraction  0.76
#define He_number_fraction 0.06
#define D_number_fraction  2.25e-5


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
extern double  T[N_T],J0_nu[N_nu],J_nu[N_nu],nu[N_nu];
extern double  k_table[k_MAX][N_T];
extern double  k_table_inverse[k_MAX][N_T];
extern double  ciHIa[N_T],ciHeIa[N_T],ciHeIIa[N_T],ciHeISa[N_T],reHIIa[N_T],brema[N_T];
extern double  ceHIa[N_T],ceHeIa[N_T],ceHeIIa[N_T],reHeII1a[N_T],reHeII2a[N_T],reHeIIIa[N_T];

/* cross-sections */
#if defined(BG_RADIATION) || defined(IONIZATION_PHOTON)
extern double sigma[s_MAX][N_nu];
extern double SIGMA_e22, SIGMA_e23, SIGMA_e24;
#endif

extern double g_weight[20];
extern double n_pop[3][20];

extern double A[3][20][3][20], B[3][20][3][20];
extern double E_vJ[3][20], E_vJ_erg[3][20];
extern double E_vJ_eV[3][20], E_vJ_K[3][20];
extern double del_E[3][20][3][20],del_E_K[3][20][3][20],del_nu[3][20][3][20];

extern double Lambda_frac[3][20][3][20], Lambda_frac_thick[3][20][3][20];
extern double alpha[3][20][3][20], n_pop_dummy[3][20][3][20],alpha_h2[3][20][3][20];
#endif


/* ----- Declare some physical constants (cgs). */
#define eV_to_K      11606.0
#define eV_to_erg    1.60184e-12
#define eV_to_Hz     2.41838e14
#define eV_to_um_inv 0.80669

#define kJ_per_mol_to_eV 1.04e-2
#define CLIGHT      2.9979e10

#define Z_SUN 0.02
#define DUST_TO_GAS_MASS_RATIO 0.00934

extern const double XH;
extern const double XHe;
extern const double XHD;
extern double Metallicity;

/* ----- Set various cutoff values in eV. */

#define  e22     13.6
#define  e23     54.4
#define  e24     24.6
#define  e25     0.755
#define  e26a    2.65
#define  e26b    11.27
#define  e26c    21.0
#define  e27a    15.42
#define  e27b    16.5
#define  e27c    17.7
#define  e28a    30.0
#define  e28b    70.0
#define  e29a    200.0
#define  e29b    500.0
#define  e29c    1.e3
#define  e29d    2.e3
#define  e29e    5.e3
#define  e29f    1.e4
#define  e30     5.12

extern double chem_e[15];

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

#endif
