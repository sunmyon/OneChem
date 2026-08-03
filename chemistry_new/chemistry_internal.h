#ifndef __CHEMISTRY
#define __CHEMISTRY

#include <math.h>
#include <stdio.h>
#include "chemistry.h"

#ifdef CHEMISTRY_SC

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

#ifndef dmax
#define dmax(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef dmin
#define dmin(a,b) ((a)<(b)?(a):(b))
#endif

extern double chem_e[15];
FILE *open_chemistry_table(const char *filename);

struct opacities {
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

extern double (*H2_cooling)(double, double, double, double, double, double, double);
#ifdef HD
extern double (*HD_cooling)(double, double, double, double);
#endif
extern double (*H2_line_escape_fraction)(double, double);
#ifdef CONT_OPACITY
extern double (*continuum_escape_fraction)(double, double, double, struct opacities *);
#endif
extern double (*partition_H2)(double);

extern int FlagPrint;
#endif
