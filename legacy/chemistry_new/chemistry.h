#pragma once

#include "species_definitions.h"

#ifdef CHEMISTRY_SC
/* ----- chemistry part ------- */
#define H_number_fraction  0.76
#define He_number_fraction 0.06
#define D_number_fraction  2.25e-5

#ifndef CLIGHT
#define CLIGHT CLIGHT_REAL
#endif

#ifndef CLIGHT_REAL
#define CLIGHT_REAL 2.99792458e10
#endif


extern const double XH;
extern const double XHe;
extern const double XHD;

extern int count_print;

#define CHEM_DECLARE_SPECIES(symbol, name, charge, h, he, d, c, o) s##symbol,
enum SPECIES{
  CHEM_SPECIES_LIST(CHEM_DECLARE_SPECIES)
  sEND
};
#undef CHEM_DECLARE_SPECIES

#if defined(CARBON_OXYGEN_CHEMISTRY) && !defined(CHEM_INTERMEDIATE_EQUILIBRIUM)
enum { sSTORED_END = sEND };
#elif defined(CARBON_OXYGEN_CHEMISTRY)
enum { sSTORED_END = sO2 + 1 };
#elif defined(HD) && defined(DMINUS)
enum { sSTORED_END = sDM + 1 };
#elif defined(HD) && defined(HD_PLUS)
enum { sSTORED_END = sHDII + 1 };
#elif defined(HD)
enum { sSTORED_END = sHDI + 1 };
#else
enum { sSTORED_END = sHM + 1 };
#endif

#define N_SP sEND
#define TGCHEM_NUM_ABUNDANCES sSTORED_END

typedef struct chemdata_in
{
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

typedef struct
{
  double H2_escape_fraction;
#ifdef HD
  double HD_escape_fraction;
#endif
  double edot_rad, edot_chem;
  double edot_H2, edot_comp, edot_ceHI, edot_radHI;
  double edot_ciHI, edot_H2form_gas, edot_HMform;
  double edot_H2form_HM, edot_H2form_dust, edot_H2form_3body;
  double n_crit;
#ifdef HELIUM
  double edot_ceHeI, edot_ceHeII;
  double edot_ciHeI, edot_ciHeII, edot_ciHeIS;
  double edot_reHeII1, edot_reHeII2, edot_reHeIII;
#endif
#ifdef HD
  double edot_HD;
#endif
#ifdef METALLICITY
  double edot_CI, edot_CII, edot_OI;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  double edot_CO, edot_OH, edot_H2O;
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
  double edot_reHI, edot_brem, edot_CIE, edot_cont;
#endif
} ChemDiagnostics;

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
#ifdef CHEMISTRY_STANDALONE_HOST
  double tau_cont;
  double continuum_reduction_factor;
  double forward_rates[CHEM_RATE_COUNT];
  double inverse_rates[CHEM_RATE_COUNT];
  ChemDiagnostics diagnostics;
#endif
} chemdata_out;

double compute_gamma_from_abund(const double *abund, double T_gas);
double chemistry_molecular_weight(const double *abund);
double chemistry_temperature(const double *abund, double gamma, double utherm_cgs);

#endif
