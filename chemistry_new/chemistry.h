#pragma once

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

double compute_gamma_from_abund(const double *abund, double T_gas);
double chemistry_molecular_weight(const double *abund);
double chemistry_temperature(const double *abund, double gamma, double utherm_cgs);

#endif
