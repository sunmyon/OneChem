#pragma once

/*
 * Legacy data-transfer structures owned by the standalone host boundary.
 * MyIDType is intentionally host-defined; include this header only after the
 * host has declared it (normally through allvars.h).
 */
#include "chemistry_build_config.h"
#include "chemistry_types.h"

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
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double rate_HD_diss;
#endif
#if CHEMISTRY_DUST_COOLING
  double G0;
#endif
#endif

#ifdef IONIZING_PHOTON
  double k22;
  double piHI;
#endif
#ifdef CHEMISTRY_SPECIES_HELIUM
  double photoionization_rate_HeI;
  double photoheating_rate_HeI;
  double photoionization_rate_HeII;
  double photoheating_rate_HeII;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  double photoionization_rate_C[3];
  double photoheating_rate_C[3];
  double photoionization_rate_O[3];
  double photoheating_rate_O[3];
  double rate_CO_diss;
#endif
} chemdata_in;

typedef struct chemdata_out
{
  double Utherm;
  double t_cool;
  double t_chem;
  double Tdust;
  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;
  char flag_implicit;
  char conv_flag;
  double kappa_dust;

  /* Standalone diagnostics are not part of the portable chemistry ABI. */
  double tau_cont;
  double continuum_reduction_factor;
  double forward_rates[CHEM_RATE_COUNT];
  double inverse_rates[CHEM_RATE_COUNT];
  ChemDiagnostics diagnostics;
} chemdata_out;
