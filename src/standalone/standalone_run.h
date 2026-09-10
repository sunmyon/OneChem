#ifndef STANDALONE_RUN_H
#define STANDALONE_RUN_H

#include "standalone_chemistry_compat.h"

enum { STANDALONE_CELL_COUNT = 1 };

typedef struct sph_particle_data
{
  double DMdensity;
  double Density;
  double DtDensity;
  double Entropy;
  double DtEntropy;
  double n_nuclei;
  double temperature;
  double temp_hold;
  double Gamma, molecular_weight;
  double pressure, ThermalEnergy;

  double elec;
  double HI, HII, HM;
#ifdef HELIUM
  double HeI, HeII, HeIII;
#endif
  double H2I;
#ifdef H2_PLUS
  double H2II;
#endif
#ifdef HD
  double DI, DII, HDI;
#ifdef HD_PLUS
  double HDII;
#endif
#ifdef DMINUS
  double DM;
#endif
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  double CI, CII, CO, CO2;
  double OI, OH, H2O, O2;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  double COII, CH, CH2;
  double OII, OHII, H2OII, H3OII, O2II;
#endif
#endif

  double t_col, t_cool, t_chem;
  double tcool_rad, tcool_chem;

  double edot_ceHI, edot_ceHeI, edot_ceHeII;
  double edot_ciHI, edot_ciHeI, edot_ciHeIS, edot_ciHeII;
  double edot_ciHI0;
  double edot_reHII, edot_reHeII1, edot_reHeII2, edot_reHeIII;
  double edot_comp, edot_brem, edot_CIE, edot_H2, edot_HD;
  double edot_H2form_3body, edot_H2form_gas, edot_chem;
  double edot_H2form_HM, edot_H2form_H2II, edot_H2form_dust;
  double edot_adia, edot_adm, edot_grav;
  double edot_radHI, edot_radHeI, edot_radHeII;
  double edot_cont, edot_2ph;
  double edot_CI, edot_CII, edot_OI;
#ifdef CARBON_OXYGEN_CHEMISTRY
  double edot_CO, edot_OH, edot_H2O;
#endif
  double edot_dust, edot_PE;
  double edot_rad;

  double edot_H2form, edot_HMform;
  double edot_H2form_HM_loss, edot_H2form_dust_loss;

  double f_reduce;
  double k[k_MAX], k_inv[k_MAX];
  double n_crit;
  double Jeans_Length;

#ifdef SHIELD_6RAY
  double ColumnDens_H2[6], ColumnDens_HI[6], ColumnDens_HD[6];
#endif
#ifdef H2_OPACITY
  double DelV[3];
#endif
  double tau_cont;
  double fesc_H2;
#ifdef HD
  double fesc_HD;
#endif

#ifdef LYMANWERNER
  FLOAT J21;
  FLOAT rate_Hminus_pd;
  FLOAT rate_H2_diss;
#ifdef HD
  FLOAT rate_HD_diss;
#endif
#ifdef DUST_COOLING
  FLOAT G0;
#endif
#endif

#ifdef THERMAL_FEEDBACK
  double stellar_radiation_flux;
#endif

#ifdef DENSITY_GRAD
  FLOAT GradRho[3];
  FLOAT l_core;
#endif

#if defined(IONIZATION_PHOTON) || defined(X_RAY)
  FLOAT d_tau_0, d_r, dist_source;
  FLOAT tau_ion, d_tau;
  double k22, piHI;
#ifdef X_RAY
  double k24, piHeI;
#ifdef X_RAY_2nd
  double k22_2nd, k24_2nd, piHI_2nd;
#endif
#endif
  FLOAT n_nuclei_ch;
  int flag_ch;
#endif

  double Hsml;
  double SHR;
  double Entropy_old;
  int flag_div;
  double Tdust;
  double Metallicity;
} StandaloneGasCell;

typedef struct particle_data
{
  FLOAT Pos[3];
  FLOAT Mass;
  MyIDType ID;
  int Type;
  int Ti_endstep;
  int Ti_begstep;
} StandaloneParticle;

extern StandaloneGasCell SphP[STANDALONE_CELL_COUNT];
extern StandaloneParticle P[STANDALONE_CELL_COUNT];

/*
 * Process-wide state owned by the standalone executable.
 *
 * This is deliberately not part of the chemistry core. A hydrodynamics
 * host supplies the corresponding values through its adapter instead.
 */
typedef struct StandaloneRunState
{
  double G;

  double UnitTime_in_s;
  double UnitMass_in_g;
  double UnitVelocity_in_cm_per_s;
  double UnitLength_in_cm;
  double UnitPressure_in_cgs;
  double UnitDensity_in_cgs;
  double UnitEnergy_in_cgs;
  /* Retained only for the historical chemistry implementation. */
  double MinEgySpec;

  double Hubble;
  double Omega0;
  double OmegaLambda;
  double OmegaBaryon;
  double HubbleParam;

  int ComovingIntegrationOn;
  char ChemistryTableDir[1024];
  double ScaleFactor;
  double Time;

#if defined(STATICNFW) || defined(PHYS_COSMO)
  double ExpansionParam;
#endif

  double Epsilon;

#ifdef LYMANWERNER
  double J21;
  double LymanWerner;
  FLOAT rate_Hminus_pd;
  FLOAT rate_H2_diss;
#ifdef HD
  FLOAT rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double G0;
#endif
#endif

#ifdef X_RAY
  double X_J21;
#endif
#ifdef SHIELD_LW
  double LenColumnDens;
#endif

  double T_CMB;
  double cf_atime;
  double cf_hubble_a;
  double cf_a3inv;
  double Metallicity;
} StandaloneRunState;

extern StandaloneRunState StandaloneRun;

/* Temporary source compatibility while the standalone driver is renamed. */
#define All StandaloneRun

#endif
