#pragma once

#ifndef CHEM_RATE_COUNT
#define CHEM_RATE_COUNT k_MAX
#endif

typedef struct {
  double rho_cgs;       /* g cm^-3 */
  double metallicity;

  double J21;
  double rate_Hminus_pd;
  double rate_H2_diss;
#ifdef HD
  double rate_HD_diss;
#endif
  double G0;
  double stellar_radiation_flux; /* erg cm^-2 s^-1 */
  
  double photoionization_rate_HI;
  double photoheating_rate_HI;

  double cmb_temperature; /* K */
  double scale_factor;
  int ID;
} ChemEnvironment;

typedef struct {
  double abund[TGCHEM_NUM_ABUNDANCES];
  double utherm_cgs; /* erg g^-1 */
  double gamma;
} ChemState;

typedef struct {
  double temperature;      /* K */
  double molecular_weight;
} ChemThermodynamicState;

typedef struct {
  double jeans_length; /* cm */
  double column_HI;    /* cm^-2 */
  double column_H2;    /* cm^-2 */
  double column_HD;    /* cm^-2 */

  double shield_HI;
  double shield_H2;
  double shield_H2self;
  double shield_HDself;

  double G0;
  double k22;
  double piHI;
} ChemStepData;

typedef struct {
  double tau;

  double escape_fraction;
  double reduction_factor;

  double dust_temperature;
  double dust_cooling_rate;
  double dust_kappa_planck;
} ChemContinuumState;

typedef struct {
  double forward[CHEM_RATE_COUNT];
  double inverse[CHEM_RATE_COUNT];
} ChemRates;

typedef enum {
  CHEMISTRY_SOLVER_EXPLICIT,
  CHEMISTRY_SOLVER_IMPLICIT,
  CHEMISTRY_SOLVER_SAHA
} ChemistrySolver;

typedef enum {
  ENERGY_SOLVER_EXPLICIT,
  ENERGY_SOLVER_IMPLICIT
} EnergySolver;

typedef struct {
  ChemistrySolver chemistry;
  EnergySolver energy;
  double timestep;
} ChemSolverPlan;

typedef struct {
  double rhoU_dot_rad;  /* erg cm^-3 s^-1 */
  double rhoU_dot_chem; /* erg cm^-3 s^-1 */

  double t_cool_rad;  /* s */
  double t_cool_chem; /* s */
  double t_chem;      /* s */

  double err_chem;
  int flag_implicit;
  int flag_convergence;
} ChemEvaluation;

typedef struct {
  double H2_escape_fraction;
#ifdef HD
  double HD_escape_fraction;
#endif
  
#if defined(DEBUG_CHEMISTRY_SC)
  double edot_H2, edot_comp, edot_ceHI, edot_radHI;
  double edot_ciHI;
  double edot_H2form_gas;

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
#endif
} ChemDiagnostics;

typedef struct {
  ChemThermodynamicState thermodynamics;
  ChemContinuumState continuum;
  ChemRates rates;
  ChemEvaluation evaluation;
  ChemDiagnostics diagnostics;
} ChemWorkspace;
