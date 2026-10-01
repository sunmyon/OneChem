#pragma once

#include "chemistry_build_config.h"
#include "species.h"

#include "include/chemistry/chemistry_status.h"

typedef struct {
  double H2_escape_fraction;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double HD_escape_fraction;
#endif
  double edot_rad, edot_chem;
  double edot_H2, edot_comp, edot_ceHI, edot_radHI;
  double edot_ciHI, edot_H2form_gas, edot_HMform;
  double edot_H2form_HM, edot_H2form_dust, edot_H2form_3body;
  double n_crit;
#ifdef CHEMISTRY_SPECIES_HELIUM
  double edot_ceHeI, edot_ceHeII;
  double edot_ciHeI, edot_ciHeII, edot_ciHeIS;
  double edot_reHeII1, edot_reHeII2, edot_reHeIII;
#endif
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double edot_HD;
#endif
#ifdef CHEMISTRY_METALLICITY
  double edot_CI, edot_CII, edot_OI;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  double edot_CO, edot_OH, edot_H2O;
#endif
#if CHEMISTRY_DUST_COOLING
  double edot_dust;
#ifdef CHEMISTRY_LW_RADIATION
  double edot_PE;
#endif
#endif
#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA
  double edot_cont;
#else
  double edot_reHI, edot_brem, edot_CIE, edot_cont;
#endif
} ChemDiagnostics;

typedef enum {
  CHEM_PROFILE_TOTAL,
  CHEM_PROFILE_ENERGY_SOLVER,
  CHEM_PROFILE_ENERGY_RESIDUAL,
  CHEM_PROFILE_IMPLICIT_CHEMISTRY,
  CHEM_PROFILE_REACTION_RATES,
  CHEM_PROFILE_INTERMEDIATES,
  CHEM_PROFILE_INTERMEDIATE_ASSEMBLY,
  CHEM_PROFILE_INTERMEDIATE_LINEAR_SOLVE,
  CHEM_PROFILE_NEWTON_RHS,
  CHEM_PROFILE_NEWTON_JACOBIAN,
  CHEM_PROFILE_JACOBIAN_SOLVER_COLUMNS,
  CHEM_PROFILE_JACOBIAN_RESPONSE,
  CHEM_PROFILE_JACOBIAN_INTERMEDIATE_SOLVE,
  CHEM_PROFILE_JACOBIAN_SCHUR,
  CHEM_PROFILE_NEWTON_LINEAR_SOLVE,
  CHEM_PROFILE_NEWTON_STEP,
  CHEM_PROFILE_CHEMICAL_COOLING,
  CHEM_PROFILE_RADIATIVE_COOLING,
  CHEM_PROFILE_MOLECULAR_COOLING,
  CHEM_PROFILE_DUST_TEMPERATURE,
  CHEM_PROFILE_COUNT
} ChemistryProfileSection;

typedef struct ChemReactionNetwork ChemReactionNetwork;
typedef struct ChemistryRateTables ChemistryRateTables;
typedef struct ChemistryCoolingTables ChemistryCoolingTables;
typedef struct ChemistryDustTables ChemistryDustTables;
typedef struct ChemistryContinuumOpacityTables ChemistryContinuumOpacityTables;

/* Copyable, non-owning inputs to evaluation code. Pointees are immutable;
 * the owning Model must outlive this view. */
typedef struct ChemRunView {
  const ChemReactionNetwork *network;
  const ChemistryRateTables *rate_tables;
  const ChemistryCoolingTables *cooling_tables;
  const ChemistryDustTables *dust_tables;
  const ChemistryContinuumOpacityTables *continuum_opacity_tables;
  double timescale_safety_factor;
#ifdef CHEMISTRY_FIXED_TEMPERATURE
  double fixed_temperature;
#endif
} ChemRunView;

typedef struct {
  double cmb_temperature; /* K */
  double scale_factor;
  double simulation_time;
  int comoving_integration;

  double background_J21;
  double background_G0;
  double background_rate_Hminus_pd;
  double background_rate_H2_diss;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double background_rate_HD_diss;
#endif
} ChemGlobalStep;

typedef struct {
  const ChemRunView *run;
  const ChemGlobalStep *global;

  double rho_cgs;       /* g cm^-3 */
  double metallicity;

  double J21;
  double rate_Hminus_pd;
  double rate_H2_diss;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double rate_HD_diss;
#endif
  double G0;
  double stellar_radiation_flux; /* erg cm^-2 s^-1 */

  double photoionization_rate_HI;
  double photoheating_rate_HI;

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

  int ID;
} ChemEnvironment;

typedef struct {
  double abund[N_SP];
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
#endif
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
  /* Non-owning reference used by solver helper paths. */
  const ChemReactionNetwork *network;
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
  /* Per-evaluation diagnostics; never shared between cells or threads. */
  ChemistryStatus status;
  int diagnostics_enabled;
  ChemThermodynamicState thermodynamics;
  ChemContinuumState continuum;
  ChemRates rates;
  ChemEvaluation evaluation;
  ChemDiagnostics diagnostics;
} ChemWorkspace;
