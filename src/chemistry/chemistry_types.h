#pragma once

#include "chemistry_build_config.h"
#include "chemistry_core.h"
#include "common/chemistry_diagnostics.h"
#include "include/chemistry/chemistry_status.h"

typedef struct ChemReactionNetwork ChemReactionNetwork;
typedef struct ChemistryRateTables ChemistryRateTables;
typedef struct ChemistryCoolingTables ChemistryCoolingTables;
typedef struct ChemistryDustTables ChemistryDustTables;
typedef struct ChemistryContinuumOpacityTables ChemistryContinuumOpacityTables;

typedef struct ChemRunContext {
  /* Owned by this run context; immutable after InitChem succeeds. */
  ChemReactionNetwork *network;
  ChemistryRateTables *rate_tables;
  ChemistryCoolingTables *cooling_tables;
  ChemistryDustTables *dust_tables;
  ChemistryContinuumOpacityTables *continuum_opacity_tables;
  double timescale_safety_factor;
#ifdef FIX_TEMPERATURE
  double fixed_temperature;
#endif
} ChemRunContext;

typedef struct {
  double cmb_temperature; /* K */
  double scale_factor;
  double simulation_time;
  int comoving_integration;

  double background_J21;
  double background_G0;
  double background_rate_Hminus_pd;
  double background_rate_H2_diss;
#ifdef HD
  double background_rate_HD_diss;
#endif
} ChemGlobalStep;

typedef struct {
  const ChemRunContext *run;
  const ChemGlobalStep *global;

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

#ifdef HELIUM
  /* Future host-radiation inputs; currently not consumed by the network. */
  double photoionization_rate_HeI;
  double photoheating_rate_HeI;
  double photoionization_rate_HeII;
  double photoheating_rate_HeII;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
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
