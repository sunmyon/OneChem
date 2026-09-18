#ifndef CHEMISTRY_STATE_H
#define CHEMISTRY_STATE_H

#include <stddef.h>
#include <stdint.h>

typedef struct ChemistryGlobalEnvironment {
  double cmb_temperature;             /* K */
  double scale_factor;
  double simulation_time;             /* host-defined diagnostic coordinate */
  int comoving_integration;
  double background_J21;
  double background_G0;
  double background_rate_Hminus_pd;   /* s^-1 */
  double background_rate_H2_diss;     /* s^-1 */
  double background_rate_HD_diss;     /* s^-1 */
} ChemistryGlobalEnvironment;

/* Prepared once per host chemistry sweep and shared read-only by all cells. */
typedef struct ChemistryStepContext {
  ChemistryGlobalEnvironment global;
} ChemistryStepContext;

typedef struct ChemistryCellEnvironment {
  double density_cgs;                  /* g cm^-3 */
  double metallicity;                  /* solar units */

  double local_J21;
  double local_G0;
  double local_rate_Hminus_pd;         /* s^-1 */
  double local_rate_H2_diss;           /* s^-1 */
  double local_rate_HD_diss;           /* s^-1 */

  double photoionization_rate_HI;      /* s^-1 */
  double photoheating_rate_HI;         /* erg s^-1 per HI */

  /* Species-specific radiation inputs supplied by the host transport. */
  double photoionization_rate_HeI;     /* s^-1 */
  double photoheating_rate_HeI;        /* erg s^-1 per HeI */
  double photoionization_rate_HeII;    /* s^-1 */
  double photoheating_rate_HeII;       /* erg s^-1 per HeII */
  double photoionization_rate_C[3];    /* CI, CII, CIII; s^-1 */
  double photoheating_rate_C[3];       /* erg s^-1 per target ion */
  double photoionization_rate_O[3];    /* OI, OII, OIII; s^-1 */
  double photoheating_rate_O[3];       /* erg s^-1 per target ion */
  double local_rate_CO_diss;           /* s^-1 */

  double stellar_radiation_flux;       /* erg cm^-2 s^-1 */

  int diagnostics_enabled;
  uint64_t diagnostic_id;
} ChemistryCellEnvironment;

/*
 * The context defines the required species ordering and count.  Pointer plus
 * count preserves the current feature-dependent layout without exposing a
 * host-controlled fixed array in the public ABI.
 */
typedef struct ChemistryState {
  double *abundances;
  size_t abundance_count;
  double specific_thermal_energy_cgs;  /* erg g^-1 */
  double gamma;
} ChemistryState;

typedef struct ChemistryResult {
  double cooling_time_s;
  double chemistry_time_s;
  double thermal_energy_change_cgs;    /* erg g^-1 */
  double dust_temperature;             /* K */
  double dust_kappa_planck;
  double continuum_optical_depth;
  double continuum_reduction_factor;
  int implicit_flag;
  int convergence_flag;
  int used_implicit_chemistry;
  int used_implicit_energy;
  int converged;
} ChemistryResult;

#endif
