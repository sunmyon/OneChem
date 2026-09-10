#pragma once

#include <stdint.h>

#include "chemistry_newton_internal.h"

/* The declaration order is the deterministic basin-search order.  Keep the
 * strategies contiguous because the Newton driver advances with reset++.
 * The two non-C/O strategies intentionally occupy the same positions as the
 * C/O-rich strategies, preserving the established build-specific ordering. */
typedef enum {
  CHEM_NEWTON_RESET_NONE = 0,
  CHEM_NEWTON_RESET_INITIAL_H2,
#ifdef CARBON_OXYGEN_CHEMISTRY
  CHEM_NEWTON_RESET_CO_RICH,
  CHEM_NEWTON_RESET_CO2_RICH,
  CHEM_NEWTON_RESET_SAVED_Z1_T3127,
  CHEM_NEWTON_RESET_SAVED_Z0P01_T32465,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T34463,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T30811,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T1261,
  CHEM_NEWTON_RESET_SAVED_Z1_T92675,
  CHEM_NEWTON_RESET_SAVED_Z1_T1501,
  CHEM_NEWTON_RESET_SAVED_Z1_T1515,
  CHEM_NEWTON_RESET_SAVED_Z1_T1308,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T1625,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T13492,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T2817,
  CHEM_NEWTON_RESET_SAVED_Z1_T3137,
  CHEM_NEWTON_RESET_SAVED_Z0P1_T2055,
  CHEM_NEWTON_RESET_MOLECULAR_RESERVOIR_CORNER,
  CHEM_NEWTON_RESET_CI_H2O_CORNER,
  CHEM_NEWTON_RESET_CO_OI_CORNER,
  CHEM_NEWTON_RESET_SAVED_SAMPLE_88,
  CHEM_NEWTON_RESET_SAVED_SAMPLE_16123,
  CHEM_NEWTON_RESET_GRID_SAMPLE_24548,
  CHEM_NEWTON_RESET_H2_DI_CO_OI_CORNER,
  CHEM_NEWTON_RESET_H2_DI_CI_OI_CORNER,
  CHEM_NEWTON_RESET_H2_CO_OI_BROAD_D_CORNER,
  CHEM_NEWTON_RESET_COLD_CH2_H2O_CORNER,
  CHEM_NEWTON_RESET_SAVED_SAMPLE_3624,
  CHEM_NEWTON_RESET_COLD_CO_O_O2_CORNER,
  CHEM_NEWTON_RESET_MINIMUM_COVER_NEUTRAL_H_DII,
  CHEM_NEWTON_RESET_MINIMUM_COVER_SECOND,
  CHEM_NEWTON_RESET_LAST = CHEM_NEWTON_RESET_MINIMUM_COVER_SECOND
#else
  CHEM_NEWTON_RESET_FULLY_IONIZED_H,
  CHEM_NEWTON_RESET_MOLECULAR_H,
  CHEM_NEWTON_RESET_LAST = CHEM_NEWTON_RESET_MOLECULAR_H
#endif
} ChemNewtonResetIndex;

void chemistry_newton_apply_deterministic_reset(ChemNewtonResetIndex reset,
                                                double metallicity,
                                                const double *initial_abundance,
                                                double *abundance);

uint64_t chemistry_newton_random_seed(double nH, double temperature, double dt);

int chemistry_find_random_newton_guess(
    double nH, double metallicity, const double *initial_abundance,
    double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation, uint64_t *random_state,
    int *iteration, int max_iterations, int *random_attempt_count);
