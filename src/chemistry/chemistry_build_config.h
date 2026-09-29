#pragma once


/* Preserve the host's historical option spelling. */
#if defined(IONIZING_PHOTON) && !defined(IONIZATION_PHOTON)
#define IONIZATION_PHOTON
#endif

#if defined(CARBON_OXYGEN_CHEMISTRY) && !defined(METALLICITY)
#error "CARBON_OXYGEN_CHEMISTRY requires METALLICITY"
#endif

/* The implicit energy residual evaluates chemistry implicitly at every trial
 * temperature.  Keep the outer timestep policy consistent even when a host
 * specifies only ENERGY_IMPLICIT. */
#if defined(ENERGY_IMPLICIT) && !defined(CHEMISTRY_IMPLICIT)
#define CHEMISTRY_IMPLICIT
#endif

#if defined(CARBON_OXYGEN_HIGH_IONIZATION) && !defined(CARBON_OXYGEN_CHEMISTRY)
#error "CARBON_OXYGEN_HIGH_IONIZATION requires CARBON_OXYGEN_CHEMISTRY"
#endif

#ifndef CHEM_CO_O2_DEPENDENT_DENSITY
#define CHEM_CO_O2_DEPENDENT_DENSITY 1.e9
#endif

#if defined(METAL) || defined(CARBON_OXYGEN_CHEMISTRY)
#define CHEM_RATE_COUNT 200
#else
#define CHEM_RATE_COUNT 60
#endif

#if defined(H2_COOLING_G15) + defined(H2_COOLING_GP12) + defined(H2_COOLING_Y06) + defined(H2_COOLING_O00) > 1
#error "Select only one H2 cooling model"
#endif
#if defined(PARTITION_H2_EBM1989) + defined(PARTITION_H2_I87) + defined(PARTITION_H2_R85) > 1
#error "Select only one H2 partition-function model"
#endif
#if defined(H2_COOLING_Y06) || defined(H2_COOLING_O00)
#error "The selected H2 cooling model has no implementation in the portable chemistry core"
#endif
#if defined(H2_OPACITY) && !defined(H2_OPACITY_FIT)
#error "H2_OPACITY requires an implemented H2 line-escape model"
#endif
