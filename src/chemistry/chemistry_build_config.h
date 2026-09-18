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
