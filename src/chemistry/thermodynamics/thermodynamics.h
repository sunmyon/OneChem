#pragma once

#include "chemistry_constants.h"
#include "chemistry_types.h"


double chemistry_temperature(const double *abund, double gamma, double utherm_cgs);
double chemistry_gamma_from_abundances(const double *abund, double temperature);
/* Legacy spelling retained for the current standalone and AREPO wrappers. */
double compute_gamma_from_abund(const double *abund, double temperature);
double chemistry_molecular_weight(const double *abund);
static inline double chemistry_hydrogen_number_density(double rho_cgs)
{
  return CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON * rho_cgs / PROTONMASS;
}
double compute_freefall_time(double rho_gas);
double chemistry_jeans_radius(double rho_cgs, double temperature,
                              double thermal_mass_in_proton_masses);
double get_core_radius(double rho_cgs, double temperature);
void update_thermodynamic_state(const ChemState *state, ChemThermodynamicState *thermodynamics);
