/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"

#include "chemistry_types.h"
#include "cooling_func.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "thermodynamics.h"


double chemistry_temperature(const double *abund, double gamma, double utherm_cgs)
{
  double molecular_weight = chemistry_molecular_weight(abund);
  double egyspec_to_temp  = (gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;

  return utherm_cgs * egyspec_to_temp;
}

double chemistry_gamma_from_abundances(const double *abund, double T_gas)
{
  double n_H= abund[sHI] + abund[sHII] + abund[sHM] + 2*abund[sH2I];

  double y_HI   = abund[sHI]   /n_H;
  double y_HII  = abund[sHII]  /n_H;
  double y_He   = XHe;

  double y_HM   = abund[sHM]   /n_H;
  double y_H2   = abund[sH2I]  /n_H;
  double y_elec = abund[sELEC] /n_H;

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  double y_metal = (abund[sCI] + abund[sCII] + abund[sCO] + abund[sCO2] +
                    abund[sOI] + abund[sOH] + abund[sH2O] + abund[sO2]) / n_H;
#else
  double y_metal = 0.0;
#endif

  double x = 6100.0/T_gas;

  double gamma_H2_minus1_inv;
#if CHEMISTRY_H2_ROTATIONAL_GAMMA
  if(T_gas < 550.){
    gamma_H2_minus1_inv
      = 1./(0.133333333333 * tanh((2.08661-log10(T_gas))*4.165) + 0.533333333333);
  }else{
    gamma_H2_minus1_inv = 0.5*(5.0 + 2.0*x*x * exp(x)/(exp(x)-1.0)/(exp(x)-1.0));
  }
#else
  if(x > 100.0)
    gamma_H2_minus1_inv = 2.5;
  else
    gamma_H2_minus1_inv = 0.5*(5.0 + 2.0*x*x * exp(x)/(exp(x)-1.0)/(exp(x)-1.0));
#endif

  double sum_y= y_HI + y_HII + y_HM + y_He + y_elec + y_H2 + y_metal;
  double denominator=1.5*(y_HI + y_HII + y_HM + y_He + y_elec + y_metal) + gamma_H2_minus1_inv * y_H2;

  double gamma = 1.0+sum_y/denominator;

  if(gamma > 1.8 || gamma < 1.2){
    printf("Gamma odd![chemistry] gamma=%g T_gas=%g\n", gamma, T_gas);
    fflush(stdout);
  }

  return gamma;
}

double compute_gamma_from_abund(const double *abund, double T_gas)
{
  return chemistry_gamma_from_abundances(abund, T_gas);
}

double chemistry_molecular_weight(const double *abund){
  double sum_fraction = 0.0;
  for(int species = 0; species < sSTORED_END; species++)
    sum_fraction += abund[species];

  /* Without an evolved He network, helium remains a fixed neutral particle
   * reservoir and therefore is not represented in abund[]. */
#ifndef CHEMISTRY_SPECIES_HELIUM
  sum_fraction += XHe;
#endif

  return 1.0 / sum_fraction;
}

double compute_freefall_time(double rho_gas){
  return sqrt(3.0*M_PI/32.0/GRAVITY/rho_gas);
}

double chemistry_jeans_radius(double rho_cgs, double T_gas,
                              double thermal_mass_in_proton_masses)
{
#ifdef CHEMISTRY_JEANS_TEMPERATURE_MAX
  if(T_gas > CHEMISTRY_JEANS_TEMPERATURE_MAX)
    T_gas = CHEMISTRY_JEANS_TEMPERATURE_MAX;
#endif

  return 0.5 * M_PI *
         sqrt(BOLTZMANN * T_gas / PROTONMASS /
              thermal_mass_in_proton_masses) /
         sqrt(GRAVITY * rho_cgs);
}

double get_core_radius(double rho_cgs, double T_gas){
  return chemistry_jeans_radius(rho_cgs, T_gas, 1.2);
}

void update_thermodynamic_state(const ChemState *state, ChemThermodynamicState *thermodynamics)
{
  thermodynamics->molecular_weight = chemistry_molecular_weight(state->abund);
  thermodynamics->temperature = state->utherm_cgs * (state->gamma - 1.0) *
      PROTONMASS * thermodynamics->molecular_weight / BOLTZMANN;
}
