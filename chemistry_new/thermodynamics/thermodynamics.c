/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"

#include "chemistry.h"
#include "chemistry_models.h"
#include "chemistry_types.h"
#include "cooling_func.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

double chemistry_temperature(const double *abund, double gamma, double utherm_cgs)
{
  double molecular_weight = chemistry_molecular_weight(abund);
  double egyspec_to_temp  = (gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;

  return utherm_cgs * egyspec_to_temp;
}

double compute_gamma_from_abund(const double *abund, double T_gas)
{
  double n_H= abund[sHI] + abund[sHII] + abund[sHM] + 2*abund[sH2I];

  double y_HI   = abund[sHI]   /n_H;
  double y_HII  = abund[sHII]  /n_H;
  double y_He   = XHe;

  double y_HM   = abund[sHM]   /n_H;
  double y_H2   = abund[sH2I]  /n_H;
  double y_elec = abund[sELEC] /n_H;

  double x = 6100.0/T_gas;

  double gamma_H2_minus1_inv;
#ifdef GAMMA_NEW
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

  double sum_y= y_HI + y_HII + y_HM + y_He + y_elec + y_H2;
  double denominator=1.5*(y_HI + y_HII + y_HM + y_He + y_elec) + gamma_H2_minus1_inv * y_H2;

  double gamma = 1.0+sum_y/denominator;

  if(gamma > 1.8 || gamma < 1.2){
    printf("Gamma odd![task%03d] gamma=%g T_gas=%g\n"
	   , ThisTask, gamma, T_gas);
    fflush(stdout);
  }

  return gamma;
}

double chemistry_molecular_weight(const double *abund){
  double sum_fraction =  abund[sELEC] + abund[sHI] + abund[sHII] + abund[sHM] + abund[sH2I]
#ifdef H2_PLUS
    + abund[sH2II]
#endif
#ifdef HD
    + abund[sDI] + abund[sDII] + abund[sHDI]
#ifdef HD_PLUS
    + abund[sHDII]
#endif
#ifdef DMINUS
    + abund[sDM]
#endif
#endif
#ifdef HELIUM
    + abund[sHeI] + abund[sHeII] + abund[sHeIII]
#else
    + XHe
#endif
    ;

  return 1.0 / sum_fraction;
}

double compute_freefall_time(double rho_gas){
  return sqrt(3.0*M_PI/32.0/GRAVITY/rho_gas);
}

double get_core_radius(double rho_cgs, double T_gas){
#ifdef JEANS_LENGTH_TMAX
  if(T_gas > JEANS_LENGTH_TMAX)
    T_gas = JEANS_LENGTH_TMAX;
#endif

  return  0.5 * M_PI * sqrt(BOLTZMANN * T_gas/ PROTONMASS / 1.2)/ sqrt(GRAVITY * rho_cgs);
}

void update_thermodynamic_state(const ChemState *state, ChemThermodynamicState *thermodynamics)
{
  thermodynamics->molecular_weight = chemistry_molecular_weight(state->abund);
  thermodynamics->temperature = state->utherm_cgs * (state->gamma - 1.0) *
      PROTONMASS * thermodynamics->molecular_weight / BOLTZMANN;
}


#endif /* CHEMISTRY_SC */
