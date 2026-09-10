#pragma once

#include "cooling_func.h"

double dust_temperature(double rho_cgs, double temperature, double dust_metallicity, double length, double stellar_flux, struct opacities *opacities, double *cooling_rate);
double H2form_on_grains(double rho_cgs, double temperature, double dust_temperature, double dust_metallicity);
double dust_PEheating(double G0, double nH, double electron_abundance,
                      double temperature, double metallicity);
#ifdef DUST_COOLING
void read_dust_opacity_tables(void);
#endif
