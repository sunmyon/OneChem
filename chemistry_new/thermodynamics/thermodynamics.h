#pragma once

#include "chemistry.h"
#include "chemistry_types.h"

double chemistry_temperature(const double *abund, double gamma, double utherm_cgs);
double compute_gamma_from_abund(const double *abund, double temperature);
double chemistry_molecular_weight(const double *abund);
double compute_freefall_time(double rho_gas);
double get_core_radius(double rho_cgs, double temperature);
void update_thermodynamic_state(const ChemState *state, ChemThermodynamicState *thermodynamics);
