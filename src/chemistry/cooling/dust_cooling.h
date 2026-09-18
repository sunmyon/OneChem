#pragma once

#include "cooling_func.h"

typedef struct ChemistryDustTables ChemistryDustTables;

double dust_temperature(const ChemistryDustTables *tables, double rho_cgs,
                        double temperature, double radiation_temperature,
                        double dust_metallicity, double length,
                        double stellar_flux, struct opacities *opacities,
                        double *cooling_rate, int diagnostics_enabled);
double H2form_on_grains(double rho_cgs, double temperature, double dust_temperature, double dust_metallicity);
double dust_PEheating(double G0, double nH, double electron_abundance,
                      double temperature, double metallicity,
                      int diagnostics_enabled);
#ifdef DUST_COOLING
int dust_opacity_tables_create(const char *table_root,
                               ChemistryDustTables **tables);
void dust_opacity_tables_destroy(ChemistryDustTables *tables);
void chemistry_dust_opacities_from_tables(
    const ChemistryDustTables *tables, double temperature,
    double *kappa_rosseland, double *kappa_planck);
#endif
