#ifndef COOLING_FUNC_H
#define COOLING_FUNC_H

typedef struct ChemistryContinuumOpacityTables ChemistryContinuumOpacityTables;

struct opacities {
  double tau;
  double tau_Rossland;
  double tau_Planck;
  double kappa_Rossland;
  double kappa_Planck;
#ifdef DUST_COOLING
  double kappa_dust_Rossland;
  double kappa_dust_Planck;
#endif
  double x;
  double f_reduce;
};


double calc_continuum_opacities(const ChemistryContinuumOpacityTables *tables,
                                double temperature, double density,
                                double length, struct opacities *opacities);
double continuum_cooling_rate_thin(double temperature, double density,
                                   const double *abundances,
                                   int diagnostics_enabled);
double continuum_cooling_rate_thin_Lambda(double temperature);
double get_PlanckOpacity(double temperature, double density);
double fesc_Lya(double nH, double temperature, double y_H, double y_e,
                double radius, double velocity, double *cooling_rate_Lya,
                double *cooling_rate_2ph, int diagnostics_enabled);

double H_criticaldensity(double nH, double temperature, double column_H, double y_H, double y_e);

#if defined(CONTINUUM_COOLING_MATSUKOBA) || defined(METAL_LINE_COOLING_SD93)


#endif

#ifdef METAL_LINE_COOLING_SD93
double get_cooling_rate_SD93(const ChemistryContinuumOpacityTables *tables, double temperature);
#endif

double OxygenI_twolevel_cooling(double nOI, double ne, double nH,
                                double temperature,
                                double radiation_temperature,
                                int diagnostics_enabled);
double CarbonI_twolevel_cooling(double nCI, double ne, double nH,
                                double temperature,
                                double radiation_temperature,
                                int diagnostics_enabled);
double CarbonII_twolevel_cooling(double nCII, double ne, double nH,
                                 double temperature,
                                 double radiation_temperature,
                                 int diagnostics_enabled);



#endif
