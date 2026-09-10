#ifndef COOLING_FUNC_H
#define COOLING_FUNC_H

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

double calc_continuum_opacities(double temperature, double density, double length, struct opacities *opacities);
double continuum_cooling_rate_thin(double temperature, double density, const double *abundances);
double continuum_cooling_rate_thin_Lambda(double temperature);
double get_PlanckOpacity(double temperature, double density);
double fesc_Lya(double nH, double temperature, double y_H, double y_e, double radius, double velocity, double *cooling_rate_Lya, double *cooling_rate_2ph);

void tabulate_opacities_in_LW_band(void);
#ifdef CONTINUUM_COOLING_MATSUKOBA
double H_criticaldensity(double nH, double temperature, double column_H, double y_H, double y_e);
void read_opacity_tables(void);
void tabulate_cooling_rate_cont_thin(void);
#endif

double OxygenI_twolevel_cooling(double nOI, double ne, double nH, double temperature);
double CarbonI_twolevel_cooling(double nCI, double ne, double nH, double temperature);
double CarbonII_twolevel_cooling(double nCII, double ne, double nH, double temperature);

#endif
