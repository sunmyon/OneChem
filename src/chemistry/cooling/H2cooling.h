#ifndef H2COOLING_H
#define H2COOLING_H


double H2_cooling(double nH, double temperature, double yHI, double yH2, double yHeI, double yHII, double yElectron);

typedef double (*H2CoolingFunction)(double nH, double temperature, double yHI, double yH2, double yHeI, double yHII, double yElectron);
typedef double (*HDCoolingFunction)(double nH, double temperature, double yHI, double yHD);
typedef double (*LineEscapeFunction)(double temperature, double columnDensity);

double H2_cooling_G15(double nH, double temperature, double yHI, double yH2, double yHeI, double yHII, double yElectron);
double H2_cooling_GP12(double nH, double temperature, double yHI, double yH2, double yHeI, double yHII, double yElectron);

#ifdef HD
double HD_cooling_GP12(double nH, double temperature,
                       double radiation_temperature, double yHI, double yHD,
                       int diagnostics_enabled);
double fesc_HD_fit(double temperature, double columnDensity);
#endif

double fesc_fit(double temperature, double columnDensity);
double fesc_unity(double temperature, double columnDensity);


#endif
