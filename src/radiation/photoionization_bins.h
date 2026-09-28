#ifndef PHOTOIONIZATION_BINS_H
#define PHOTOIONIZATION_BINS_H

#define PHOTOIONIZATION_GROUP_COUNT 17
#define PHOTOIONIZATION_ABSORBER_COUNT 3
#define PHOTOIONIZATION_TRANSITION_COUNT 9
#define PHOTOIONIZATION_TEMPERATURE_COUNT 300

enum PhotoionizationTransition {
  PHOTOIONIZATION_HI,
  PHOTOIONIZATION_HEI,
  PHOTOIONIZATION_HEII,
  PHOTOIONIZATION_CI,
  PHOTOIONIZATION_CII,
  PHOTOIONIZATION_CIII,
  PHOTOIONIZATION_OI,
  PHOTOIONIZATION_OII,
  PHOTOIONIZATION_OIII
};

typedef struct {
  double temperature[PHOTOIONIZATION_TEMPERATURE_COUNT];
  double rate[PHOTOIONIZATION_TEMPERATURE_COUNT]
             [PHOTOIONIZATION_GROUP_COUNT]
             [PHOTOIONIZATION_TRANSITION_COUNT];
  double heating[PHOTOIONIZATION_TEMPERATURE_COUNT]
                [PHOTOIONIZATION_GROUP_COUNT]
                [PHOTOIONIZATION_TRANSITION_COUNT];
  double opacity_ratio[PHOTOIONIZATION_TEMPERATURE_COUNT]
                [PHOTOIONIZATION_GROUP_COUNT]
                [PHOTOIONIZATION_ABSORBER_COUNT];
  int initialized;
} PhotoionizationBinTables;

/* Temperature-interpolated coefficients, shared by all rays to one source. */
typedef struct {
  double rate[PHOTOIONIZATION_GROUP_COUNT][PHOTOIONIZATION_TRANSITION_COUNT];
  double heating[PHOTOIONIZATION_GROUP_COUNT][PHOTOIONIZATION_TRANSITION_COUNT];
  double opacity_ratio[PHOTOIONIZATION_GROUP_COUNT][PHOTOIONIZATION_ABSORBER_COUNT];
  int initialized;
} PhotoionizationBinCoefficients;

void photoionization_bins_prepare(const PhotoionizationBinTables *tables,
                                 double effective_temperature,
                                 PhotoionizationBinCoefficients *coefficients);
void photoionization_bins_evaluate_prepared(
    const PhotoionizationBinCoefficients *coefficients,
    double radius_squared_cm2, int absorber_count,
    const double tau_near[PHOTOIONIZATION_ABSORBER_COUNT],
    const double delta_tau[PHOTOIONIZATION_ABSORBER_COUNT],
    double rate[PHOTOIONIZATION_TRANSITION_COUNT],
    double heating[PHOTOIONIZATION_TRANSITION_COUNT]);

double photoionization_threshold_energy(enum PhotoionizationTransition transition);
double photoionization_cross_section(enum PhotoionizationTransition transition,
                                     double energy_eV);
double photoionization_blackbody_integral(double effective_temperature);
int photoionization_bins_init(PhotoionizationBinTables *tables);
void photoionization_bins_evaluate(
    const PhotoionizationBinTables *tables, double effective_temperature,
    double radius_squared_cm2, int absorber_count,
    const double tau_near[PHOTOIONIZATION_ABSORBER_COUNT],
    const double delta_tau[PHOTOIONIZATION_ABSORBER_COUNT],
    double rate[PHOTOIONIZATION_TRANSITION_COUNT],
    double heating[PHOTOIONIZATION_TRANSITION_COUNT]);

#endif
