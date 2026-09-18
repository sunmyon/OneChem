#include "photoionization_bins.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#define PHOTO_QUADRATURE_COUNT 96
#define PHOTO_TEMPERATURE_MIN 1.0e3
#define PHOTO_TEMPERATURE_MAX 1.0e6
#define PHOTO_REFERENCE_ENERGY 13.6
#define PHOTO_EV_TO_ERG 1.602176634e-12
#define PHOTO_EV_TO_HZ 2.417989242084918e14
#define PHOTO_BOLTZMANN 1.380649e-16
#define PHOTO_PLANCK 6.62607015e-27
#define PHOTO_LIGHT_SPEED 2.99792458e10
#define PHOTO_PI 3.14159265358979323846
#define PHOTO_FOUR_PI 12.56637061435917295385

typedef struct {
  double threshold, maximum, energy_scale, sigma_scale;
  double ya, power, yw, y0, y1;
} PhotoCrossSectionFit;

/* Verner et al. (1996), sigma_scale in Mb. */
static const PhotoCrossSectionFit photo_fit[PHOTOIONIZATION_TRANSITION_COUNT] = {
  {13.60, 50000.0, 0.4298, 5.475e4, 32.88, 2.963, 0.0, 0.0, 0.0},
  {24.59, 50000.0, 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434, 2.136},
  {54.42, 50000.0, 1.720, 1.369e4, 32.88, 2.963, 0.0, 0.0, 0.0},
  {11.26, 291.0, 2.144, 502.7, 62.16, 5.101, 0.09157, 1.133, 1.607},
  {24.38, 307.6, 0.4058, 8.709, 126.1, 8.578, 2.093, 49.29, 3.234},
  {47.89, 328.9, 4.614, 1.539e4, 1.737, 15.93, 5.922, 0.004378, 0.02528},
  {13.62, 538.0, 1.240, 1745.0, 3.784, 17.64, 0.07589, 8.698, 0.1271},
  {35.12, 558.1, 1.386, 59.67, 31.75, 8.943, 0.01934, 21.31, 0.01503},
  {54.94, 584.0, 0.1723, 675.3, 385.2, 6.822, 0.1191, 0.003839, 0.4569}
};

/* Threshold-aligned groups shared with the C/O extension planned next. */
static const double photo_group_edge[PHOTOIONIZATION_GROUP_COUNT + 1] = {
  11.26, 13.60, 13.62, 14.53, 24.38, 24.59, 29.60, 35.12, 47.45,
  47.89, 54.42, 54.94, 70.0, 100.0, 150.0, 250.0, 400.0, 740.0
};

static double log_expm1_positive(double value)
{
  if(value > 50.0)
    return value + log1p(-exp(-value));
  return log(expm1(value));
}

static double normalized_blackbody(double energy_eV, double temperature)
{
  const double x = energy_eV * PHOTO_EV_TO_ERG /
                   (PHOTO_BOLTZMANN * temperature);
  const double x_reference = PHOTO_REFERENCE_ENERGY * PHOTO_EV_TO_ERG /
                             (PHOTO_BOLTZMANN * temperature);
  const double log_value = 3.0 * log(energy_eV / PHOTO_REFERENCE_ENERGY) +
      log_expm1_positive(x_reference) - log_expm1_positive(x);

  if(log_value < log(DBL_MIN))
    return 0.0;
  if(log_value > log(DBL_MAX))
    return DBL_MAX;
  return exp(log_value);
}

double photoionization_threshold_energy(enum PhotoionizationTransition transition)
{
  if(transition < 0 || transition >= PHOTOIONIZATION_TRANSITION_COUNT)
    return 0.0;
  return photo_fit[transition].threshold;
}

double photoionization_cross_section(enum PhotoionizationTransition transition,
                                     double energy_eV)
{
  const PhotoCrossSectionFit *fit;
  double x, y, value;

  if(transition < 0 || transition >= PHOTOIONIZATION_TRANSITION_COUNT)
    return 0.0;
  fit = &photo_fit[transition];
  if(energy_eV < fit->threshold || energy_eV > fit->maximum)
    return 0.0;

  x = energy_eV / fit->energy_scale - fit->y0;
  y = sqrt(x * x + fit->y1 * fit->y1);
  if(y <= 0.0)
    return 0.0;
  value = fit->sigma_scale * ((x - 1.0) * (x - 1.0) + fit->yw * fit->yw);
  value *= pow(y, 0.5 * fit->power - 5.5) *
           pow(1.0 + sqrt(y / fit->ya), -fit->power);
  return value * 1.0e-18;
}

double photoionization_blackbody_integral(double temperature)
{
  double hc, reference_energy, x_reference, log_norm, log_integral;
  if(temperature <= 0.0)
    return 0.0;
  hc = PHOTO_LIGHT_SPEED * PHOTO_PLANCK;
  reference_energy = PHOTO_REFERENCE_ENERGY * PHOTO_EV_TO_ERG;
  x_reference = reference_energy / (PHOTO_BOLTZMANN * temperature);
  log_norm = log(2.0) + 3.0 * log(reference_energy) - 2.0 * log(hc) -
             log_expm1_positive(x_reference);
  log_integral = log(2.0 * PHOTO_LIGHT_SPEED / 15.0) +
                 4.0 * log(PHOTO_PI * PHOTO_BOLTZMANN * temperature) -
                 3.0 * log(hc) - log_norm;
  return log_integral > log(DBL_MAX) ? DBL_MAX : exp(log_integral);
}

static void integrate_group(double temperature, int group,
                            double rate[PHOTOIONIZATION_TRANSITION_COUNT],
                            double heating[PHOTOIONIZATION_TRANSITION_COUNT],
                            double opacity[PHOTOIONIZATION_ABSORBER_COUNT])
{
  double photon_integral = 0.0;
  double opacity_sum[PHOTOIONIZATION_ABSORBER_COUNT] = {0.0, 0.0, 0.0};
  const double energy_min = photo_group_edge[group];
  const double energy_step =
      (photo_group_edge[group + 1] - energy_min) / PHOTO_QUADRATURE_COUNT;
  int transition, absorber, sample;

  for(transition = 0; transition < PHOTOIONIZATION_TRANSITION_COUNT; transition++) {
    rate[transition] = 0.0;
    heating[transition] = 0.0;
  }
  for(sample = 0; sample < PHOTO_QUADRATURE_COUNT; sample++) {
    const double energy = energy_min + (sample + 0.5) * energy_step;
    const double spectrum = normalized_blackbody(energy, temperature);
    const double photon_weight = spectrum / (energy * PHOTO_EV_TO_ERG) *
                                 energy_step * PHOTO_EV_TO_HZ;
    photon_integral += photon_weight;
    for(absorber = 0; absorber < PHOTOIONIZATION_ABSORBER_COUNT; absorber++)
      opacity_sum[absorber] += photon_weight *
          photoionization_cross_section(absorber, energy);
    for(transition = 0; transition < PHOTOIONIZATION_TRANSITION_COUNT; transition++) {
      const double sigma = photoionization_cross_section(transition, energy);
      const double excess = (energy - photo_fit[transition].threshold) *
                            PHOTO_EV_TO_ERG;
      rate[transition] += photon_weight * sigma / PHOTO_FOUR_PI;
      if(excess > 0.0)
        heating[transition] += photon_weight * sigma * excess / PHOTO_FOUR_PI;
    }
  }
  for(absorber = 0; absorber < PHOTOIONIZATION_ABSORBER_COUNT; absorber++)
    opacity[absorber] = photon_integral > 0.0 ?
        opacity_sum[absorber] / photon_integral : 0.0;
}

int photoionization_bins_init(PhotoionizationBinTables *tables)
{
  const double step = log(PHOTO_TEMPERATURE_MAX / PHOTO_TEMPERATURE_MIN) /
                      (PHOTOIONIZATION_TEMPERATURE_COUNT - 1);
  int index, group;
  if(tables == NULL)
    return 0;
  memset(tables, 0, sizeof(*tables));
  for(index = 0; index < PHOTOIONIZATION_TEMPERATURE_COUNT; index++) {
    tables->temperature[index] = PHOTO_TEMPERATURE_MIN * exp(index * step);
    for(group = 0; group < PHOTOIONIZATION_GROUP_COUNT; group++)
      integrate_group(tables->temperature[index], group,
                      tables->rate[index][group],
                      tables->heating[index][group],
                      tables->opacity[index][group]);
  }
  tables->initialized = 1;
  return 1;
}

static void temperature_coordinate(double temperature, int *index,
                                   double *fraction)
{
  const double step = log(PHOTO_TEMPERATURE_MAX / PHOTO_TEMPERATURE_MIN) /
                      (PHOTOIONIZATION_TEMPERATURE_COUNT - 1);
  double position;
  if(temperature <= PHOTO_TEMPERATURE_MIN) {
    *index = 0; *fraction = 0.0; return;
  }
  if(temperature >= PHOTO_TEMPERATURE_MAX) {
    *index = PHOTOIONIZATION_TEMPERATURE_COUNT - 2;
    *fraction = 1.0; return;
  }
  position = log(temperature / PHOTO_TEMPERATURE_MIN) / step;
  *index = (int)position;
  *fraction = position - *index;
}

static double interpolate(double first, double second, double fraction)
{
  if(first <= 0.0 || second <= 0.0)
    return (1.0 - fraction) * first + fraction * second;
  return exp((1.0 - fraction) * log(first) + fraction * log(second));
}

static double cell_average_attenuation(double tau_near, double delta_tau)
{
  double average;
  if(tau_near >= 745.0)
    return 0.0;
  if(delta_tau < 0.0)
    delta_tau = 0.0;
  average = delta_tau < 1.0e-7 ?
      1.0 - 0.5 * delta_tau + delta_tau * delta_tau / 6.0 :
      -expm1(-delta_tau) / delta_tau;
  return exp(-fmax(0.0, tau_near)) * average;
}

void photoionization_bins_evaluate(
    const PhotoionizationBinTables *tables, double temperature,
    double radius_squared_cm2, int absorber_count,
    const double tau_near[PHOTOIONIZATION_ABSORBER_COUNT],
    const double delta_tau[PHOTOIONIZATION_ABSORBER_COUNT],
    double rate[PHOTOIONIZATION_TRANSITION_COUNT],
    double heating[PHOTOIONIZATION_TRANSITION_COUNT])
{
  int index, group, transition, opacity_absorber;
  double fraction;
  for(transition = 0; transition < PHOTOIONIZATION_TRANSITION_COUNT; transition++) {
    rate[transition] = 0.0;
    heating[transition] = 0.0;
  }
  if(tables == NULL || !tables->initialized || radius_squared_cm2 <= 0.0)
    return;
  if(absorber_count < 1)
    absorber_count = 1;
  if(absorber_count > PHOTOIONIZATION_ABSORBER_COUNT)
    absorber_count = PHOTOIONIZATION_ABSORBER_COUNT;
  temperature_coordinate(temperature, &index, &fraction);

  for(group = 0; group < PHOTOIONIZATION_GROUP_COUNT; group++) {
    double group_tau = 0.0, group_delta_tau = 0.0;
    for(opacity_absorber = 0; opacity_absorber < absorber_count;
        opacity_absorber++) {
      const double sigma = (1.0 - fraction) *
          tables->opacity[index][group][opacity_absorber] + fraction *
          tables->opacity[index + 1][group][opacity_absorber];
      const double reference = photoionization_cross_section(
          opacity_absorber,
          photoionization_threshold_energy(opacity_absorber));
      if(reference > 0.0 && sigma > 0.0) {
        group_tau += fmax(0.0, tau_near[opacity_absorber]) * sigma / reference;
        group_delta_tau += fmax(0.0, delta_tau[opacity_absorber]) *
                           sigma / reference;
      }
    }
    const double attenuation =
        cell_average_attenuation(group_tau, group_delta_tau);
    for(transition = 0; transition < PHOTOIONIZATION_TRANSITION_COUNT; transition++) {
      rate[transition] += attenuation * interpolate(
          tables->rate[index][group][transition],
          tables->rate[index + 1][group][transition], fraction) /
          radius_squared_cm2;
      heating[transition] += attenuation * interpolate(
          tables->heating[index][group][transition],
          tables->heating[index + 1][group][transition], fraction) /
          radius_squared_cm2;
    }
  }
}
