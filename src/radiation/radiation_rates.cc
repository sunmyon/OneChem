#include "radiation_rates.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#define RADIATION_PLANCK 6.6262e-27
#define RADIATION_LIGHT_SPEED 2.99792458e10
/* Preserve the historical normalization, which used C_LIGHT while spectral
 * integrals used CLIGHT_REAL. */
#define RADIATION_REFERENCE_LIGHT_SPEED 2.9979e10
#define RADIATION_BOLTZMANN 1.3806e-16
#define RADIATION_EV_TO_ERG 1.60184e-12
#define RADIATION_EV_TO_HZ 2.41838e14
#define RADIATION_PI 3.14159265358979323846

double radiation_hminus_photodetachment_cross_section(
    double photon_energy_ev)
{
  static const double coefficient[6] = {
      152.519, 49.534, -118.858, 92.536, -34.194, 4.982};
  const double lambda_0 = 1.6419;
  double lambda;
  double f_lambda = 0.0;

  if(!(photon_energy_ev > 0.0) || !isfinite(photon_energy_ev))
    return 0.0;

  lambda = RADIATION_PLANCK * RADIATION_LIGHT_SPEED /
           (photon_energy_ev * RADIATION_EV_TO_ERG) * 1.e4;
  if(lambda > lambda_0)
    return 0.0;

  for(int i = 0; i < 6; ++i)
    f_lambda += coefficient[i] *
                pow(1.0 / lambda - 1.0 / lambda_0, 0.5 * i);

  return 1.e-18 * pow(lambda, 3) *
         pow(1.0 / lambda - 1.0 / lambda_0, 1.5) * f_lambda;
}

double Hminus_detachment_John1988(double photon_energy_ev)
{
  return radiation_hminus_photodetachment_cross_section(photon_energy_ev);
}

static double blackbody_intensity(double photon_energy_ev,
                                  double temperature)
{
  const double energy = photon_energy_ev * RADIATION_EV_TO_ERG;
  const double hc = RADIATION_PLANCK * RADIATION_LIGHT_SPEED;
  return 2.0 * pow(energy, 3) / (hc * hc) /
         (exp(energy / (RADIATION_BOLTZMANN * temperature)) - 1.0);
}

static double reference_blackbody_intensity(double photon_energy_ev,
                                            double temperature)
{
  const double energy = photon_energy_ev * RADIATION_EV_TO_ERG;
  const double hc = RADIATION_PLANCK * RADIATION_REFERENCE_LIGHT_SPEED;
  return 2.0 * pow(energy, 3) / hc / hc /
         (exp(energy / (RADIATION_BOLTZMANN * temperature)) - 1.0);
}

static double integrate_hminus_rate(double temperature)
{
  const double energy_min = 0.76;
  const double energy_max = 13.6;
  const int bin_count = 500;
  const double dlog_energy = log(energy_max / energy_min) / bin_count;
  double energy_old = energy_min;
  double integrand_old = 0.0;
  double sum = 0.0;

  for(int i = 0; i < bin_count; ++i) {
    const double energy = energy_min * exp(i * dlog_energy);
    const double dnu = (energy - energy_old) * RADIATION_EV_TO_HZ;
    const double intensity = blackbody_intensity(energy, temperature);
    const double sigma =
        radiation_hminus_photodetachment_cross_section(energy);
    const double integrand =
        4.0 * RADIATION_PI * intensity * sigma /
        (energy * RADIATION_EV_TO_ERG);

    if(i > 0)
      sum += 0.5 * (integrand_old + integrand) * dnu;

    energy_old = energy;
    integrand_old = integrand;
  }

  return sum;
}

static double integrate_g0(double temperature)
{
  const double energy_min = 5.12;
  const double energy_max = 13.6;
  const int bin_count = 500;
  const double dlog_energy =
      log(energy_max / energy_min) / (double)bin_count;
  double energy_old = energy_min;
  double intensity_old = 0.0;
  double sum = 0.0;

  for(int i = 0; i < bin_count; ++i) {
    const double energy = energy_min * exp(i * dlog_energy);
    const double dnu = (energy - energy_old) * RADIATION_EV_TO_HZ;
    const double intensity = blackbody_intensity(energy, temperature);

    if(i > 0)
      sum += 0.5 * (intensity_old + intensity) * dnu;

    energy_old = energy;
    intensity_old = intensity;
  }

  return 4.0 * RADIATION_PI * sum / 1.6e-3;
}

int radiation_rate_context_init(RadiationRateContext *context)
{
  const double temperature_min = 3000.0;
  const double temperature_max = 1.e6;

  if(context == NULL)
    return 0;

  memset(context, 0, sizeof(*context));
  context->effective_temperature_count =
      RADIATION_RATE_TEMPERATURE_BINS;
  context->effective_temperature_min = temperature_min;
  context->effective_temperature_max = temperature_max;
  context->log_temperature_step =
      log(temperature_max / temperature_min) /
      RADIATION_RATE_TEMPERATURE_BINS;

  for(int i = 0; i < context->effective_temperature_count; ++i) {
    const double temperature =
        temperature_min * exp(context->log_temperature_step * i);
    context->effective_temperature[i] = temperature;
    context->hminus_photodetachment_rate[i] =
        integrate_hminus_rate(temperature);
    context->g0[i] = integrate_g0(temperature);
  }

  context->initialized = 1;
  return 1;
}

void radiation_rate_context_destroy(RadiationRateContext *context)
{
  if(context != NULL)
    memset(context, 0, sizeof(*context));
}

int radiation_rate_blackbody_spectrum(
    const RadiationRateContext *context,
    double effective_temperature,
    RadiationBlackbodySpectrumRates *rates)
{
  int index;
  double input_temperature;
  double temperature_1;
  double temperature_2;
  double weight_1;
  double weight_2;
  int count;

  if(context == NULL || rates == NULL || !context->initialized ||
     context->effective_temperature_count < 2 ||
     !isfinite(effective_temperature) ||
     effective_temperature < context->effective_temperature_min ||
     effective_temperature > context->effective_temperature_max)
    return 0;

  count = context->effective_temperature_count;
  input_temperature = effective_temperature;
  if(input_temperature < context->effective_temperature[0])
    input_temperature = context->effective_temperature[0];
  if(input_temperature > context->effective_temperature[count - 1])
    input_temperature = context->effective_temperature[count - 1];

  index = (int)(log(input_temperature /
                    context->effective_temperature_min) /
                context->log_temperature_step);
  if(index < 0)
    index = 0;
  if(index > count - 2)
    index = count - 2;

  temperature_1 = context->effective_temperature[index];
  temperature_2 = context->effective_temperature[index + 1];
  weight_1 = (temperature_2 - input_temperature) /
             (temperature_2 - temperature_1);
  weight_2 = (input_temperature - temperature_1) /
             (temperature_2 - temperature_1);

  rates->hminus_photodetachment_rate =
      weight_1 * context->hminus_photodetachment_rate[index] +
      weight_2 * context->hminus_photodetachment_rate[index + 1];
  rates->g0 = weight_1 * context->g0[index] +
              weight_2 * context->g0[index + 1];
  return 1;
}

int radiation_rate_uniform_background(
    const RadiationRateContext *context,
    const RadiationBlackbodyInput *input,
    RadiationRates *rates)
{
  const double reference_energy_ev = 13.6;
  RadiationBlackbodySpectrumRates spectrum;
  double reference_intensity;
  double normalization;

  if(context == NULL || input == NULL || rates == NULL ||
     !isfinite(input->J21) || input->J21 < 0.0 ||
     !radiation_rate_blackbody_spectrum(
         context, input->effective_temperature, &spectrum))
    return 0;

  reference_intensity =
      reference_blackbody_intensity(reference_energy_ev,
                                    input->effective_temperature);
  if(!(reference_intensity > 0.0) || !isfinite(reference_intensity))
    return 0;

  normalization = input->J21 * 1.e-21 / reference_intensity;
  rates->hminus_photodetachment_rate =
      normalization * spectrum.hminus_photodetachment_rate;
  rates->h2_photodissociation_rate =
      1.39e9 * (input->J21 * 1.e-21);
  rates->hd_photodissociation_rate =
      (1.55 / 1.39) * rates->h2_photodissociation_rate;
  rates->g0 = normalization * spectrum.g0;
  return 1;
}

void radiation_background_state_reset(RadiationBackgroundState *state)
{
  if(state != NULL)
    memset(state, 0, sizeof(*state));
}

int radiation_background_state_update(
    const RadiationRateContext *context,
    const RadiationBlackbodyInput *input,
    RadiationBackgroundState *state)
{
  RadiationRates rates;

  if(context == NULL || input == NULL || state == NULL)
    return 0;

  if(state->initialized && state->input.J21 == input->J21 &&
     state->input.effective_temperature == input->effective_temperature)
    return 1;

  if(!radiation_rate_uniform_background(context, input, &rates))
    return 0;

  state->input = *input;
  state->rates = rates;
  state->initialized = 1;
  state->revision++;
  return 1;
}
