#include <math.h>
#include <stdio.h>

#include "../../radiation/radiation_rates.h"

static int nearly_equal(double actual, double expected, double relative_error)
{
  const double scale = fabs(expected) > 1.e-300 ? fabs(expected) : 1.0;
  return fabs(actual - expected) <= relative_error * scale;
}

int main(void)
{
  RadiationRateContext context;
  RadiationBlackbodyInput background = {1.0, 1.e4};
  RadiationRates rates;
  RadiationRates doubled;
  RadiationRates zero;
  RadiationRates hot_rates;
  RadiationBlackbodySpectrumRates spectrum;
  RadiationRateContext independent_context;
  RadiationBackgroundState cached_background = {0};

  if(!radiation_rate_context_init(&context) ||
     !radiation_rate_uniform_background(&context, &background, &rates))
    return 1;

  if(!radiation_background_state_update(
         &context, &background, &cached_background) ||
     cached_background.revision != 1 ||
     !radiation_background_state_update(
         &context, &background, &cached_background) ||
     cached_background.revision != 1) {
    fprintf(stderr, "unchanged background was recomputed\n");
    return 1;
  }

  if(!radiation_rate_blackbody_spectrum(&context, 1.e4, &spectrum) ||
     !(spectrum.hminus_photodetachment_rate > 0.0) ||
     !(spectrum.g0 > 0.0) ||
     radiation_hminus_photodetachment_cross_section(0.5) != 0.0 ||
     !(radiation_hminus_photodetachment_cross_section(12.4) > 0.0))
    return 1;

  background.J21 = 2.0;
  if(!radiation_background_state_update(
         &context, &background, &cached_background) ||
     cached_background.revision != 2) {
    fprintf(stderr, "changed background was not recomputed\n");
    return 1;
  }
  if(!radiation_rate_uniform_background(&context, &background, &doubled))
    return 1;
  background.J21 = 0.0;
  if(!radiation_rate_uniform_background(&context, &background, &zero))
    return 1;

  if(!nearly_equal(doubled.hminus_photodetachment_rate,
                   2.0 * rates.hminus_photodetachment_rate, 1.e-14) ||
     !nearly_equal(doubled.h2_photodissociation_rate,
                   2.0 * rates.h2_photodissociation_rate, 1.e-14) ||
     !nearly_equal(doubled.hd_photodissociation_rate,
                   2.0 * rates.hd_photodissociation_rate, 1.e-14) ||
     !nearly_equal(doubled.g0, 2.0 * rates.g0, 1.e-14) ||
     zero.hminus_photodetachment_rate != 0.0 ||
     zero.h2_photodissociation_rate != 0.0 ||
     zero.hd_photodissociation_rate != 0.0 || zero.g0 != 0.0) {
    fprintf(stderr, "background radiation rates are not linear in J21\n");
    return 1;
  }

  /* Fixed regression values from the pre-migration 200-bin lookup. */
  if(!nearly_equal(rates.hminus_photodetachment_rate,
                   1.9699751793582524e-7, 2.e-12) ||
     !nearly_equal(rates.h2_photodissociation_rate, 1.39e-12, 1.e-14) ||
     !nearly_equal(rates.hd_photodissociation_rate, 1.55e-12, 1.e-14) ||
     !nearly_equal(rates.g0, 2.803495007223459, 2.e-12)) {
    fprintf(stderr, "J21=1 Teff=1e4: Hminus=%.17g H2=%.17g G0=%.17g\n",
            rates.hminus_photodetachment_rate,
            rates.h2_photodissociation_rate, rates.g0);
    return 1;
  }

  background.J21 = 1.0;
  background.effective_temperature = 1.e5;
  if(!radiation_rate_uniform_background(&context, &background, &hot_rates))
    return 1;
  if(!nearly_equal(hot_rates.hminus_photodetachment_rate,
                   1.255908313016088e-11, 2.e-12) ||
     !nearly_equal(hot_rates.h2_photodissociation_rate, 1.39e-12, 1.e-14) ||
     !nearly_equal(hot_rates.hd_photodissociation_rate, 1.55e-12, 1.e-14) ||
     !nearly_equal(hot_rates.g0, 0.010164195519360708, 2.e-12)) {
    fprintf(stderr, "J21=1 Teff=1e5: Hminus=%.17g H2=%.17g G0=%.17g\n",
            hot_rates.hminus_photodetachment_rate,
            hot_rates.h2_photodissociation_rate, hot_rates.g0);
    return 1;
  }

  if(!radiation_rate_context_init(&independent_context))
    return 1;
  radiation_rate_context_destroy(&context);
  if(radiation_rate_uniform_background(&context, &background, &rates) ||
     !radiation_rate_uniform_background(
         &independent_context, &background, &rates))
    return 1;
  radiation_rate_context_destroy(&independent_context);
  radiation_background_state_reset(&cached_background);
  return 0;
}
