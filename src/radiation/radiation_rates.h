#pragma once

#define RADIATION_RATE_TEMPERATURE_BINS 200

/* Rank/task-owned, run-lifetime lookup data.  Initialize explicitly before
 * use; the context is immutable after successful initialization. */
typedef struct RadiationRateContext {
  int initialized;
  int effective_temperature_count;
  double effective_temperature_min;
  double effective_temperature_max;
  double log_temperature_step;
  double effective_temperature[RADIATION_RATE_TEMPERATURE_BINS];
  double hminus_photodetachment_rate[RADIATION_RATE_TEMPERATURE_BINS];
  double g0[RADIATION_RATE_TEMPERATURE_BINS];
} RadiationRateContext;

typedef struct RadiationBlackbodyInput {
  double J21;
  double effective_temperature;
} RadiationBlackbodyInput;

/* Undiluted blackbody quantities used by point-source radiation transport.
 * A host applies its own geometric dilution or transfer solution. */
typedef struct RadiationBlackbodySpectrumRates {
  double hminus_photodetachment_rate;
  double g0;
} RadiationBlackbodySpectrumRates;

/* Final spatially uniform rates ready to be passed to chemistry. */
typedef struct RadiationRates {
  double hminus_photodetachment_rate; /* s^-1 */
  double h2_photodissociation_rate;   /* s^-1 */
  double hd_photodissociation_rate;   /* s^-1 */
  double g0;
} RadiationRates;

/* Host-owned cache for a spatially uniform background. The spectral
 * conversion is repeated only when J21 or effective temperature changes. */
typedef struct RadiationBackgroundState {
  int initialized;
  RadiationBlackbodyInput input;
  RadiationRates rates;
  unsigned long revision;
} RadiationBackgroundState;

#ifdef __cplusplus
extern "C" {
#endif

int radiation_rate_context_init(RadiationRateContext *context);
void radiation_rate_context_destroy(RadiationRateContext *context);

int radiation_rate_blackbody_spectrum(
    const RadiationRateContext *context,
    double effective_temperature,
    RadiationBlackbodySpectrumRates *rates);

int radiation_rate_uniform_background(
    const RadiationRateContext *context,
    const RadiationBlackbodyInput *input,
    RadiationRates *rates);

void radiation_background_state_reset(RadiationBackgroundState *state);
int radiation_background_state_update(
    const RadiationRateContext *context,
    const RadiationBlackbodyInput *input,
    RadiationBackgroundState *state);

/* Public because LW opacity/transport needs the same physical cross section;
 * chemistry reaction and cooling modules do not call this function. */
double radiation_hminus_photodetachment_cross_section(
    double photon_energy_ev);

/* Temporary source-compatible entry point for the selected AREPO RayTrace. */
double Hminus_detachment_John1988(double photon_energy_ev);

#ifdef __cplusplus
}
#endif
