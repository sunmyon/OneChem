#include <string.h>

#include "allvars.h"

#include "chemistry_arepo_host.h"
#include "radiation/radiation_rates.h"

/*
 * Existing AREPO calls cooling_only() without an argument, so the thin host
 * layer must retain one run-lifetime handle. Tables and numerical state remain
 * inside the explicitly passed context; the chemistry core has no accessor for
 * this object.
 */
static ChemArepoContext arepo_chemistry_context;

#ifdef LYMANWERNER
static RadiationRateContext arepo_radiation_context;
#ifdef BACKGROUND_LW_RADIATION
static RadiationBackgroundState arepo_background_state;

static double background_effective_temperature(void)
{
#ifdef TYPE_B
  return 1.e4;
#elif defined(TYPE_3E4)
  return 3.e4;
#else
  return 1.e5;
#endif
}
#endif
#endif

int chemistry_arepo_host_init(void)
{
  const ChemArepoConfig config = {
      All.ChemistryTableDir,
      All.Epsilon};

  if(chemistry_arepo_context_init(&config, &arepo_chemistry_context) !=
     CHEMISTRY_STATUS_OK)
    return 0;

#ifdef LYMANWERNER
  if(!radiation_rate_context_init(&arepo_radiation_context)) {
    chemistry_arepo_context_destroy(&arepo_chemistry_context);
    return 0;
  }
#ifdef BACKGROUND_LW_RADIATION
  RadiationBlackbodyInput background = {
      All.J21,
      background_effective_temperature()};
  if(!radiation_background_state_update(
         &arepo_radiation_context, &background,
         &arepo_background_state)) {
    radiation_rate_context_destroy(&arepo_radiation_context);
    radiation_background_state_reset(&arepo_background_state);
    chemistry_arepo_context_destroy(&arepo_chemistry_context);
    return 0;
  }
#endif
#endif
  return 1;
}

void chemistry_arepo_host_destroy(void)
{
#ifdef LYMANWERNER
#ifdef BACKGROUND_LW_RADIATION
  radiation_background_state_reset(&arepo_background_state);
#endif
  radiation_rate_context_destroy(&arepo_radiation_context);
#endif
  chemistry_arepo_context_destroy(&arepo_chemistry_context);
}

const ChemArepoContext *chemistry_arepo_host_context(void)
{
  return arepo_chemistry_context.chemistry != NULL ?
      &arepo_chemistry_context : NULL;
}

int chemistry_arepo_host_prepare_step(ChemistryStepContext *step)
{
  if(step == NULL)
    return 0;

  ChemistryGlobalEnvironment global = {0};
  global.cmb_temperature = All.T_CMB;
  global.scale_factor = All.cf_atime;
  global.simulation_time = All.Time;
  global.comoving_integration = All.ComovingIntegrationOn;

#if defined(LYMANWERNER) && defined(BACKGROUND_LW_RADIATION)
  const RadiationBlackbodyInput background = {
      All.J21,
      background_effective_temperature()};
  if(!radiation_background_state_update(
         &arepo_radiation_context, &background,
         &arepo_background_state))
    return 0;

  global.background_J21 = arepo_background_state.input.J21;
  global.background_rate_Hminus_pd =
      arepo_background_state.rates.hminus_photodetachment_rate;
  global.background_rate_H2_diss =
      arepo_background_state.rates.h2_photodissociation_rate;
#ifdef HD
  global.background_rate_HD_diss =
      arepo_background_state.rates.hd_photodissociation_rate;
#endif
#ifdef METALLICITY
  global.background_G0 = arepo_background_state.rates.g0;
#endif
#endif
  return chemistry_step_prepare(
      arepo_chemistry_context.chemistry, &global, step) == CHEMISTRY_STATUS_OK;
}

int chemistry_arepo_blackbody_rates(double effective_temperature,
                                    double *photodetachment_rate, double *g0)
{
  RadiationBlackbodySpectrumRates rates = {0};
  if(photodetachment_rate == NULL || g0 == NULL)
    return 0;

#ifdef LYMANWERNER
  if(effective_temperature < arepo_radiation_context.effective_temperature_min)
    effective_temperature = arepo_radiation_context.effective_temperature_min;
  if(effective_temperature > arepo_radiation_context.effective_temperature_max)
    effective_temperature = arepo_radiation_context.effective_temperature_max;
  if(!radiation_rate_blackbody_spectrum(&arepo_radiation_context,
                                        effective_temperature, &rates))
    return 0;
#else
  return 0;
#endif
  *photodetachment_rate = rates.hminus_photodetachment_rate;
  *g0 = rates.g0;
  return 1;
}

int chemistry_arepo_dust_opacities(double temperature,
                                   double *kappa_rosseland,
                                   double *kappa_planck)
{
  if(kappa_rosseland == NULL || kappa_planck == NULL)
    return 0;

  *kappa_rosseland = 0.0;
  *kappa_planck = 0.0;
#ifdef DUST_COOLING
  return chemistry_dust_opacities(arepo_chemistry_context.chemistry,
                                  temperature, kappa_rosseland,
                                  kappa_planck) == CHEMISTRY_STATUS_OK;
#else
  (void)temperature;
  return 0;
#endif
}
