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

#if defined(LYMANWERNER) && defined(BACKGROUND_LW_RADIATION)
static RadiationRateContext arepo_radiation_context;
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

int chemistry_arepo_host_init(void)
{
  const ChemArepoConfig config = {
      All.ChemistryTableDir,
      All.Epsilon};

  if(chemistry_arepo_context_init(&config, &arepo_chemistry_context) !=
     CHEMISTRY_STATUS_OK)
    return 0;

#if defined(LYMANWERNER) && defined(BACKGROUND_LW_RADIATION)
  RadiationBlackbodyInput background = {
      All.J21,
      background_effective_temperature()};
  if(!radiation_rate_context_init(&arepo_radiation_context) ||
     !radiation_background_state_update(
         &arepo_radiation_context, &background,
         &arepo_background_state)) {
    radiation_rate_context_destroy(&arepo_radiation_context);
    radiation_background_state_reset(&arepo_background_state);
    chemistry_arepo_context_destroy(&arepo_chemistry_context);
    return 0;
  }
#endif
  return 1;
}

void chemistry_arepo_host_destroy(void)
{
#if defined(LYMANWERNER) && defined(BACKGROUND_LW_RADIATION)
  radiation_background_state_reset(&arepo_background_state);
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
