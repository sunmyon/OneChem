#include "chemistry_standalone_adapter.h"

#include <stddef.h>

int chemistry_standalone_prepare_step_with_background(
    const ChemGlobalStep *step,
    const RadiationRateContext *radiation_context,
    const RadiationBlackbodyInput *background)
{
  ChemGlobalStep prepared_step;
  RadiationRates rates;

  if(step == NULL ||
     !radiation_rate_uniform_background(
         radiation_context, background, &rates))
    return 0;

  prepared_step = *step;
  prepared_step.background_J21 = background->J21;
  prepared_step.background_G0 = rates.g0;
  prepared_step.background_rate_Hminus_pd =
      rates.hminus_photodetachment_rate;
  prepared_step.background_rate_H2_diss =
      rates.h2_photodissociation_rate;
#ifdef HD
  prepared_step.background_rate_HD_diss =
      rates.hd_photodissociation_rate;
#endif
  chemistry_standalone_prepare_step(&prepared_step);
  return 1;
}
