#include <math.h>
#include <string.h>

#include "chemistry_standalone_adapter.h"

static ChemGlobalStep captured_step;
static int captured;

void chemistry_standalone_prepare_step(const ChemGlobalStep *step)
{
  captured_step = *step;
  captured = 1;
}

static int nearly_equal(double actual, double expected)
{
  return fabs(actual - expected) <= 1.e-12 * fabs(expected);
}

int main(void)
{
  RadiationRateContext context;
  RadiationBlackbodyInput background = {3.0, 1.e4};
  ChemGlobalStep step;

  memset(&step, 0, sizeof(step));
  step.cmb_temperature = 20.0;
  step.scale_factor = 0.125;
  step.simulation_time = 42.0;

  if(!radiation_rate_context_init(&context) ||
     !chemistry_standalone_prepare_step_with_background(
         &step, &context, &background) || !captured)
    return 1;

  if(captured_step.cmb_temperature != step.cmb_temperature ||
     captured_step.scale_factor != step.scale_factor ||
     captured_step.simulation_time != step.simulation_time ||
     captured_step.background_J21 != background.J21 ||
     !nearly_equal(captured_step.background_rate_Hminus_pd,
                   3.0 * 1.9699751793582524e-7) ||
     !nearly_equal(captured_step.background_rate_H2_diss,
                   3.0 * 1.39e-12) ||
#ifdef HD
     !nearly_equal(captured_step.background_rate_HD_diss,
                   3.0 * 1.55e-12) ||
#endif
     !nearly_equal(captured_step.background_G0,
                   3.0 * 2.803495007223459))
    return 1;

  radiation_rate_context_destroy(&context);
  return 0;
}
