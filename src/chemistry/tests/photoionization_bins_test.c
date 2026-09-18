#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "../../radiation/photoionization_bins.h"

int main(void)
{
  static PhotoionizationBinTables tables;
  double clear_tau[PHOTOIONIZATION_ABSORBER_COUNT] = {0.0, 0.0, 0.0};
  double thick_tau[PHOTOIONIZATION_ABSORBER_COUNT] = {10.0, 10.0, 10.0};
  double rate_clear[PHOTOIONIZATION_TRANSITION_COUNT];
  double heat_clear[PHOTOIONIZATION_TRANSITION_COUNT];
  double rate_thick[PHOTOIONIZATION_TRANSITION_COUNT];
  double heat_thick[PHOTOIONIZATION_TRANSITION_COUNT];

  assert(photoionization_bins_init(&tables));
  assert(photoionization_cross_section(PHOTOIONIZATION_HI, 13.59) == 0.0);
  assert(photoionization_cross_section(PHOTOIONIZATION_HEI, 24.58) == 0.0);
  assert(photoionization_cross_section(PHOTOIONIZATION_HEII, 54.41) == 0.0);
  assert(photoionization_cross_section(PHOTOIONIZATION_CI, 11.25) == 0.0);
  assert(photoionization_cross_section(PHOTOIONIZATION_CI, 20.0) > 0.0);
  assert(photoionization_cross_section(PHOTOIONIZATION_OIII, 70.0) > 0.0);
  assert(photoionization_blackbody_integral(1.0e5) > 0.0);

  photoionization_bins_evaluate(&tables, 1.0e5, 1.0e36, 3,
                                clear_tau, clear_tau,
                                rate_clear, heat_clear);
  photoionization_bins_evaluate(&tables, 1.0e5, 1.0e36, 3,
                                thick_tau, thick_tau,
                                rate_thick, heat_thick);
  for(int absorber = 0; absorber < PHOTOIONIZATION_ABSORBER_COUNT; absorber++) {
    assert(isfinite(rate_clear[absorber]) && rate_clear[absorber] > 0.0);
    assert(isfinite(heat_clear[absorber]) && heat_clear[absorber] > 0.0);
    assert(rate_thick[absorber] >= 0.0);
    assert(heat_thick[absorber] >= 0.0);
    assert(rate_thick[absorber] < rate_clear[absorber]);
    assert(heat_thick[absorber] < heat_clear[absorber]);
  }
  for(int transition = PHOTOIONIZATION_CI;
      transition < PHOTOIONIZATION_TRANSITION_COUNT; transition++) {
    assert(isfinite(rate_clear[transition]) && rate_clear[transition] > 0.0);
    assert(isfinite(heat_clear[transition]) && heat_clear[transition] > 0.0);
    assert(rate_thick[transition] >= 0.0);
    assert(heat_thick[transition] >= 0.0);
    assert(rate_thick[transition] < rate_clear[transition]);
    assert(heat_thick[transition] < heat_clear[transition]);
  }
  puts("photoionization bin kernel check passed");
  return 0;
}
