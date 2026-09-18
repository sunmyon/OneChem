#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "co_photodissociation.h"

static int nearly_equal(double actual, double expected)
{
  return fabs(actual - expected) <= 1.e-12 * fmax(fabs(expected), 1.e-30);
}

int main(void)
{
  assert(nearly_equal(radiation_co_line_shielding_lee1996(0.0, 0.0), 1.0));
  assert(nearly_equal(radiation_co_line_shielding_lee1996(pow(10.0, 20.20),
                                                          pow(10.0, 22.01)),
                      1.619e-5 * 4.666e-3));
  assert(nearly_equal(radiation_co_photodissociation_rate(2.0, 0.0, 0.0, 0.0),
                      2.06e-10));
  assert(nearly_equal(radiation_co_photodissociation_rate(2.0, 0.0, 0.0, 1.0),
                      2.06e-10 * exp(-1.0)));
  assert(radiation_co_photodissociation_rate(-1.0, 0.0, 0.0, 0.0) == 0.0);
  puts("CO photodissociation check passed");
  return 0;
}
