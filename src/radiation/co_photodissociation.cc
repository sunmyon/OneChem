#include "co_photodissociation.h"

#include <math.h>
#include <stddef.h>

/* log10(column / cm^-2) tables from the Lee et al. (1996) data already
 * distributed with AREPO's SGChem OLD_CO_SELF_SHIELDING implementation. */
static const double co_columns[] = {
  12.00, 12.22, 12.48, 12.78, 13.12, 13.50, 13.93, 14.39, 14.90,
  15.22, 15.41, 15.65, 15.78, 15.95, 16.13, 16.22, 16.36, 16.49,
  16.63, 16.77, 16.90, 17.04, 17.18, 17.31, 17.45, 17.63, 17.81,
  17.98, 18.16, 18.34, 18.52, 18.62, 18.75, 18.89, 19.03, 19.11,
  19.23, 19.34, 19.40, 19.51, 19.61, 19.72, 19.77, 19.87, 19.92,
  19.97, 20.01, 20.06, 20.11, 20.15, 20.20
};

static const double co_factors[] = {
  9.990e-1, 9.981e-1, 9.961e-1, 9.912e-1, 9.815e-1, 9.601e-1,
  9.113e-1, 8.094e-1, 6.284e-1, 4.808e-1, 3.889e-1, 2.827e-1,
  2.293e-1, 1.695e-1, 1.224e-1, 1.017e-1, 7.764e-2, 5.931e-2,
  4.546e-2, 3.506e-2, 2.728e-2, 2.143e-2, 1.700e-2, 1.360e-2,
  1.094e-2, 8.273e-3, 6.283e-3, 4.773e-3, 3.611e-3, 2.704e-3,
  1.986e-3, 1.657e-3, 1.258e-3, 9.332e-4, 6.745e-4, 5.596e-4,
  4.123e-4, 2.982e-4, 2.490e-4, 1.827e-4, 1.324e-4, 9.473e-5,
  7.891e-5, 5.668e-5, 4.732e-5, 3.967e-5, 3.327e-5, 2.788e-5,
  2.331e-5, 1.944e-5, 1.619e-5
};

static const double h2_columns[] = {
  13.00, 13.43, 14.58, 15.82, 16.95, 17.97, 18.00, 18.31, 18.48,
  18.61, 18.70, 18.78, 18.85, 18.91, 18.97, 19.01, 19.31, 19.48,
  19.61, 19.71, 19.78, 19.85, 19.91, 19.97, 20.01, 20.31, 20.48,
  20.61, 20.71, 20.78, 20.85, 20.91, 20.97, 21.01, 21.31, 21.48,
  21.61, 21.71, 21.78, 21.85, 21.91, 21.97, 22.01
};

static const double h2_factors[] = {
  1.000e0, 9.999e-1, 9.893e-1, 9.678e-1, 9.465e-1, 9.137e-1,
  9.121e-1, 8.966e-1, 8.862e-1, 8.781e-1, 8.716e-1, 8.660e-1,
  8.612e-1, 8.569e-1, 8.524e-1, 8.497e-1, 8.262e-1, 8.118e-1,
  8.011e-1, 7.921e-1, 7.841e-1, 7.769e-1, 7.702e-1, 7.626e-1,
  7.579e-1, 7.094e-1, 6.712e-1, 6.378e-1, 6.074e-1, 5.791e-1,
  5.524e-1, 5.271e-1, 4.977e-1, 4.793e-1, 2.837e-1, 1.526e-1,
  7.774e-2, 3.952e-2, 2.093e-2, 1.199e-2, 7.666e-3, 5.333e-3,
  4.666e-3
};

static double interpolate_factor(double column, const double *log_columns,
                                 const double *factors, size_t count)
{
  if(!(column > 0.0) || !isfinite(column))
    return column == INFINITY ? factors[count - 1] : 1.0;

  const double log_column = log10(column);
  if(log_column <= log_columns[0])
    return 1.0;
  if(log_column >= log_columns[count - 1])
    return factors[count - 1];

  size_t upper = 1;
  while(upper < count && log_column > log_columns[upper])
    upper++;

  const double fraction = (log_column - log_columns[upper - 1]) /
      (log_columns[upper] - log_columns[upper - 1]);
  return factors[upper - 1] +
      fraction * (factors[upper] - factors[upper - 1]);
}

double radiation_co_line_shielding_lee1996(double column_co,
                                           double column_h2)
{
  const double self = interpolate_factor(
      column_co, co_columns, co_factors,
      sizeof(co_columns) / sizeof(co_columns[0]));
  const double overlap = interpolate_factor(
      column_h2, h2_columns, h2_factors,
      sizeof(h2_columns) / sizeof(h2_columns[0]));
  return self * overlap;
}

double radiation_co_photodissociation_rate(double g_fuv,
                                           double column_co,
                                           double column_h2,
                                           double dust_optical_depth)
{
  if(!(g_fuv > 0.0) || !isfinite(g_fuv))
    return 0.0;

  double tau = dust_optical_depth;
  if(!(tau > 0.0) || !isfinite(tau))
    tau = tau == INFINITY ? INFINITY : 0.0;

  return 1.03e-10 * g_fuv *
      radiation_co_line_shielding_lee1996(column_co, column_h2) * exp(-tau);
}
