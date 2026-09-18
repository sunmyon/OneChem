#include <math.h>

#include "chemistry_constants.h"
#include "forward_rates.h"
#include "chemistry_tables.h"

#ifdef CHEM_HIRANO
#ifndef O05_CHEM
#define O05_CHEM
#endif
#endif

#if defined(TBR_FLOW07) + defined(TBR_ABEL02) + defined(TBR_FORR13) > 1
#error "Select only one three-body forward-rate prescription"
#endif

static double rate_zero(double T)
{
  (void)T;
  return 0.;
}

#ifdef CARBON_OXYGEN_CHEMISTRY
/* Omukai (2005), Table 1.  Metal-network rates use table slots 100 + Z. */
static double rate_Z1(double T)  { (void)T; return 6.80e-10; }
static double rate_Z2(double T)  { return 7.00e-10 * exp(-232.0 / T); }
static double rate_Z3(double T)  { (void)T; return 1.70e-9; }
static double rate_Z4(double T)  { (void)T; return 1.01e-9; }
static double rate_Z5(double T)  { (void)T; return 8.30e-10; }
static double rate_Z6(double T)  { return 2.00e-7 * pow(T / 300.0, -0.50); }
static double rate_Z7(double T)  { return 1.60e-7 * pow(T / 300.0, -0.50); }
static double rate_Z8(double T)  { return 6.50e-7 * pow(T / 300.0, -0.50); }
static double rate_Z9(double T)  { return 3.50e-7 * pow(T / 300.0, -0.50); }
static double rate_Z10(double T) { return 9.90e-19 * pow(T / 300.0, -0.38); }
static double rate_Z11(double T) { return 3.43e-13 * pow(T / 300.0, 2.67) * exp(-3160.0 / T); }
static double rate_Z12(double T) { return 1.55e-12 * pow(T / 300.0, 1.60) * exp(-1660.0 / T); }
static double rate_Z13(double T) { return 1.65e-12 * pow(T / 300.0, 1.14) * exp(-50.0 / T); }
static double rate_Z14(double T) { (void)T; return 2.10e-9; }
static double rate_Z15(double T) { (void)T; return 6.90e-9; }
static double rate_Z16(double T) { return 7.00e-14 * pow(T / 300.0, 2.80) * exp(-1950.0 / T); }
static double rate_Z17(double T) { return 6.83e-12 * pow(T / 300.0, 1.60) * exp(-9720.0 / T); }
static double rate_Z18(double T) { return 4.90e-20 * pow(T / 300.0, 1.58); }
static double rate_Z19(double T) { return 4.34e-11 * pow(T / 300.0, -0.50) * exp(-30.0 / T); }
static double rate_Z20(double T) { return 3.30e-10 * exp(-8460.0 / T); }
static double rate_Z21(double T) { (void)T; return 2.00e-9; }
static double rate_Z22(double T) { return 1.95e-7 * pow(T / 300.0, -0.70); }
static double rate_Z23(double T) { (void)T; return 6.60e-11; }
static double rate_Z24(double T) { return 1.10e-10 * pow(T / 300.0, 0.50); }
static double rate_Z25(double T) { (void)T; return 3.30e-11; }
static double rate_Z26(double T) { (void)T; return 6.20e-10; }
static double rate_Z27(double T) { (void)T; return 1.00e-13; }
static double rate_Z28(double T) { return 4.40e-12 * pow(T / 300.0, -0.61); }
static double rate_Z29(double T) { (void)T; return 7.70e-10; }
static double rate_Z30(double T) { (void)T; return 7.50e-10; }
static double rate_Z31(double T) { (void)T; return 1.00e-17; }
static double rate_Z32(double T) { return 6.64e-10 * exp(-11700.0 / T); }
static double rate_Z33(double T) { (void)T; return 4.98e-11; }
static double rate_Z34(double T) { return 2.38e-10 * exp(-1760.0 / T); }
static double rate_Z35(double T) { (void)T; return 2.70e-10; }
static double rate_Z40(double T) { (void)T; return 1.00e-17; }

#ifdef CARBON_OXYGEN_HIGH_IONIZATION
#define METAL_ION_KB_EV 8.617385e-5

typedef struct { double threshold, p, a, x, k; } MetalCollisionalFit;
typedef struct { double a, b, t0, t1; } MetalRadiativeFit;
typedef struct { double a, b, c, d, tmin, tmax; } MetalChargeTransferFit;
typedef struct {
  double a, b, c, d, tmin, tmax, delta_t;
} MetalChargeTransferIonizationFit;
typedef struct { int count; double c[6], energy[6]; } MetalDielectronicFit;

/* C I--III followed by O I--III. Fits match the Gadget3 reference module. */
static const MetalCollisionalFit metal_collisional_fit[6] = {
  {11.3, 0., 6.85e-8, .1930, .25}, {24.4, 1., 1.86e-8, .2860, .24},
  {47.9, 1., 6.35e-9, .4270, .21}, {13.6, 0., 3.59e-8, .0730, .34},
  {35.1, 1., 1.39e-8, .2120, .22}, {54.9, 1., 9.31e-9, .2700, .27}
};
static const MetalRadiativeFit metal_radiative_fit[6] = {
  {7.651e-9, .8027, .001193, 9.334e12},
  {8.577e-10, .7837, .7286, 1.140e7},
  {2.020e-9, .7798, .6690, 2.425e6},
  {1.341e-10, .6159, 1.673, 6.366e16},
  {1.816e-8, .7170, .006717, 3.286e12},
  {3.672e-9, .7676, .2900, 8.521e7}
};
static const MetalDielectronicFit metal_dielectronic_fit[6] = {
  {5, {6.346e-9, 9.793e-9, 1.634e-6, 8.369e-4, 3.355e-4},
      {12.17, 73.8, 1.523e4, 1.207e5, 2.144e5}},
  {6, {3.489e-6, 2.222e-7, 1.954e-5, 4.212e-3, 2.037e-4, 2.936e-4},
      {2660., 3756., 2.566e4, 1.400e5, 1.801e6, 4.307e6}},
  {6, {4.673e-7, 1.887e-5, 1.305e-5, 3.099e-3, 3.001e-4, 2.553e-3},
      {723.3, 2847., 1.054e4, 8.915e4, 2.812e5, 3.254e6}},
  {4, {5.629e-8, 2.550e-7, 6.173e-4, 1.627e-4},
      {5395., 1.770e4, 1.671e5, 2.687e5}},
  {6, {1.627e-7, 1.262e-7, 6.663e-7, 3.925e-6, 2.406e-3, 1.146e-3},
      {45.35, 284.7, 4166., 2.877e4, 1.953e5, 3.646e5}},
  {5, {3.932e-7, 2.523e-7, 3.447e-5, 5.776e-3, 5.101e-3},
      {150.9, 621.1, 1.562e4, 1.936e5, 4.700e5}}
};
static const MetalChargeTransferFit metal_ct_recombination_fit[6] = {
  {1.60e-4, 2.80, 350., -4.20, 5.e3, 5.e4},
  {3.65, .25, 0., 0., 5.e3, 5.e4},
  {.76, 1.40, 0., 0., 1.e3, 1.e4},
  {.80, .43, 0., 0., 5.e3, 5.e4},
  {.85, .44, 0., 0., 5.e3, 5.e4},
  {.31, .44, 0., 0., 1.e3, 3.e4}
};
static const MetalChargeTransferIonizationFit metal_ct_ionization_fit[6] = {
  {1.07e-6, 3.15, 176.43, -4.29, 1.e3, 1.e5, 0.},
  {0., 0., 0., 0., 0., 0., 0.},
  {0., 0., 0., 0., 0., 0., 0.},
  {7.40e-2, .47, 24.37, -.74, 1.e1, 1.e4, .023},
  {0., 0., 0., 0., 0., 0., 0.},
  {0., 0., 0., 0., 0., 0., 0.}
};

static double metal_collisional_rate(int transition, double T)
{
  const MetalCollisionalFit *fit = &metal_collisional_fit[transition];
  if(T <= 0.0)
    return 0.0;
  const double u = fit->threshold / (METAL_ION_KB_EV * T);
  if(u > 80.0)
    return 0.0;
  return fit->a * (1.0 + fit->p * sqrt(u)) / (fit->x + u) *
      pow(u, fit->k) * exp(-u);
}

static double metal_recombination_rate(int transition, double T)
{
  const MetalRadiativeFit *radiative = &metal_radiative_fit[transition];
  const MetalDielectronicFit *dielectronic =
      &metal_dielectronic_fit[transition];
  if(T <= 0.0)
    return 0.0;
  const double x = sqrt(T / radiative->t0);
  double rate = radiative->a /
      (x * pow(1.0 + x, 1.0 - radiative->b) *
       pow(1.0 + sqrt(T / radiative->t1), 1.0 + radiative->b));
  double sum = 0.0;
  for(int term = 0; term < dielectronic->count; term++)
    sum += dielectronic->c[term] * exp(-dielectronic->energy[term] / T);
  return rate + fmax(0.0, sum * pow(T, -1.5));
}

static double metal_ct_recombination_rate(int transition, double T)
{
  const MetalChargeTransferFit *fit =
      &metal_ct_recombination_fit[transition];
  if(T <= 0.0 || fit->a <= 0.0)
    return 0.0;
  const double limited_T = fmin(fmax(T, fit->tmin), fit->tmax);
  return 1.0e-9 * fit->a * pow(limited_T / 1.0e4, fit->b) *
      (1.0 + fit->c * exp(fit->d * limited_T / 1.0e4));
}

static double metal_ct_ionization_rate(int transition, double T)
{
  const MetalChargeTransferIonizationFit *fit =
      &metal_ct_ionization_fit[transition];
  if(T <= 0.0 || fit->a <= 0.0)
    return 0.0;
  const double limited_T = fmin(fmax(T, fit->tmin), fit->tmax);
  return 1.0e-9 * fit->a * pow(limited_T / 1.0e4, fit->b) *
      (1.0 + fit->c * exp(fit->d * limited_T / 1.0e4)) *
      exp(-fit->delta_t * 1.0e4 / T);
}

#define DEFINE_METAL_RATE_WRAPPERS(index) \
  static double rate_metal_ci_##index(double T) \
    { return metal_collisional_rate(index, T); } \
  static double rate_metal_rr_##index(double T) \
    { return metal_recombination_rate(index, T); } \
  static double rate_metal_ctr_##index(double T) \
    { return metal_ct_recombination_rate(index, T); } \
  static double rate_metal_cti_##index(double T) \
    { return metal_ct_ionization_rate(index, T); }
DEFINE_METAL_RATE_WRAPPERS(0)
DEFINE_METAL_RATE_WRAPPERS(1)
DEFINE_METAL_RATE_WRAPPERS(2)
DEFINE_METAL_RATE_WRAPPERS(3)
DEFINE_METAL_RATE_WRAPPERS(4)
DEFINE_METAL_RATE_WRAPPERS(5)
#undef DEFINE_METAL_RATE_WRAPPERS
#endif
#endif

static double rate_1(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8;

  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;

  if(T/eV_to_K < 0.01)
    return 0.;

  return exp(-32.71396786375
             + 13.53655609057      * x
             - 5.739328757388      * x2
             + 1.563154982022      * x3
             - 0.2877056004391     * x4
             + 0.03482559773736999 * x5
             - 0.00263197617559    * x6
             + 0.0001119543953861  * x7
             - 2.039149852002e-6   * x8);
}

static double rate_2_case_a(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8, x9;

#ifdef CHEM_HIRANO
  if(T < 5500.)
    return 0.;
#endif

  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;
  x9 = x*x8;

  return exp(-28.61303380689232
             - 0.7241125657826851     * x
             - 0.02026044731984691    * x2
             - 0.002380861877349834   * x3
             - 0.0003212605213188796  * x4
             - 0.00001421502914054107 * x5
             + 4.989108920299513e-6   * x6
             + 5.755614137575758e-7   * x7
             - 1.856767039775261e-8   * x8
             - 3.071135243196595e-9   * x9);
}

static double rate_2_case_b(double T)
{
#ifdef CHEM_HIRANO
  if(T < 5500.)
    return 0.;
#endif
  return 2.753e-14 * pow(315614.0/T, 1.5)/pow(1.0 + pow(115188.0/T, 0.406), 2.242);
}

static double rate_3(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8;

  if(T/eV_to_K < 0.01)
    return 0.;

  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;

  return exp(-44.09864886561001
             + 23.91596563469       * x
             - 10.75323019821       * x2
             + 3.058038757198       * x3
             - 0.5685118909884001   * x4
             + 0.06795391233790001  * x5
             - 0.005009056101857001 * x6
             + 0.0002067236157507   * x7
             - 3.649161410833e-6    * x8);
}

static double rate_4_case_a(double T)
{
  double T_eV = T/eV_to_K;

  if(T_eV < 0.01)
    return 3.92e-13/pow(T_eV, 0.6353);

  return 1.54e-9*(1. + 0.3/exp(8.099328789667/T_eV))
    /(exp(40.49664394833662/T_eV)*pow(T_eV, 1.5))
    + 3.92e-13/pow(T_eV, 0.6353);
}

static double rate_4_case_b(double T)
{
  double radiative = 1.26e-14*pow(2.*285335.0/T, 0.75);

  if(T/eV_to_K < 0.01)
    return radiative;

  return radiative + 1.9e-3*pow(T, -1.5)*exp(-470000./T)*(1. + 0.3*exp(-94000./T));
}

static double rate_5(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8;

  if(T/eV_to_K < 0.8)
    return 0.;

  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;

  return exp(-68.71040990212001
             + 43.93347632635       * x
             - 18.48066993568       * x2
             + 4.701626486759002    * x3
             - 0.7692466334492      * x4
             + 0.08113042097303     * x5
             - 0.005324020628287001 * x6
             + 0.0001975705312221   * x7
             - 3.165581065665e-6    * x8);
}

static double rate_6_case_a(double T)
{
  return 3.36e-10/sqrt(T)/pow(T/1.e3, 0.2)/(1. + pow(T/4.e6, 0.7));
}

static double rate_6_case_b(double T)
{
  return 2.0*2.753e-14*pow(315614.0/(T/4.0), 1.5)
    /pow(1.0 + pow(115188.0/(T/4.0), 0.406), 2.242);
}

static double rate_7_enzo(double T)
{
  return 6.775e-15*pow(T/eV_to_K, 0.8779);
}

static double rate_7_gp98(double T)
{
  return 1.4e-18*pow(T, 0.928)*exp(-T/16200.0);
}

static double rate_7_matsukoba(double T)
{
  double x = log10(T);

  if(T < 6000.)
    return pow(10., -17.845 + 0.762*x + 0.1523*x*x - 0.03274*x*x*x);

  return pow(10., -16.4199 + 0.1998*x*x - 5.447e-3*x*x*x*x + 4.0415e-5*x*x*x*x*x*x);
}

static double rate_8_enzo(double T)
{
  double x, x2, x3, x4, x5, x6, x7;

  if(T/eV_to_K <= 0.1)
    return 1.43e-9;

  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;

  return exp(-20.06913897587003
             + 0.2289800603272916    * x
             + 0.03599837721023835   * x2
             - 0.004555120027032095  * x3
             - 0.0003105115447124016 * x4
             + 0.0001073294010367247 * x5
             - 8.36671960467864e-6   * x6
             + 2.238306228891639e-7  * x7);
}

static double rate_8_gp98(double T)
{
  if(T > 300.)
    return 4.0e-9*pow(T, -0.17);
  return 1.516e-9;
}

static double rate_8_matsukoba(double T)
{
  return 1.35e-9*(pow(T, 0.098493) + 0.32852*pow(T, 0.5561) + 2.771e-7*pow(T, 2.1826))
    /(1. + 6.191e-3*pow(T, 1.0461) + 8.9712e-11*pow(T, 3.0424) + 3.2576e-14*pow(T, 3.7741));
}

static double rate_9(double T)
{
  if(T > 6.7e3)
    return 5.81e-16*pow(T/56200., -0.6657*log10(T/56200.));
  return 1.85e-23*pow(T, 1.8);
}

static double rate_10(double T)
{
  (void)T;
  return 6.0e-10;
}

static double rate_11_default(double T)
{
  double T_eV = T/eV_to_K;
  if(T_eV <= 0.3)
    return 0.;
  return 1.0670825e-10*pow(T_eV, 2.012)/(exp(4.463/T_eV)*pow(1. + 0.2472*T_eV, 3.512));
}

static double rate_11_enzo(double T)
{
  double T_eV = T/eV_to_K;
  return 1.0670825e-10*pow(T_eV, 2.012)/(exp(4.463/T_eV)*pow(1. + 0.2472*T_eV, 3.512));
}

static double log_rate_12_default(double T)
{
  if(T < 10000.)
    return log(3.0e-10) - 21050.0/T;
  return log(1.5e-10) - 14000.0/T;
}

static double rate_12_default(double T)
{
  if(T < 10000.)
    return 3.0e-10*exp(-21050.0/T);
  return 1.5e-10*exp(-14000.0/T);
}

static double log_rate_12_enzo(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8;
  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;
  return -24.249146877095
             + 3.400824447095291     * x
             - 3.898003964650152     * x2
             + 2.045587822403071     * x3
             - 0.5416182856220388    * x4
             + 0.0841077503763412    * x5
             - 0.007879026154483455  * x6
             + 0.0004138398421504563 * x7
             - 9.36345888928611e-6   * x8;
}

static double rate_12_enzo(double T)
{
  return exp(log_rate_12_enzo(T));
}

double chemistry_log_forward_rate_12(double temperature)
{
#ifdef ENZO_CHEM
  return log_rate_12_enzo(temperature);
#else
  return log_rate_12_default(temperature);
#endif
}

static double rate_13_default(double T)
{
  if(T/eV_to_K <= 0.3)
    return 0.;
  return 4.38e-10*exp(-102000.0/T)*pow(T, 0.35);
}

static double rate_13_enzo(double T)
{
  return 5.6e-11*exp(-102124.0/T)*pow(T, 0.5);
}

static double rate_13_matsukoba(double T)
{
  return 1.91e-9*pow(T, 0.136)*exp(-53407.1/T);
}

static double rate_14_default(double T) { return 1.2e-15*pow(T/300.0, 0.25); }

static double rate_15_default(double T)
{
  if(T > 10000.)
    return 4.0e-37*pow(T, 4.74);
  return 1.26e-9*pow(T, -0.75)*exp(-127500.0/T);
}

static double rate_16_default(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8;
  if(T/eV_to_K <= 0.04)
    return 0.;
  x = log(T/eV_to_K);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;
  return exp(-18.01849334273
             + 2.360852208681      * x
             - 0.2827443061704     * x2
             + 0.01623316639567    * x3
             - 0.03365012031362999 * x4
             + 0.01178329782711    * x5
             - 0.001656194699504   * x6
             + 0.0001068275202678  * x7
             - 2.631285809207e-6   * x8);
}

static double rate_17_default(double T)
{
  double x, x2, x3, x4, x5, x6, x7, x8, x9;
  double T_eV = T/eV_to_K;
  if(T_eV <= 0.1)
    return 2.56e-9*pow(T_eV, 1.78186);
  x = log(T_eV);
  x2 = x*x;
  x3 = x*x2;
  x4 = x*x3;
  x5 = x*x4;
  x6 = x*x5;
  x7 = x*x6;
  x8 = x*x7;
  x9 = x*x8;
  return exp(-20.37260896533324
             + 1.139449335841631   * x
             - 0.1421013521554148  * x2
             + 0.00846445538663    * x3
             - 0.0014327641212992  * x4
             + 0.0002012250284791  * x5
             + 0.0000866396324309  * x6
             - 0.00002585009680264 * x7
             + 2.4555011970392e-6  * x8
             - 8.06838246118e-8    * x9);
}

static double rate_18_default(double T)
{
  return 6.92e-14*pow(T, 0.5) + 8.00e-13*pow(T, 0.83) + 6.18e-7*pow(T, -0.27);
}

static double rate_18_enzo(double T) { return 6.5e-9/sqrt(T/eV_to_K); }
static double rate_18_hirano(double T) { return 6.3e-8 + 5.7e-6/sqrt(T) - 9.2e-11*sqrt(T) + 4.4e-13*T; }
static double rate_18_matsukoba(double T) { return 2.4e-6*pow(T, -0.5)*(1. + T/2.e4); }

static double rate_19(double T)
{
  if(T > 1.e4)
    return 4.0e-4*pow(T, -1.4)*exp(-15100.0/T);
  return 1.0e-8*pow(T, -0.4);
}

static double rate_20_default(double T) { return 2.0e-7/pow(T, 0.5); }
static double rate_20_enzo(double T) { return T > 617. ? 1.32e-6*pow(T, -0.76) : 1.e-8; }
static double rate_21(double T) { return 4.64e-8/sqrt(T/eV_to_K); }
static double rate_30_matsukoba(double T) { return 1.2e-17*pow(T, 1.2)*exp(-1.578e5/T); }

#ifdef HD
static double rate_31_hirano(double T) { return 3.6e-12*pow(T/300.0, -0.75); }
static double rate_32_hirano(double T) { return 2.0e-10*pow(T, 0.402)*exp(-37.1/T) - 3.31e-17*pow(T, 1.48); }
static double rate_32_default(double T) { return 3.7e-10*pow(T, 0.28)*exp(-43.0/T); }
static double rate_33_hirano(double T) { return 2.06e-10*pow(T, 0.396)*exp(-33.0/T) + 2.03e-9*pow(T, -0.332); }
static double rate_33_default(double T) { return 3.7e-10*pow(T, 0.28); }
static double rate_34_hirano(double T) { (void)T; return 1.0e-25; }
static double rate_35(double T) { return 9.0e-11*exp(-3876.0/T); }
static double rate_37_hirano(double T) { (void)T; return 1.6e-9; }
static double rate_37_default(double T) { (void)T; return 2.1e-9; }
static double rate_38(double T) { return 3.2e-11*exp(-3624.0/T); }
static double rate_39(double T) { return 1.0e-9*exp(-464.0/T); }
#ifdef HD_PLUS
static double rate_36(double T) { (void)T; return 6.4e-10; }
static double rate_40(double T)
{
  double x;
  if(T <= 10.)
    return 0.;
  x = log10(T);
  return pow(10.0, -19.38 - 1.523*x + 1.118*x*x - 0.1269*x*x*x);
}
static double rate_42(double T) { return 7.2e-8/sqrt(T); }
#endif
#ifdef DMINUS
static double rate_43(double T) { return 3.0e-16*pow(T/300.0, 0.95)*exp(-T/9320.0); }
static double rate_44(double T) { return 5.7e-8*pow(T/300.0, -0.5); }
static double rate_45(double T) { return 4.6e-8*pow(T/300.0, -0.5); }
static double rate_46(double T) { return 6.4e-9*pow(T/300.0, 0.41); }
static double rate_48(double T) { return 1.5e-9*pow(T/300.0, -0.1); }
#endif
#endif

static double rate_55_flower07(double T) { return 1.44e-26/pow(T, 1.54); }
static double rate_55_abel02(double T) { return T < 300. ? 1.14e-31/pow(T, 0.38) : 3.9e-30/T; }
static double rate_55_forrey13(double T) { return 6.e-32/pow(T, 0.25) + 2.e-31/pow(T, 0.5); }
static double rate_56(double T) { return 6.5e-7/sqrt(T)*exp(-52000.0/T)*(1.0 - exp(-6000.0/T)); }
static double rate_57(double T) { return (4.5e-29/T)/40.0; }
static double rate_58(double T) { return 5.0e-8/sqrt(T)*exp(-52000.0/T)*(1.0 - exp(-6000.0/T)); }

#ifdef CASE_B
#define FORWARD_RATE_2 rate_2_case_b
#define FORWARD_RATE_4 rate_4_case_b
#define FORWARD_RATE_6 rate_6_case_b
#else
#define FORWARD_RATE_2 rate_2_case_a
#define FORWARD_RATE_4 rate_4_case_a
#define FORWARD_RATE_6 rate_6_case_a
#endif

#ifdef ENZO_CHEM
#define FORWARD_RATE_7 rate_7_enzo
#define FORWARD_RATE_8 rate_8_enzo
#elif defined(O05_CHEM)
#define FORWARD_RATE_7 rate_7_gp98
#define FORWARD_RATE_8 rate_8_gp98
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
#define FORWARD_RATE_7 rate_7_matsukoba
#define FORWARD_RATE_8 rate_8_matsukoba
#else
#define FORWARD_RATE_7 rate_7_gp98
#define FORWARD_RATE_8 rate_8_gp98
#endif

#ifdef ENZO_CHEM
#define FORWARD_RATE_11 rate_11_enzo
#define FORWARD_RATE_12 rate_12_enzo
#define FORWARD_RATE_14 rate_zero
#define FORWARD_RATE_15 rate_zero
#define FORWARD_RATE_20 rate_20_enzo
#else
#define FORWARD_RATE_11 rate_11_default
#define FORWARD_RATE_12 rate_12_default
#define FORWARD_RATE_14 rate_14_default
#define FORWARD_RATE_15 rate_15_default
#define FORWARD_RATE_20 rate_20_default
#endif

#ifdef ENZO_CHEM
#define FORWARD_RATE_13 rate_13_enzo
#elif defined(O05_CHEM)
#define FORWARD_RATE_13 rate_13_default
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
#define FORWARD_RATE_13 rate_13_matsukoba
#else
#define FORWARD_RATE_13 rate_13_default
#endif

#ifdef ENZO_CHEM
#define FORWARD_RATE_18 rate_18_enzo
#elif defined(CHEM_HIRANO)
#define FORWARD_RATE_18 rate_18_hirano
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
#define FORWARD_RATE_18 rate_18_matsukoba
#else
#define FORWARD_RATE_18 rate_18_default
#endif

#if !defined(ENZO_CHEM) && !defined(CHEM_HIRANO) && defined(CONTINUUM_COOLING_MATSUKOBA)
#define FORWARD_RATE_16 rate_zero
#define FORWARD_RATE_17 rate_zero
#else
#define FORWARD_RATE_16 rate_16_default
#define FORWARD_RATE_17 rate_17_default
#endif

#ifdef O05_CHEM
#define FORWARD_RATE_30 rate_zero
#elif defined(CONTINUUM_COOLING_MATSUKOBA)
#define FORWARD_RATE_30 rate_30_matsukoba
#else
#define FORWARD_RATE_30 rate_zero
#endif

#if defined(TBR_FLOW07)
#define FORWARD_RATE_55 rate_55_flower07
#elif defined(TBR_ABEL02)
#define FORWARD_RATE_55 rate_55_abel02
#else
#define FORWARD_RATE_55 rate_55_forrey13
#endif

const ForwardRateFunction forward_rate_table[CHEM_RATE_COUNT] = {
  [1]  = rate_1,
  [2]  = FORWARD_RATE_2,
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  [3]  = rate_3,
  [4]  = FORWARD_RATE_4,
  [5]  = rate_5,
  [6]  = FORWARD_RATE_6,
#endif
  [7]  = FORWARD_RATE_7,
  [8]  = FORWARD_RATE_8,
#ifdef H2_PLUS
  [9]  = rate_9,
  [10] = rate_10,
#endif
  [11] = FORWARD_RATE_11,
#ifdef H2_PLUS
  [12] = FORWARD_RATE_12,
#endif
  [13] = FORWARD_RATE_13,
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  [14] = FORWARD_RATE_14,
  [15] = FORWARD_RATE_15,
#endif
  [16] = FORWARD_RATE_16,
  [17] = FORWARD_RATE_17,
  [18] = FORWARD_RATE_18,
#ifdef H2_PLUS
  [19] = rate_19,
  [20] = FORWARD_RATE_20,
  [21] = rate_21,
#endif
  [30] = FORWARD_RATE_30,
#ifdef HD
#ifdef CHEM_HIRANO
  [31] = rate_31_hirano,
  [32] = rate_32_hirano,
  [33] = rate_33_hirano,
  [34] = rate_34_hirano,
  [37] = rate_37_hirano,
#else
  [31] = rate_zero,
  [32] = rate_32_default,
  [33] = rate_33_default,
  [34] = rate_zero,
  [37] = rate_37_default,
#endif
  [35] = rate_35,
  [38] = rate_38,
  [39] = rate_39,
#ifdef HD_PLUS
  [36] = rate_36,
  [40] = rate_40,
  [41] = rate_40,
  [42] = rate_42,
#endif
#ifdef DMINUS
  [43] = rate_43,
  [44] = rate_44,
  [45] = rate_45,
  [46] = rate_46,
  [47] = rate_46,
  [48] = rate_48,
#endif
#endif
  [55] = FORWARD_RATE_55,
  [56] = rate_56,
  [57] = rate_57,
  [58] = rate_58,
#ifdef CARBON_OXYGEN_CHEMISTRY
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
  [101] = rate_metal_ctr_3,
  [102] = rate_metal_cti_3,
#else
  [101] = rate_Z1,
  [102] = rate_Z2,
#endif
  [103] = rate_Z3,
  [104] = rate_Z4,
  [105] = rate_Z5,
  [106] = rate_Z6,
  [107] = rate_Z7,
  [108] = rate_Z8,
  [109] = rate_Z9,
  [110] = rate_Z10,
  [111] = rate_Z11,
  [112] = rate_Z12,
  [113] = rate_Z13,
  [114] = rate_Z14,
  [115] = rate_Z15,
  [116] = rate_Z16,
  [117] = rate_Z17,
  [118] = rate_Z18,
  [119] = rate_Z19,
  [120] = rate_Z20,
  [121] = rate_Z21,
  [122] = rate_Z22,
  [123] = rate_Z23,
  [124] = rate_Z24,
  [125] = rate_Z25,
  [126] = rate_Z26,
  [127] = rate_Z27,
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
  [128] = rate_metal_rr_0,
#else
  [128] = rate_Z28,
#endif
  [129] = rate_Z29,
  [130] = rate_Z30,
  [131] = rate_Z31,
  [132] = rate_Z32,
  [133] = rate_Z33,
  [134] = rate_Z34,
  [135] = rate_Z35,
  [140] = rate_Z40,
#ifdef CARBON_OXYGEN_HIGH_IONIZATION
  [141] = rate_metal_ci_0,
  [142] = rate_metal_ci_1,
  [143] = rate_metal_ci_2,
  [144] = rate_metal_ci_3,
  [145] = rate_metal_ci_4,
  [146] = rate_metal_ci_5,
  [148] = rate_metal_rr_1,
  [149] = rate_metal_rr_2,
  [150] = rate_metal_rr_3,
  [151] = rate_metal_rr_4,
  [152] = rate_metal_rr_5,
  [159] = rate_metal_ctr_0,
  [160] = rate_metal_ctr_1,
  [161] = rate_metal_ctr_2,
  [163] = rate_metal_ctr_4,
  [164] = rate_metal_ctr_5,
  [165] = rate_metal_cti_0,
  [166] = rate_metal_cti_1,
  [167] = rate_metal_cti_2,
  [169] = rate_metal_cti_4,
  [170] = rate_metal_cti_5,
#endif
#endif
};
