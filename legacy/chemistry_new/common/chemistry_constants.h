#pragma once

#define N_THREEBODY 1.0e14
#define N_EQCHEM 1.0e18

#ifndef T_CMB0
#define T_CMB0 2.728
#endif

#define eV_to_K 11606.0
#define eV_to_erg 1.60184e-12
#define eV_to_Hz 2.41838e14
#define eV_to_um_inv 0.80669
#define kJ_per_mol_to_eV 1.04e-2

#define Z_SUN 0.02
#define DUST_TO_GAS_MASS_RATIO 0.00934
#define STEFAN 5.6692e-5

#define PROTONMASS 1.6726e-24
#define BOLTZMANN 1.3806e-16
#define GRAVITY 6.672e-8
#define PLANCK 6.6262e-27
#define ELECTRONMASS 9.10953e-28
#define MAX_REAL_NUMBER 1.e30

#ifdef METALLICITY
#define Y_C_sun 0.927e-4
#define Y_O_sun 3.568e-4
#define A_C_sun 3.16e-4
#define A_O_sun 5.75e-4
#define A_Mg_sun 3.47e-5
#define A_Si_sun 3.31e-5
#define A_Fe_sun 3.31e-5
#define Zmet_sun 0.0153
#endif

static inline double dmax(double left, double right)
{
  return left > right ? left : right;
}

static inline double dmin(double left, double right)
{
  return left < right ? left : right;
}
