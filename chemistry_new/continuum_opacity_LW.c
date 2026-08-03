#if defined(CHEMISTRY_SC) && defined(LYMANWERNER)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>

#include "../allvars.h"
#include "../proto.h"
#include "chemistry.h"

double partition_H(double T_gas);
double Hminus_detachment_John1988(double hv);

enum dust_species{
  dWater,
  dVolorg,
  dReforg,
  dTroilite,
  dIron,
  dOlivine,
  dOrthopy,
  dEND
};

static double kappa_dust_abs[7];
static double kappa_dust_scat[7];
static double mass_fraction_dust[7];
static double kappa_dust_T[5];

const double E_LW = 11.2;
const double E_Lyl = 13.6;
const double E_2s2p = 10.2;

static double Sigma_HM_bf_at_ELW;
static double Sigma_HI_bf_2s2p_at_ELW;
static double Sigma_H_Rayleigh_at_ELW;
static double Sigma_Thomson;

void tabulate_opacities_in_LW_band(void){
  kappa_dust_abs[dWater] = 4.66888e+01;
  kappa_dust_scat[dWater] = 4.95918e+01;
  mass_fraction_dust[dWater] = 1.19e-3 / DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dVolorg] = 5.93716e+01;
  kappa_dust_scat[dVolorg] = 3.36187e+01;
  mass_fraction_dust[dVolorg] = 6.02e-4 / DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dReforg] = 2.32095e+02;
  kappa_dust_scat[dReforg] = 1.31422e+02;
  mass_fraction_dust[dReforg] = 3.53e-3/ DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dTroilite] = 9.60442e+00;
  kappa_dust_scat[dTroilite] = 1.06814e+01;
  mass_fraction_dust[dTroilite] = 5.69e-4/ DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dIron] = 3.21235e+00;
  kappa_dust_scat[dIron] = 1.81988e+00;
  mass_fraction_dust[dIron] = 2.53e-4/ DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dOlivine] = 1.06432e+01;
  kappa_dust_scat[dOlivine] = 6.92057e+01;
  mass_fraction_dust[dOlivine] = 2.51e-3/ DUST_TO_GAS_MASS_RATIO;

  kappa_dust_abs[dOrthopy] = 1.49784e+00;
  kappa_dust_scat[dOrthopy] = 2.13093e+01;
  mass_fraction_dust[dOrthopy] = 7.33e-4/ DUST_TO_GAS_MASS_RATIO;

  double sum = 0.;
  sum += kappa_dust_abs[dWater] * mass_fraction_dust[dWater];
  sum += kappa_dust_scat[dWater] * mass_fraction_dust[dWater];
  kappa_dust_T[0] = sum;

  sum += kappa_dust_abs[dVolorg] * mass_fraction_dust[dVolorg];
  sum += kappa_dust_scat[dVolorg] * mass_fraction_dust[dVolorg];
  kappa_dust_T[1] = sum;

  sum += kappa_dust_abs[dReforg] * mass_fraction_dust[dReforg];
  sum += kappa_dust_scat[dReforg] * mass_fraction_dust[dReforg];
  kappa_dust_T[2] = sum;

  sum += kappa_dust_abs[dTroilite] * mass_fraction_dust[dTroilite];
  sum += kappa_dust_scat[dTroilite] * mass_fraction_dust[dTroilite];
  kappa_dust_T[3] = sum;

  sum += kappa_dust_abs[dIron] * mass_fraction_dust[dIron];
  sum += kappa_dust_scat[dIron] * mass_fraction_dust[dIron];

  sum += kappa_dust_abs[dOlivine] * mass_fraction_dust[dOlivine];
  sum += kappa_dust_scat[dOlivine] * mass_fraction_dust[dOlivine];

  sum += kappa_dust_abs[dOrthopy] * mass_fraction_dust[dOrthopy];
  sum += kappa_dust_scat[dOrthopy] * mass_fraction_dust[dOrthopy];
  kappa_dust_T[4] = sum;

  Sigma_HM_bf_at_ELW = Hminus_detachment_John1988(E_LW);
  Sigma_HI_bf_2s2p_at_ELW = 1.58e-17 * pow(4.*E_LW/E_Lyl, -3);

  double lambda_i = 1. / (eV_to_um_inv * E_LW);
  Sigma_H_Rayleigh_at_ELW = 5.799e-29*pow(lambda_i, -4) + 1.422e-30*pow(lambda_i, -6) + 2.784e-32*pow(lambda_i,-8);;
  Sigma_Thomson = 6.6525e-25;
}


double calc_continuum_opacities_LW(double *abund, double Tdust, double metallicity){
  double T_gas = compute_temp(ithis);
  double beta_LW = E_LW * eV_to_erg/ BOLTZMANN / T_gas;
  double beta_2s2p = E_2s2p * eV_to_erg/ BOLTZMANN / T_gas;

  const double g2 = 8.;

  double f_HM = abund[sHM];
  double f_HI = abund[sHI];
  double f_HI_2s2p = f_HI * g2 * exp(-beta_2s2p)/partition_H(T_gas);
  double f_elec = abund[sELEC];

  double kappa_nu = (Sigma_HM_bf_at_ELW * f_HM + Sigma_HI_bf_2s2p_at_ELW * f_HI_2s2p) * (1. - exp(-beta_LW)) 
    + Sigma_H_Rayleigh_at_ELW * f_HI + Sigma_Thomson * f_elec;

  double kappa_gas[4];
  kappa_gas[0] = Sigma_HM_bf_at_ELW * f_HM * (1. - exp(-beta_LW));
  kappa_gas[1] = Sigma_HI_bf_2s2p_at_ELW * f_HI_2s2p * (1. - exp(-beta_LW));
  kappa_gas[2] = Sigma_H_Rayleigh_at_ELW * f_HI;
  kappa_gas[3] = Sigma_Thomson * f_elec;
  

  double kappa_dust = 0.;
#ifdef METALLICITY
  if(T_dust < 159.)
    kappa_dust = kappa_dust_T[0];
  else if(T_dust < 275.)
    kappa_dust = kappa_dust_T[1];
  else if(T_dust < 425.)
    kappa_dust = kappa_dust_T[2];
  else if(T_dust < 680.)
    kappa_dust = kappa_dust_T[3];
  else if(T_dust < 1230.)
    kappa_dust = kappa_dust_T[4];

  kappa_dust *= metallicity;
#endif

  if(FlagPrint)
    printf("kappa_gas=%g %g %g %g kappa_dust=%g beta_LW=%g, f_HM=%g f_HI=%g, f_HI_2s2p=%g f_elec=%g\n"
	   , kappa_gas[0]/PROTONMASS, kappa_gas[1]/PROTONMASS, kappa_gas[2]/PROTONMASS, kappa_gas[3]/PROTONMASS, kappa_dust
	   , beta_LW, f_HM, f_HI, f_HI_2s2p, f_elec);

  double kappa = kappa_nu/PROTONMASS + kappa_dust;

  return kappa;
}


#endif

 
