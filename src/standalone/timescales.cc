#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "allvars.h"
#include "standalone_constants.h"
#include "standalone_math.h"
#ifndef CHEMISTRY_OLD
#include "thermodynamics.h"
#endif

double compute_freefall_time(double rho_gas);
int compute_molecular_weight(int ithis);
int compute_Jeans_length(int ithis, double T_gas);
double compute_effective_gamma(double rho, double e, double Gamma, double L_rad, double dt, double *gamma_eff0);

double timestep_flag(double dt_phys_pre, int flag_chem){
  double dt_cool  = fabs(SphP[0].t_cool);       /* in seconds */
  dt_cool *= All.Epsilon;

  double dt_phys=dt_cool;
  /* another criterion given by the electron number density change */

  if(flag_chem){
    double dt_chem  = fabs(SphP[0].t_chem);       /* in seconds */
    dt_chem *= All.Epsilon;

    if(dt_chem > 0 && dt_chem < dt_phys)
      dt_phys=dt_chem;
  }

  if(dt_phys > 1.1 * dt_phys_pre && dt_phys_pre !=0)
    dt_phys = 1.1 * dt_phys_pre;

  return dt_phys;
}



double timestep(double dt_phys_pre){
  const double Courant = 0.05;

  double dt_col  = fabs(SphP[0].t_col);
  double dt_cool = fabs(SphP[0].t_cool);
  double dt_chem = fabs(SphP[0].t_chem);

#if defined(CHEMISTRY_ENERGY_IMPLICIT) && !defined(LIMIT_DT_IMPLICIT)
  dt_cool = MAX_REAL_NUMBER;
#endif
#ifdef CHEMISTRY_IMPLICIT
  dt_chem = MAX_REAL_NUMBER;
#endif

  double dt_phys = Courant * standalone_min(
      dt_col, standalone_min(dt_cool, dt_chem));

  if(SphP[0].n_nuclei > 0.99* 1.e18)
    dt_phys *= 0.5;
  dt_phys *= All.HubbleParam/All.UnitTime_in_s;

  if(dt_phys > 1.1 * dt_phys_pre)
    dt_phys = 1.1 * dt_phys_pre;

  return dt_phys;
}




int density_update(int mode, int ithis, double theta, double a_now, double a_end, double *da_require){
  double density_old  = SphP[ithis].Density;

  if(mode == 0){
    SphP[ithis].DMdensity=(9.0*M_PI*M_PI/2.0) * All.UnitDensity_in_cgs * RHOCRITICAL * All.Omega0
      * pow((1.0 + StandaloneCollapseRedshift) / (1.0 - cos(theta)), 3.0);

    if(SphP[ithis].DMdensity > StandaloneVirialDarkMatterDensity)
      SphP[ithis].DMdensity = StandaloneVirialDarkMatterDensity;
  }

  double da = a_end - a_now;
  double hubble_a= All.Hubble*sqrt(All.Omega0/(a_now*a_now*a_now)
				   + (1-All.Omega0-All.OmegaLambda)/(a_now*a_now) + All.OmegaLambda);

  /* dt = da/ ( a * H ) */

  double dt = da / (a_now * hubble_a);
  double dt_in_sec = dt * All.UnitTime_in_s/All.HubbleParam;
  double t_freefall = compute_freefall_time(SphP[ithis].Density);
  double density_new = density_old +  (density_old/t_freefall) * dt_in_sec;

  SphP[ithis].Density = density_new;

  double drho = density_new/density_old;
  if(drho > 1.2 || drho < 0.8)
    printf("WARNING TOO LARGE TIMESTEP! (%g)",da);

  /* timestep constraint */
  dt = All.Epsilon * t_freefall;

  dt *= All.HubbleParam;
  dt /= All.UnitTime_in_s;

  *da_require = dt * a_now * hubble_a;

  return ithis;
}



int density_update_phys(int ithis, double dt_in_sec, double *dt_require){
  if(SphP[ithis].DMdensity > StandaloneVirialDarkMatterDensity)
    SphP[ithis].DMdensity = StandaloneVirialDarkMatterDensity;

  double density_old  = SphP[ithis].Density;

  double t_freefall = compute_freefall_time(SphP[ithis].Density);

  double density_new = density_old +  (density_old/t_freefall) * dt_in_sec;

  SphP[ithis].Density = density_new;

  double drho = density_new/density_old;
  if(drho > 1.2 || drho < 0.8)
    printf("WARNING TOO LARGE TIMESTEP!");

  /* timestep constraint */
  *dt_require = All.Epsilon * t_freefall;

  return ithis;
}



int density_update_noncosmo(int mode, int ithis, double theta, double a_now, double a_end, double *da_require){
  double density_now  = SphP[ithis].Density;

  double da = a_end - a_now;
  double dt_in_sec = da;

  double t_freefall =  compute_freefall_time(SphP[ithis].Density);

  double density_new = density_now +  (density_now/t_freefall) * dt_in_sec;

  SphP[ithis].Density = density_new;


  double drho = density_new/density_now;
  if(drho > 1.3 || drho < 0.8)
    printf("WARNING TOO LARGE TIMESTEP!");

  /* timestep constraint */
  double dt = 0.02 * t_freefall;
  *da_require = dt;

  return ithis;
}


int density_isobaric(int ithis, double rho_init, double T_init, double *da_require){
  double density_now  = SphP[ithis].Density;
  double density_new = rho_init * T_init/SphP[ithis].temperature;

  SphP[ithis].Density = density_new;
  SphP[ithis].temp_hold = SphP[ithis].temperature;

  double drho = density_new/density_now;
  if(drho > 1.3 || drho < 0.8)
    printf("WARNING TOO LARGE TIMESTEP!");

  *da_require = 1.0e20;

  return ithis;
}


double compute_collapse_time(double SHR, double rho_gas){
  double t_ff = sqrt(3.0*M_PI/32.0/GRAVITY/rho_gas); //unit in sec
  double factor = 0.;
  if(rho_gas > 1.0e2*PROTONMASS/XH) {
    if(SHR<0.83)
      factor = 0.0;
    else if(SHR<1.0)
      factor = 0.6 + 2.5*(SHR-1.0) - 6.0*pow(SHR-1.0,2.0);
    else if(SHR<4.0/3.0)
      factor = 1.0 + 0.2*(SHR-4.0/3.0) - 2.9*pow(SHR-4.0/3.0,2.0);
    else
      factor = 1.0;

    if(factor>0.95) factor = 0.95;
  }

  return t_ff / sqrt(1.0-factor);
}

double compute_collapse_factor(double SHR){
  double factor = 1.;
  if(SHR<0.83)
    factor = 0.0;
  else if(SHR<1.0)
    factor = 0.6 + 2.5*(SHR-1.0) - 6.0*pow(SHR-1.0,2.0);
  else if(SHR<4.0/3.0)
    factor = 1.0 + 0.2*(SHR-4.0/3.0) - 2.9*pow(SHR-4.0/3.0,2.0);

  if(factor>0.95) factor = 0.95;

  return factor;
}


double update_phys(int ithis,double dt,double SHR0){
  double SHR = SHR0;
  double dt_in_gadget = dt/All.UnitTime_in_s * All.HubbleParam;   //dt must be unit in sec
  int ifunc = compute_molecular_weight(ithis);

  double Gamma  = SphP[ithis].Gamma;
  double mu     = SphP[ithis].molecular_weight;
  double rho0   = SphP[ithis].Density  * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs; //in cgs
  double p0     = SphP[ithis].pressure                                     * All.UnitMass_in_g     / All.UnitEnergy_in_cgs;
  double e0     = SphP[ithis].ThermalEnergy                                * All.UnitEnergy_in_cgs / All.UnitMass_in_g; // in cgs
  double ent0   = SphP[ithis].Entropy;
  double t_col  = compute_freefall_time(rho0); // in sec (if rho0 is in cgs)

#ifdef SLOW_COLLAPSE
  t_col /= SLOW_COLLAPSE;
#endif

#ifdef SHR_ITER
  double gamma_eff;
  double L_net = (SphP[ithis].edot_rad + SphP[ithis].edot_chem)/rho0; //in cgs
  if(rho0 / PROTONMASS > 100.)
    t_col = compute_effective_gamma(rho0, e0, Gamma, L_net, dt, &gamma_eff);

  SphP[ithis].SHR = gamma_eff;
  SHR = gamma_eff;

  printf("gamma_eff=SHR=%g\n", SHR);
#endif

  double G_comp = (Gamma - 1.0) * e0/t_col;
  double e = e0 + G_comp * dt;

  double drho = rho0/t_col*dt;
  double rho = rho0 + drho;
  double T = e *  (Gamma-1.0) * (mu*PROTONMASS)/BOLTZMANN;
  double p = rho * BOLTZMANN * T / (mu*PROTONMASS);

  //printf("%g %g %g %g %g %g\n",drho,SphP[ithis].edot_rad,SphP[ithis].edot_chem[0],mu,Gamma,t_col);

  SphP[ithis].edot_adia = G_comp*rho;
  SphP[ithis].t_col     = t_col;
  SphP[ithis].Density   = rho  / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
  SphP[ithis].DtDensity = drho / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs) / dt_in_gadget;
  SphP[ithis].pressure  = p    / (All.UnitMass_in_g/All.UnitEnergy_in_cgs);
  //SphP[ithis].ThermalEnergy = e / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);
  SphP[ithis].temperature = T;
  SphP[ithis].Entropy     = e / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) * (SphP[ithis].Gamma-1.0) / pow(SphP[ithis].Density, SphP[ithis].Gamma-1.0);
  SphP[ithis].DtEntropy   = (SphP[ithis].Entropy - ent0)/dt_in_gadget;

  ifunc = compute_Jeans_length(ithis, T);

  double column_dens_HI  = SphP[ithis].Jeans_Length * SphP[ithis].HI  * SphP[ithis].n_nuclei ;
#ifdef CHEMISTRY_SPECIES_HELIUM
  double column_dens_HeI = SphP[ithis].Jeans_Length * SphP[ithis].HeI * SphP[ithis].n_nuclei ;
#endif




  return SHR;
}

int compute_molecular_weight(int ithis){
#ifdef CHEMISTRY_OLD
  double sum_fraction = SphP[ithis].elec + SphP[ithis].HI +
    SphP[ithis].HII + SphP[ithis].HM + SphP[ithis].H2I
#ifdef CHEMISTRY_SPECIES_H2_PLUS
    + SphP[ithis].H2II
#endif
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    + SphP[ithis].DI + SphP[ithis].DII + SphP[ithis].HDI
#ifdef CHEMISTRY_SPECIES_HD_PLUS
    + SphP[ithis].HDII
#endif
#ifdef CHEMISTRY_SPECIES_D_MINUS
    + SphP[ithis].DM
#endif
#endif
#ifdef CHEMISTRY_SPECIES_HELIUM
    + SphP[ithis].HeI + SphP[ithis].HeII + SphP[ithis].HeIII
#else
    + XHe
#endif
    ;
  SphP[ithis].molecular_weight = 1.0 / sum_fraction;
#else
  double abund[TGCHEM_NUM_ABUNDANCES] = {0};
  abund[sELEC] = SphP[ithis].elec;
  abund[sHI] = SphP[ithis].HI;
  abund[sHII] = SphP[ithis].HII;
  abund[sHM] = SphP[ithis].HM;
  abund[sH2I] = SphP[ithis].H2I;
#ifdef CHEMISTRY_SPECIES_H2_PLUS
  abund[sH2II] = SphP[ithis].H2II;
#endif
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  abund[sDI] = SphP[ithis].DI;
  abund[sDII] = SphP[ithis].DII;
  abund[sHDI] = SphP[ithis].HDI;
#ifdef CHEMISTRY_SPECIES_HD_PLUS
  abund[sHDII] = SphP[ithis].HDII;
#endif
#ifdef CHEMISTRY_SPECIES_D_MINUS
  abund[sDM] = SphP[ithis].DM;
#endif
#endif
#ifdef CHEMISTRY_SPECIES_HELIUM
  abund[sHeI] = SphP[ithis].HeI;
  abund[sHeII] = SphP[ithis].HeII;
  abund[sHeIII] = SphP[ithis].HeIII;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  abund[sCI] = SphP[ithis].CI;
  abund[sCII] = SphP[ithis].CII;
  abund[sCO] = SphP[ithis].CO;
  abund[sCO2] = SphP[ithis].CO2;
  abund[sOI] = SphP[ithis].OI;
  abund[sOH] = SphP[ithis].OH;
  abund[sH2O] = SphP[ithis].H2O;
  abund[sO2] = SphP[ithis].O2;
#endif

  SphP[ithis].molecular_weight = chemistry_molecular_weight(abund);
#endif

  return 1;
}

int compute_Jeans_length(int ithis, double T_gas){
#ifdef CHEMISTRY_JEANS_TEMPERATURE_MAX
  if(T_gas > CHEMISTRY_JEANS_TEMPERATURE_MAX)
    T_gas = CHEMISTRY_JEANS_TEMPERATURE_MAX;
#endif
#ifdef CHEMISTRY_OLD
  SphP[ithis].Jeans_Length =
      0.5 * M_PI * sqrt(BOLTZMANN * T_gas / PROTONMASS) /
      sqrt(GRAVITY * SphP[ithis].Density * All.HubbleParam *
           All.HubbleParam * All.UnitDensity_in_cgs);
#else
  const double rho_cgs =
      SphP[ithis].Density * All.HubbleParam * All.HubbleParam *
      All.UnitDensity_in_cgs;
  SphP[ithis].Jeans_Length =
      chemistry_jeans_radius(rho_cgs, T_gas, 1.0);
#endif
  //factor 0.5 because we define the length scale (radius) as 0.5 Ljeans
  return 1;
}

double compute_effective_gamma(double rho, double e, double Gamma, double L_rad, double dt, double *gamma_eff0){
  double err, err0;
  double gamma_min, gamma_max;
  double t_col_min, t_col_max;
  double resid_min, resid_max;
  double f_new, resid_new, t_col_new;
  double t_col0, gamma_eff;
  double d_logE;

  gamma_min = 0.83;
  gamma_max = 4./3.;

  t_col0 = compute_freefall_time(rho);
  d_logE = L_rad * dt / e;

  gamma_eff = gamma_min;
  t_col_min = t_col0; //unit in second
  resid_min = gamma_eff - Gamma - d_logE * t_col_min / dt;

  printf("gamma_min=%g resid=%g\n", gamma_min, resid_min);

  gamma_eff = gamma_max;
  f_new = compute_collapse_factor(gamma_eff);
  t_col_max = 1./sqrt(1 - f_new) * t_col0; //unit in second
  resid_max = gamma_eff - Gamma - d_logE * t_col_max / dt;

  printf("gamma_max=%g resid=%g\n", gamma_max, resid_max);

  if(resid_min < 0. && resid_max < 0.){
    *gamma_eff0 = gamma_max;
    return t_col_max;
  }

  if(resid_min > 0. && resid_max > 0.){
    *gamma_eff0 = gamma_min;
    return t_col_min;
  }

  if(resid_min > 0. && resid_max < 0.){
    printf("really this happens?\n");
    exit(1);
  }

  /**** OK, effective gamma lies in the range of [0.83:4./3.] ****/
  /********** Let's find the value by bisection method ***********/
  const int itr_max=20;
  for(int itr=0;itr<itr_max;itr++){
    gamma_eff = 0.5 * (gamma_min + gamma_max);

    double f_new = compute_collapse_factor(gamma_eff);
    double t_col_new = 1./sqrt(1 - f_new) * t_col0; //unit in second
    double resid_new = gamma_eff - Gamma - d_logE * t_col_new / dt;

    if(resid_new < 0.){
      gamma_min = gamma_eff;
      t_col_min = t_col_new;
    }else{
      gamma_max = gamma_eff;
      t_col_max = t_col_new;
    }

    printf("gamma_eff=%g resid=%g f=%g\n", gamma_eff, resid_new, f_new);
  }

  *gamma_eff0 = gamma_eff;
  return t_col_new;
}
