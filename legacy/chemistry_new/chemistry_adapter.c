#include <math.h>

#include "allvars.h"
#include "standalone_constants.h"
#include "chemistry_adapter.h"
#include "chemistry_types.h"

#ifndef CHEMISTRY_OLD
double compute_gamma_from_abund(double *abund, double T_gas);
double get_core_radius(double rho_cgs, double T_gas);
double calc_continuum_opacities(double T_gas, double rho_cgs, double l_core, struct opacities *opac);
#if defined(DUST_COOLING) && defined(METALLICITY)
double dust_temperature(double rho_cgs, double T_gas, double metallicity, double l_Jeans,
  double stellar_flux, struct opacities *opac, double *edot_dust);
#endif

static double particle_n_nuclei_1zone(int i)
{
  return All.UnitDensity_in_cgs * SphP[i].Density * All.HubbleParam * All.HubbleParam *
    All.cf_a3inv / PROTONMASS;
}

static double particle_utherm_from_entropy_1zone(int i)
{
  return SphP[i].Entropy / (SphP[i].Gamma - 1.0) *
    pow(SphP[i].Density * All.cf_a3inv, SphP[i].Gamma - 1.0);
}

static double abundances_molecular_weight_1zone(const double *abund)
{
  double sum_fraction = abund[sELEC] + abund[sHI] + abund[sHII] + abund[sHM] + abund[sH2I]
#ifdef H2_PLUS
    + abund[sH2II]
#endif
#ifdef HD
    + abund[sDI] + abund[sDII] + abund[sHDI]
#ifdef HD_PLUS
    + abund[sHDII]
#endif
#ifdef DMINUS
    + abund[sDM]
#endif
#endif
#ifdef HELIUM
    + abund[sHeI] + abund[sHeII] + abund[sHeIII]
#else
    + XHe
#endif
    ;

  return 1.0 / sum_fraction;
}

static void update_particle_derived_1zone(int i, const double *abund)
{
  SphP[i].n_nuclei = particle_n_nuclei_1zone(i);
  SphP[i].molecular_weight = abundances_molecular_weight_1zone(abund);
}

static void update_particle_opacity_1zone(int i, const double *abund)
{
#ifdef CONT_OPACITY
  double rho_cgs = SphP[i].n_nuclei * PROTONMASS;
  double T_gas = (All.UnitEnergy_in_cgs / All.UnitMass_in_g * SphP[i].ThermalEnergy) *
    (SphP[i].Gamma - 1.0) / BOLTZMANN * PROTONMASS * abundances_molecular_weight_1zone(abund);
  double l_core = get_core_radius(rho_cgs, T_gas);
  SphP[i].Jeans_Length = l_core;
  struct opacities opac = {0};

  calc_continuum_opacities(T_gas, rho_cgs, l_core, &opac);
#if defined(DUST_COOLING) && defined(METALLICITY)
  double edot_dust = 0.0;
  SphP[i].Tdust = dust_temperature(rho_cgs, T_gas, SphP[i].Metallicity, l_core, 0.0,
    &opac, &edot_dust);
  SphP[i].edot_dust = edot_dust;
#endif
  SphP[i].tau_cont = opac.tau;
  SphP[i].f_reduce = 1.0 / (1.0 + 0.75 * opac.tau * opac.tau);
#endif
}

void particle_to_abundances_1zone(int i, double *abund)
{
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    abund[j] = 0.0;

  abund[sELEC] = SphP[i].elec;
  abund[sHI] = SphP[i].HI;
  abund[sHII] = SphP[i].HII;
  abund[sH2I] = SphP[i].H2I;
  abund[sHM] = SphP[i].HM;
#ifdef HELIUM
  abund[sHeI] = SphP[i].HeI;
  abund[sHeII] = SphP[i].HeII;
  abund[sHeIII] = SphP[i].HeIII;
#endif
#ifdef H2_PLUS
  abund[sH2II] = SphP[i].H2II;
#endif
#ifdef HD
  abund[sDI] = SphP[i].DI;
  abund[sDII] = SphP[i].DII;
  abund[sHDI] = SphP[i].HDI;
#ifdef HD_PLUS
  abund[sHDII] = SphP[i].HDII;
#endif
#ifdef DMINUS
  abund[sDM] = SphP[i].DM;
#endif
#endif
}

void abundances_to_particle_1zone(int i, const double *abund)
{
  SphP[i].elec = abund[sELEC];
  SphP[i].HI = abund[sHI];
  SphP[i].HII = abund[sHII];
  SphP[i].H2I = abund[sH2I];
  SphP[i].HM = abund[sHM];
#ifdef HELIUM
  SphP[i].HeI = abund[sHeI];
  SphP[i].HeII = abund[sHeII];
  SphP[i].HeIII = abund[sHeIII];
#endif
#ifdef H2_PLUS
  SphP[i].H2II = abund[sH2II];
#endif
#ifdef HD
  SphP[i].DI = abund[sDI];
  SphP[i].DII = abund[sDII];
  SphP[i].HDI = abund[sHDI];
#ifdef HD_PLUS
  SphP[i].HDII = abund[sHDII];
#endif
#ifdef DMINUS
  SphP[i].DM = abund[sDM];
#endif
#endif
}

void particle2in_chemdata_1zone(int i, struct chemdata_in *in, struct chemdata_out *out)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  particle_to_abundances_1zone(i, abund);
  update_particle_derived_1zone(i, abund);

  in->ID = P[i].ID;
  in->flag_print = 0;
#ifdef ID_CHEM_DEBUG
  if(P[i].ID == ID_CHEM_DEBUG)
    in->flag_print = 1;
#endif
#ifdef FORCE_CHEM_PRINT
  in->flag_print = 1;
#endif
  in->Density = SphP[i].Density;
  in->rho_cgs = SphP[i].n_nuclei * PROTONMASS;
#ifdef METALLICITY
  in->Metallicity = SphP[i].Metallicity;
#endif
#if defined(LYMANWERNER) && defined(LW_RADIATION_LOCAL)
  in->rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  in->rate_H2_diss = SphP[i].rate_H2_diss;
  in->J21 = SphP[i].J21;
#ifdef DUST_COOLING
  in->G0 = SphP[i].G0;
#endif
#endif
#ifdef IONIZING_PHOTON
  in->k22 = SphP[i].k22;
  in->piHI = SphP[i].piHI;
#endif

  SphP[i].ThermalEnergy = particle_utherm_from_entropy_1zone(i);
  out->Utherm = SphP[i].ThermalEnergy;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    out->Abund[j] = abund[j];
  out->Gamma = SphP[i].Gamma;
  out->t_chem = SphP[i].t_chem;
  out->t_cool = SphP[i].t_cool;
  out->Tdust = SphP[i].Tdust;
  out->kappa_dust = 0.0;
  out->flag_implicit = 0;
  out->conv_flag = 0;
}

void out2particle_chemdata_1zone(int i, const struct chemdata_out *out)
{
  SphP[i].ThermalEnergy = out->Utherm;
  abundances_to_particle_1zone(i, out->Abund);
  update_particle_derived_1zone(i, out->Abund);
  SphP[i].Gamma = out->Gamma;
  SphP[i].t_chem = out->t_chem;
  SphP[i].t_cool = out->t_cool;
  SphP[i].Tdust = out->Tdust;
  update_particle_opacity_1zone(i, out->Abund);
  SphP[i].Entropy = SphP[i].ThermalEnergy * (SphP[i].Gamma - 1.0) /
    pow(SphP[i].Density * All.cf_a3inv, SphP[i].Gamma - 1.0);
}

int compute_gamma(int ithis)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  particle_to_abundances_1zone(ithis, abund);
  SphP[ithis].Gamma = compute_gamma_from_abund(abund, SphP[ithis].temperature);
  return 1;
}

double compute_temp(int ithis)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  particle_to_abundances_1zone(ithis, abund);
  double molecular_weight = abundances_molecular_weight_1zone(abund);
  double egyspec_to_temp = (SphP[ithis].Gamma - 1.0) / BOLTZMANN * PROTONMASS * molecular_weight;
  double utherm = particle_utherm_from_entropy_1zone(ithis);

  return (All.UnitEnergy_in_cgs / All.UnitMass_in_g * utherm) * egyspec_to_temp;
}
#endif
