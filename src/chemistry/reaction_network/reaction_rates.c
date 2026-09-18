/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"
#include "chemistry_diagnostics.h"
#include "chemistry_models.h"
#include "chemistry_profile.h"
#include "chemistry_types.h"
#include "cooling_func.h"
#include "dust_cooling.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "reaction_rates.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

static void evaluate_chemical_energy_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *work)
{
  const double nH = XH * environment->rho_cgs / PROTONMASS;
  const double temperature = work->thermodynamics.temperature;

  const double *a = state->abund;
  const double n_elec = nH * a[sELEC] / XH;
  const double n_HI = nH * a[sHI] / XH;
  const double n_HM = nH * a[sHM] / XH;
  const double n_H2I = nH * a[sH2I] / XH;
#ifdef H2_PLUS
  const double n_H2II = nH * a[sH2II] / XH;
#endif
#ifdef HELIUM
  const double n_HeI = nH * a[sHeI] / XH;
  const double n_HeII = nH * a[sHeII] / XH;
#endif

  const double n_crit = (1.e6 / sqrt(temperature)) /
      (1.6 * a[sHI] / XH * exp(-pow(400.0 / temperature, 2.0)) +
       1.4 * a[sH2I] / XH * exp(-12000.0 / (temperature + 1200.0)));

  double edot_H2form = eV_to_erg * work->rates.forward[8] * n_HM * n_HI *
                       3.73 / (1.0 + n_crit / n_HI);
#ifdef H2_PLUS
  edot_H2form += eV_to_erg * work->rates.forward[10] * n_H2II * n_HI *
                 1.83 / (1.0 + n_crit / n_HI);
#endif
#ifdef H2_FORMATION_ON_DUST
  edot_H2form += eV_to_erg * work->rates.forward[51] * n_HI *
                 (0.2 + 4.28 / (1.0 + n_crit / n_HI));
#endif
#ifdef THREEBODY
  edot_H2form += 4.48 * eV_to_erg *
      (((work->rates.forward[55] * n_HI + work->rates.forward[57] * n_H2I) * n_HI * n_HI) -
       ((work->rates.forward[11] * n_HI + work->rates.forward[58] * n_H2I) * n_H2I));
#endif

  double energy_rate = edot_H2form -
                       13.598 * eV_to_erg * work->rates.forward[1] * n_HI * n_elec;
#ifdef HELIUM
  energy_rate -= 24.587 * eV_to_erg * work->rates.forward[3] * n_HeI * n_elec;
  energy_rate -= 54.416 * eV_to_erg * work->rates.forward[5] * n_HeII * n_elec;
#endif

  work->evaluation.rhoU_dot_chem = energy_rate;
  if(fabs(energy_rate) < 1.e-40) {
    work->evaluation.t_cool_chem = MAX_REAL_NUMBER;
  } else {
    work->evaluation.t_cool_chem = state->utherm_cgs * environment->rho_cgs / energy_rate;
  }
}

static double chemistry_species_timescale(double abundance, double production, double destruction)
{
  const double derivative = production - destruction * abundance;
  if(derivative == 0.0)
    return MAX_REAL_NUMBER;
  return fabs(abundance / derivative);
}

int evaluate_chemistry_timescale(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *workspace)
{
  const double n_nuclei = environment->rho_cgs / PROTONMASS;
  const double nH = XH * n_nuclei;
  double abundance[N_SP];
  double production[N_SP];
  double destruction[N_SP];

  for(int species = 0; species < N_SP; species++)
    abundance[species] = state->abund[species] / XH;

  chemistry_network_production_destruction(environment->run->network, nH,
      abundance, workspace->rates.forward, workspace->rates.inverse,
      production, destruction);

  double timescale = dmin(chemistry_species_timescale(abundance[sELEC], production[sELEC], destruction[sELEC]),
                          chemistry_species_timescale(abundance[sHI], production[sHI], destruction[sHI]));
  if(n_nuclei > 1.0 || n_nuclei > N_THREEBODY)
    timescale = dmin(timescale, chemistry_species_timescale(abundance[sH2I], production[sH2I], destruction[sH2I]));

  workspace->evaluation.t_chem = timescale;
  evaluate_chemical_energy_timescale(environment, state, workspace);
  return 0;
}

int evaluate_reaction_rates(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *workspace){
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  workspace->rates.network = environment->run->network;
  const ChemistryRateTables *tables = environment->run->rate_tables;
  const double rho_cgs = environment->rho_cgs;
  const double n_nuclei = rho_cgs / PROTONMASS;
  const double nH = XH * n_nuclei;

  const double T_gas = workspace->thermodynamics.temperature;
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_cgs, T_gas);
#else
  const double l_core = step->jeans_length;
#endif
#ifdef METALLICITY
  const double metallicity = environment->metallicity;
#endif

  //update_radiation_state(sph_in, sph_out, sph_state, T_gas);

  /*********** initialize ************/
  for(int j=0; j<CHEM_RATE_COUNT; j++)
    workspace->rates.forward[j] = 0.0;
  /***********************************/

  /* ------- Compute temperature index. */
  const double ln_T0 = log(T_start);
  const double ln_T9 = log(T_end);

  const double log_T  = log10(T_gas);
  const double log_T2 = log_T * log_T;
  const double log_T3 = log_T * log_T * log_T;

  double ln_T  = log(T_gas);
  if(ln_T < ln_T0) ln_T = ln_T0;
  if(ln_T > ln_T9) ln_T = ln_T9;
  int i_T = (int) ((ln_T-ln_T0)/tables->dln_temperature);

  if(i_T<0) i_T=0;
  if(i_T>N_T-2) i_T=N_T-2;

  if (i_T > N_T-2) i_T = N_T-2;  /* [i_T+1] appears below */
  if (i_T < 0    ) i_T = 0;

  const double c_T = (T_gas - tables->temperature[i_T]) /
      (tables->temperature[i_T + 1] - tables->temperature[i_T]);

  /* Every pure forward rate registered in forward_rate_table is tabulated.
     Density-, radiation-, and dust-dependent rates are overwritten below. */
  for(int j = 0; j < CHEM_RATE_COUNT; j++) {
    if(forward_rate_table[j] != NULL)
      workspace->rates.forward[j] = tables->forward[j][i_T] +
          c_T * (tables->forward[j][i_T + 1] - tables->forward[j][i_T]);
  }

  workspace->rates.forward[22] = workspace->rates.forward[23] = workspace->rates.forward[24] = 0.;

#ifdef IONIZATION_PHOTON
  workspace->rates.forward[22] = step->k22;
#ifdef HELIUM
  workspace->rates.forward[23] = step->photoionization_rate_HeI;
  workspace->rates.forward[24] = step->photoionization_rate_HeII;
#endif
#endif

#ifdef RYTR_CARBON_OXYGEN_IONIZATION
  for(int stage = 0; stage < 3; stage++) {
    workspace->rates.forward[153 + stage] =
        step->photoionization_rate_C[stage];
    workspace->rates.forward[156 + stage] =
        step->photoionization_rate_O[stage];
  }
#endif

#ifdef RYTR_NCO
  workspace->rates.forward[171] = environment->rate_CO_diss;
#endif

#ifdef LYMANWERNER    /* with Photo-Dissociation */
  workspace->rates.forward[25] = 0.;
  workspace->rates.forward[29] = 0.;
  workspace->rates.forward[53] = 0.;

  double J21 = 0.;
#ifdef LW_RADIATION_LOCAL
  workspace->rates.forward[25] = environment->rate_Hminus_pd;
  workspace->rates.forward[29] = environment->rate_H2_diss;
#ifdef RYTR_NHD
  workspace->rates.forward[53] = environment->rate_HD_diss;
#else
  workspace->rates.forward[53] = 1.55/ 1.39 * workspace->rates.forward[29];      /* optically-thin limit found in Wolcott-Green & Haiman 2011 */
#endif
#endif // LW_RADIATION_LOCAL

#ifdef BACKGROUND_LW_RADIATION
  workspace->rates.forward[25] += environment->global->background_rate_Hminus_pd;

  double k29_ext = environment->global->background_rate_H2_diss;
#ifdef HD
  double k53_ext = environment->global->background_rate_HD_diss;
#endif
#ifdef SHIELD_LW
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double column_H2 = l_core * n_nuclei * state->abund[sH2I];
#ifdef HD
  const double column_HI = l_core * n_nuclei * state->abund[sHI];
  const double column_HD = l_core * n_nuclei * state->abund[sHDI];
  const double shield_HI = EachShield(column_HI, 1.49e-1, 1.62e+0, 2.85e23);
  const double shield_H2 = EachShield(column_H2, 5.20e-3, 2.38e-1, 2.34e19);
#endif
  const double shield_H2self = SelfShield(column_H2, T_gas,
                                          9.12 * sqrt(T_gas / 1.e4));
#ifdef HD
  const double shield_HDself = SelfShield(column_HD, T_gas,
                                          7.45 * sqrt(T_gas / 1.e4));
#endif
#else
  const double shield_H2self = step->shield_H2self;
#ifdef HD
  const double shield_HI = step->shield_HI;
  const double shield_H2 = step->shield_H2;
  const double shield_HDself = step->shield_HDself;
#endif
#endif
  const double f_print_H2self = shield_H2self;
#else
  const double f_print_H2self = 1.0;
#endif
#ifdef SHIELD_LW
  k29_ext *= shield_H2self /* SphP[ithis].f_byHI*/;
#ifdef HD
  k53_ext *= shield_HDself * shield_HI * shield_H2;
#endif
#endif //SHIELD_LW

  workspace->rates.forward[29] += k29_ext;
#ifdef HD
  workspace->rates.forward[53] += k53_ext;
#endif

  if(workspace->diagnostics_enabled)
    printf("k25=%g %g k29=%g %g f=%g\n", workspace->rates.forward[25], environment->global->background_rate_Hminus_pd, workspace->rates.forward[29], k29_ext, f_print_H2self);

#endif //BACKGROUND_LW_RADIATION
#endif // LYMANWERNER

  /****************************************************************************************************************/
  /******************************** We evaluate continuum opacities here ******************************************/
  /****************************************************************************************************************/

#ifdef DUST_COOLING
  double edot_dust = 0.;
#endif

  double tau_cont=0.;
  workspace->continuum.tau = 0.;
  workspace->continuum.reduction_factor = 1.;
#ifdef CONT_OPACITY
  struct opacities opac;

  double dummy = calc_continuum_opacities(
      environment->run->continuum_opacity_tables, T_gas, rho_cgs, l_core,
      &opac);
  tau_cont = opac.tau;

#ifdef DUST_COOLING
  double stellar_flux = 0.;
#ifdef THERMAL_FEEDBACK
  stellar_flux = environment->stellar_radiation_flux;
#endif
  workspace->continuum.dust_temperature = dust_temperature(
      environment->run->dust_tables, rho_cgs, T_gas,
      environment->global->cmb_temperature, metallicity, l_core,
      stellar_flux, &opac, &edot_dust, workspace->diagnostics_enabled);
  workspace->continuum.dust_cooling_rate = edot_dust;

#ifdef RYTR_TAU_DUST
  workspace->continuum.dust_kappa_planck = opac.kappa_dust_Planck;
  if(workspace->continuum.dust_kappa_planck > 2.06 * metallicity)
    workspace->continuum.dust_kappa_planck = 2.06 * metallicity;
#endif

  tau_cont = opac.tau;
#endif

  workspace->continuum.tau = tau_cont;
  workspace->continuum.reduction_factor = 1./(1. + 3./4.*tau_cont*tau_cont);
#endif
  /****************************************************************************************************************/
  /****************************************************************************************************************/
  /****************************************************************************************************************/

#ifdef INVERSE_REACTIONS
  const double (*k_table_inverse)[N_T] = tables->inverse;
  double n_crit = H_criticaldensity(nH, T_gas, nH*state->abund[sHI]*l_core, state->abund[sHI], state->abund[sELEC]);

  workspace->rates.inverse[1] = k_table_inverse[1][i_T] + c_T*(k_table_inverse[1][i_T+1] - k_table_inverse[1][i_T]);
  workspace->rates.inverse[2] = k_table_inverse[2][i_T] + c_T*(k_table_inverse[2][i_T+1] - k_table_inverse[2][i_T]);
  workspace->rates.inverse[7] = k_table_inverse[7][i_T] + c_T*(k_table_inverse[7][i_T+1] - k_table_inverse[7][i_T]);
  workspace->rates.inverse[8] = k_table_inverse[8][i_T] + c_T*(k_table_inverse[8][i_T+1] - k_table_inverse[8][i_T]);
  workspace->rates.inverse[13] = k_table_inverse[13][i_T] + c_T*(k_table_inverse[13][i_T+1] - k_table_inverse[13][i_T]);
  workspace->rates.inverse[16] = k_table_inverse[16][i_T] + c_T*(k_table_inverse[16][i_T+1] - k_table_inverse[16][i_T]);
  workspace->rates.inverse[17] = k_table_inverse[17][i_T] + c_T*(k_table_inverse[17][i_T+1] - k_table_inverse[17][i_T]);
  workspace->rates.inverse[18] = k_table_inverse[18][i_T] + c_T*(k_table_inverse[18][i_T+1] - k_table_inverse[18][i_T]);

  workspace->rates.inverse[30] = k_table_inverse[30][i_T] + c_T*(k_table_inverse[30][i_T+1] - k_table_inverse[30][i_T]);
  workspace->rates.inverse[55] = k_table_inverse[55][i_T] + c_T*(k_table_inverse[55][i_T+1] - k_table_inverse[55][i_T]);

  //workspace->rates.inverse[1] *= nH / (n_crit + nH); //temporally
  workspace->rates.inverse[18] *= nH / (n_crit + nH);

  //printf("tau=%g n_crit=%g\n", tau_cont, n_crit);

  workspace->rates.inverse[2] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  workspace->rates.inverse[7] *= (1 - exp(-tau_cont)) * n_nuclei / (n_crit + n_nuclei);
  //workspace->rates.inverse[2] *= (1 - exp(-tau_cont));
  //workspace->rates.inverse[7] *= (1 - exp(-tau_cont));
  workspace->rates.forward[25]        *= exp(-tau_cont);

  if(workspace->diagnostics_enabled){
    printf("k2=%g k_inv=%g %g %g %g %g %g %g %g %g %g tau_cont=%g n=%g n_c=%g l_core=%g\n"
	   , workspace->rates.forward[2]
	   , workspace->rates.inverse[1], workspace->rates.inverse[2], workspace->rates.inverse[7], workspace->rates.inverse[8]
	   , workspace->rates.inverse[13], workspace->rates.inverse[16], workspace->rates.inverse[17], workspace->rates.inverse[18]
	   , workspace->rates.inverse[30], workspace->rates.inverse[55], tau_cont, n_nuclei, n_crit, l_core);
  }


#ifdef ID_CHEM_DEBUG
  if(workspace->diagnostics_enabled){
    printf("k_inv=%g %g tau_cont=%g n=%g n_c=%g l_core=%g i_T=%d c_T=%g\n"
	   , workspace->rates.inverse[2], workspace->rates.inverse[7], tau_cont, n_nuclei, n_crit, l_core, i_T, c_T);
  }
#endif
#endif

#ifdef H2_FORMATION_ON_DUST
  workspace->rates.forward[51] = H2form_on_grains(rho_cgs, T_gas, workspace->continuum.dust_temperature, metallicity);
#endif

  if(T_gas > 10.){
#ifdef MARTIN97
    double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21;
    double p, n_c1, n_c2, log_g_h1, log_g_l1, log_g_h2, log_g_l2, log_n_c1, log_n_c2;
    double log_gamma_cd, log_gamma_dt;

    a1  = -1.784239e2;
    a2  = -6.842243e1;
    a3  =  4.320243e1;
    a4  = -4.633167;
    a5  =  6.970086e1;
    a6  =  4.087038e4;
    a7  = -2.370570e4;
    a8  =  1.288953e2;
    a9  = -5.391334e1;
    a10 =  5.315517;
    a11 = -1.973427e1;
    a12 =  1.678905e4;
    a13 = -2.578611e4;
    a14 =  1.482123e1;
    a15 = -4.890915;
    a16 =  4.749030e-1;
    a17 = -1.338283e2;
    a18 = -1.164408;
    a19 =  8.227443e-1;
    a20 =  5.864073e-1;
    a21 = -2.056313;


    log_g_h1 = a1 + a2 * log_T + a3 * log_T2 + a4 * log_T3 + a5 * log10(1 + a6/T_gas);
    log_g_h2 = a7/T_gas;

    log_g_l1 = a8 + a9 * log_T + a10 * log_T2 + a11 * log10(1 + a12/T_gas);
    log_g_l2 = a13/T_gas;

    log_n_c1 = a14 + a15*log_T + a16*log_T2 + a17/T_gas;
    n_c1 = pow(10.0, log_n_c1);

    log_n_c2 = a18 + log_n_c1;
    n_c2 = pow(10.0, log_n_c2);

    p = a19 + a20 *exp(-T_gas/1850.0) + a21 * exp(-T_gas/440.0);

    log_gamma_cd = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+pow(nH/n_c1, p)) + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+pow(nH/n_c2, p));

    a1  = -1.427664e2;
    a2  =  4.270741e1;
    a3  = -2.027365;
    a4  = -2.582097e-1;
    a5  =  2.136094e1;
    a6  =  2.753531e4;
    a7  = -2.146779e4;
    a8  =  6.034928e1;
    a9  = -2.743096e1;
    a10  = 2.676150;
    a11 = -1.128215e1;
    a12 =  1.425455e4;
    a13 = -2.312520e4;
    a14 =  9.305564;
    a15 = -2.464009;
    a16 =  1.985955e-1;
    a17 =  7.430600e2;
    a18 = -1.174242;
    a19 =  7.502286e-1;
    a20 =  2.358848e-1;
    a21 =  2.937507;


    log_g_h1 = a1 + a2*log_T + a3*log_T2 + a4*log_T3 + a5*log10(1+a6/T_gas);
    log_g_h2 = a7/T_gas;

    log_g_l1 = a8 + a9*log_T + a10*log_T2 + a11*log10(1+a12/T_gas);
    log_g_l2 = a13/T_gas;

    log_n_c1 = a14 + a15*log_T + a16*log_T2 + a17/T_gas;
    n_c1 = pow(10.0, log_n_c1);

    log_n_c2 = a18 + log_n_c1;
    n_c2 = pow(10.0, log_n_c2);

    p = a19 + a20*exp(-T_gas/1850.0) + a21*exp(-T_gas/440.0);

    log_gamma_dt = log_g_h1 - (log_g_h1 - log_g_l1)/(1.0+pow(nH/n_c1, p)) + log_g_h2 - (log_g_h2 - log_g_l2)/(1.0+pow(nH/n_c2, p));

    workspace->rates.forward[11] = pow(10.0, log_gamma_cd) + pow(10.0, log_gamma_dt);

    if(workspace->diagnostics_enabled)
      printf("k11=%g gamma=%g %g\n", workspace->rates.forward[11], log_gamma_cd, log_gamma_dt);


#ifdef INVERSE_REACTIONS
    double xk_L = workspace->rates.forward[11];
    double xk_H  = workspace->rates.inverse[55];

    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);

    double a = 1.0/(1.0 + nH/rn_cr);

    if(a == 1.0)       workspace->rates.forward[11] = xk_L;
    else if(a == 0.0)  workspace->rates.forward[11] = xk_H;
    else               workspace->rates.forward[11] = pow(xk_H,1.0-a)*pow(xk_L,a);
#endif
#else
    double xk_L  = 1.12e-10*exp(-7.035e+4/T_gas);
    double xk_H  = 6.5e-7/sqrt(T_gas)*exp(-5.2e+4/T_gas)*(1.0-exp(-6.e+3/T_gas));
    double rlgT4 = log10(T_gas*1.e-4);
    double rn_cr = pow(10.0, 4.0 - 0.416*rlgT4 - 0.327*rlgT4*rlgT4);

    double a = 1.0/(1.0 + nH/rn_cr);

    if(a == 1.0)       workspace->rates.forward[11] = xk_L;
    else if(a == 0.0)  workspace->rates.forward[11] = xk_H;
    else               workspace->rates.forward[11] = pow(xk_H,1.0-a)*pow(xk_L,a);
#endif
  }else
    workspace->rates.forward[11] = 0.;


#if defined(STROMGREN_TEST)
  for(int j=1; j<=21; j++){
#ifdef IONIZATION_PHOTON
    if(j == 2)  continue;
#endif
    workspace->rates.forward[j] = 0.;
  }

  workspace->rates.forward[30] = 0.;
#endif
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_REACTION_RATES, profile_timer);
#endif
  return 0;
}

#endif /* CHEMISTRY_SC */
