/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"

#include "chemistry.h"
#include "chemistry_constants.h"
#include "chemistry_cooling.h"
#include "chemistry_diagnostics.h"
#include "chemistry_models.h"
#include "chemistry_types.h"
#include "cooling_func.h"
#include "dust_cooling.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "chemistry_tables.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

int evaluate_radiative_cooling(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *work)
{
  const double T_gas = work->thermodynamics.temperature;
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas/PROTONMASS;
#ifdef METALLICITY
  const double metallicity = environment->metallicity;
#endif
#ifdef CONT_OPACITY
  const double tau_cont = work->continuum.tau;
  const double continuum_reduction_factor = work->continuum.reduction_factor;
#endif

  const double ln_T0 = log(T_start);
  const double ln_T9 = log(T_end);

  /* --------- Compute temperature index.  */
  double ln_T = log(T_gas);
  if (ln_T < ln_T0) ln_T = ln_T0;
  if (ln_T > ln_T9) ln_T = ln_T9;

  volatile int i_T = (int)( (ln_T-ln_T0)/dln_T ) ;

  if (i_T > N_T-2) i_T = N_T - 2;     /* [i_T+1] appears below */
  if (i_T < 0    ) i_T = 0;

  const double c_T = (T_gas - T_table[i_T])/(T_table[i_T+1]-T_table[i_T]);

  /* --------- Interpolate temperature dependent coefficients. */
#ifdef HLINE_COOLING_THIN
  double ceHI  = ceHIa[i_T]   + c_T*(ceHIa[i_T+1] -ceHIa[i_T]);
#endif

#ifndef CONTINUUM_COOLING_MATSUKOBA
  double reHII = reHIIa[i_T]  + c_T*(reHIIa[i_T+1]-reHIIa[i_T]);
  double brem  = brema[i_T]   + c_T*(brema[i_T+1] -brema[i_T]);
  double cie   = H2_CIE[i_T]  + c_T*(H2_CIE[i_T+1]-H2_CIE[i_T]);
#endif

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  double ceHeI  = ceHeIa[i_T]  + c_T*(ceHeIa[i_T+1] - ceHeIa[i_T]);
  double ceHeII = ceHeIIa[i_T] + c_T*(ceHeIIa[i_T+1]- ceHeIIa[i_T]);

#ifdef HELIUM
  double ciHeI  = ciHeIa[i_T]  + c_T*(ciHeIa[i_T+1]  -ciHeIa[i_T]);
  double ciHeII = ciHeIIa[i_T] + c_T*(ciHeIIa[i_T+1] -ciHeIIa[i_T]);
  double ciHeIS = ciHeISa[i_T] + c_T*(ciHeISa[i_T+1] -ciHeISa[i_T]);

  double reHeII1= reHeII1a[i_T] + c_T*(reHeII1a[i_T+1]-reHeII1a[i_T]);
  double reHeII2= reHeII2a[i_T] + c_T*(reHeII2a[i_T+1]-reHeII2a[i_T]);
  double reHeIII= reHeIIIa[i_T] + c_T*(reHeIIIa[i_T+1]-reHeIIIa[i_T]);
#endif
#endif

  /* ========= Compton heating =========================================================*/
  double comp1 = 5.65e-36*pow(environment->scale_factor,-4);
  double edot_comp = -comp1*(T_gas-environment->cmb_temperature)*state->abund[sELEC]*n_nuclei;

  /**********************************************************************************************/
  /************************ line cooling by H, He, H2, HD, CI, OII ******************************/
  /**********************************************************************************************/

#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_gas, T_gas);
#else
  const double l_core = step->jeans_length;
#endif

  /******************************* Hydrogen lines ***********************************************/
  double edot_ceHI = 0.;
#ifdef HLINE_COOLING_THIN
  edot_ceHI   = -ceHI  *state->abund[sHI]  *state->abund[sELEC]     *n_nuclei*n_nuclei;
#else
  if(n_nuclei < 1.e8){
    double edot_Lya, edot_2ph;
    double v_thermal = sqrt(BOLTZMANN * T_gas/PROTONMASS);
    double fesc = fesc_Lya(n_nuclei, T_gas, state->abund[sHI], state->abund[sELEC], l_core, v_thermal, &edot_Lya, &edot_2ph);

    edot_ceHI = edot_Lya;
  }
#endif

  /******************************* Helium lines **************************************************/
#ifdef HELIUM
  double edot_ceHeI  = -ceHeI *state->abund[sHeII]*state->abund[sELEC]*state->abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
  double edot_ceHeII = -ceHeII*state->abund[sHeII]*state->abund[sELEC]     *n_nuclei*n_nuclei;

  double edot_ciHeI  = -ciHeI *state->abund[sHeI] *state->abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeII = -ciHeII*state->abund[sHeII]*state->abund[sELEC]     *n_nuclei*n_nuclei;
  double edot_ciHeIS  = -ciHeIS *state->abund[sHeII]*state->abund[sELEC]*state->abund[sELEC]*n_nuclei*n_nuclei*n_nuclei;
#elif defined(HELIUM_EQ_COOLING)
  double denom = work->rates.forward[4]*work->rates.forward[6] + work->rates.forward[3]*work->rates.forward[6] + work->rates.forward[3]*work->rates.forward[5];
  if(denom < 1.e-100)
    denom = 1.e-100;

  double fHeI_eq   = XHe * work->rates.forward[4]*work->rates.forward[6] / denom;
  double fHeII_eq  = XHe * work->rates.forward[3]*work->rates.forward[6] / denom;
  double fHeIII_eq = XHe * work->rates.forward[3]*work->rates.forward[5] / denom;
  double felec_new = state->abund[sELEC] + fHeII_eq + 2 * fHeIII_eq;

  double edot_ceHeI  = -ceHeI * fHeI_eq  * felec_new * n_nuclei * n_nuclei;
  double edot_ceHeII = -ceHeII* fHeII_eq * felec_new * n_nuclei * n_nuclei;
#endif


  /******************************* H2 & HD lines ***********************************************/
  double nH  = n_nuclei * XH;
  double yHI = state->abund[sHI]    / XH;
  double yHII = state->abund[sHII]  / XH;
  double yH2I = state->abund[sH2I]  / XH;
  double y_e  = state->abund[sELEC] / XH;
#ifdef HELIUM
  double yHeI = state->abund[sHeI]  / XH;
#else
  double yHeI = XHe / XH;
#endif
#ifdef HD
  double yHDI = state->abund[sHDI]/XH;
#endif

  double edot_H2 = chemistry_models.h2_cooling(nH, T_gas, yHI, yH2I, yHeI, yHII, y_e);
#ifdef HD
  double edot_HD = HD_cooling_GP12(nH, T_gas, yHI, yHDI);
#endif


#ifdef H2_OPACITY  /* accounting for the local optical depth to H2 lines   */
  if(n_nuclei > 1.0e4){
    double N_H2 = l_core * n_nuclei * state->abund[sH2I];
    double f_reduce_H2 = chemistry_models.h2_line_escape(T_gas, N_H2);
    work->diagnostics.H2_escape_fraction = f_reduce_H2;

    if(n_nuclei >= N_EQCHEM) f_reduce_H2 = 0.;
    if(f_reduce_H2 > 1.) f_reduce_H2 = 1.;

    if(f_reduce_H2 < 0){ /* something wrong */
      fflush(stdout);
      endrun();
    }

    edot_H2 *= f_reduce_H2;
    if(FlagPrint)
      printf("f_reduce=%g\n", f_reduce_H2);

#ifdef HD
    double N_HD = l_core * n_nuclei * state->abund[sHDI];
    double f_reduce_HD = chemistry_models.hd_line_escape(T_gas, N_HD);
    work->diagnostics.HD_escape_fraction = f_reduce_HD;

    if(n_nuclei >= N_EQCHEM) f_reduce_HD = 0.;
    if(f_reduce_HD > 1.) f_reduce_HD = 1.;

    if(f_reduce_HD < 0){ /* something wrong */
      fflush(stdout);
      endrun();
    }

    edot_HD *= f_reduce_HD;
#endif
  }
#endif

  /******************************* CII & OI lines ***********************************************/
#ifdef METALLICITY
  if(FlagPrint)
    printf("f_metal=%g\n", metallicity);

  double nOI  = n_nuclei * XH * Y_O_sun * metallicity;
  double nCI  = n_nuclei * XH * Y_C_sun * metallicity;
  double nCII = n_nuclei * XH * Y_C_sun * metallicity;
  double ne   = state->abund[sELEC] * n_nuclei;

  double cool_rate_OI  = OxygenI_twolevel_cooling(nOI, ne, nH, T_gas);
  double cool_rate_CI  = CarbonI_twolevel_cooling(nCI, ne, nH, T_gas);
  double cool_rate_CII = CarbonII_twolevel_cooling(nCII, ne, nH, T_gas);

  double edot_OI  = -1.0 * cool_rate_OI;
  double edot_CI  = -1.0 * cool_rate_CI;
  double edot_CII = -1.0 * cool_rate_CII;
#endif

  /******************************************************************************************/
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/
  /******************************************************************************************/

#ifdef CONTINUUM_COOLING_MATSUKOBA
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, state->abund);
#else
  double edot_reHI    = -reHII*  state->abund[sHII] * state->abund[sELEC] * n_nuclei*n_nuclei;
  double edot_CIE     = -cie  *  state->abund[sH2I] * state->abund[sH2I]  * n_nuclei*n_nuclei;
  double edot_brem    = -brem * state->abund[sHII]  *  state->abund[sELEC]*n_nuclei*n_nuclei;

#ifdef HELIUM
  edot_brem    -= brem*(state->abund[sHeII]+4.0*state->abund[sHeIII]) * state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII1 = -reHeII1*state->abund[sHeII] *state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeII2 = -reHeII2*state->abund[sHeII] *state->abund[sELEC]*n_nuclei*n_nuclei;
  double edot_reHeIII = -reHeIII*state->abund[sHeIII]*state->abund[sELEC]*n_nuclei*n_nuclei;
#endif
#endif

  /********************** opacity for the continuum process *******************************/

  double edot_dust = 0.;
#ifdef DUST_COOLING
  edot_dust = work->continuum.dust_cooling_rate;
#endif
#if defined(DUST_COOLING) && defined(LYMANWERNER)
  double edot_PE = 0.;
  if(work->continuum.dust_temperature < 1000.){
    double G0 = environment->G0;
    edot_PE = dust_PEheating(G0, nH, y_e, T_gas, metallicity);
    if(FlagPrint)
      printf("edot_dust=%g cont=%g T_dust=%g\n", edot_dust, edot_cont, work->continuum.dust_temperature);
  }
#endif

#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * metallicity * n_nuclei * n_nuclei;
#endif

#ifdef CONT_OPACITY
  edot_ceHI *= exp(-tau_cont);
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
  edot_ceHeI *= exp(-tau_cont);
  edot_ceHeII *= exp(-tau_cont);
#endif

  edot_H2   *= exp(-tau_cont);
#ifdef HD
  edot_HD   *= exp(-tau_cont);
#endif
#ifdef METALLICITY
  edot_CI   *= exp(-tau_cont);
  edot_CII  *= exp(-tau_cont);
  edot_OI   *= exp(-tau_cont);
#endif

#ifdef METAL_LINE_COOLING_SD93
  edot_metal_SD93 *= exp(-tau_cont);
#endif

  edot_cont *= continuum_reduction_factor;
#endif

  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = step->piHI  *state->abund[sHI]  *n_nuclei;
#ifdef HELIUM
  double edot_radHeI  = step->piHeI *state->abund[sHeI] *n_nuclei;
  double edot_radHeII = step->piHeII*state->abund[sHeII]*n_nuclei;
#endif

  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  double edot = edot_ceHI + edot_comp + edot_H2
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    + edot_ceHeI + edot_ceHeII
#ifdef HELIUM
    + edot_ciHeI + edot_ciHeII + edot_ciHeIS
    + edot_reHeII1 + edot_reHeII2 + edot_reHeIII
#endif
#endif
#ifdef HD
    + edot_HD
#endif
#ifdef METALLICITY
    + edot_CI + edot_CII + edot_OI
#endif
#ifdef METAL_LINE_COOLING_SD93
    + edot_metal_SD93
#endif
#ifdef DUST_COOLING
    + edot_dust
#ifdef LYMANWERNER
    + edot_PE
#endif
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
    + edot_cont
#else
    + edot_reHI + edot_brem + edot_CIE
#endif
    + edot_radHI
#ifdef HELIUM
    + edot_radHeI + edot_radHeII
#endif
    ;


  if(FlagPrint){
#ifdef METALLICITY
    printf("cool_rate:tot=%g HI=%g H2=%g HD=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g T_gas=%g radHI=%g\n",
	   edot, edot_ceHI, edot_H2, edot_HD, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont,
	   work->continuum.reduction_factor, work->continuum.dust_temperature, T_gas, edot_radHI);
#else
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g cont=%g f_reduce=%g T_gas=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_cont, f_reduce, T_gas);
#endif
  }

#ifdef ID_CHEM_DEBUG
  if(FlagPrint)
    printf("cool_rate:tot=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g Tdust=%g\n",
	   edot, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust,  edot_cont, f_reduce, work->continuum.dust_temperature);
#else
#if defined(THERMAL_FEEDBACK) && defined(METALLICITY)
  if(FlagPrint)
    printf("cool_rate:T_gas=%g n=%g tot=%g rad=%g HI=%g H2=%g comp=%g CI=%g CII=%g OI=%g dust=%g cont=%g f_reduce=%g f_line=%g tau=%g Tdust=%g Trad=%g\n",
	   T_gas, n_nuclei, edot, edot_radHI, edot_ceHI, edot_H2, edot_comp, edot_CI, edot_CII, edot_OI, edot_dust, edot_cont,
	   f_reduce, work->diagnostics.H2_escape_fraction->EscapeProb, tau_cont,
	   work->continuum.dust_temperature, pow(environment->stellar_rad_flux / STEFAN / 4., 0.25));
#endif
#endif

#if defined(DEBUG_CHEMISTRY_SC)
  work->diagnostics.edot_ceHI = edot_ceHI;
  work->diagnostics.edot_comp = edot_comp;
  work->diagnostics.edot_H2 = edot_H2;

#ifdef HELIUM
  work->diagnostics.edot_ceHeI = edot_ceHeI;
  work->diagnostics.edot_ceHeII = edot_ceHeII;
  work->diagnostics.edot_ciHeI = edot_ciHeI;
  work->diagnostics.edot_ciHeII = edot_ciHeII;
  work->diagnostics.edot_ciHeIS = edot_ciHeIS;

  work->diagnostics.edot_reHeII1 = edot_reHeII1;
  work->diagnostics.edot_reHeII2 = edot_reHeII2;
  work->diagnostics.edot_reHeIII = edot_reHeIII;
#endif

#ifdef CONTINUUM_COOLING_MATSUKOBA
  work->diagnostics.edot_cont = edot_cont;
#else
  work->diagnostics.edot_reHI = edot_reHI;
  work->diagnostics.edot_brem = edot_brem;
  work->diagnostics.edot_CIE = edot_CIE;
#endif

#ifdef HD
  work->diagnostics.edot_HD = edot_HD;
#endif
#ifdef METALLICITY
  work->diagnostics.edot_CI  = edot_CI;
  work->diagnostics.edot_CII = edot_CII;
  work->diagnostics.edot_OI  = edot_OI;
#endif
#ifdef DUST_COOLING
  work->diagnostics.edot_dust = edot_dust;
#ifdef LYMANWERNER
  work->diagnostics.edot_PE = edot_PE;
#endif
#endif
#endif

  work->evaluation.rhoU_dot_rad = edot;

  double rhoU = state->utherm_cgs * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.

  work->evaluation.t_cool_rad = (rhoU/edot); /* in seconds */
  if(fabs(edot) < 1.0e-40)
    work->evaluation.t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}


int radiative_cooling_high_T(const ChemEnvironment *environment, const ChemStepData *step, const ChemState *state, ChemWorkspace *work)
{
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei  = rho_gas/PROTONMASS;
  const double T_gas = work->thermodynamics.temperature;
#ifdef UPDATE_L_CORE_IN_SOLVER
  const double l_core = get_core_radius(rho_gas, T_gas);
#else
  const double l_core = step->jeans_length;
#endif

  /******************************************************************************************/
  /********************** continuu cooling (brems, H-fb, CIE, etc) **************************/
  /******************************************************************************************/

#ifdef CONTINUUM_COOLING_MATSUKOBA
  double edot_cont = -continuum_cooling_rate_thin(T_gas, rho_gas, state->abund);
#endif
  double edot_reHI = -8.418e-26*pow(T_gas,0.11) * n_nuclei*n_nuclei;
  double edot_brem = -1.43e-27*sqrt(T_gas) * (1. + 4.*XHe/XH) * n_nuclei * n_nuclei;

  if(FlagPrint)
    printf("T_gas=%g edot_reHI=%g brem=%g cont=%g\n", T_gas, edot_reHI, edot_brem, edot_cont);

  /********************** opacity for the continuum process *******************************/

#ifdef METAL_LINE_COOLING_SD93
  double edot_metal_SD93 = -1. * get_cooling_rate_SD93(T_gas) * environment->metallicity * n_nuclei * n_nuclei;
#endif

  /* ========= Heating functions ============================================================== */
  /* --------- Photoinization heating */
  double edot_radHI   = step->piHI  *state->abund[sHI] *n_nuclei;

  /* ========= net cooling rate in erg/sec/cm3 ================================================= */
  double edot = edot_reHI + edot_brem
#ifdef METAL_LINE_COOLING_SD93
    + edot_metal_SD93
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
    + edot_cont
#endif
    + edot_radHI;

  work->evaluation.rhoU_dot_rad = edot;

  double rhoU = state->utherm_cgs * rho_gas; //modified[0814] SphP[ithis].Entropy is not updated.

  work->evaluation.t_cool_rad = (rhoU/edot); /* in seconds */

  if(FlagPrint){
#ifdef METAL_LINE_COOLING_SD93
    printf("cool_rate:tot=%g reHI=%g brem=%g cont_metal=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, edot_metal_SD93, T_gas, rhoU, work->evaluation.t_cool_rad);
#else
    printf("cool_rate:tot=%g reHI=%g brem=%g T_gas=%g rhoU=%g t_cool=%g\n"
	   , edot, edot_reHI, edot_brem, T_gas, rhoU, work->evaluation.t_cool_rad);
#endif
  }

  if(fabs(edot) < 1.0e-40)
    work->evaluation.t_cool_rad = MAX_REAL_NUMBER;

  return 0;
}

int chemical_cooling(const ChemEnvironment *environment, const ChemState *state, ChemWorkspace *work, double dt, double *abund_old){
  const double T_gas = work->thermodynamics.temperature;
  const double rho_gas = environment->rho_cgs;
  const double n_nuclei = rho_gas / PROTONMASS;

  double edot_H2form,edot_H2form_3body;
  double edot_ciHI, edot_ciHeI, edot_ciHeII;
  double edot_HMform;

  enum SPECIES species;
  double d_abundance[N_SP];

  /* calculate chemical heating/cooling */
  double n_elec = n_nuclei * state->abund[sELEC];
  double n_HI   = n_nuclei * state->abund[sHI];
  double n_HM   = n_nuclei * state->abund[sHM];
  double n_H2I  = n_nuclei * state->abund[sH2I];
#ifdef H2_PLUS
  double n_H2II = n_nuclei * state->abund[sH2II];
#endif
#ifdef HELIUM
  double n_HeI  = n_nuclei * state->abund[sHeI];
  double n_HeII = n_nuclei * state->abund[sHeII];
#endif
#ifdef HD
  double n_DI = n_nuclei * state->abund[sDI];
  double n_DII = n_nuclei * state->abund[sDII];
  double n_HDI = n_nuclei * state->abund[sHDI];
#endif

  for(int j=0;j<N_SP;j++)
    d_abundance[j] = state->abund[j] - abund_old[j];

  double n_crit = (1.e6 / sqrt(T_gas))
    / (1.6 * state->abund[sHI] * exp(-(400.0/T_gas)*(400.0/T_gas)) + 1.4 * state->abund[sH2I] * exp(-12000.0/(T_gas + 1200.0)));

  double edot_H2form_HM, edot_H2form_dust, edot_ciHI_EUV;
  double edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV, edot_H2form_3body_loss;

  if(dt == 0.){
    edot_H2form_HM = eV_to_erg * ( work->rates.forward[8] * n_HM * n_HI * 3.73/ (1. + n_crit/n_HI)
#ifdef H2_PLUS
				   + work->rates.forward[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/n_HI))
#endif
				   ); //1.38?? modified to 1.83


    edot_H2form_dust = 0.;
#ifdef H2_FORMATION_ON_DUST
    edot_H2form_dust = eV_to_erg * work->rates.forward[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/n_HI));
#endif

    edot_H2form_3body = 0.;
#ifdef THREEBODY
    edot_H2form_3body = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI
					     - ( work->rates.forward[11]*n_HI + work->rates.forward[58]*n_H2I) * n_H2I );
#endif
    edot_H2form = edot_H2form_HM + edot_H2form_3body + edot_H2form_dust;

    edot_ciHI   = -13.598 * eV_to_erg * work->rates.forward[1] * n_HI * n_elec;
#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * work->rates.forward[3] * n_HeI  * n_elec;
    edot_ciHeII = -54.416 * eV_to_erg * work->rates.forward[5] * n_HeII * n_elec;
#endif

    edot_HMform = 0.;
  }else{
    if(XH * n_nuclei < N_THREEBODY){
      edot_H2form_HM_loss = eV_to_erg * ( work->rates.forward[8]  * n_HM   * n_HI * 3.73 / (1. + n_nuclei/n_crit)
#ifdef H2_PLUS
					  + work->rates.forward[10] * n_H2II * n_HI * 1.83 / (1. + n_nuclei/n_crit)
#endif
					  );

      edot_H2form_dust_loss = 0.;
#ifdef H2_FORMATION_ON_DUST
      edot_H2form_dust_loss = eV_to_erg * work->rates.forward[51] * n_HI * 4.28 / (1. + n_nuclei/n_crit);
#endif
      edot_H2diss_FUV = eV_to_erg * work->rates.forward[29] * n_H2I * 4.48;

      edot_H2form_3body = 4.48 * eV_to_erg * d_abundance[sH2I] * n_nuclei / dt;
      edot_H2form_3body_loss = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI ) / (1. + n_nuclei/n_crit);


      edot_H2form = edot_H2form_3body - edot_H2form_HM_loss - edot_H2form_dust_loss + edot_H2diss_FUV;

      edot_H2form_HM = eV_to_erg * ( work->rates.forward[8] * n_HM * n_HI * 3.73/ (1. + n_crit/n_nuclei)
#ifdef H2_PLUS
				   + work->rates.forward[10] * n_H2II * n_HI * (1.83 / (1. + n_crit/n_nuclei))
#endif
				   ); //1.38?? modified to 1.83
      edot_H2form_dust = eV_to_erg * work->rates.forward[51] * n_HI * (0.2 + 4.28/ (1. + n_crit/n_nuclei));
      edot_H2form_3body = 4.48 * eV_to_erg * ( ( work->rates.forward[55]*n_HI + work->rates.forward[57]*n_H2I ) * n_HI * n_HI
					       - ( work->rates.forward[11]*n_HI + work->rates.forward[58]*n_H2I) * n_H2I );
      edot_HMform = 0.78 * eV_to_erg * d_abundance[sHM]  * n_nuclei / dt;
    }else{
      edot_H2form = edot_H2form_3body = 4.48 * eV_to_erg * d_abundance[sH2I] * n_nuclei / dt;
      edot_HMform = 0.78 * eV_to_erg * d_abundance[sHM]  * n_nuclei / dt;
    }

    edot_ciHI   = -13.598 * eV_to_erg * d_abundance[sHII]    * n_nuclei / dt;
#ifdef IONIZATION_PHOTON
    edot_ciHI_EUV = 13.598 * eV_to_erg * work->rates.forward[22] * n_HI;
    edot_ciHI  += edot_ciHI_EUV;
#endif

#ifdef HELIUM
    edot_ciHeI  = -24.587 * eV_to_erg * d_abundance[sHeII]   * n_nuclei / dt;
    edot_ciHeII = -54.416 * eV_to_erg * d_abundance[sHeIII]  * n_nuclei / dt;
#endif
  }

  if(FlagPrint){
    printf("edot_chem: n_nuclei=%g H2=%g 3b=%g 3b_l=%g HM_l=%g dust_l=%g FUV=%g ciHI=%g EUV=%g HM=%g ncrit=%g dy_H2=%g\n"
	   , n_nuclei, edot_H2form, edot_H2form_3body, edot_H2form_3body_loss, edot_H2form_HM_loss, edot_H2form_dust_loss, edot_H2diss_FUV
	   , edot_ciHI, edot_ciHI_EUV, edot_HMform, n_crit, d_abundance[sH2I]);
  }

  work->evaluation.rhoU_dot_chem = edot_H2form + edot_ciHI + edot_HMform;
  if(dt < 1.e-30){
    work->evaluation.rhoU_dot_chem = 0.;
    edot_H2form = edot_ciHI = edot_HMform = 0.;
  }

  double rhoU = state->utherm_cgs * environment->rho_cgs;
  //modified[0814] in case for SphP[ithis].Entropy is not updated.

  work->evaluation.t_cool_chem = (rhoU/work->evaluation.rhoU_dot_chem); /* in seconds */

  return 0;
}


#endif /* CHEMISTRY_SC */
