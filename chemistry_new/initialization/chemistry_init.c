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
#include "chemistry_diagnostics.h"
#include "chemistry_init.h"
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

double dln_T;
double chem_e[15];
int FlagPrint = 0;
int N_SP, N_RE;
int count_print;

const double XH = 0.76;
const double XHe = 0.0625;
const double XHD = 3.e-5;

double T_table[N_T];
double k_table[k_MAX][N_T];
double k_table_inverse[k_MAX][N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T], H2_CIE[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

FILE *open_chemistry_table(const char *filename)
{
  char path[2048];

  snprintf(path, sizeof(path), "%s/tables/%s", All.ChemistryTableDir, filename);

  FILE *fp = fopen(path, "r");
  if(fp == NULL)
    printf("Can't open chemistry table %s (All.ChemistryTableDir=%s).\n", path, All.ChemistryTableDir);

  return fp;
}


#ifdef BACKGROUND_LW_RADIATION
void update_background_J21(void){
#ifdef TYPE_B
  double Teff = 1.e4;
#else
  double Teff = 1.e5;
#endif
  double J21 = All.J21 * 1.e-21;

  double hc = PLANCK * CLIGHT_REAL;
  const double E_LW = 13.6;
  double J_13eV = 2. * pow(E_LW * eV_to_erg, 3) / hc / hc / (exp(E_LW*eV_to_erg/BOLTZMANN/Teff) - 1);
  double Norm = J21 / J_13eV;

  double PD_rate, g0;
  get_Hminus_PD_rate(Teff, &PD_rate, &g0);

  All.rate_Hminus_pd = Norm * PD_rate;
#ifdef DUST_COOLING
  All.G0           = Norm * g0 / 1.6e-3;
#endif
  All.rate_H2_diss = 1.39e9 * J21;
}
#endif

int InitChem(void)
{
  enum SPECIES species;

  /* ----- cross-sections, etc ... */
  int ifunc = chem_tabulate();

#ifdef CONTINUUM_COOLING_MATSUKOBA
  read_opacity_tables();
  tabulate_cooling_rate_cont_thin();
#endif
#ifdef DUST_COOLING
  read_dust_opacity_tables();
#endif

#ifdef RYTR_TAU_LW
  tabulate_opacities_in_LW_band();
#endif

#ifdef LYMANWERNER
  make_rate_tables_LW_radiation();
#ifdef BACKGROUND_LW_RADIATION
  update_background_J21();
#endif
#endif

#ifdef METAL_LINE_COOLING_SD93
  read_continuum_cooling_rate_table_SD93();
#endif

#if defined(INIT_TEMPERATURE) || defined(FIX_TEMPERATURE)
  compute_init_gamma();
#endif

#ifdef CHEMISTRY_CORRECT_NEGATIVE_ABUNDANCE
  All.AbundMinimum[sELEC] = 1.e-10;
  All.AbundMinimum[sHI]  = 1.e-10;
  All.AbundMinimum[sHII] = 1.e-10;
  All.AbundMinimum[sH2I] = 1.e-10;
  All.AbundMinimum[sHM]  = 1.e-20;
#ifdef HD
  All.AbundMinimum[sDI]  = 1.e-20;
  All.AbundMinimum[sDII] = 1.e-20;
  All.AbundMinimum[sHDI] = 1.e-20;
#endif
#endif

  N_SP = sEND - sELEC;
  N_RE = (int)chemistry_network_reaction_count();

  return(ifunc);
}

  /*--- [List of Chemistry] ---*/
  /*
    ------- 1:	HI    + e   -> HII   + 2e    hydrogen collisional ionization
    ------- 2:	HII   + e   -> HI    + p     hydrogen radiative recombination
    ------- 3:	HeI   + e   -> HeII  + 2e    He collisional ionization
    ------- 4:	HeII  + e   -> HeI   + p     He recombination
    ------- 5:	HeII  + e   -> HeIII + 2e    HeII collisional ionization
    ------- 6:	HeIII + e   -> HeII  + p     HeII recombination
    ------- 7:	HI    + e   -> HM    + p     photo-attachment
    ------- 8:	HM    + HI  -> H2I*  + e     H2 formation path (H-)
    ------- 9:  HI    + HII -> H2II  + p     H2+ channel
    ------- 10:	H2II  + HI  -> H2I*  + HII   H2+ channel
    ------- 11:	H2I   + HI  -> 3HI
    ------- 12:	H2I   + HII -> H2II  + HI    charge transfer
    ------- 13:	H2I   + e   -> 2HI   + e     collisional dissociation
    ------- 14: HeII  + HI  -> HeI   + HII + p  charge transfer
    ------- 15: HeI   + HII -> HeII  + HI
    ------- 16:	HM    + e   -> HI    + 2e
    ------- 17:	HM    + HI  -> 2H    + e
    ------- 18:	HM    + HII -> 2HI
    ------- 19:	HM    + HII -> H2II  + e
    ------- 20:	H2II  + e   -> 2HI
    ------- 21:	H2II  + HM  -> HI    + H2I

    ------- 22:	HI    + p   -> HII   + e     photo-ionization H
    ------- 23:	HeII  + p   -> HeIII + e     HeII ionization
    ------- 24:	HeI   + p   -> HeII  + e     photo-ionization He
    ------- 25:	HM    + p   -> HI    + e     photo-detachment
    ------- 26:	H2II  + p   -> HI    + HII
    ------- 27:	H2I   + p   -> H2II  + e
    ------- 28:	H2II  + p   -> 2HII  + e
    ------- 29:	H2I   + p   -> 2HI           photo-dissociation

    ------- 60,61: H2    + e   -> HM    + HI
    ------- 62,63: HI    + HI  -> H2II  + e

    HD (Nakamura-Umemura complete model)

    ------- 31:  DII  + e    -> DI   + p      radiative recombination
    ------- 32:  DI   + HII  -> DII  + HI     charge transfer
    ------- 33:  DII  + HI   -> DI   + HII    reverse reaction of 32
    ------- 34:  DI   + HI   -> HDI  + p      * neutral formation path
    ------- 35:  DI   + H2I  -> HI   + HDI    * neutral path by H2
    ------- 36:  HDII + HI   -> HII  + HDI    charge transfer
    ------- 37:  DII  + H2I  -> HII  + HDI    * D+ formation path
    ------- 38:  HDI  + HI   -> H2I  + DI
    ------- 39:  HDI  + HII  -> H2I  + DII
    ------- 40:  DI   + HII  -> HDII + p
    ------- 41:  DII  + HI   -> HDII + p
    ------- 42:  HDII + e    -> HI   + DI     collisional dissociation
    ------- 43:  DI   + e    -> DM   + p      D- formation/photo-attachment
    ------- 44:  DII  + DM   -> 2DI           mutual neutrization
    ------- 45:  HII  + DM   -> DI   + HI     mutual neutrization
    ------- 46:  HM   + DI   -> HI   + DM
    ------- 47:  DM   + HI   -> DI   + HM
    ------- 48:  DM   + HI   -> HDI  + e      * D- formation path

    ------- 50:  DI   + p    -> DII  + e      photo-ionization D
    ------- 51:  HDII + p    -> HI   + DII    photo-dissociation HD+
    ------- 52:  HDII + p    -> HII  + DI     photo-dissociation HD+
    ------- 53:  HDI  + p    -> HI   + DI     photo-dissociation HD

    3 body reactions

    ------- 55:  HI + HI + HI -> H2I + HI
    ------- 56:  H2I + HI     -> HI + HI + HI  same as reaction 11
    ------- 57:  HI + HI + H2I-> H2I + H2I
    ------- 58:  H2I + H2I    -> HI + HI + H2I  reverse reaction
  */


  /* probably this part should be parallelized in the near future */
int chem_tabulate(void)
{
  enum SPECIES species;

  chem_e[sELEC] = 0.;
  chem_e[sHI]   = 216.034;
  chem_e[sHII]  = 1528.084;
  chem_e[sH2I]  = 0.;
#ifdef H2_PLUS
  chem_e[sH2II] = 1488.365;
#endif
  chem_e[sHM]   = 143.246;
#ifdef HD
  chem_e[sDI]   = 219.804;
  chem_e[sDII]  = 1532.211;
  chem_e[sHDI]  = 0.328;
#endif

#ifdef CHECK_PARTITION_FUNCTIONS
  check_partition_functions_H2();
#endif

  /* ----- Compute temperature grid. */
  dln_T= (log(T_end) - log(T_start))/N_T;
  for(int i = 0;i< N_T;i++){
    double ln_T = log(T_start) + i*dln_T;
    T_table[i] = exp(ln_T);
  }

  /* ----- Initialize arrays. */
  for(int i=0;i<N_T;i++){
    for(int j=0;j<k_MAX;j++) k_table[j][i] = 0;
    for(int j=0;j<k_MAX;j++) k_table_inverse[j][i] = 0;

    ceHIa[i]    = 0.;
    ceHeIa[i]   = 0.;
    ceHeIIa[i]  = 0.;
    ciHIa[i]    = 0.;
    ciHeIa[i]   = 0.;
    ciHeISa[i]  = 0.;
    ciHeIIa[i]  = 0.;
#ifndef CONTINUUM_COOLING_MATSUKOBA
    reHIIa[i]   = 0.;
#ifdef HELIUM
    reHeII1a[i] = 0.;
    reHeII2a[i] = 0.;
    reHeIIIa[i] = 0.;
#endif
    brema[i]    = 0.;
    H2_CIE[i]   = 0.;
#endif
  }


  /* ----- Fill in tables over the range T_start to T_end. */
  for(int i=0 ; i<N_T ; i++){
    double z_H2I = chemistry_models.h2_partition(T_table[i]);
    double z_HI  = partition_H(T_table[i]);
#ifdef H2_PLUS
    double z_H2II= pf_z_H2p(T_table[i]);
#endif
    double z_HM  = 1.;
    double z_elec = 2.;
    double z_HII = 1.;

    /* Pure forward reaction rates are selected in forward_rate_table. */
    for(int j=0; j<k_MAX; j++) {
      ForwardRateFunction evaluate = forward_rate_table[j];
      if(evaluate != NULL)
        k_table[j][i] = evaluate(T_table[i]);
    }

#ifdef INVERSE_REACTIONS
    double K_T;
    double Ene, dE, Z_pf, Mass;
    double m_HI, m_HII, m_HM, m_elec, m_H2I, m_H2II;

    m_HI = m_HII = m_HM = PROTONMASS;
    m_H2I = m_H2II = 2. * m_HI;
    m_elec = ELECTRONMASS;

    Ene = 2 * M_PI * BOLTZMANN * T_table[i] / PLANCK / PLANCK;

    /**** HI + e- -> Hp + 2e-****/
    if(T_table[i] > 1.e3){
      dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHII] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_elec/ z_HII / z_elec / z_elec;
      Mass = m_HI * m_elec/ m_HII / m_elec / m_elec;

      K_T = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[1][i] = k_table[1][i] * K_T;
    }

    /**** Hp + e- -> HI  ****/
    dE   = (chem_e[sHII] + chem_e[sELEC] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HII * z_elec / z_HI;
    Mass = m_HII * m_elec / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[2][i] = k_table[2][i] * K_T; //note this value is optically thick limit.

    /**** HI + e- -> HM  ****/
    dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHM]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_elec / z_HM;
    Mass = m_HI * m_elec / m_HM;
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[7][i] = k_table[7][i] * K_T; //note this value is optically thick limit.

    /**** HI + HM -> H2I + e-  ****/
    dE   = (chem_e[sHI] + chem_e[sHM] - chem_e[sH2I] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_HM / z_H2I / z_elec;
    Mass = m_HI * m_HM / m_H2I / m_elec;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[8][i] = k_table[8][i] * K_T;

#ifdef H2_PLUS
    /**** H2I + Hp -> H2II + HI  ****/
    dE   = (chem_e[sH2I] + chem_e[sHII] - chem_e[sH2II] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_H2I * z_HII / z_HI / z_H2II;
    Mass = m_H2I * m_HII / m_HI / m_H2II;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    //k_table_inverse[10][i] = k_table[10][i] * K_T;
    k_table[10][i] = k_table[12][i] * K_T;
#endif

    if(T_table[i]>200.){
      /**** H2I + e- -> 2HI + e-  ****/
      dE   = (chem_e[sH2I] + chem_e[sELEC] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_elec / z_HI / z_HI / z_elec;
      Mass = m_H2I * m_elec / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[13][i] = k_table[13][i] * K_T;
    }

    if(T_table[i]>200.){
      /**** HM + e- -> HI + 2e-  ****/
      dE   = (chem_e[sHM] + chem_e[sELEC] - chem_e[sHI] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_elec / z_HI / z_elec / z_elec;
      Mass = m_HM * m_elec / m_HI / m_elec / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
#ifdef AIU
      k_table_inverse[16][i] = k_table[16][i] * K_T;
#endif
    }

    /**** HM + Hp -> HI + HI  ****/
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[18][i] = k_table[18][i] * K_T;

#ifdef CONTINUUM_COOLING_MATSUKOBA
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    //k_table_HI_inverse[1][i] = k_table_HI[1][i] * K_T;
#endif

    if(T_table[i] > 1.e3){
      /***** HI + HI -> Hp + HI + e- *****/
      dE   = (2*chem_e[sHI] - chem_e[sHII] - chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI / z_HII / z_elec / z_HI;
      Mass = m_HI * m_HI / m_HII / m_elec / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[30][i] = k_table[30][i] * K_T;
    }

    if(T_table[i] > 500.){
      /****  HI + HI + HI -> HI + H2I  ****/
      dE   = (3*chem_e[sHI] - chem_e[sH2I] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI * z_HI / z_H2I / z_HI;
      Mass = m_HI * m_HI * m_HI / m_H2I / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
      k_table_inverse[55][i] = k_table[55][i] * K_T;


      /****  H2I + H2I -> H2I + HI + HI  ****/
      dE   = (2*chem_e[sH2I] - 2*chem_e[sHI] - chem_e[sH2I]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_H2I / z_HI / z_HI / z_H2I;
      Mass = m_H2I * m_H2I / m_HI / m_HI / m_H2I;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table[57][i] = k_table[58][i] * K_T;
    }
#endif
  }


  /* coefficients for cooling and heating */
  for(int i=0; i<N_T; i++){
#ifdef CASE_B
    double lambda_HI   = 2 * 157807.0/T_table[i];
#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    double lambda_HeI  = 2 * 285335.0/T_table[i];
    double lambda_HeII = 2 * 631515.0/T_table[i];
#endif
#endif

    /* ------- Collisional excitations (Black 1981; Cen 1992) */
    ceHIa[i]   = 7.5e-19*exp( -118348./T_table[i] )
      /(1.+sqrt(T_table[i]/1.e5));

    ceHeIa[i]   = 9.1e-27*exp( -13179./T_table[i] )
      *pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5));

    ceHeIIa[i]  = 5.54e-17*exp( -473638./T_table[i] )
      *pow(T_table[i],-0.397)/(1.+sqrt(T_table[i]/1.e5));


    /* ------- Collisional ionizations (Tom's polynomial fits) */
    ciHIa[i]    = 2.18e-11*k_table[1][i];
    ciHeIa[i]   = 3.94e-11*k_table[3][i];
    ciHeIIa[i]  = 8.72e-11*k_table[5][i];

    ciHeISa[i]  = 5.01e-27*pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5))
      *exp( -55338./T_table[i] );

    /* ------- Recombinations */
#ifdef CASE_B
    reHIIa[i]   = 3.435e-30 * T_table[i] * pow(lambda_HI, 1.97)/pow(1.0+pow(lambda_HI/2.25, 0.376), 3.72);

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    reHeII1a[i] = BOLTZMANN * T_table[i] * 1.26e-14 * pow(lambda_HeI, 0.75);

    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */

    reHeIIIa[i] = 8 * 3.435e-30 * T_table[i] * pow(lambda_HeII, 1.97)/pow(1.0+pow(lambda_HeII/2.25, 0.376), 3.72);
#endif
#else
    reHIIa[i]   = 6.50e-27*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/1.e6,0.7));

#if defined(HELIUM) || defined(HELIUM_EQ_COOLING)
    reHeII1a[i] = 1.55e-26*pow(T_table[i],0.3647);

    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */

    reHeIIIa[i] = 3.48e-26*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/4.e6,0.7));
#endif

    /* ------- Bremsstrahlung (Black 1981)(Spitzer & Hart 1979) */
    brema[i]    = 1.43e-27*sqrt(T_table[i])
      *(1.1+0.34*exp(-pow(5.5-log10(T_table[i]),2.0)/3.0));

    log10_T = log10(T_table[i]);
    H2_CIE[i]   = pow(10.0, - 116.6
		      + 96.3   *log10_T
		      - 47.15  *log10_T*log10_T
		      + 10.744 *log10_T*log10_T*log10_T
		      -  0.9155*log10_T*log10_T*log10_T*log10_T);
#endif
  }

  return N_T;
}

#endif /* CHEMISTRY_SC */
