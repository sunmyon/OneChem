#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../../src/standalone/allvars.h"
#include "chemistry_old_api.h"
#include "chemistry.h"

//#define ThisTask 0
#define  CMB_BG 0

double compute_b_esc(double tau_line, int ithis);
double compute_Q_bg(double T_nu);

int H2line_tabulate(void);
double compute_freefall_time(double rho_gas);

void gaussj5(double a[5][5], int n, double b[5][2], int m);


double H2line_LTE(double T);
double H2cool_LTEfit(double T);

double H2cool_G15(double xnH,double T_K,double y_H,double y_H2,double y_e,double y_Hp);
double H2cool_LowD_G15(double xnH,double T_K,double y_H,double y_H2,double y_e,double y_Hp,double y_He);
double fesc_fit(int ithis, double T_K, double NH2);

double g_weight[20];
double n_pop[3][20];

double A[3][20][3][20], B[3][20][3][20];
double E_vJ[3][20], E_vJ_erg[3][20];
double E_vJ_eV[3][20], E_vJ_K[3][20];
double del_E[3][20][3][20],del_E_K[3][20][3][20],del_nu[3][20][3][20];

double Lambda_frac[3][20][3][20], Lambda_frac_thick[3][20][3][20];
double alpha[3][20][3][20], n_pop_dummy[3][20][3][20], alpha_h2[3][20][3][20];


double H2_cooling_O00(int ithis, double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e)
{
  double rho,n_HI,n_H2;
  double mu,T3;
  int vib,vi,vf, Jrot,Ji,Jf;

  rho = nH*PROTONMASS / XH;
  //mu = SphP[ithis].molecular_weight;
  mu = 1.;

  n_HI = yHI * nH;
  n_H2 = yH2 * nH;

  T3 = T_gas/1000.0;

// begin computing population, f_vJ ////////////////////////////////////////////////////////////////
  double g_0,g_1,g_2;
  double A_10,A_20,A_21;
  double gammaV_H_10,gammaV_H_20,gammaV_H_21,
        gammaV_H2_10,gammaV_H2_20,gammaV_H2_21;
  double Erg10,Erg20,Erg21;
  double DT_10,DT_20,DT_21;
  double C_10,C_20,C_21, C_01,C_02,C_12;
  double R_10,R_20,R_21, R_01,R_02,R_12;
  double Q_10,Q_20,Q_21;
  double f_v[3];


  g_0 = 1.0;
  g_1 = 1.0;
  g_2 = 1.0;

  // vib spon.
  A_10 = 8.3e-7;
  A_20 = 4.1e-7;
  A_21 = 1.1e-6;

  // vib coll.
  gammaV_H_10 = 1.0e-12 * sqrt(T_gas) * exp(-1000.0/T_gas);
  gammaV_H_20 = 1.6e-12 * sqrt(T_gas) * exp(-pow(400.0/T_gas,2.0));
  gammaV_H_21 = 4.5e-12 * sqrt(T_gas) * exp(-pow(500.0/T_gas,2.0));

  gammaV_H2_10 = 1.4e-12 * sqrt(T_gas) * exp(-(18100.0/(T_gas+1200.0)));
  gammaV_H2_20 = 0.0;
  gammaV_H2_21 = gammaV_H2_10;

  Erg10 = E_vJ_erg[1][0]-E_vJ_erg[0][0];
  Erg20 = E_vJ_erg[2][0]-E_vJ_erg[0][0];
  Erg21 = E_vJ_erg[2][0]-E_vJ_erg[1][0];

  DT_10 = Erg10/BOLTZMANN;
  DT_20 = Erg20/BOLTZMANN;
  DT_21 = Erg21/BOLTZMANN;

  // total coll. rate
  C_10 = gammaV_H_10*n_HI + gammaV_H2_10*n_H2;
  C_20 = gammaV_H_20*n_HI + gammaV_H2_20*n_H2;
  C_21 = gammaV_H_21*n_HI + gammaV_H2_21*n_H2;

  C_01 = C_10 * exp(-DT_10/T_gas);
  C_02 = C_20 * exp(-DT_20/T_gas);
  C_12 = C_21 * exp(-DT_21/T_gas);

  Q_10 = compute_Q_bg(DT_10);
  Q_20 = compute_Q_bg(DT_20);
  Q_21 = compute_Q_bg(DT_21);

  // total transition rate
  R_10 = A_10*(1.0+Q_10) + C_10;
  R_20 = A_20*(1.0+Q_20) + C_20;
  R_21 = A_21*(1.0+Q_21) + C_21;

  R_01 = (g_1/g_0)*A_10*Q_10 + C_01;
  R_02 = (g_2/g_0)*A_20*Q_20 + C_02;
  R_12 = (g_2/g_1)*A_21*Q_21 + C_12;

  // frac. of population (vib)
  f_v[0] = ( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
       /( (R_01+R_02+R_20)*(R_10+R_12+R_21) - (R_01-R_21)*(R_10-R_20) );
  f_v[1] = (f_v[0]*(R_01-R_21)+R_21)/(R_10+R_12+R_21);
  f_v[2] = (f_v[0]*R_02+f_v[1]*R_12)/(R_21+R_20);



  double g_u,g_l, C_ul,C_lu;
  double gammaJ_H[20],gammaJ_H2[20];
  double DE,DT,Q_ul;
  double ratio;
  double fJ[20];
  double f_para,f_ortho, z_para,z_ortho;
  double f_vJ[3][20];

  // rotational HM89
  for(Jrot=2;Jrot<20;Jrot++) {
    gammaJ_H[Jrot]  = 4.6e-12 * (2*Jrot-3) * sqrt(T_gas) * sqrt(1.0 + 2.0*85.25*(2*Jrot-1)/T_gas)
       * exp(-(10.0*85.25*(2*Jrot-1))/(T_gas + 85.25*Jrot*(Jrot+1)) + 0.1187*(4*Jrot-2));

    gammaJ_H2[Jrot] = (3.3e-12 + 6.6e-12*T3)
                    * (0.276*Jrot*Jrot*exp(-pow(Jrot/3.18,1.7)));
  }

  // rotational GP98
  gammaJ_H[2]  = 2.93e-14 + 1.21e-15*T_gas + 2.16e-19*T_gas*T_gas + 1.32e-21*T_gas*T_gas*T_gas;
  gammaJ_H[3]  = 8.34e-14 + 5.97e-16*T_gas + 7.76e-19*T_gas*T_gas + 1.72e-21*T_gas*T_gas*T_gas;
  gammaJ_H[4]  = 7.54e-14 + 2.65e-16*T_gas + 2.14e-19*T_gas*T_gas + 2.65e-21*T_gas*T_gas*T_gas;
  gammaJ_H[5]  = 2.95e-14 + 1.21e-16*T_gas + 5.53e-19*T_gas*T_gas + 2.51e-21*T_gas*T_gas*T_gas;


  f_para = 0.25;
  f_ortho = 0.75;
  for(vib=0;vib<3;vib++) {

  // frac. of population (rot)
    z_para = 0.0;
    z_ortho = 0.0;

    fJ[0] = 1.0;
    fJ[1] = 1.0;

    for(Jrot=2;Jrot<20;Jrot++) {

      if(Jrot%2==0) z_para = z_para + fJ[Jrot-2];
      if(Jrot%2==1) z_ortho = z_ortho + fJ[Jrot-2];

      C_ul = gammaJ_H[Jrot]*n_HI + gammaJ_H2[Jrot]*n_H2;
      DE = E_vJ_erg[vib][Jrot]-E_vJ_erg[vib][Jrot-2];
      DT = DE/BOLTZMANN;
      Q_ul = compute_Q_bg(DT);
      g_u = 2*Jrot+1;
      g_l = 2*Jrot-3;
      ratio = ( g_u/g_l*A[vib][Jrot][vib][Jrot-2]*Q_ul +  g_u/g_l*C_ul*exp(-DE/BOLTZMANN/T_gas) )
            / ( A[vib][Jrot][vib][Jrot-2]*(1.0+Q_ul) + C_ul );
//          printf("%3d %3d %3d %3d %9.3e\n",vib,Jrot,vib,Jrot-2,A[vib][Jrot][vib][Jrot-2]);
      fJ[Jrot] = fJ[Jrot-2] * ratio;
    }

    for(Jrot=0;Jrot<20;Jrot++) {
      if(Jrot%2==0) f_vJ[vib][Jrot] = (fJ[Jrot]/z_para) * f_para * f_v[vib];
      if(Jrot%2==1) f_vJ[vib][Jrot] = (fJ[Jrot]/z_ortho) * f_ortho * f_v[vib];
    }
  }


// end computing population, f_vJ //////////////////////////////////////////////////////////////////


// begin computing escape fraction, tau_ul & cooling rate //////////////////////////////////////////
  double nu_if, tau_if, b_esc_if;
  double g_i,g_f,Q_if;
  double cool_rate_H2;
  double l_sh,Dv_D,t_col,Dv_bulk,dvdz;
  double cool_rate_H2_thin;


  l_sh = sqrt(M_PI*BOLTZMANN*T_gas/(GRAVITY*mu*PROTONMASS*rho));
  t_col = compute_freefall_time(rho);
  Dv_bulk = l_sh/3.0/t_col;
  Dv_D = sqrt(2.0*BOLTZMANN*T_gas/(2.0*PROTONMASS));
  if(Dv_D>Dv_bulk) dvdz = Dv_D/l_sh;
  else             dvdz = Dv_bulk/l_sh;

  cool_rate_H2 = cool_rate_H2_thin = 0.0;

  for(Ji=0;Ji<18;Ji++) {
    for(vi=1;vi<3;vi++) {
      for(vf=0;vf<=vi-1;vf++) {
        Jf=Ji+2;
          DE = E_vJ_erg[vi][Ji]-E_vJ_erg[vf][Jf];
          DT = DE / BOLTZMANN;
          nu_if = DE / PLANCK;
          g_i = 2*Ji+1;
          g_f = 2*Jf+1;
          Q_if = compute_Q_bg(DT);
          tau_if = A[vi][Ji][vf][Jf]*pow(CLIGHT,3)/(8.0*M_PI*pow(nu_if,3))
                 * (g_i/g_f*f_vJ[vf][Jf] - f_vJ[vi][Ji])*n_H2/dvdz;
          b_esc_if = compute_b_esc(tau_if, ithis);
          if(f_vJ[vi][Ji]!=0){
            cool_rate_H2 += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE * b_esc_if
                          * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
            cool_rate_H2_thin += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE
                          * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
	  }
      }
    }
  }
  for(Ji=1;Ji<20;Ji++) {
    for(vi=1;vi<3;vi++) {
      for(vf=0;vf<=vi-1;vf++) {
        Jf=Ji;
          DE = E_vJ_erg[vi][Ji]-E_vJ_erg[vf][Jf];
          DT = DE / BOLTZMANN;
          nu_if = DE / PLANCK;
          g_i = 2*Ji+1;
          g_f = 2*Jf+1;
          Q_if = compute_Q_bg(DT);
          tau_if = A[vi][Ji][vf][Jf]*pow(CLIGHT,3)/(8.0*M_PI*pow(nu_if,3))
                 * (g_i/g_f*f_vJ[vf][Jf] - f_vJ[vi][Ji])*n_H2/dvdz;
          b_esc_if = compute_b_esc(tau_if, ithis);
          if(f_vJ[vi][Ji]!=0){
            cool_rate_H2 += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE * b_esc_if
                          * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
            cool_rate_H2_thin += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE
                          * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
	  }
      }
    }
  }
  for(Ji=2;Ji<20;Ji++) {
    for(vi=0;vi<3;vi++) {
      for(vf=0;vf<=vi;vf++) {
        Jf=Ji-2;
          DE = E_vJ_erg[vi][Ji]-E_vJ_erg[vf][Jf];
          DT = DE / BOLTZMANN;
          nu_if = DE / PLANCK;
          g_i = 2*Ji+1;
          g_f = 2*Jf+1;
          Q_if = compute_Q_bg(DT);
          tau_if = A[vi][Ji][vf][Jf]*pow(CLIGHT,3)/(8.0*M_PI*pow(nu_if,3))
                 * (g_i/g_f*f_vJ[vf][Jf] - f_vJ[vi][Ji])*n_H2/dvdz;
          b_esc_if = compute_b_esc(tau_if, ithis);
          if(f_vJ[vi][Ji]!=0){
            cool_rate_H2 += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE * b_esc_if
	      * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
            cool_rate_H2_thin += f_vJ[vi][Ji] * A[vi][Ji][vf][Jf] * DE
	      * (1.0-Q_if*(g_i/g_f*f_vJ[vf][Jf]/f_vJ[vi][Ji] - 1.0));
	  }

      }
    }
  }

  return -n_H2 * cool_rate_H2;
}



double H2_cooling_Y06(int ithis, double nH, double T_gas, double yHI, double yH2, double yHeI, double yHII, double y_e)
{
  double  T_CMB,TM,LT,T3,cool_rate,heat_rate; 
  double  cool_rate_H_H2, cool_rate_H2_H2; 
  double  E2_E0, E3_E1, E10, E20;
  double  gamma2, gamma3, gamma21, gamma10, gamma20; 
  double  HighDensLimit_H_rot, HighDensLimit_H_vib, HighDensLimit_H, LowDensLimit_H;
  double  LowDensLimit_H2_rot, HighDensLimit_H2_rot, LowDensLimit_H2_vib, HighDensLimit_H2_vib;
  double  Lrot_H2, Lvib_H2;
  double  cool_rate_H2;

  /* --------- Cooling from H2 (H-H2 collisions) */
  
  if(T_gas<13.0) TM=13.0;
  else if(T_gas>1.0e5) TM = 1.0e5;
  else TM=T_gas;

  LT = log10(TM);
  
  /*    low density limit from Galli and Palla   */  
  LowDensLimit_H  = pow(10.0, -103.0+97.59*LT-48.05*LT*LT+10.80*LT*LT*LT-0.9032*LT*LT*LT*LT);
    
  /* high density limit from HM79   */ 
  T3 = TM/1000.;
  HighDensLimit_H_rot = ((9.5e-22*pow(T3,3.76))/(1.+0.12*pow(T3,2.1))*exp(-pow(0.13/T3,3))+3.e-24*exp(-0.51/T3));      
  HighDensLimit_H_vib = (6.7e-19*exp(-5.86/T3) + 1.6e-18*exp(-11.7/T3));
  
  HighDensLimit_H  = (HighDensLimit_H_rot + HighDensLimit_H_vib)/(yHI*nH);
  
  /*    cooling rate per molecule [erg/s]   */  
  cool_rate_H_H2 = yHI*nH * HighDensLimit_H/(1.+(HighDensLimit_H/LowDensLimit_H));       
  
  /*  H2-H2 collisions Hollenbach-McKee eqs (6.39)-(6.41)  */
  
  /* high density limit from HM79    */
  T3 = TM/1000.;
  HighDensLimit_H2_rot = ((9.5e-22*pow(T3, 3.76))/(1.+0.12*pow(T3, 2.1))*exp(-pow(0.13/T3, 3))+3.e-24*exp(-0.51/T3));
  HighDensLimit_H2_rot /= (yH2 * nH);
  
  HighDensLimit_H2_vib = (6.7e-19*exp(-5.86/T3) + 1.6e-18*exp(-11.7/T3));
  HighDensLimit_H2_vib /= (yH2 * nH);
  
  
  gamma2=(3.3e-12 + 6.6e-12 * T3)*(0.276*2*2*exp(-pow(2/3.18, 1.7)));
  gamma3=(3.3e-12 + 6.6e-12 * T3)*(0.276*3*3*exp(-pow(3/3.18, 1.7)));
  gamma10 = 1.4e-12*sqrt(TM)*exp(-(12000.0/(TM+1200.0)));
  gamma20 = 0;
  gamma21 = 1.4e-12*sqrt(TM)*exp(-(12000.0/(TM+1200.0)));
  
  
  E2_E0 = 7.0366909e-14; /* = 509.7 K. The commoly used 512K neglects higher order terms of J, J+1, see my code energy_H2.pro  */ 
  E3_E1 = 1.1655254e-13; /* = 844.2 K  */
    
  E10 = 8.2635748e-13;
  E20 = 2.0367850e-12;
  
  LowDensLimit_H2_rot = 0.25 * 5        * gamma2 * exp(-(E2_E0)/(BOLTZMANN * TM))*E2_E0
    +0.75 * (7.0/3.0)* gamma3 * exp(-(E3_E1)/(BOLTZMANN * TM))*E3_E1;
  
  LowDensLimit_H2_vib = gamma10 * exp(-E10/(BOLTZMANN * TM))*E10 + gamma20 * exp(-E20/(BOLTZMANN * TM))*E20;
  
    
  Lrot_H2 = HighDensLimit_H2_rot/(1.+(HighDensLimit_H2_rot/LowDensLimit_H2_rot));
  Lvib_H2 = HighDensLimit_H2_vib/(1.+(HighDensLimit_H2_vib/LowDensLimit_H2_vib));
  
  /*    cooling rate per molecule [erg/s]  */
  cool_rate_H2_H2 = yH2*nH * (Lrot_H2 + Lvib_H2); 
  
  cool_rate_H2 = (cool_rate_H_H2 + cool_rate_H2_H2);   
  /* cooling rate per volume [erg/s/cm^3]  */

  return cool_rate_H2;
}





double H2_cooling_GP98(int i, double n_nuclei, double T_gas, double SphP_HI, double SphP_H2I)
{
  /* new H2 cooling rate in Gall & Palla (1998) */

  double T, logT, logT_2, logT_3, logT_4;
  double cool_rate_GP98_H2H;
  double cool_rate_LOW_H2, cool_rate_LTE_H2, cool_rate_H2;

  T      = dmin(T_gas, 1.e4);
  logT   = log10(T);
  logT_2 = logT   * logT;
  logT_3 = logT_2 * logT;
  logT_4 = logT_3 * logT;

  /* ----- Glover&Abel'08 (Table 8) : H2 cooling rate in low-density limit, assumes 3:1 o:p ratio */
  /* ----- H2-H : based on Wrathmall&Flower'07 */
  cool_rate_GP98_H2H = pow(10.0, -103.
			   + 97.59 * logT
			   - 48.05 * logT_2
			   + 10.80 * logT_3
			   - .9032 * logT_4);
			          
  /* cooling rate per molecule [erg sec^-1] */
  cool_rate_LOW_H2 = cool_rate_GP98_H2H  * SphP_HI   * n_nuclei;
			      
  /* H2 : LTE(high-density) cooling function */
  cool_rate_LTE_H2 = H2line_LTE(T_gas);

  cool_rate_H2 = cool_rate_LTE_H2/(1.+(cool_rate_LTE_H2/cool_rate_LOW_H2));  //=[erg sec^-1]

  /* ----- cooling rate per volume [erg/s/cm^3]  */
  return - (SphP_H2I * n_nuclei) * cool_rate_H2;
}


double H2_cooling_GP12(int i, double nH, double T_gas, double yHI, double yH2I, double yHeI, double yHII, double y_e)
{
  /* new H2 cooling rate in Gall & Palla (2012) */

  double T3, log10T;
  double logT3, logT3_p2, logT3_p3, logT3_p4, logT3_p5;
  double logT3_MAX, logT3_p2_MAX, logT3_p3_MAX, logT3_p4_MAX, logT3_p5_MAX;
  double cool_rate_G08_H2H, cool_rate_G08_H2H2, cool_rate_G08_H2He, cool_rate_G08_H2p, cool_rate_G08_H2e;
  double cool_rate_LOW_H2, cool_rate_LTE_H2, cool_rate_H2;

  T3 = T_gas/1.e3;
  logT3    = log10(T3);
  logT3_p2 = logT3    * logT3;
  logT3_p3 = logT3_p2 * logT3;
  logT3_p4 = logT3_p3 * logT3;
  logT3_p5 = logT3_p4 * logT3;

  /* At T > 10000K, assume rate same as at T=10000K */
  logT3_MAX    = log10(1.e4/1.e3);
  logT3_p2_MAX = logT3_MAX    * logT3;
  logT3_p3_MAX = logT3_p2_MAX * logT3;
  logT3_p4_MAX = logT3_p3_MAX * logT3;
  logT3_p5_MAX = logT3_p4_MAX * logT3;

  log10T = log10(T_gas);


  /* ----- Glover&Abel'08 (Table 8) : H2 cooling rate in low-density limit, assumes 3:1 o:p ratio */
  
  /* ----- H2-H : based on Wrathmall&Flower'07 */
  if(T_gas > 1.e4){
    cool_rate_G08_H2H = pow(10.0, -24.311209  
			          + 4.6450521 * logT3_MAX
			          - 3.7209846 * logT3_p2_MAX
			          + 5.9369081 * logT3_p3_MAX
			          - 5.5108047 * logT3_p4_MAX
			          + 1.5538288 * logT3_p5_MAX); 
  }else if(T_gas > 1.e3){
    cool_rate_G08_H2H = pow(10.0, -24.311209  
			          + 4.6450521 * logT3 
			          - 3.7209846 * logT3_p2 
			          + 5.9369081 * logT3_p3 
			          - 5.5108047 * logT3_p4 
			          + 1.5538288 * logT3_p5); 
  }else if(T_gas > 1.e2){
    cool_rate_G08_H2H = pow(10.0, -24.311209 
		     	          + 3.5692468 * logT3 
		     	          -11.332860  * logT3_p2 
		     	          -27.850082  * logT3_p3 
		     	          -21.328264  * logT3_p4 
		     	          - 4.2519023 * logT3_p5);
  }else{
    cool_rate_G08_H2H = pow(10.0, -16.818342 
			          +37.383713  * logT3 
			          +58.145166  * logT3_p2 
			          +48.656103  * logT3_p3 
			          +20.159831  * logT3_p4 
			          + 3.8479610 * logT3_p5);
  } 
  

  /* ----- H2-H2 : based on Le Bourlot+'99 */  
  if(T_gas > 1.e4){ 
    cool_rate_G08_H2H2 = pow(10.0, -23.962112 
			           + 2.09433740  * logT3_MAX
			           - 0.77151436  * logT3_p2_MAX
			           + 0.43693353  * logT3_p3_MAX
			           - 0.14913216  * logT3_p4_MAX
			           - 0.033638326 * logT3_p5_MAX);
  }else{
    cool_rate_G08_H2H2 = pow(10.0, -23.962112   
			           + 2.09433740  * logT3  
			           - 0.77151436  * logT3_p2 
			           + 0.43693353  * logT3_p3 
			           - 0.14913216  * logT3_p4 
			           - 0.033638326 * logT3_p5);
  }
   
  
  /* ----- H2-He : based on Flower+'98, Balakrishnan+'99 */
  if(T_gas > 1.e4){
    cool_rate_G08_H2He = pow(10.0, -23.703432 
			           + 2.2601756  * logT3_MAX
			           - 1.0736440  * logT3_p2_MAX
			           + 0.66060709 * logT3_p3_MAX
			           - 0.26963286 * logT3_p4_MAX
			           + 0.12261185 * logT3_p5_MAX);
  }else{
    cool_rate_G08_H2He = pow(10.0, -23.689237   
			           + 2.1892372  * logT3  
			           - 0.81520438 * logT3_p2 
			           + 0.29036281 * logT3_p3 
			           - 0.16596184 * logT3_p4 
			           + 0.19191375 * logT3_p5);
  }
  
  
  /* ----- H2-H+ : based on Gerlich'90 and Krstic'02; Fit accurate to within 5%. */
  if (T_gas > 1.e4){
    cool_rate_G08_H2p = pow(10.0, -21.716699   
			          + 1.3865783   * logT3_MAX  
			          - 0.37915285  * logT3_p2_MAX  
			          + 0.11453688  * logT3_p3_MAX  
			          - 0.23214154  * logT3_p4_MAX  
			          + 0.058538864 * logT3_p5_MAX);
  }else{
    cool_rate_G08_H2p = pow(10.0, -21.716699  
			          + 1.3865783   * logT3  
			          - 0.37915285  * logT3_p2  
			          + 0.11453688  * logT3_p3  
			          - 0.23214154  * logT3_p4  
			          + 0.058538864 * logT3_p5);
  }

  
  /* ----- H2-e : based on Draine+'83*/  
  if(T_gas > 1.e4){
    cool_rate_G08_H2e = pow(10.0, -22.190316
      			          + 1.5728955  * logT3_MAX
			          - 0.21335100 * logT3_p2_MAX
			          + 0.96149759 * logT3_p3_MAX
			          - 0.91023195 * logT3_p4_MAX
			          + 0.13749749 * logT3_p5_MAX);
  }else if(T_gas > 2.e2){
    cool_rate_G08_H2e = pow(10.0, -22.190316  
			          + 1.5728955  * logT3  
			          - 0.21335100 * logT3_p2 
			          + 0.96149759 * logT3_p3 
			          - 0.91023195 * logT3_p4 
			          + 0.13749749 * logT3_p5);
  }else{
    cool_rate_G08_H2e = pow(10.0, -34.286155 
			          -48.537163   * logT3  
			          -77.121176   * logT3_p2 
			          -51.352459   * logT3_p3 
			          -15.169160   * logT3_p4 
			          - 0.98120322 * logT3_p5);
  }


  /* cooling rate per molecule [erg sec^-1] */
  cool_rate_LOW_H2 = (cool_rate_G08_H2H  * yHI   * nH
		      + cool_rate_G08_H2H2 * yH2I * nH
		      + cool_rate_G08_H2He * yHeI * nH
		      + cool_rate_G08_H2p  * yHII * nH
		      + cool_rate_G08_H2e  * y_e  * nH);
			      

  cool_rate_LTE_H2 = H2line_LTE(T_gas);
  cool_rate_H2 = cool_rate_LTE_H2/(1.+(cool_rate_LTE_H2/cool_rate_LOW_H2));  //=[erg sec^-1]
  
  /* ----- cooling rate per volume [erg/s/cm^3]  */
  return - (yH2I * nH) * cool_rate_H2;
}

double H2line_LTE(double T)
{  
  int v,J,vI,JI,vF,JF;
  double Z, cool_rate_thin;

  if(T < 75.)
    return 0.;

  /* Partition function */
  Z=0;
  for(v=0;v<=2;v++)
    for(J=0;J<=19;J++)
      Z += g_weight[J] * (2*J + 1) * exp(-E_vJ_K[v][J]/T);        
  
  /* Boltzmann distribution */
  for(v=0;v<=2;v++)
    for(J=0;J<=19;J++)
      n_pop[v][J]= g_weight[J] * (2*J + 1) * exp(-E_vJ_K[v][J]/T)/Z;
  
  /* initialization */
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++)
	  Lambda_frac[vI][JI][vF][JF] = Lambda_frac_thick [vI][JI][vF][JF] = 0;
  

  /*  line emissivity
   computing the fractional contribution */
  
  /* J=0 */
  Lambda_frac[1][0][0][0] = A[1][0][0][0] * n_pop[1][0] * (E_vJ_erg[1][0] - E_vJ_erg[0][0]);
  Lambda_frac[2][0][0][0] = A[2][0][0][0] * n_pop[2][0] * (E_vJ_erg[2][0] - E_vJ_erg[0][0]);
  Lambda_frac[2][0][1][0] = A[2][0][1][0] * n_pop[2][0] * (E_vJ_erg[2][0] - E_vJ_erg[1][0]);

  /* J=1 */
  Lambda_frac[1][1][0][1] = A[1][1][0][1] * n_pop[1][1] * (E_vJ_erg[1][1] - E_vJ_erg[0][1]);
  Lambda_frac[1][1][0][3] = A[1][1][0][3] * n_pop[1][1] * (E_vJ_erg[1][1] - E_vJ_erg[0][3]);
  
  Lambda_frac[2][1][0][1] = A[2][1][0][1] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[0][1]);
  Lambda_frac[2][1][0][3] = A[2][1][0][3] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[0][3]);
  Lambda_frac[2][1][1][1] = A[2][1][1][1] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[1][1]);
  Lambda_frac[2][1][1][3] = A[2][1][1][3] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[1][3]);
  
  
  for(J=2;J<=17;J++){
    
    /* v=0 */
    Lambda_frac[0][J][0][J-2] = A[0][J][0][J-2] * n_pop[0][J] * (E_vJ_erg[0][J] - E_vJ_erg[0][J-2]);
    
    /* v=1 */
    Lambda_frac[1][J][0][J-2] = A[1][J][0][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[1][J][0][J]   = A[1][J][0][J]   * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J]);
    Lambda_frac[1][J][0][J+2] = A[1][J][0][J+2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J+2]);
    Lambda_frac[1][J][1][J-2] = A[1][J][1][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[1][J-2]);

    /* v=2 */
    Lambda_frac[2][J][0][J-2] = A[2][J][0][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[2][J][0][J]   = A[2][J][0][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J]);
    Lambda_frac[2][J][0][J+2] = A[2][J][0][J+2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J+2]);
    Lambda_frac[2][J][1][J-2] = A[2][J][1][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J-2]);
    Lambda_frac[2][J][1][J]   = A[2][J][1][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J]);
    Lambda_frac[2][J][1][J+2] = A[2][J][1][J+2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J+2]);

    Lambda_frac[2][J][2][J-2] = A[2][J][2][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[2][J-2]);
  }

  
  for(J=18;J<=19;J++){
    
    /* v=0 */
    Lambda_frac[0][J][0][J-2] = A[0][J][0][J-2] * n_pop[0][J] * (E_vJ_erg[0][J] - E_vJ_erg[0][J-2]);
    
    /* v=1 */
    Lambda_frac[1][J][0][J-2] = A[1][J][0][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[1][J][0][J]   = A[1][J][0][J]   * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J]);
    
    Lambda_frac[1][J][1][J-2] = A[1][J][1][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[1][J-2]);

    /* v=2 */
    Lambda_frac[2][J][0][J-2] = A[2][J][0][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[2][J][0][J]   = A[2][J][0][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J]);
    
    Lambda_frac[2][J][1][J-2] = A[2][J][1][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J-2]);
    Lambda_frac[2][J][1][J]   = A[2][J][1][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J]);
  
    Lambda_frac[2][J][2][J-2] = A[2][J][2][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[2][J-2]);
  }


  cool_rate_thin = 0;  
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++){
	  cool_rate_thin += Lambda_frac[vI][JI][vF][JF];
	}
  
  return cool_rate_thin;
}


/*-----------------------------------------------------------------------
  ! LTE H2 cooling fit by Glover(2015) and fesc fit by Omukai-san
  !
  !     xLmbd_H2: H2 cooling rate per unit volume
  !     skeMeMo:I use this for whole region of density
  !             Of course, it is mainly for n > 1.d8
  !-----------------------------------------------------------------------*/
double H2_cooling_G15(int i, double nH, double T_K,double yHI,double yH2,double yHeI,double yHII, double y_e){
  double xLmbd_H2;
  double xLdH2_n0,xLdH2_LTE,y_He,T_rad, xLdH2_dbg,xLdH2_dbg2;

  //y_He = yHe;
  T_rad = 0.;

  if(T_K > 2.e4)
    T_K = 2.e4;

  //lineがthickになった後に、opt-thin low-density cooling rateが現れないように場合分け
  if(nH < 1.e8){
    xLdH2_n0  = H2cool_LowD_G15(nH, T_K, yHI, yH2, y_e, yHII, yHeI);
    xLdH2_LTE = H2cool_LTEfit(T_K) * nH * yH2;
    xLmbd_H2  = xLdH2_LTE / (1. + xLdH2_LTE/xLdH2_n0);

    if(FlagPrint)
      printf("[H2cool]net=%g xLdH2_lowd=%g LTE=%g (%g)\n"
	     , xLmbd_H2, xLdH2_n0, xLdH2_LTE, H2line_LTE(T_K) * nH * yH2);
  }else{
    xLdH2_LTE = H2cool_LTEfit(T_K) * nH * yH2;
    xLmbd_H2  = xLdH2_LTE;
  }

  return -1 * xLmbd_H2;
}




double H2cool_LTEfit(double T_K){
  double xLmbd_H2;
  double xLmbd,xLmbdr=0.;
  double xlgT3,a,ar;
  
  if(T_K > 100.){
    xlgT3=log10(T_K*1.e-3);
    a=-20.584225
      +5.0194035*xlgT3
      -1.5738805*pow(xlgT3,2) 
      -4.7155769*pow(xlgT3,3)
      +2.4714161*pow(xlgT3,4)
      +5.4710750*pow(xlgT3,5)
      -3.9467356*pow(xlgT3,6)
      -2.2148338*pow(xlgT3,7)
      +1.8161874*pow(xlgT3,8);
	
    xLmbd=pow(10.,a);
  }else
    xLmbd=0.;
    
  /*if(T_rad > 100.){
    xlgT3=log10(T_rad*1.e-3);
    ar=-20.584225
      +5.0194035*xlgT3
      -1.5738805*pow(xlgT3,2)
      -4.7155769*pow(xlgT3,3)
      +2.4714161*pow(xlgT3,4)
      +5.4710750*pow(xlgT3,5)
      -3.9467356*pow(xlgT3,6)
      -2.2148338*pow(xlgT3,7)
      +1.8161874*pow(xlgT3,8);
    xLmbdr=pow(10.,ar);
  }else
    xLmbdr=0.;*/

  xLmbd_H2 = xLmbd-xLmbdr;
  //xLmbd_H2 = xLmbd_H2 * xnH * y_H2;

  //fesc =  fesc_fit(xNc_H2,T_K);
  //xLmbd_H2=xLmbd_H2*fesc;

  return xLmbd_H2;
}

/*  gives fitting value for fesc  */
double fesc_fit(int ithis, double T_K, double xNc_H2){
  double fesc;
  double a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5,alpha,beta,xlnf;
  double T_K2,T_K3,T_K4,T_K5;

  T_K2= T_K  * T_K;
  T_K3= T_K2 * T_K;
  T_K4= T_K3 * T_K;
  T_K5= T_K4 * T_K;

  if(T_K < 1.e3){
    a0= 0.978382;
    a1= 0.000399572;
    a2= -1.73108e-6;
    a3= 1.15363e-9;
    a4= 8.24607e-13;
    a5= -7.65975e-16;
  }else if(T_K <4.e3){
    a0=0.827472;       
    a1= 0.000107697;
    a2= -8.25123e-8;
    a3= 1.92812e-12;
    a4= 5.7192e-15;
    a5= -7.869e-19;
  }else{
    a0= 1.11569;
    a1= -0.000329302;
    a2= 1.01846e-7;     
    a3= -1.46666e-11;
    a4= 1.00764e-15;
    a5= -2.68873e-20;
  }
  
  alpha=a0+a1*T_K+a2*T_K2+a3*T_K3+a4*T_K4+a5*T_K5;

  b0 = 24.0561;
  b1 = 0.00110043;
  b2 = -2.87224e-7;
  b3 = 6.11525e-11;
  b4 = -6.55034e-15;
  b5 = 2.54997e-19;
  beta=b0+b1*T_K+b2*T_K2+b3*T_K3+b4*T_K4+b5*T_K5;

  xlnf=-alpha*log10(1.+xNc_H2/pow(10.,beta));
  fesc=pow(10., xlnf);

  return fesc;
}


double fesc_unity(int ithis, double T_K, double xNc_H2){
  return 1.;
}


/*-----------------------------------------------------------------------
  ! low-density H2 cooling fit by Glover(2015) 
  !
  !     xLmbd_H2 : H2 cooling rate per unit volume
  !     !skeMeMo:I use this subroutine at n > 1.d8
  !-----------------------------------------------------------------------*/
double H2cool_LowD_G15(double xnH,double T_K,double y_H,double y_H2,double y_e,double y_Hp,double y_He){
  double xLmbd_H2;
  double x_o,x_p,a0,a1,a2,a3,a4,a5,a6,a7,a8;
  double T3,xlgT3,xlgT3_2,xlgT3_3,xlgT3_4,xlgT3_5,xlgT3_6,xlgT3_7,xlgT3_8;
  double xLd_pH2H,xlgLd_pH2H,xLd_oH2H,xlgLd_oH2H,xLd_H2H,xlgLd_pH2pH2,xLd_pH2pH2;
  double xlgLd_pH2oH2,xLd_pH2oH2,xlgLd_oH2pH2,xLd_oH2pH2;
  double xlgLd_oH2oH2,xLd_oH2oH2,xLd_H2H2;
  double xlgLd_pH2He,xLd_pH2He,xlgLd_oH2He,xLd_oH2He,xLd_H2He;
  double xlgLd_H2Hp,xLd_H2Hp,xlgLd_H2e,xLd_H2e;

  x_o=0.75;
  x_p=0.25;

  /*****     ortho-H2 excited by collisions with H  ****/
  if(T_K < 1000.){
    a0=-24.330855;
    a1=4.4404496;
    a2=-4.0460989;
    a3=-1.1390725;
    a4=9.8094223;
    a5=8.6273872;
    a6=0.;
    a7=0.;
  }else{
    a0=-24.329086;
    a1=4.6105087;
    a2=-3.9505350;
    a3=12.363818;
    a4=-32.403165;
    a5=48.853562;
    a6=-38.542008;
    a7=12.066770;
  }

  T3=T_K*1.e-3;
  xlgT3   = log10(T3);
  xlgT3_2 = xlgT3   * xlgT3;
  xlgT3_3 = xlgT3_2 * xlgT3;
  xlgT3_4 = xlgT3_3 * xlgT3;
  xlgT3_5 = xlgT3_4 * xlgT3;
  xlgT3_6 = xlgT3_5 * xlgT3;
  xlgT3_7 = xlgT3_6 * xlgT3;
  xlgT3_8 = xlgT3_7 * xlgT3;
  
  if(T_K < 100.)
    xLd_oH2H=5.09e-27*sqrt(T3)*exp(-852.5/T_K);
  else{
    xlgLd_oH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3
      +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7;
    xLd_oH2H=pow(10.,xlgLd_oH2H);
  }

  /*******     para-H2 excited by collisions with H   ******/
  if(T_K < 1000.){
    a0=-24.216387;
    a1=3.3237480;
    a2=-11.642384;
    a3=-35.553366;
    a4=-35.105689;
    a5=-10.922078;
    a6=0.0;
    a7=0.0;
  }else{
    a0=-24.216387;
    a1=4.2046488;
    a2=-1.3155285;
    a3=-1.6552763;
    a4=4.1780102;
    a5=-0.56949697;
    a6=-3.3824407;
    a7=1.0904027;
  }

  if(T_K < 100.){
    xLd_pH2H=8.16e-26*sqrt(T3)*exp(-509.85/T_K);
  }else{
       xlgLd_pH2H=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
	 +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6+a7*xlgT3_7;
       xLd_pH2H=pow(10.,xlgLd_pH2H);
  }


  xLd_H2H=(x_o/(x_o+x_p))*xLd_oH2H +(x_p/(x_o+x_p))*xLd_pH2H;


  /***********     para-H2 with para-H2     ***************/
  a0=-23.889798;
  a1=1.8550774;
  a2=-0.55593388;
  a3=0.28429361;
  a4=-0.20581113;
  a5=0.13112378;
  xlgLd_pH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_pH2pH2=pow(10.,xlgLd_pH2pH2);

  /***********     para-H2 with ortho-H2   ****************/
  a0=-23.748534;
  a1=1.76676480;
  a2=-0.58634325;
  a3=0.31074159;
  a4=-0.17455629;
  a5=0.18530758;
  xlgLd_pH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_pH2oH2=pow(10.,xlgLd_pH2oH2);

  /************     ortho-H2 with para-H2   **************/
  a0=-24.126177;
  a1=2.3258217;
  a2=-1.0082491;
  a3=0.54823768;
  a4=-0.33679759;
  a5=0.20771406;
  xlgLd_oH2pH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_oH2pH2=pow(10.,xlgLd_oH2pH2);

  /************     ortho-H2 with ortho-H2   **************/
  a0=-24.020047;
  a1=2.2687566;
  a2=-1.0200304;
  a3=0.83561432;
  a4=-0.40772247;
  a5=0.096025713;
  xlgLd_oH2oH2=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_oH2oH2=pow(10.,xlgLd_oH2oH2);


  xLd_H2H2=x_p*x_p*xLd_pH2pH2+x_p*x_o*xLd_pH2oH2 
    +x_o*x_p*xLd_oH2pH2+x_o*x_o*xLd_oH2oH2;

  /***************     para-H2 with He      ****************/
  a0=-23.489029;
  a1=1.8210825;
  a2=-0.59110559;
  a3=0.42280623;
  a4=-0.30171138;
  a5=0.12872839;
  xlgLd_pH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_pH2He=pow(10.,xlgLd_pH2He);

  /****************     ortho-H2 with He    ****************/
  a0=-23.7749;
  a1=2.40654;
  a2=-1.23449;
  a3=0.739874;
  a4=-0.258940;
  a5=0.120573;
  xlgLd_oH2He=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
    +a4*xlgT3_4+a5*xlgT3_5;
  xLd_oH2He=pow(10.,xlgLd_oH2He);


  xLd_H2He=(x_o/(x_o+x_p))*xLd_oH2He+(x_p/(x_o+x_p))*xLd_pH2He;


  /*****************     H2 with proton     ****************/
  if(10.<T_K){
    a0=-22.089523;
    a1=1.5714711;
    a2=0.015391166;
    a3=-0.23619985;
    a4=-0.51002221;
    a5=0.32168730;
       xlgLd_H2Hp=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
	 +a4*xlgT3_4+a5*xlgT3_5;
       xLd_H2Hp=pow(10.,xlgLd_H2Hp);
  }else{
      xLd_H2Hp=0.;
  }

   
  /******************     H2 with e (KS MODIFIED) **********/
  if(100.<T_K && T_K<500.){
    a0=-21.928796;
    a1=16.815730;
    a2=96.743155;
    a3=343.19180;
    a4=734.71651;
    a5=983.67576;
    a6=801.81247;
    a7=364.14446;
    a8=70.609154;
    xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
      +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 
      +a7*xlgT3_7+a8*xlgT3_8;
    xLd_H2e=pow(10.,xlgLd_H2e);
  }else if(500.<T_K){
    a0=-22.92118;
    a1=1.6802758;
    a2=0.93310622;
    a3=4.0406627;
    a4=-4.7274036;
    a5=-8.8077017;
    a6=8.9167183;
    a7=6.4380698;
    a8=-6.3701156;
    xlgLd_H2e=a0+a1*xlgT3+a2*xlgT3_2+a3*xlgT3_3 
      +a4*xlgT3_4+a5*xlgT3_5+a6*xlgT3_6 
      +a7*xlgT3_7+a8*xlgT3_8;
    xLd_H2e=pow(10.,xlgLd_H2e);
  }else{
    xLd_H2e=0.;
  }


  xLmbd_H2=xLd_H2H*y_H+xLd_H2H2*y_H2+xLd_H2He*y_He 
    +xLd_H2Hp*y_Hp+xLd_H2e*y_e;

  xLmbd_H2 = xLmbd_H2 * xnH * xnH * y_H2;
  return xLmbd_H2;
}



#ifdef HD

double HD_cooling_O00(int ithis, double nH, double T_gas, double yHI, double yHD)
{
  double cool_rate_HD;
  double rho,n_HI,n_HD;
  double mu;
  
  rho = nH*PROTONMASS / XH;
  //mu = SphP[ithis].molecular_weight;
  mu = 1.;

  n_HI = yHI * nH;
  n_HD = yHD * nH;

// tabulate coefs. /////////////////////////////////////////////////////////////////////////////////
  double g0,g1,g2,g3;
  double A10d,A21d,A32d,A20d,A31d;
  double nu10d,nu21d,nu32d,nu20d,nu31d;
  double DE10d,DE21d,DE32d,DE20d,DE31d;
  double DT10d,DT21d,DT32d,DT20d,DT31d;
  double C10d,C21d,C32d,C20d,C31d;
  double C01d,C12d,C23d,C02d,C13d;
  double Q10d,Q21d,Q32d,Q20d,Q31d;

  g0 = 1.0;
  g1 = 3.0;
  g2 = 5.0;
  g3 = 7.0;

  nu10d = 2.67e12;
  nu21d = 5.32e12;
  nu32d = 7.93e12;
  nu20d = nu21d+nu10d;
  nu31d = nu32d+nu21d;
  
  DE10d = PLANCK*nu10d;
  DE21d = PLANCK*nu21d;
  DE32d = PLANCK*nu32d;
  DE20d = PLANCK*nu20d;
  DE31d = PLANCK*nu31d;

  DT10d = DE10d/BOLTZMANN;
  DT21d = DE21d/BOLTZMANN;
  DT32d = DE32d/BOLTZMANN;
  DT20d = DE20d/BOLTZMANN;
  DT31d = DE31d/BOLTZMANN;

  A10d = 5.12e-8;
  A21d = 4.86e-7;
  A32d = 1.72e-6;
  A20d = 7.05e-12;
  A31d = 1.15e-10;

// collision with HI
  C10d = (4.4e-12+3.6e-13*pow(T_gas,0.77))*n_HI/1.27;
  C21d = (4.1e-12+2.0e-13*pow(T_gas,0.92))*n_HI/1.27;
  C32d = (2.4e-12+8.7e-14*pow(T_gas,1.03))*n_HI/1.27;
  C20d = (3.4e-13+1.1e-14*pow(T_gas,1.12))*n_HI/1.27;
  C31d = (3.2e-13+1.3e-15*pow(T_gas,1.47))*n_HI/1.27;

  C01d = (g1/g0)*C10d*exp(-DT10d/T_gas);
  C12d = (g2/g1)*C21d*exp(-DT21d/T_gas);
  C23d = (g3/g2)*C32d*exp(-DT32d/T_gas);
  C02d = (g2/g0)*C20d*exp(-DT20d/T_gas);
  C13d = (g3/g1)*C31d*exp(-DT31d/T_gas);

  Q10d = compute_Q_bg(DT10d);
  Q21d = compute_Q_bg(DT21d);
  Q32d = compute_Q_bg(DT32d);
  Q20d = compute_Q_bg(DT20d);
  Q31d = compute_Q_bg(DT31d);

// compute population //////////////////////////////////////////////////////////////////////////////
  double R10d,R21d,R32d,R20d,R31d,R30d;
  double R01d,R12d,R23d,R02d,R13d,R03d;
  double R0d,R1d,R2d;
  double an0d,an1d,an2d,an3d,antotd;
  double f0d,f1d,f2d,f3d;

  R10d = A10d*(1.0 + Q10d) + C10d;
  R21d = A21d*(1.0 + Q21d) + C21d;
  R32d = A32d*(1.0 + Q32d) + C32d;
  R20d = A20d*(1.0 + Q20d) + C20d;
  R31d = A31d*(1.0 + Q31d) + C31d;
  R30d = 0.0;

  R01d = (g1/g0)*A10d*Q10d + C01d;
  R12d = (g2/g1)*A21d*Q21d + C12d;
  R23d = (g3/g2)*A32d*Q32d + C23d;
  R02d = (g2/g0)*A20d*Q20d + C02d;
  R13d = (g3/g1)*A31d*Q31d + C13d;
  R03d = 0.0;

  R0d = R01d+R02d+R03d;
  R1d = R10d+R12d+R13d;
  R2d = R20d+R21d+R23d;

  an0d = -(R10d*R21d*R32d + R30d*R2d *R1d  + R20d*R31d*R12d -
           R12d*R21d*R30d + R2d *R31d*R10d + R32d*R20d*R1d );
  an1d = - R0d *R21d*R32d + R20d*R31d*R02d - R30d*R2d *R01d -
           R02d*R21d*R30d - R2d *R31d*R0d  - R32d*R20d*R01d;
  an2d = -(R0d *R1d *R32d + R10d*R31d*R02d + R30d*R12d*R01d +
           R02d*R1d *R30d + R12d*R31d*R0d  - R32d*R10d*R01d);
  an3d = - R0d *R1d * R2d + R10d*R21d*R02d + R20d*R12d*R01d +
           R02d*R1d *R20d + R12d*R21d*R0d  + R2d *R10d*R01d;
  antotd = an0d+an1d+an2d+an3d;

  f0d = an0d/antotd;
  f1d = an1d/antotd;
  f2d = an2d/antotd;
  f3d = an3d/antotd;

// radiative transfer //////////////////////////////////////////////////////////////////////////////
  double tau10d,tau21d,tau32d,tau20d,tau31d,tau30d;
  double esc10d,esc21d,esc32d,esc20d,esc31d,esc30d;
  double l_sh,Dv_D,t_col,Dv_bulk,dvdz;

  l_sh = sqrt(M_PI*BOLTZMANN*T_gas/(GRAVITY*mu*PROTONMASS*rho));
  t_col = compute_freefall_time(rho);
  Dv_bulk = l_sh/3.0/t_col;
  Dv_D = sqrt(2.0*BOLTZMANN*T_gas/(3.0*PROTONMASS));
  if(Dv_D>Dv_bulk) dvdz = Dv_D/l_sh;
  else            dvdz = Dv_bulk/l_sh;

  tau10d = (A10d/8.0/M_PI)*pow(CLIGHT/nu10d,3.0)
           *(f0d*g1/g0-f1d)*n_HD/dvdz;
  tau21d = (A21d/8.0/M_PI)*pow(CLIGHT/nu21d,3.0)
           *(f1d*g2/g1-f2d)*n_HD/dvdz;
  tau32d = (A32d/8.0/M_PI)*pow(CLIGHT/nu32d,3.0)
           *(f2d*g3/g2-f3d)*n_HD/dvdz;
  tau20d = (A20d/8.0/M_PI)*pow(CLIGHT/nu20d,3.0)
           *(f0d*g2/g0-f2d)*n_HD/dvdz;
  tau31d = (A31d/8.0/M_PI)*pow(CLIGHT/nu31d,3.0)
           *(f1d*g3/g1-f3d)*n_HD/dvdz;

  esc10d=compute_b_esc(tau10d, ithis);
  esc21d=compute_b_esc(tau21d, ithis);
  esc32d=compute_b_esc(tau32d, ithis);
  esc20d=compute_b_esc(tau20d, ithis);
  esc31d=compute_b_esc(tau31d, ithis);

// cooling rate ////////////////////////////////////////////////////////////////////////////////////
  double L10d,L21d,L32d,L20d,L31d,L30d;
  double S10d,S21d,S32d,S20d,S31d,S30d;
  double X10d,X21d,X32d,X20d,X31d,X30d;

  if(f1d!=0.0) {
    S10d = 1.0/(g1/g0*f0d/f1d-1.0);
    X10d = esc10d*(1.0-Q10d/S10d);
  } else {
    X10d = 0.0;
  }
  if(f2d!=0.0) {
    S21d = 1.0/(g2/g1*f1d/f2d-1.0);
    X21d = esc21d*(1.0-Q21d/S21d);
    S20d = 1.0/(g2/g0*f0d/f2d-1.0);
    X20d = esc20d*(1.0-Q20d/S20d);
  } else {
    X21d = 0.0;
    X20d = 0.0;
  }
  if(f3d!=0.0) {
    S32d = 1.0/(g3/g2*f2d/f3d-1.0);
    X32d = esc32d*(1.0-Q32d/S32d);
    S31d = 1.0/(g3/g1*f1d/f3d-1.0);
    X31d = esc31d*(1.0-Q31d/S31d);
  } else {
    X32d = 0.0;
    X31d = 0.0;
  }
   

  L10d = X10d * f1d * A10d * DE10d; 
  L21d = X21d * f2d * A21d * DE21d;
  L32d = X32d * f3d * A32d * DE32d;
  L20d = X20d * f2d * A20d * DE20d;
  L31d = X31d * f3d * A31d * DE31d;    

  cool_rate_HD = L10d+L21d+L32d+L20d+L31d;

  return cool_rate_HD;
}

double HD_cooling_L05(int ithis, double nH, double T_gas, double yHI, double yHDI)
{
  /* ------- HD cooling by Lipovka 2005 */
  double cool_rate_HD;
  int lll, mmm;
  double n_HI,log_nH,log_nH_ref, log_Tk, logW;
  double DD[5][5];

  n_HI = yHI * nH;
  log_nH = log10(n_HI);
  log_Tk = log10(T_gas);

  DD[0][0]=-42.57688; DD[0][1]=0.92433; DD[0][2]=0.54962; DD[0][3]=-0.07676; DD[0][4]=0.00275;
  DD[1][0]=21.93385; DD[1][1]=0.77952; DD[1][2]=-1.06447; DD[1][3]=0.11864; DD[1][4]=-0.00366; 
  DD[2][0]=-10.19097; DD[2][1]=-0.54263; DD[2][2]=0.62343; DD[2][3]=-0.07366; DD[2][4]=0.002514; 
  DD[3][0]=2.19906; DD[3][1]=0.11711; DD[3][2]=-0.13768; DD[3][3]=0.01759; DD[3][4]=-0.000666317; 
  DD[4][0]=-0.17334; DD[4][1]=-0.00835; DD[4][2]=0.0106; DD[4][3]=-0.001482; DD[4][4]=0.000061926; 

  if(n_HI<1.0)
    log_nH_ref = 0.0;
  else if(n_HI>1.0e8)
    log_nH_ref = 8.0;
  else 
    log_nH_ref = log_nH;

  logW=0.0;
  for(lll=0;lll<5;lll++)
    for(mmm=0;mmm<5;mmm++)
      logW += DD[lll][mmm] * pow(log_Tk,lll) * pow(log_nH_ref,mmm);

  if(n_HI<1.0)
    cool_rate_HD = pow(10.0, logW)*(n_HI/1.0);
  else
    cool_rate_HD = pow(10.0, logW);

  return cool_rate_HD;

}


double HD_cooling_GP12(int i, double nH, double T_gas, double yHI, double yHDI)
{
  double log10T = log10(T_gas);
  double cool_rate_LOW_HD, cool_rate_LTE_HD, cool_rate_HD;
  double rate;

  /* HD : Low-density cooling function (Lipovka et al. 2005) */
  /* ----- HD-H collision cooling : 1<n_H[cm^-3]<1E8 & 1E2<T_gas[K]<1E4 */
  /*       >> eq(4) & Table(1) in Lipovka et al. (2005) */

  double n_HI; 
  double log10n_HI1, log10n_HI2, log10n_HI3, log10n_HI4;
  n_HI  = yHI * nH;
  log10n_HI1 = log10(n_HI);
  log10n_HI2 = log10n_HI1*log10n_HI1;
  log10n_HI3 = log10n_HI2*log10n_HI1;
  log10n_HI4 = log10n_HI3*log10n_HI1;
  
  if(nH > 1. && nH < 1.e7){
    cool_rate_LOW_HD 
      = pow(10.0,
	    ( -42.57688 +0.92433*log10n_HI1 +0.54962*log10n_HI2 -0.07676 *log10n_HI3 +0.00275    *log10n_HI4 )
	    + ( +21.93385 +0.77952*log10n_HI1 -1.06447*log10n_HI2 +0.11864 *log10n_HI3 -0.00366    *log10n_HI4 ) * log10T
	    + ( -10.19097 -0.54263*log10n_HI1 +0.62343*log10n_HI2 -0.07366 *log10n_HI3 +0.002514   *log10n_HI4 ) * pow(log10T, 2.)
	    + ( + 2.19906 +0.11711*log10n_HI1 -0.13768*log10n_HI2 +0.01759 *log10n_HI3 -0.000666317*log10n_HI4 ) * pow(log10T, 3.)
	    + ( - 0.17334 -0.00835*log10n_HI1 +0.0106 *log10n_HI2 -0.001482*log10n_HI3 +0.000061926*log10n_HI4 ) * pow(log10T, 4.)
	    );  //=[erg sec^-1]
  }else
    cool_rate_LOW_HD = 0.;

  /* HD : LTE(High-density) cooling function */
  /*       >> eq.(7) & Tanle(1) in Coppola et al. (2011) */
  
  cool_rate_LTE_HD 
    = pow(10.0, 
	  (- 55.5725
	   + 56.649   * log10T
	   - 37.9102  * pow(log10T, 2.)
	   + 12.698   * pow(log10T, 3.)
	   - 2.02424  * pow(log10T, 4.)
	   + 0.122393 * pow(log10T, 5.)
	   ));  //=[erg sec^-1]


  rate = dmin(nH/1.e7, 1.);
  cool_rate_HD = dmin(cool_rate_LTE_HD, cool_rate_LOW_HD*(1-rate) + cool_rate_LTE_HD*rate);  //=[erg sec^-1] 


  /* ---------- Heating by CMB */  
  /* --------- Cooling by HD : Flower et al. 2000 */
  double aa,bb,omeg,phi,c1,c2,d1,d2,dens,y,x,w,dwdx;
  double T_CMB, heat_rate, edot_HD;
  aa   = -26.2982e0;
  bb   = - 0.215807e0;
  omeg =   2.03657e0;
  phi  =   4.63258e0;
  c1   =   0.283978e0;
  c2   = - 1.27333e0;
  d1   = - 2.08189e0;
  d2   =   4.66288e0;

  dens = yHI * nH;

  /* Limit of this fitting (Flower et al 2000)*/
  if (dens < 1.e0) dens=1.e0;
  if (dens > 1.e8) dens=1.e8;
  y = log10(dens);

  heat_rate = 0;
  
  x = log10(All.T_CMB);
  w = 0.5 * y + aa * pow(x,bb) - sqrt( 0.25 * y*y + (c1*y+c2) * sin(omeg*x+phi) + (d1*y+d2) );

  heat_rate = pow(10.0,w);  //=[erg sec^-1]

  /* ---------- Cooling Rate by HD (-CMB) */
  edot_HD = - yHDI * nH * (cool_rate_HD - heat_rate);
  if (T_gas < 10.0 || T_gas > 3000.0)
    edot_HD = 0.0;

#ifdef ID_CHEM_DEBUG
  if(FlagPrint){
    printf("yHDI=%g nH=%g cool_HD=%g low=%g LTE=%g rate=%g heat=%g\n",
	   yHDI, nH, cool_rate_HD, cool_rate_LOW_HD, cool_rate_LTE_HD, rate, heat_rate);
  }
#endif
  
  /* ----- HD cooling rate per volume [erg/s/cm^3] */
  return edot_HD;
}

double fesc_HD_fit(double T_K, double xNc_HD){
  double fesc;
  double alpha,Ncr,lgNc;
  double a0,a1,a2,a3,a4,a5;
  double b0,b1,b2,b3,b4,b5;
  double lgT1,lgT2,lgT3,lgT4,lgT5;
  double T_K2,T_K3,T_K4,T_K5;

  if(T_K < 1.e3){
    a0 =  1.02339;
    a1 = -9.04450e-4;
    a2 =  3.26141e-6;
    a3 = -5.73674e-9;
    a4 =  4.98247e-12;
    a5 = -1.68078e-15;
  }else{
    a0 =  8.94334e-1;
    a1 =  6.21583e-5;
    a2 = -1.47443e-9;
    a3 = -9.07612e-12;
    a4 =  3.31028e-15;
    a5 = -3.79804e-19;
  }

  T_K2= T_K  * T_K;
  T_K3= T_K2 * T_K;
  T_K4= T_K3 * T_K;
  T_K5= T_K4 * T_K;

  alpha = a0 + a1*T_K + a2*T_K2 + a3*T_K3  + a4*T_K4 + a5*T_K5;

  lgT1 = log10(T_K);
  lgT2 = lgT1*lgT1;
  lgT3 = lgT1*lgT2;
  lgT4 = lgT1*lgT3;
  lgT5 = lgT1*lgT4;

  b0 =  1.88924e1;
  b1 = -1.48020e0;
  b2 =  8.59050e-1;
  b3 =  4.38401e-1;
  b4 = -3.05691e-1;
  b5 =  4.64115e-2;

  lgNc = b0 + b1*lgT1 + b2*lgT2 + b3*lgT3 + b4*lgT4 + b5*lgT5;
  Ncr = pow(10., lgNc);
            
  fesc = 1. * pow(1.+xNc_HD/Ncr,-alpha);
            
  return fesc;
}

#endif



double compute_b_esc(double tau_line, int ithis)
{
  double b_esc;

  if(tau_line<0.0) 
    b_esc = 1.0;
  else if (tau_line<1.0e-5) 
    b_esc = exp(-SphP[ithis].tau_cont);
  else
    b_esc = (1.0 - exp(-tau_line)) / tau_line * exp(-SphP[ithis].tau_cont);

  return b_esc;
}


double compute_Q_bg(double T_nu)
{
  double Q_bg,Q_bg_CMB;
  double xx;

  if(!CMB_BG)
    Q_bg_CMB = 0.0;
  else {
    xx = T_nu/All.T_CMB;
    if(xx>100.0)
      Q_bg_CMB = 0.0;
    else
      Q_bg_CMB = 1.0/(exp(xx)-1.0);
  }

  Q_bg = Q_bg_CMB;

  return Q_bg;
}


int H2line_tabulate(void)
{

  int v,J, JI, JF; 
  
  double  we, wexe, weye, weze;
  double  Be, De, alpha_e, beta_e;
//  double Bv, Dv;
  
  double hc, hc2, cminverse_to_eV, eV_to_lambda;
  
  
  hc = PLANCK * CLIGHT ;// cm erg
  hc2 = hc*hc;

  cminverse_to_eV = 1.23985e-4; // 1.0/8065.47
  eV_to_lambda = 12401.0 ;// 911.5069 (Ly-cont) * 13.605 (Rydberg) 
  
  
  //Einstein A-coefficients for v=0, 1, 2
  // A(vI,JI,vF,JF)
  
  //JI = 0
  A[1][0][0][0]=8.54e-7; //?
  A[2][0][0][0]=3.47e-7; //?
  A[2][0][1][0]=1.29e-6; //?

  A[1][0][0][2]=8.54e-7;
  A[2][0][0][2]=3.47e-7;
  A[2][0][1][2]=1.29e-6;


  //JI = 1
  //vI=1
  A[1][1][0][1]=4.29e-7;
  A[1][1][0][3]=4.23e-7;
  
  //vI=2
  A[2][1][0][1]=1.94e-7;
// A[2][1][0][3] = ***;
  A[2][1][0][3]=1.61e-7;
//  A[2][1][1][1]=6.37e-6;
  A[2][1][1][1]=6.37e-7;
//  A[2][1][1][3]=6.40e-6;
  A[2][1][1][3]=6.40e-7;



  //JI = 2
  //vI=0
  A[0][2][0][0]=2.94e-11;
    
  //vI=1
  A[1][2][0][0]=2.53e-7;
  A[1][2][0][2]=3.03e-7;
  A[1][2][0][4]=2.90e-7;
  A[1][2][1][0]=2.79e-11;
  
  //vI=2
  A[2][2][0][0]=1.27e-7;
  A[2][2][0][2]=1.38e-7;
  A[2][2][0][4]=1.03e-7;
  A[2][2][1][0]=3.68e-7;
  A[2][2][1][2]=4.50e-7;
  A[2][2][1][4]=4.41e-7;

//  A[2][2][2][0]=4.41e-7;
  A[2][2][2][0]=2.56e-11;
  

  //JI = 3
  //vI=0
  A[0][3][0][1]=4.76e-10;
 
  //vI=1
  A[1][3][0][1]=3.47e-7;
  A[1][3][0][3]=2.78e-7;
  A[1][3][0][5]=2.09e-7;
  A[1][3][1][1]=4.50e-10;
 
  //vI=2
  A[2][3][0][1]=1.90e-7;
  A[2][3][0][3]=1.29e-7;
  A[2][3][0][5]=6.98e-8;
  A[2][3][1][1]=4.98e-7;
  A[2][3][1][3]=4.12e-7;
  A[2][3][1][5]=3.18e-7;
 
  A[2][3][2][1]=4.12e-10;


  //JI = 4
  //vI=0
  A[0][4][0][2]=2.76e-9;

  //vI=1
  A[1][4][0][2]=3.98e-7;
  A[1][4][0][4]=2.65e-7;
  A[1][4][0][6]=1.50e-7;
  A[1][4][1][2]=2.59e-9;

  //vI=2
  A[2][4][0][2]=2.38e-7;
  A[2][4][0][4]=1.25e-7;
  A[2][4][0][6]=4.72e-8;
  A[2][4][1][2]=5.60e-7;
  A[2][4][1][4]=3.91e-7;
  A[2][4][1][6]=2.28e-7;

//  A[2][4][2][2]=2.37e-10;
  A[2][4][2][2]=2.37e-9;


  //JI = 5
  //vI=0
  A[0][5][0][3]=9.84e-9;

  //vI=1
  A[1][5][0][3]=4.21e-7;
  A[1][5][0][5]=2.55e-7;
  A[1][5][0][7]=1.06e-7;
  A[1][5][1][3]=9.21e-9;

  //vI=2
  A[2][5][0][3]=2.77e-7;
  A[2][5][0][5]=1.23e-7;
  A[2][5][0][7]=3.15e-8;
  A[2][5][1][3]=5.77e-7;
  A[2][5][1][5]=3.74e-7;
  A[2][5][1][7]=1.62e-7;

//  A[2][5][2][3]=8.37e-10;
  A[2][5][2][3]=8.37e-9;


  //JI = 6
  //vI=0
  A[0][6][0][4]=2.64e-8;

  //vI=1
  A[1][6][0][4]=4.19e-7;
  A[1][6][0][6]=2.45e-7;
  A[1][6][0][8]=7.38e-8;
  A[1][6][1][4]=2.46e-8;

  //vI=2
  A[2][6][0][4]=3.07e-7;
  A[2][6][0][6]=1.21e-7;
  A[2][6][0][8]=2.05e-8;
  A[2][6][1][4]=5.57e-7;
  A[2][6][1][6]=3.58e-7;
  A[2][6][1][8]=1.12e-7;

  A[2][6][2][4]=2.22e-8;


  //JI = 7
  //vI=0
  A[0][7][0][5]=5.88e-8;

  //vI=1
  A[1][7][0][5]=3.96e-7;
  A[1][7][0][7]=2.34e-7;
//  A[1][7][0][9]=4.98e-7;
  A[1][7][0][9]=4.98e-8;
  A[1][7][1][5]=5.44e-8;

  //vI=2
  A[2][7][0][5]=3.28e-7;
  A[2][7][0][7]=1.20e-7;
  A[2][7][0][9]=1.31e-8;
  A[2][7][1][5]=5.05e-7;
  A[2][7][1][7]=3.40e-7;
//  A[2][7][1][9]=7.50e-7;
  A[2][7][1][9]=7.50e-8;

  //CHECK THIS VALUE!
  A[2][7][2][5]=4.88e-8;



  //JI = 8
  //vI=0
  A[0][8][0][6]=1.14e-7;

  //vI=1
  A[1][8][0][6]=3.54e-7;
  A[1][8][0][8]=2.23e-7;
//  A[1][8][0][10]=3.27e-7;
  A[1][8][0][10]=3.27e-8;
  A[1][8][1][6]=1.05e-7;

  //vI=2
  A[2][8][0][6]=3.39e-7;
  A[2][8][0][8]=1.18e-7;
  A[2][8][0][10]=8.07e-9;
  A[2][8][1][6]=4.30e-7;
  A[2][8][1][8]=3.22e-7;
//  A[2][8][1][10]=4.88e-7;
  A[2][8][1][10]=4.88e-8;

  A[2][8][2][6]=9.33e-8;


  //JI = 9
  //vI=0
  A[0][9][0][7]=2.00e-7;

  //vI=1
  A[1][9][0][7]=2.98e-7;
  A[1][9][0][9]=2.12e-7;
  A[1][9][0][11]=2.07e-8;
  A[1][9][1][7]=1.82e-7;

  //vI=2
  A[2][9][0][7]=3.40e-7;
  A[2][9][0][9]=1.17e-7;
  A[2][9][0][11]=4.83e-9;
  A[2][9][1][7]=3.38e-7;
  A[2][9][1][9]=3.03e-7;
  A[2][9][1][11]=3.06e-8;

  A[2][9][2][7]=1.61e-7;


  //JI = 10
  //vI=0
  A[0][10][0][8]=3.24e-7;

  //vI=1
  A[1][10][0][8]=2.34e-7;
  A[1][10][0][10]=1.99e-7;
  A[1][10][0][12]=1.27e-8;
  A[1][10][1][8]=2.92e-7;

  //vI=2
  A[2][10][0][8]=3.30e-7;
  A[2][10][0][10]=1.15e-7;
  A[2][10][0][12]=2.79e-9;
  A[2][10][1][8]=2.41e-7;
  A[2][10][1][10]=2.84e-7;
  A[2][10][1][12]=1.85e-8;

  A[2][10][2][8]=2.55e-7;




  //JI = 11
  //vI=0
  A[0][11][0][9]=4.90e-7;

  //vI=1
  A[1][11][0][9]=1.68e-7;
  A[1][11][0][11]=1.87e-7;
  A[1][11][0][13]=7.44e-9;
  A[1][11][1][9]=4.38e-7;

  //vI=2
  A[2][11][0][9]=3.12e-7;
  A[2][11][0][11]=1.13e-7;
  A[2][11][0][13]=1.55e-9;
  A[2][11][1][9]=1.49e-7;
  A[2][11][1][11]=2.63e-7;
  A[2][11][1][13]=1.07e-8;

  A[2][11][2][9]=3.79e-7;



  //JI = 12
  //vI=0
  A[0][12][0][10]=7.03e-7;

  //vI=1
  A[1][12][0][10]=1.05e-7;
  A[1][12][0][12]=1.74e-7;
  A[1][12][0][14]=4.16e-9;
  A[1][12][1][10]=6.21e-7;

  //vI=2
  A[2][12][0][10]=2.85e-7;
  A[2][12][0][12]=1.11e-7;
  A[2][12][0][14]=8.27e-10;
  A[2][12][1][10]=7.20e-8;
  A[2][12][1][12]=2.42e-7;
  A[2][12][1][14]=5.84e-9;

  A[2][12][2][10]=5.32e-7;


  //JI = 13
  //vI=0
  A[0][13][0][11]=9.64e-7;

  //vI=1
  A[1][13][0][11]=5.30e-8;
  A[1][13][0][13]=1.61e-7;
  A[1][13][0][15]=2.19e-9;
  A[1][13][1][11]=8.42e-7;

  //vI=2
  A[2][13][0][11]=2.53e-7;
  A[2][13][0][13]=1.08e-7;
  A[2][13][0][15]=4.21e-10;
  A[2][13][1][11]=1.96e-8;
  A[2][13][1][13]=2.21e-7;
  A[2][13][1][15]=3.00e-9;

  A[2][13][2][11]=7.13e-7;




  //JI = 14
  //vI=0
  A[0][14][0][12]=1.27e-6;

  //vI=1
  A[1][14][0][12]=1.65e-8;
  A[1][14][0][14]=1.47e-7;
  A[1][14][0][16]=1.08e-9;
  A[1][14][1][12]=1.10e-6;

  //vI=2
  A[2][14][0][12]=2.16e-7;
  A[2][14][0][14]=1.05e-7;
  A[2][14][0][16]=2.04e-10;
  A[2][14][1][12]=1.49e-11;
  A[2][14][1][14]=2.01e-7;
  A[2][14][1][16]=1.43e-9;

  A[2][14][2][12]=9.18e-7;


  //JI = 15
  //vI=0
  A[0][15][0][13]=1.62e-6;

  //vI=1
  A[1][15][0][13]=4.26e-10;
  A[1][15][0][15]=1.34e-7;
  A[1][15][0][17]=4.87e-10;
  A[1][15][1][13]=1.38e-6;

  //vI=2
  A[2][15][0][13]=1.78e-7;
  A[2][15][0][15]=1.02e-7;
  A[2][15][0][17]=9.45e-11;
  A[2][15][1][13]=1.91e-8;
  A[2][15][1][15]=1.80e-7;
  A[2][15][1][17]=6.17e-10;

  A[2][15][2][13]=1.14e-6;


  //JI = 16
  //vI=0
  A[0][16][0][14]=2.00e-6;

  //vI=1
//  A[1][16][0][14]=8.28e-9;
  A[1][16][0][14]=8.38e-9;
  A[1][16][0][16]=1.21e-7;
  A[1][16][0][18]=1.96e-10;
  A[1][16][1][14]=1.69e-6;

  //vI=2
  A[2][16][0][14]=1.39e-7;
  A[2][16][0][16]=9.88e-8;
  A[2][16][0][18]=4.23e-11;
  A[2][16][1][14]=8.07e-8;
  A[2][16][1][16]=1.60e-7;
  A[2][16][1][18]=2.33e-10;
  
  A[2][16][2][14]=1.37e-6;
  
  
/*  
  //JI = 16
  //vI=0
  A[0][16][0][14]=2.41e-6;
  
  //vI=1
  A[1][16][0][14]=1.04e-7;
  A[1][16][0][16]=9.61e-8;
  A[1][16][0][18]=1.82e-11;
  A[1][16][1][14]=2.32e-6;
  
  //vI=2
  A[2][16][0][14]=1.78e-7;
  A[2][16][0][16]=1.02e-7;
  A[2][16][0][18]=9.45e-11;
  A[2][16][1][14]=1.91e-8;
  A[2][16][1][16]=1.80e-7;
  A[2][16][1][18]=6.17e-10;

  A[2][16][2][14]=1.14e-6;
*/



  //JI = 17
  //vI=0
  A[0][17][0][15]=2.41e-6;

  //vI=1
  A[1][17][0][15]=4.27e-8;
  A[1][17][0][17]=1.08e-7;
  A[1][17][0][19]=6.71e-11;
    A[1][17][1][15]=2.00e-6;

  //vI=2
  A[2][17][0][15]=1.01e-7;
  A[2][17][0][17]=9.50e-8;
  A[2][17][0][19]=1.91e-11;
  A[2][17][1][15]=1.86e-7;
  A[2][17][1][17]=1.41e-7;
  A[2][17][1][19]=7.32e-11;

  A[2][17][2][15]=1.61e-6;





  //JI = 18
  //vI=0
  A[0][18][0][16]=2.83e-6;

  //vI=1
  A[1][18][0][16]=1.04e-7;
  A[1][18][0][18]=9.61e-8;

  A[1][18][1][16]=2.32e-6;

  //vI=2
  A[2][18][0][16]=6.80e-8;
  A[2][18][0][18]=9.07e-8;

  A[2][18][1][16]=3.35e-7;
  A[2][18][1][18]=1.22e-7;

  A[2][18][2][16]=1.84e-6;



  //JI = 19
  //vI=0
  A[0][19][0][17]=3.26e-6;

  //vI=1
  A[1][19][0][17]=1.93e-7;
  A[1][19][0][19]=8.43e-8;

  A[1][19][1][17]=2.64e-6;

  //vI=2
  A[2][19][0][17]=3.98e-8;
  A[2][19][0][19]=8.62e-8;

  A[2][19][1][17]=5.24e-7;
  A[2][19][1][19]=1.05e-7;

  A[2][19][2][17]=2.05e-6;

 
  //statistical weight
  for(J=0;J<=18;J=J+2) //Para H2
    g_weight[J]=1;

  for(J=1;J<=19;J=J+2) //Ortho H2
    g_weight[J]=3;

/*
  //energy spacing
  we=4395.2;
  wexe=117.99;
  weye=0.29;
  weze=0.046;
  Be=60.809;
  De=0.04656;
  alpha_e=2.993;
  beta_e=-1.799e-3;

  for(v=0;v<=2;v++)
    for(J=0;J<=19;J++){
      Bv = Be - alpha_e * (v+1.0/2.0);
      Dv = De +  beta_e * (v+1.0/2.0);
      
      E_vJ[v][J]=we*(v+1.0/2.0) - wexe*pow(v+1.0/2.0, 2) + weye*pow(v+1.0/2.0, 3) - weze*pow(v+1.0/2.0, 4) 
	+ Bv * J*(J+1) - Dv * J*J*(J+1)*(J+1);  // in cm^{-1}
      
      E_vJ_erg[v][J] = E_vJ[v][J] * hc;  //erg
      E_vJ_eV[v][J]  = E_vJ[v][J] * cminverse_to_eV; //eV
      E_vJ_K[v][J]   = E_vJ_eV[v][J] * eV_to_K; //K
    }
*/

  //energy spacing (Borysow, Frommhold & Moraldi 1989)
  double Ev0,Ev1,Ev2,Ev3,Ev4;
  double Bv0,Bv1,Bv2,Bv3,Bv4;
  double Dv0,Dv1,Dv2,Dv3,Dv4;
  double Fv0,Fv1,Fv2,Fv3,Fv4;
  double Gv0,Gv1,Gv2,Gv3,Gv4;
  double Hv0,Hv1,Hv2,Hv3,Hv4;
  double Ov0,Ov1,Ov2,Ov3,Ov4;
  double Pv0,Pv1,Pv2,Pv3,Pv4;
  double Ev,Bv,Dv,Fv,Gv,Hv,Ov,Pv;
  double vv,vv2,vv3,vv4;
  double rrot,rrot2,rrot3,rrot4,rrot5,rrot6,rrot7;
  
      Ev0=0.38496;
      Ev1=-0.04609;
      Ev2=0.00178;
      Ev3=-7.0e-5;
      Ev4=2.9511e-6;

      Bv0=54.438;
      Bv1=4.6063;
      Bv2=-2.0050;
      Bv3=0.19260;
      Bv4=-6.2953e-3;
      
      Dv0=-0.72593;
      Dv1=5.9990;
      Dv2=-1.5187;
      Dv3=0.12721;
      Dv4=-3.3391e-3;
      
      Fv0=-12.662;
      Fv1=20.047;
      Fv2=-4.7873;
      Fv3=0.36900;
      Fv4=-9.003e-3;

      Gv0=-24.006;
      Gv1=33.989;
      Gv2=-7.6841;
      Gv3=0.52413;
      Gv4=-1.1297e-2;

      Hv0=-22.384;
      Hv1=31.150;
      Hv2=-6.6139;
      Hv3=0.39226;
      Hv4=-9.496e-3;

      Ov0=-10.541;
      Ov1=14.746;
      Ov2=-2.9476;
      Ov3=0.16016;
      Ov4=-6.5005e-3;

      Pv0=-2.0021;
      Pv1=2.8400;
      Pv2=-0.54654;
      Pv3=0.031636;
      Pv4=-2.4398e-3;

  for(v=0;v<=2;v++) {
      vv  = (double)v + 0.5;
      vv2 = vv*vv;
      vv3 = vv*vv*vv;
      vv4 = vv*vv*vv*vv;
  
      Ev = Ev0 + Ev1*vv + Ev2*vv2 + Ev3*vv3 + Ev4*vv4;
      Bv = Bv0 + Bv1*vv + Bv2*vv2 + Bv3*vv3 + Bv4*vv4;
      Dv = Dv0 + Dv1*vv + Dv2*vv2 + Dv3*vv3 + Dv4*vv4;
      Fv = Fv0 + Fv1*vv + Fv2*vv2 + Fv3*vv3 + Fv4*vv4;
      Gv = Gv0 + Gv1*vv + Gv2*vv2 + Gv3*vv3 + Gv4*vv4;
      Hv = Hv0 + Hv1*vv + Hv2*vv2 + Hv3*vv3 + Hv4*vv4;
      Ov = Ov0 + Ov1*vv + Ov2*vv2 + Ov3*vv3 + Ov4*vv4;
      Pv = Pv0 + Pv1*vv + Pv2*vv2 + Pv3*vv3 + Pv4*vv4;

    for(J=0;J<=19;J++){
      rrot = (double)J*((double)J+1.0);
      rrot2 = rrot*rrot;
      rrot3 = rrot*rrot2;
      rrot4 = rrot*rrot3;
      rrot5 = rrot*rrot4;
      rrot6 = rrot*rrot5;
      rrot7 = rrot*rrot6;

      E_vJ_K[v][J] = 1.43879 *
        ( -       1.0e5*Ev +         rrot *Bv -  1.0e-2*rrot2*Dv +  1.0e-5*rrot3*Fv
          -1.0e-8*rrot4*Gv + 1.0e-11*rrot5*Hv - 1.0e-14*rrot6*Ov + 1.0e-17*rrot7*Pv);
      E_vJ_erg[v][J] = BOLTZMANN*E_vJ_K[v][J];
      E_vJ_eV[v][J]  = E_vJ_erg[v][J] / eV_to_erg;
      E_vJ[v][J]     = E_vJ_erg[v][J] / hc;
    }
  }

  
  //JI=0
  del_E[1][0][0][0]   = E_vJ_erg[1][0] - E_vJ_erg[0][0];
  del_E[2][0][0][0]   = E_vJ_erg[2][0] - E_vJ_erg[0][0];
  del_E[2][0][1][0]   = E_vJ_erg[2][0] - E_vJ_erg[1][0];

  //JI=1
  del_E[1][1][0][1]   = E_vJ_erg[1][1] - E_vJ_erg[0][1];
  del_E[1][1][0][3]   = E_vJ_erg[1][1] - E_vJ_erg[0][3];

  del_E[2][1][0][1]   = E_vJ_erg[2][1] - E_vJ_erg[0][1];
  del_E[2][1][0][3]   = E_vJ_erg[2][1] - E_vJ_erg[0][3];
  del_E[2][1][1][1]   = E_vJ_erg[2][1] - E_vJ_erg[1][1];
  del_E[2][1][1][3]   = E_vJ_erg[2][1] - E_vJ_erg[1][3];



  for(J=2;J<=17;J++){
    //vI=0
    del_E[0][J][0][J-2]   = E_vJ_erg[0][J] - E_vJ_erg[0][J-2];
      
    //vI=0
    del_E[1][J][0][J-2]   = E_vJ_erg[1][J] - E_vJ_erg[0][J-2];
    del_E[1][J][0][J]     = E_vJ_erg[1][J] - E_vJ_erg[0][J];
    del_E[1][J][0][J+2]   = E_vJ_erg[1][J] - E_vJ_erg[0][J+2];
    del_E[1][J][1][J-2]   = E_vJ_erg[1][J] - E_vJ_erg[1][J-2];

    //vI=2
    del_E[2][J][0][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[0][J-2];
    del_E[2][J][0][J]     = E_vJ_erg[2][J] - E_vJ_erg[0][J];
    del_E[2][J][0][J+2]   = E_vJ_erg[2][J] - E_vJ_erg[0][J+2];
    del_E[2][J][1][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[1][J-2];
    del_E[2][J][1][J]     = E_vJ_erg[2][J] - E_vJ_erg[1][J];
    del_E[2][J][1][J+2]   = E_vJ_erg[2][J] - E_vJ_erg[1][J+2];

    del_E[2][J][2][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[2][J-2];
  }

  for(J=18;J<=19;J++){
    //vI=0
    del_E[0][J][0][J-2]   = E_vJ_erg[0][J] - E_vJ_erg[0][J-2];
      
    //vI=0
    del_E[1][J][0][J-2]   = E_vJ_erg[1][J] - E_vJ_erg[0][J-2];
    del_E[1][J][0][J]     = E_vJ_erg[1][J] - E_vJ_erg[0][J];
  
    del_E[1][J][1][J-2]   = E_vJ_erg[1][J] - E_vJ_erg[1][J-2];

    //vI=2
    del_E[2][J][0][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[0][J-2];
    del_E[2][J][0][J]     = E_vJ_erg[2][J] - E_vJ_erg[0][J];
  
    del_E[2][J][1][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[1][J-2];
    del_E[2][J][1][J]     = E_vJ_erg[2][J] - E_vJ_erg[1][J];
  

    del_E[2][J][2][J-2]   = E_vJ_erg[2][J] - E_vJ_erg[2][J-2];
  }




  //B coefficients
  // Indices are ordered in the reverse way such that B(*,6,*,4) =  B_4,6(in standard notation)

  //JI=0
  JI=0; JF=0;
  B[1][0][0][0] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[1][0][0][0] / (2.0 * pow(del_E[1][0][0][0], 3));
    
  B[2][0][0][0] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][0][0][0] / (2.0 * pow(del_E[2][0][0][0], 3));

  B[2][0][1][0] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][0][1][0] / (2.0 * pow(del_E[2][0][1][0], 3));

  //JI=1
  JI=1; JF=1;
  B[1][1][0][1] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[1][1][0][1] / (2.0 * pow(del_E[1][1][0][1], 3));

  JI=1; JF=3;
  B[1][1][0][3] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[1][1][0][3] / (2.0 * pow(del_E[1][1][0][3], 3));

  JI=1; JF=1;
  B[2][1][0][1] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][1][0][1] / (2.0 * pow(del_E[2][1][0][1], 3));

  JI=1; JF=3;
  B[2][1][0][3] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][1][0][3] / (2.0 * pow(del_E[2][1][0][3], 3));

  JI=1; JF=1;
  B[2][1][1][1] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][1][1][1] / (2.0 * pow(del_E[2][1][1][1], 3));

  JI=1; JF=3;
  B[2][1][1][3] = (2.0*JI+1.0)/(2.0*JF + 1.0) * hc2 * A[2][1][1][3] / (2.0 * pow(del_E[2][1][1][3], 3));



  for(J=2;J<=17;J++){
    //vI=0
    B[0][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[0][J][0][J-2]  / (2.0 * pow(del_E[0][J][0][J-2], 3));
    

    //vI=1
    B[1][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[1][J][0][J-2]  / (2.0 * pow(del_E[1][J][0][J-2], 3));
    B[1][J][0][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[1][J][0][J]    / (2.0 * pow(del_E[1][J][0][J],   3));
    B[1][J][0][J+2] = (2.0*J+1.0)/(2.0*(J+2) + 1.0) * hc2 * A[1][J][0][J+2]  / (2.0 * pow(del_E[1][J][0][J+2], 3));
    B[1][J][1][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[1][J][1][J-2]  / (2.0 * pow(del_E[1][J][1][J-2], 3));


    //vI=2
    B[2][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][0][J-2]  / (2.0 * pow(del_E[2][J][0][J-2], 3));
    B[2][J][0][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[2][J][0][J]    / (2.0 * pow(del_E[2][J][0][J],   3));
    B[2][J][0][J+2] = (2.0*J+1.0)/(2.0*(J+2) + 1.0) * hc2 * A[2][J][0][J+2]  / (2.0 * pow(del_E[2][J][0][J+2], 3));   
    B[2][J][1][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][1][J-2]  / (2.0 * pow(del_E[2][J][1][J-2], 3));
    B[2][J][1][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[2][J][1][J]    / (2.0 * pow(del_E[2][J][1][J],   3));
    B[2][J][1][J+2] = (2.0*J+1.0)/(2.0*(J+2) + 1.0) * hc2 * A[2][J][1][J+2]  / (2.0 * pow(del_E[2][J][1][J+2], 3));
    B[2][J][2][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][2][J-2]  / (2.0 * pow(del_E[2][J][2][J-2], 3));
  }


  for(J=18;J<=19;J++){
    //vI=0
    B[0][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[0][J][0][J-2]  / (2.0 * pow(del_E[0][J][0][J-2], 3));
    

    //vI=1
    B[1][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[1][J][0][J-2]  / (2.0 * pow(del_E[1][J][0][J-2], 3));
    B[1][J][0][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[1][J][0][J]    / (2.0 * pow(del_E[1][J][0][J],   3));
    B[1][J][1][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[1][J][1][J-2]  / (2.0 * pow(del_E[1][J][1][J-2], 3));

    
    //vI=2
    B[2][J][0][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][0][J-2]  / (2.0 * pow(del_E[2][J][0][J-2], 3));
    B[2][J][0][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[2][J][0][J]    / (2.0 * pow(del_E[2][J][0][J],   3));
    B[2][J][1][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][1][J-2]  / (2.0 * pow(del_E[2][J][1][J-2], 3));
    B[2][J][1][J]   = (2.0*J+1.0)/(2.0*J + 1.0)     * hc2 * A[2][J][1][J]    / (2.0 * pow(del_E[2][J][1][J],   3));
    B[2][J][2][J-2] = (2.0*J+1.0)/(2.0*(J-2) + 1.0) * hc2 * A[2][J][2][J-2]  / (2.0 * pow(del_E[2][J][2][J-2], 3));
  }
 
/*
      for(J=0;J<20;J++) {
        printf("%3d",J);
        for(v=0;v<3;v++) {
          printf("%10.3e ",(E_vJ_erg[v][J])/BOLTZMANN);
        }
        printf("\n");
      }
        printf("\n");
	*/
	return 1; 
}


double H2line_reduction(int ithis, double T, double NH2)
{
  
  int v,J,vI,JI,vF,JF;
  double Z, hc;
  double f_reduce, tau, H2column;
  double cool_rate_thin, cool_rate_thick;
  double v_thermal, t_freefall, Sobolev_Length;
  
  hc = PLANCK * CLIGHT;

  //Partition function
  Z=0;
  for(v=0;v<=2;v++)
    for(J=0;J<=19;J++)
      Z += g_weight[J] * (2*J + 1) * exp(-E_vJ_K[v][J]/T);        
  
  //Boltzmann distribution
  for(v=0;v<=2;v++)
    for(J=0;J<=19;J++)
      n_pop[v][J] = g_weight[J] * (2*J + 1) * exp(-E_vJ_K[v][J]/T)/Z;
  
  if(T < 70.0){
    for(v=0;v<=2;v++)
      for(J=0;J<=19;J++)
	if(v==0 && J==0)
	  n_pop[v][J] = 1.;
	else
	  n_pop[v][J] = 0.;      

    return 1.;
  }
  
  

  // line emissivity
  //computing the fractional contribution
  
  //J=0
  Lambda_frac[1][0][0][0] = A[1][0][0][0] * n_pop[1][0] * (E_vJ_erg[1][0] - E_vJ_erg[0][0]);
  Lambda_frac[2][0][0][0] = A[2][0][0][0] * n_pop[2][0] * (E_vJ_erg[2][0] - E_vJ_erg[0][0]);
  Lambda_frac[2][0][1][0] = A[2][0][1][0] * n_pop[2][0] * (E_vJ_erg[2][0] - E_vJ_erg[1][0]);

  //J=1
  Lambda_frac[1][1][0][1] = A[1][1][0][1] * n_pop[1][1] * (E_vJ_erg[1][1] - E_vJ_erg[0][1]);
  Lambda_frac[1][1][0][3] = A[1][1][0][3] * n_pop[1][1] * (E_vJ_erg[1][1] - E_vJ_erg[0][3]);
  
  Lambda_frac[2][1][0][1] = A[2][1][0][1] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[0][1]);
  Lambda_frac[2][1][0][3] = A[2][1][0][3] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[0][3]);
  Lambda_frac[2][1][1][1] = A[2][1][1][1] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[1][1]);
  Lambda_frac[2][1][1][3] = A[2][1][1][3] * n_pop[2][1] * (E_vJ_erg[2][1] - E_vJ_erg[1][3]);
  
  
  for(J=2;J<=17;J++){
    
    //v=0
    Lambda_frac[0][J][0][J-2] = A[0][J][0][J-2] * n_pop[0][J] * (E_vJ_erg[0][J] - E_vJ_erg[0][J-2]);
    
    //v=1
    Lambda_frac[1][J][0][J-2] = A[1][J][0][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[1][J][0][J]   = A[1][J][0][J]   * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J]);
    Lambda_frac[1][J][0][J+2] = A[1][J][0][J+2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J+2]);
    Lambda_frac[1][J][1][J-2] = A[1][J][1][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[1][J-2]);

    //v=2
    Lambda_frac[2][J][0][J-2] = A[2][J][0][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[2][J][0][J]   = A[2][J][0][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J]);
    Lambda_frac[2][J][0][J+2] = A[2][J][0][J+2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J+2]);
    Lambda_frac[2][J][1][J-2] = A[2][J][1][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J-2]);
    Lambda_frac[2][J][1][J]   = A[2][J][1][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J]);
    Lambda_frac[2][J][1][J+2] = A[2][J][1][J+2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J+2]);

    Lambda_frac[2][J][2][J-2] = A[2][J][2][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[2][J-2]);
  }

  
  for(J=18;J<=19;J++){
    
    //v=0
    Lambda_frac[0][J][0][J-2] = A[0][J][0][J-2] * n_pop[0][J] * (E_vJ_erg[0][J] - E_vJ_erg[0][J-2]);
    
    //v=1
    Lambda_frac[1][J][0][J-2] = A[1][J][0][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[1][J][0][J]   = A[1][J][0][J]   * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[0][J]);
    
    Lambda_frac[1][J][1][J-2] = A[1][J][1][J-2] * n_pop[1][J] * (E_vJ_erg[1][J] - E_vJ_erg[1][J-2]);

    //v=2
    Lambda_frac[2][J][0][J-2] = A[2][J][0][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J-2]);
    Lambda_frac[2][J][0][J]   = A[2][J][0][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[0][J]);
    
    Lambda_frac[2][J][1][J-2] = A[2][J][1][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J-2]);
    Lambda_frac[2][J][1][J]   = A[2][J][1][J]   * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[1][J]);
  
    Lambda_frac[2][J][2][J-2] = A[2][J][2][J-2] * n_pop[2][J] * (E_vJ_erg[2][J] - E_vJ_erg[2][J-2]);
  }
  

  
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++){
	  
	  del_nu[vI][JI][vF][JF]  = del_E[vI][JI][vF][JF] / hc * sqrt(BOLTZMANN * T/PROTONMASS); //factor 2 canceled
	  del_E_K[vI][JI][vF][JF] = del_E[vI][JI][vF][JF] / BOLTZMANN;
	  
	}
  
  
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++){
	  alpha[vI][JI][vF][JF] = del_E[vI][JI][vF][JF] * n_pop[vF][JF]
	    * B[vI][JI][vF][JF] * (1.0 - exp(-del_E_K[vI][JI][vF][JF]/T)); 
	  
	  if(del_nu[vI][JI][vF][JF] > 0)
	    alpha[vI][JI][vF][JF] = alpha[vI][JI][vF][JF]/(sqrt(M_PI) * del_nu[vI][JI][vF][JF] * 4.0 * M_PI);
	  else
	    alpha[vI][JI][vF][JF] = 0;
	   
	}

         
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++){
	  tau = alpha[vI][JI][vF][JF] * NH2;

	  if(tau > 0 && tau < 30.0)
	    //Lambda_frac_thick[vI][JI][vF][JF] = exp(-tau) * Lambda_frac[vI][JI][vF][JF];
	    Lambda_frac_thick[vI][JI][vF][JF] = ((1.0-exp(-tau))/tau) * Lambda_frac[vI][JI][vF][JF];
	  //if velocity is proportional to the radius
	  else
	    Lambda_frac_thick[vI][JI][vF][JF] = 0;	  	 
	}
  

  
  //optically-thin/thick cooling rate
  cool_rate_thin =0; 
  cool_rate_thick=0; 
  
  for(vI=0;vI<=2;vI++)
    for(JI=0;JI<=19;JI++)
      for(vF=0;vF<=2;vF++)
	for(JF=0;JF<=19;JF++){
	  cool_rate_thin  += Lambda_frac[vI][JI][vF][JF];
	  cool_rate_thick += Lambda_frac_thick[vI][JI][vF][JF];
	}
  
  f_reduce = cool_rate_thick/cool_rate_thin;
  
  if(isnan(f_reduce)){
    f_reduce = 1.;
    //raise(11);
  }
  		
  return f_reduce;
}



#ifdef METAL_COOLING
double OxygenI_cooling(double nOI, double ne, double nH, double T_gas)
{

  double g1, g2, g3, g4, g5;

  double gamma21H, gamma31H, gamma41H, gamma51H;
  double gamma21e, gamma31e, gamma41e, gamma51e;

  double gamma32H, gamma42H, gamma52H;
  double gamma32e, gamma42e, gamma52e;

  double gamma43H, gamma53H;
  double gamma43e, gamma53e;

  double gamma54H;
  double gamma54e;


  double gamma12H, gamma13H, gamma14H, gamma15H;
  double gamma12e, gamma13e, gamma14e, gamma15e;

  double gamma23H, gamma24H, gamma25H;
  double gamma23e, gamma24e, gamma25e;

  double gamma34H, gamma35H;
  double gamma34e, gamma35e;

  double gamma45H;
  double gamma45e;

  double R12, R13, R14, R15, R21, R31, R41, R51;
  double R23, R24, R25, R32, R42, R52;
  double R34, R35, R43, R53;
  double R45, R54;

  double A21, A31, A41, A51, A32, A42, A52, A43, A53, A54;
  double delE_21, delE_31, delE_41, delE_51, delE_32, delE_42, delE_52, delE_43, delE_53, delE_54;


  double T100, Lambda;
  
  double Coeff[5][5], b[5][2];

  int dim, column;



  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 5;  // 3P2
  g2 = 3;  // 3P1
  g3 = 1;  // 3P0
  g4 = 5;  // 1D2
  g5 = 1;  // 1S0

  gamma21H = 9.2e-11 * pow(T100, 0.67);

  gamma31H = 4.3e-11 * pow(T100, 0.8);
  gamma32H = 1.1e-10 * pow(T100, 0.44);

  gamma41H = 1.0e-12;
  gamma42H = 1.0e-12;
  gamma43H = 1.0e-12;

  gamma51H = 1.0e-12;
  gamma52H = 1.0e-12;
  gamma53H = 1.0e-12;
  gamma54H = 1.0e-12;


  gamma21e = 1.4e-8;

  gamma31e = 1.4e-8;
  gamma32e = 5.0e-9;

  gamma41e = 1.0e-10;
  gamma42e = 1.0e-10;
  gamma43e = 1.0e-10;

  gamma51e = 1.0e-10;
  gamma52e = 1.0e-10;
  gamma53e = 1.0e-10;
  gamma54e = 1.0e-10;


  A21 = 8.9e-5;
  A31 = 1.0e-10;
  A32 = 1.7e-5;
  A41 = 6.3e-3;
  A42 = 2.1e-3;
  A43 = 7.3e-7;
  A51 = 2.9e-4;
  A52 = 7.3e-2;
  A53 = 0; // no J=0 -> 0 radiative transition
  A54 = 1.2;

  delE_21 = 3.144e-14; // erg Texc = 227.7 K
  delE_32 = 1.365e-14; // erg Texc = 98.8 K
  delE_43 = 3.14e-12; // erg Texc = 2.283e4 K
  delE_54 = 3.56e-12; // erg Texc = 2.578e4 K

  delE_31 = delE_32 + delE_21;
  delE_41 = delE_43 + delE_32 + delE_21;
  delE_51 = delE_54 + delE_43 + delE_32 + delE_21;

  delE_42 = delE_43 + delE_32;
  delE_52 = delE_54 + delE_43 + delE_32;

  delE_53 = delE_54 + delE_43;


  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  gamma13H = (g3/g1) * gamma31H * exp(-delE_31/(BOLTZMANN * T_gas));
  gamma13e = (g3/g1) * gamma31e * exp(-delE_31/(BOLTZMANN * T_gas));
  
  gamma23H = (g3/g2) * gamma32H * exp(-delE_32/(BOLTZMANN * T_gas));
  gamma23e = (g3/g2) * gamma32e * exp(-delE_32/(BOLTZMANN * T_gas));

  gamma14H = (g4/g1) * gamma41H * exp(-delE_41/(BOLTZMANN * T_gas));
  gamma14e = (g4/g1) * gamma41e * exp(-delE_41/(BOLTZMANN * T_gas));

  gamma24H = (g4/g2) * gamma42H * exp(-delE_42/(BOLTZMANN * T_gas));
  gamma24e = (g4/g2) * gamma42e * exp(-delE_42/(BOLTZMANN * T_gas));

  gamma34H = (g4/g3) * gamma43H * exp(-delE_43/(BOLTZMANN * T_gas));
  gamma34e = (g4/g3) * gamma43e * exp(-delE_43/(BOLTZMANN * T_gas));

  gamma15H = (g5/g1) * gamma51H * exp(-delE_51/(BOLTZMANN * T_gas));
  gamma15e = (g5/g1) * gamma51e * exp(-delE_51/(BOLTZMANN * T_gas));

  gamma25H = (g5/g2) * gamma52H * exp(-delE_52/(BOLTZMANN * T_gas));
  gamma25e = (g5/g2) * gamma52e * exp(-delE_52/(BOLTZMANN * T_gas));

  gamma35H = (g5/g3) * gamma53H * exp(-delE_53/(BOLTZMANN * T_gas));
  gamma35e = (g5/g3) * gamma53e * exp(-delE_53/(BOLTZMANN * T_gas));

  gamma45H = (g5/g4) * gamma54H * exp(-delE_54/(BOLTZMANN * T_gas));
  gamma45e = (g5/g4) * gamma54e * exp(-delE_54/(BOLTZMANN * T_gas));
 



  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R31 = A21 + nH * gamma31H  + gamma31e * ne;
  R41 = A41 + nH * gamma41H  + gamma41e * ne;
  R51 = A51 + nH * gamma51H  + gamma51e * ne;

  R32 = A32 + nH * gamma32H  + gamma32e * ne;
  R42 = A42 + nH * gamma42H  + gamma42e * ne;
  R52 = A52 + nH * gamma52H  + gamma52e * ne;

  R43 = A43 + nH * gamma43H  + gamma43e * ne;
  R53 = A53 + nH * gamma53H  + gamma53e * ne;

  R54 = A54 + nH * gamma54H * nH + gamma54e * ne;


  R12 = nH * gamma12H  + gamma12e * ne;
  R13 = nH * gamma13H  + gamma13e * ne;
  R14 = nH * gamma14H  + gamma14e * ne;
  R15 = nH * gamma15H  + gamma15e * ne;

  R23 = nH * gamma23H  + gamma23e * ne;
  R24 = nH * gamma24H  + gamma24e * ne;
  R25 = nH * gamma25H  + gamma25e * ne;

  R34 = nH * gamma34H  + gamma34e * ne;
  R35 = nH * gamma35H  + gamma35e * ne;
  
  R45 = nH * gamma45H  + gamma45e * ne;



  Coeff[0][0] = 1.0;
  Coeff[0][1] = 1.0;
  Coeff[0][2] = 1.0;
  Coeff[0][3] = 1.0;
  Coeff[0][4] = 1.0;

  Coeff[1][0] = R12;
  Coeff[1][1] = -(R21 + R23 + R24 + R25);
  Coeff[1][2] = R32;
  Coeff[1][3] = R42;
  Coeff[1][4] = R52;

  Coeff[2][0] = R13;
  Coeff[2][1] = R23;
  Coeff[2][2] = -(R31 + R32 + R34 + R35);
  Coeff[2][3] = R43;
  Coeff[2][4] = R53;

  Coeff[3][0] = R14;
  Coeff[3][1] = R24;
  Coeff[3][2] = R34;
  Coeff[3][3] = -(R41 + R42 + R43 + R45);
  Coeff[3][4] = R54;

  Coeff[4][0] = R15;
  Coeff[4][1] = R25;
  Coeff[4][2] = R35;
  Coeff[4][3] = R45;
  Coeff[4][4] = -(R51 + R52 + R53 + R54);


  b[0][0] = nOI;
  b[1][0] = 0;
  b[2][0] = 0;
  b[3][0] = 0;
  b[4][0] = 0;

  b[0][1] = nOI;
  b[1][1] = 0;
  b[2][1] = 0;
  b[3][1] = 0;
  b[4][1] = 0;

  dim=5; column=2;

  gaussj5(Coeff, dim, b, column);

  
  double n1, n2, n3, n4, n5;
  
  n1 = b[0][0];
  n2 = b[1][0];
  n3 = b[2][0];
  n4 = b[3][0];
  n5 = b[4][0];

  Lambda = n5 * (A51 * delE_51 + A52 * delE_52 + A53 * delE_53 + A54 * delE_54)
    + n4 * (A41 * delE_41 + A42 * delE_42 + A43 * delE_43)
    + n3 * (A31 * delE_31 + A32 * delE_32)
    + n2 * A21 * delE_21;



  return Lambda; // erg/sec per OI atom

}


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
/* modified from original gaussj so that array[0..n] can be used */
void gaussj5(double a[5][5], int n, double b[5][2], int m)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	int indxc[5],indxr[5],ipiv[5];

	if(n != 5 || m!= 2)
	  exit(3348);
	
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0){
		  printf("gaussj: Singular Matrix");
		  exit(38276);
		}
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}
#undef SWAP

#endif /* METAL_LINE_COOLING */




