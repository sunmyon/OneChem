#include "chemistry_model.hpp"
#include <new>
#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include "chemistry_constants.h"
#include "chemistry_types.h"
#include "chemistry_tables.h"
#include "cooling_func.h"

#ifdef CHEMISTRY_SC

double C0[6]={152.519  ,     49.534,   -118.858 ,     92.536 ,    -34.194 ,    4.982  };

/*********************************************************************/
/****************** Cooling by OxygenI *******************************/
/*********************************************************************/

static int gaussj5(double a[5][5], int n, double b[5][2], int m);
static int gaussj3(double a[3][3], int n, double b[3][2], int m);

double OxygenI_fivelevel_cooling(double nOI, double ne, double nH, double T_gas)
{
  double g_OI[6], n_OI[6];
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

  double A_OI[6][6], delE_OI[6][6], tau_OI[6][6], tau;
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

  g_OI[1] = g1;
  g_OI[2] = g2;
  g_OI[3] = g3;
  g_OI[4] = g4;
  g_OI[5] = g5;

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

  A_OI[2][1] = A21;
  A_OI[3][1] = A31;
  A_OI[4][1] = A41;
  A_OI[5][1] = A51;
  A_OI[3][2] = A32;
  A_OI[4][2] = A42;
  A_OI[5][2] = A52;
  A_OI[4][3] = A43;
  A_OI[5][3] = A53;
  A_OI[5][4] = A54;

  delE_OI[2][1] = delE_21;
  delE_OI[3][1] = delE_31;
  delE_OI[4][1] = delE_41;
  delE_OI[5][1] = delE_51;
  delE_OI[3][2] = delE_32;
  delE_OI[4][2] = delE_42;
  delE_OI[5][2] = delE_52;
  delE_OI[4][3] = delE_43;
  delE_OI[5][3] = delE_53;
  delE_OI[5][4] = delE_54;


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



  // cout<<"Oxygen coefficients "<<endl;
  //cout<<Coeff[1][0]<<" "<<Coeff[1][1]<<" "<<Coeff[1][2]<<" "<<Coeff[1][3]<<" "<<Coeff[1][4]<<endl;
  //cout<<Coeff[2][0]<<" "<<Coeff[2][1]<<" "<<Coeff[2][2]<<" "<<Coeff[2][3]<<" "<<Coeff[2][4]<<endl;
  //cout<<Coeff[3][0]<<" "<<Coeff[3][1]<<" "<<Coeff[3][2]<<" "<<Coeff[3][3]<<" "<<Coeff[3][4]<<endl;
  //cout<<Coeff[4][0]<<" "<<Coeff[4][1]<<" "<<Coeff[4][2]<<" "<<Coeff[4][3]<<" "<<Coeff[4][4]<<endl;


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

  if(!gaussj5(Coeff, dim, b, column))
    return NAN;


  //cout<<"Oxygen population "<<endl;
  //cout<<b[0][0]<<" "<<b[1][0]<<" "<<b[2][0]<<" "<<b[3][0]<<" "<<b[4][0]<<" and "
  //    <<b[0][0]+b[1][0]+b[2][0]+b[3][0]+b[4][0]<<endl;


  double n1, n2, n3, n4, n5;

  n1 = b[0][0];
  n2 = b[1][0];
  n3 = b[2][0];
  n4 = b[3][0];
  n5 = b[4][0];

  double Lambda_thin;
  Lambda_thin = n5 * (A51 * delE_51 + A52 * delE_52 + A53 * delE_53 + A54 * delE_54)
    + n4 * (A41 * delE_41 + A42 * delE_42 + A43 * delE_43)
    + n3 * (A31 * delE_31 + A32 * delE_32)
    + n2 * A21 * delE_21;
  //Lambda_thin = n2 * A21 * delE_21;

  n_OI[1] = n1;
  n_OI[2] = n2;
  n_OI[3] = n3;
  n_OI[4] = n4;
  n_OI[5] = n5;

  double t_freefall;
  t_freefall = sqrt(3.0*M_PI/32.0/GRAVITY/(1.2 * PROTONMASS * nH));  /*  in second */

  double hc;
  hc  = PLANCK*CLIGHT;

  int iu, il;
  for(iu=5, Lambda=0.;iu>1;iu--){
    for(il=iu-1;il>=1;il--){
      tau =
	g_OI[iu]/g_OI[il] * A_OI[iu][il] * hc * hc * hc / (8 * M_PI * pow(delE_OI[iu][il],3))
	* (n_OI[il] - n_OI[iu] * g_OI[il]/g_OI[iu]) * t_freefall;

      if(tau > 1.e-5)
	Lambda += n_OI[iu] * A_OI[iu][il] * delE_OI[iu][il] * (1. - exp(-tau))/tau;
      else
	Lambda += n_OI[iu] * A_OI[iu][il] * delE_OI[iu][il];

      tau_OI[iu][il] = tau;
    }
  }

  double fac;
  if(tau_OI[2][1] > 1.e-5)
    fac = (1. - exp(-tau_OI[2][1]))/tau_OI[2][1];
  else
    fac = 1.;

  double Lambda_two_level;
  Lambda_two_level = n_OI[2] * A_OI[2][1] * delE_OI[2][1] * fac;

  printf("[OI] tau[5]=%g %g %g %g\n", tau_OI[5][4], tau_OI[5][3], tau_OI[5][2], tau_OI[5][1]);
  printf("[OI] tau[4]=%g %g %g\n", tau_OI[4][3], tau_OI[4][2], tau_OI[4][1]);
  printf("[OI] tau[3]=%g %g\n", tau_OI[3][2], tau_OI[3][1]);
  printf("[OI] tau[2]=%g\n", tau_OI[2][1]);

  printf("[OI] y[]=%g %g %g %g %g\n", n1/nOI, n2/nOI, n3/nOI, n4/nOI, n5/nOI);

  printf("[OI] Lambda=%g thin=%g twolevel=%g %g %g\n"
	 , Lambda/nOI, Lambda_thin/nOI
	 , n2/nOI * A21 * delE_21 * fac, Lambda/nOI, Lambda/nH);

  return Lambda; // erg/sec per OI atom
}


double OxygenI_twolevel_cooling(double nOI, double ne, double nH,
                                double T_gas, double radiation_temperature,
                                int diagnostics_enabled)
{
  double g1, g2;

  double gamma21H;
  double gamma21e;
  double gamma12H;
  double gamma12e;

  double R12, R21;

  double A21;
  double delE_21;

  double T100, Lambda;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 5;  // 3P2
  g2 = 3;  // 3P1

  gamma21H = 9.2e-11 * pow(T100, 0.67);
  gamma21e = 1.4e-8;
  A21 = 8.9e-5;

  delE_21 = 3.144e-14; // erg Texc = 227.7 K

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R12 = nH * gamma12H  + gamma12e * ne;

#ifdef CMB_HEATING
  double Texc, Q21;
  Texc = 227.7;
  Q21  = 1. / (exp(Texc/radiation_temperature) - 1.);

  R21 = A21 * (1. + Q21) + nH * gamma21H  + gamma21e * ne;
  R12 = g2/g1 * A21 *Q21 + nH * gamma12H  + gamma12e * ne;
#endif


  double n1, n2;
  n1 = nOI * R21 / (R12 + R21);
  n2 = nOI * R12 / (R12 + R21);

  Lambda = n2 * A21 * delE_21;

  double t_freefall;
  t_freefall = sqrt(3.0*M_PI/32.0/GRAVITY/(1.2 * PROTONMASS * nH));  /*  in second */

  double hc;
  hc  = PLANCK*CLIGHT;

  double tau21;
  tau21 =
    g2/g1 * A21 * hc * hc * hc / (8 * M_PI * pow(delE_21,3))
    * (n1 - n2 * g1/g2) * t_freefall;

  if(tau21 > 1.e-5)
    Lambda = n2 * A21 * delE_21 * (1. - exp(-tau21))/tau21;
  else
    Lambda = n2 * A21 * delE_21;

#ifdef CMB_HEATING
  double fac;
  //fac = 1. - Q21 * (g2/g1 * n1/n2 - 1);
  fac = 1. - Q21 * (g2/g1 * R21/R12 - 1);
  Lambda *= fac;

  //printf("[OI]nH=%g Lambda=%g fac=%g f_red=%g\n", nH, Lambda, fac, (1. - exp(-tau21))/tau21);
#endif

  if(diagnostics_enabled)
    printf("[OI]nH=%g nOI=%g Lambda=%g tau21=%g n2=%g R12=%g R21=%g\n", nH, nOI, Lambda, tau21, n2, R12, R21);

  return Lambda; // erg/sec per OI atom
}

double OxygenI_twolevel_cooling_thin(double nOI, double ne, double nH, double T_gas)
{
  double g1, g2;

  double gamma21H, gamma12H;
  double gamma21e, gamma12e;

  double R12, R21;

  double A21;
  double delE_21;

  double T100, Lambda;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 5;  // 3P2
  g2 = 3;  // 3P1

  gamma21H = 9.2e-11 * pow(T100, 0.67);
  gamma21e = 1.4e-8;
  A21 = 8.9e-5;

  delE_21 = 3.144e-14; // erg Texc = 227.7 K

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R12 = nH * gamma12H  + gamma12e * ne;


  double n2;
  n2 = nOI * R12 / (R12 + R21);

  Lambda = n2 * A21 * delE_21;

  return Lambda; // erg/sec per OI atom
}


/*********************************************************************/
/****************** Cooling by CarbonII *******************************/
/*********************************************************************/


double CarbonII_twolevel_cooling(double nCII, double ne, double nH,
                                 double T_gas, double radiation_temperature,
                                 int diagnostics_enabled)
{
  double g1, g2;

  double gamma21H, gamma12H;
  double gamma21e, gamma12e;

  double R12, R21;

  double A21;
  double delE_21;

  double T100, Lambda;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 2.0;
  g2 = 4.0;

  gamma21H = 8.0e-10 * pow(T100, 0.07);
  gamma21e = 2.8e-7  * pow(T100, -0.5);
  A21 = 2.4e-6;
  delE_21 = 1.259e-14;

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R12 = nH * gamma12H  + gamma12e * ne;

#ifdef CMB_HEATING
  double Texc, Q21;
  Texc = 92.;
  Q21  = 1. / (exp(Texc/radiation_temperature) - 1.);

  R21 = A21 * (1. + Q21) + nH * gamma21H  + gamma21e * ne;
  R12 = g2/g1 * A21 *Q21 + nH * gamma12H  + gamma12e * ne;
#endif


  double n1, n2;
  n1 = nCII * R21 / (R12 + R21);
  n2 = nCII * R12 / (R12 + R21);

  Lambda = n2 * A21 * delE_21;

  double t_freefall;
  t_freefall = sqrt(3.0*M_PI/32.0/GRAVITY/(1.2 * PROTONMASS * nH));  /*  in second */

  double hc;
  hc  = PLANCK*CLIGHT;

  double tau21;
  tau21 =
    g2/g1 * A21 * hc * hc * hc / (8 * M_PI * pow(delE_21,3))
    * (n1 - n2 * g1/g2) * t_freefall;

  if(tau21 > 1.e-5)
    Lambda = n2 * A21 * delE_21 * (1. - exp(-tau21))/tau21;
  else
    Lambda = n2 * A21 * delE_21;


#ifdef CMB_HEATING
  double fac;
  //  fac = 1. - Q21 * (g2/g1 * n1/n2 - 1);
  fac = 1. - Q21 * (g2/g1 * R21/R12 - 1);
  Lambda *= fac;

  //printf("[CII]nH=%g Lambda=%g fac=%g f_red=%g\n", nH, Lambda, fac, (1. - exp(-tau21))/tau21);
#endif

  if(diagnostics_enabled)
    printf("[CII]nH=%g nCII=%g Lambda=%g tau21=%g n2=%g R12=%g R21=%g\n", nH, nCII, Lambda,  tau21, n2, R12, R21);

  return Lambda; // erg/sec per CII atom
}


double CarbonII_twolevel_cooling_thin(double nCII, double ne, double nH, double T_gas)
{
  double g1, g2;

  double gamma21H, gamma12H;
  double gamma21e, gamma12e;

  double R12, R21;

  double A21;
  double delE_21;

  double T100, Lambda;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 2.0;
  g2 = 4.0;

  gamma21H = 8.0e-10 * pow(T100, 0.07);
  gamma21e = 2.8e-7  * pow(T100, -0.5);
  A21 = 2.4e-6;
  delE_21 = 1.259e-14;

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R12 = nH * gamma12H  + gamma12e * ne;

  double n2;
  n2 = nCII * R12 / (R12 + R21);

  Lambda = n2 * A21 * delE_21;

  return Lambda; // erg/sec per OI atom
}

/*********************************************************************/
/****************** Cooling by CarbonI *******************************/
/*********************************************************************/

double CarbonI_threelevel_cooling(double nCI, double ne, double nH, double T_gas)
{
  double g_CI[4], n_CI[4];
  double g1, g2, g3;

  double gamma21H, gamma31H;
  double gamma21e, gamma31e;

  double gamma32H;
  double gamma32e;

  double gamma12H, gamma13H;
  double gamma12e, gamma13e;

  double gamma23H;
  double gamma23e;

  double R12, R13, R21, R31;
  double R23, R32;

  double A_CI[4][4], delE_CI[4][4], tau_CI[4][4], tau;
  double A21, A31, A32;
  double delE_21, delE_31, delE_32;

  double T100, Lambda;

  double Coeff[3][3], b[3][2];

  int dim, column;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 1;  // 3P0
  g2 = 3;  // 3P1
  g3 = 5;  // 3P2

  g_CI[1] = g1;
  g_CI[2] = g2;
  g_CI[3] = g3;

  gamma21H = 1.6e-10 * pow(T100, 0.14);
  gamma31H = 9.2e-11 * pow(T100, 0.26);
  gamma32H = 2.9e-10 * pow(T100, 0.26);

  gamma21e = 3.0e-9;
  gamma31e = 5.0e-9;
  gamma32e = 1.5e-8;

  A21 = 7.9e-8;
  A31 = 2.0e-14;
  A32 = 2.7e-7;

  delE_21 = 3.312e-15; // erg Texc = 24 K
  delE_32 = 5.382e-15; // erg Texc = 39 K
  delE_31 = delE_32 + delE_21;

  A_CI[2][1] = A21;
  A_CI[3][1] = A31;
  A_CI[3][2] = A32;

  delE_CI[2][1] = delE_21;
  delE_CI[3][1] = delE_31;
  delE_CI[3][2] = delE_32;

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  gamma13H = (g3/g1) * gamma31H * exp(-delE_31/(BOLTZMANN * T_gas));
  gamma13e = (g3/g1) * gamma31e * exp(-delE_31/(BOLTZMANN * T_gas));

  gamma23H = (g3/g2) * gamma32H * exp(-delE_32/(BOLTZMANN * T_gas));
  gamma23e = (g3/g2) * gamma32e * exp(-delE_32/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R31 = A21 + nH * gamma31H  + gamma31e * ne;
  R32 = A32 + nH * gamma32H  + gamma32e * ne;

  R12 = nH * gamma12H  + gamma12e * ne;
  R13 = nH * gamma13H  + gamma13e * ne;
  R23 = nH * gamma23H  + gamma23e * ne;

  Coeff[0][0] = 1.0;
  Coeff[0][1] = 1.0;
  Coeff[0][2] = 1.0;

  Coeff[1][0] = R12;
  Coeff[1][1] = -(R21 + R23);
  Coeff[1][2] = R32;

  Coeff[2][0] = R13;
  Coeff[2][1] = R23;
  Coeff[2][2] = -(R31 + R32);

  b[0][0] = nCI;
  b[1][0] = 0;
  b[2][0] = 0;

  b[0][1] = nCI;
  b[1][1] = 0;
  b[2][1] = 0;

  dim=3; column=2;

  if(!gaussj3(Coeff, dim, b, column))
    return NAN;

  double n1, n2, n3;
  n1 = b[0][0];
  n2 = b[1][0];
  n3 = b[2][0];

  /*Lambda = n5 * (A51 * delE_51 + A52 * delE_52 + A53 * delE_53 + A54 * delE_54)
    + n4 * (A41 * delE_41 + A42 * delE_42 + A43 * delE_43)
    + n3 * (A31 * delE_31 + A32 * delE_32)
    + n2 * A21 * delE_21;*/

  Lambda = n2 * A21 * delE_21;

  n_CI[1] = n1;
  n_CI[2] = n2;
  n_CI[3] = n3;

  double t_freefall;
  t_freefall = sqrt(3.0*M_PI/32.0/GRAVITY/(1.2 * PROTONMASS * nH));  /*  in second */

  double hc;
  hc  = PLANCK*CLIGHT;

  int iu, il;
  for(iu=3, Lambda=0.;iu>1;iu--){
    for(il=iu-1;il>=1;il--){
      tau =
	g_CI[iu]/g_CI[il] * A_CI[iu][il] * hc * hc * hc / (8 * M_PI * pow(delE_CI[iu][il],3))
	* (n_CI[il] - n_CI[iu] * g_CI[il]/g_CI[iu]) * t_freefall;


      if(tau > 1.e-5)
	Lambda += n_CI[iu] * A_CI[iu][il] * delE_CI[iu][il] * (1. - exp(-tau))/tau;
      else
	Lambda += n_CI[iu] * A_CI[iu][il] * delE_CI[iu][il];

      tau_CI[iu][il] = tau;
    }
  }

  double fac;
  if(tau_CI[2][1] > 1.e-5)
    fac = (1. - exp(-tau_CI[2][1]))/tau_CI[2][1];
  else
    fac = 1.;

  Lambda = n_CI[2] * A_CI[2][1] * delE_CI[2][1] * fac;

  printf("[CI] tau[3]=%g %g\n", tau_CI[3][2], tau_CI[3][1]);
  printf("[CI] tau[2]=%g\n", tau_CI[2][1]);
  printf("[CI] y[]=%g %g %g\n", n1/nCI, n2/nCI, n3/nCI);
  printf("[CI] Lambda=%g %g %g\n", n2/nCI * A21 * delE_21 * fac, Lambda/nCI, Lambda/nH);

  return Lambda; // erg/sec per OI atom
}


double CarbonI_twolevel_cooling(double nCI, double ne, double nH,
                                double T_gas, double radiation_temperature,
                                int diagnostics_enabled)
{
  double g1, g2;

  double gamma21H;
  double gamma21e;
  double gamma12H;
  double gamma12e;

  double R12, R21;

  double A21;
  double delE_21;

  double T100, Lambda;

  T100 = T_gas/100.0;

  /* statistical weights */
  g1 = 1;  // 3P0
  g2 = 3;  // 3P1

  gamma21H = 1.6e-10 * pow(T100, 0.14);
  gamma21e = 3.0e-9;
  A21 = 7.9e-8;

  delE_21 = 3.312e-15; // erg Texc = 24 K

  gamma12H = (g2/g1) * gamma21H * exp(-delE_21/(BOLTZMANN * T_gas));
  gamma12e = (g2/g1) * gamma21e * exp(-delE_21/(BOLTZMANN * T_gas));

  R21 = A21 + nH * gamma21H  + gamma21e * ne;
  R12 = nH * gamma12H  + gamma12e * ne;

#ifdef CMB_HEATING
  double Texc, Q21;
  Texc = 23.6;
  Q21  = 1. / (exp(Texc/radiation_temperature) - 1.);

  R21 = A21 * (1. + Q21) + nH * gamma21H  + gamma21e * ne;
  R12 = g2/g1 * A21 *Q21 + nH * gamma12H  + gamma12e * ne;

  /*printf("R21=%g R21=%g A21=%g Q21=%g Q12=%g C21=%g C12=%g\n"
	 , R21, R12, A21, Q21, g2/g1*Q21
	 , nH * gamma21H  + gamma21e * ne, nH * gamma12H  + gamma12e * ne );*/
#endif


  double n1, n2;
  n1 = nCI * R21 / (R12 + R21);
  n2 = nCI * R12 / (R12 + R21);

  Lambda = n2 * A21 * delE_21;

  double t_freefall;
  t_freefall = sqrt(3.0*M_PI/32.0/GRAVITY/(1.2 * PROTONMASS * nH));  /*  in second */

  double hc;
  hc  = PLANCK*CLIGHT;

  double tau21;
  tau21 =
    g2/g1 * A21 * hc * hc * hc / (8 * M_PI * pow(delE_21,3))
    * (n1 - n2 * g1/g2) * t_freefall;

  if(tau21 > 1.e-5)
    Lambda = n2 * A21 * delE_21 * (1. - exp(-tau21))/tau21;
  else
    Lambda = n2 * A21 * delE_21;

#ifdef CMB_HEATING
  double fac;
  //fac = 1. - Q21 * (g2/g1 * n1/n2 - 1);
  fac = 1. - Q21 * (g2/g1 * R21/R12 - 1);

  Lambda *= fac;

  //printf("[CI]nH=%g Lambda=%g fac=%g f_red=%g\n", nH, Lambda, fac, (1. - exp(-tau21))/tau21);
#endif

  if(diagnostics_enabled)
    printf("[CI]nH=%g nCI=%g Lambda=%g tau21=%g n2=%g R12=%g R21=%g\n", nH, nCI, Lambda, tau21, n2, R12, R21);

  return Lambda; // erg/sec per OI atom
}


/***************************************************************************************/
/***************************************************************************************/

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
/* modified from original gaussj so that array[0..n] can be used */
static int gaussj5(double a[5][5], int n, double b[5][2], int m)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	int indxc[5],indxr[5],ipiv[5];

	if(n != 5 || m!= 2)
	  return 0;

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
		  return 0;
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
	return 1;
}

static int gaussj3(double a[3][3], int n, double b[3][2], int m)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	int indxc[3],indxr[3],ipiv[3];

	if(n != 3 || m!= 2)
	  return 0;

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
		  return 0;
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
	return 1;
}

#undef SWAP


#if defined(CONTINUUM_COOLING_MATSUKOBA) || defined(METAL_LINE_COOLING_SD93)
struct ChemistryContinuumOpacityTables {
public:
  ChemistryContinuumOpacityTables() = default;
  ChemistryContinuumOpacityTables(const ChemistryContinuumOpacityTables &) = delete;
  ChemistryContinuumOpacityTables &operator=(const ChemistryContinuumOpacityTables &) = delete;
  ~ChemistryContinuumOpacityTables();
  double calc_continuum_opacities(double T_gas,
    double rho_cgs, double l_Jeans, struct opacities *opac) const;

  static chemistry::ResourceOwner<ChemistryContinuumOpacityTables> create(const char *table_root) noexcept;
#ifdef METAL_LINE_COOLING_SD93
  double metal_cooling_rate(double temperature) const;
#endif
private:
#ifdef METAL_LINE_COOLING_SD93
  int read_metal_cooling(const char *table_root);
  double *SD93_Tgas = nullptr, *SD93_CoolRate = nullptr;
  double TMax_SD93 = 0., TMin_SD93 = 0.;
  int nT_SD93 = 0;
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
  int read_opacities(const char *table_root);
  static void lntrim(char *str);
  static int opacity_upper_index(const double *grid, int size, double value);
  void get_opacities(double T_gas,
    double rho_cgs, double *kappa_Rossland, double *kappa_Planck) const;
  void get_opacities_log(double T_gas,
    double rho_cgs, double *kappa_Rossland, double *kappa_Planck) const;

  int temperature_count{};
  int density_count{};
  double *temperature{};
  double *density{};
  double *rosseland{};
  double *planck{};
  double density_min{};
  double density_max{};
  double temperature_min{};
  double temperature_max{};

#endif
};

chemistry::ResourceOwner<ChemistryContinuumOpacityTables> ChemistryContinuumOpacityTables::create(const char *table_root) noexcept {
  chemistry::ResourceOwner<ChemistryContinuumOpacityTables> tables{new(std::nothrow) ChemistryContinuumOpacityTables{}};
  if(!tables) return {};
#ifdef CONTINUUM_COOLING_MATSUKOBA
  if(!tables->read_opacities(table_root)) {  return {}; }
#endif
#ifdef METAL_LINE_COOLING_SD93
  if(!tables->read_metal_cooling(table_root)) {  return {}; }
#endif
  return tables;
}
ChemistryContinuumOpacityTables::~ChemistryContinuumOpacityTables() {
#ifdef CONTINUUM_COOLING_MATSUKOBA
  free(temperature);
  free(density);
  free(rosseland);
  free(planck);
#endif
#ifdef METAL_LINE_COOLING_SD93
  free(SD93_Tgas);
  free(SD93_CoolRate);
#endif
}
void chemistry::ModelResourceDeleter::operator()(ChemistryContinuumOpacityTables *p) const noexcept { delete p; }

chemistry::ResourceOwner<ChemistryContinuumOpacityTables> chemistry::load_continuum_opacity(const char *root) noexcept {
  return ChemistryContinuumOpacityTables::create(root);
}


#endif

#ifdef METAL_LINE_COOLING_SD93
static const char fname_SD93[]="continuum_metal.dat";

int ChemistryContinuumOpacityTables::read_metal_cooling(const char *table_root){
  {
    FILE *fp;

    if(!(fp = open_chemistry_table(table_root, fname_SD93))){
      return 0;
    }

    char buf[200];
    char *tp;

    int count = 0;
    while(fgets(buf, sizeof(buf), fp) != NULL){
      double Tgas, cool_rate;
      int n=sscanf(buf, "%le %le", &Tgas, &cool_rate);
      if(n != 2) continue;

      count++;
    }

    rewind(fp);
    nT_SD93   = count;

    if(!(SD93_Tgas = static_cast<double *>(malloc(sizeof(double)*nT_SD93))) ||
       !(SD93_CoolRate = static_cast<double *>(malloc(sizeof(double)*nT_SD93)))){
      printf("Can't allocate memory.\n");
      free(SD93_Tgas);
      free(SD93_CoolRate);
      SD93_Tgas = SD93_CoolRate = NULL;
      fclose(fp);
      return 0;
    }

    count = 0;
    TMax_SD93= -1.e30;
    TMin_SD93 = 1.e30;

    while(fgets(buf, sizeof(buf), fp) != NULL){
      double Tgas, cool_rate;
      int n=sscanf(buf, "%le %le", &Tgas, &cool_rate);
      if(n != 2) continue;

      if(TMax_SD93 < Tgas) TMax_SD93 = Tgas;
      if(TMin_SD93 > Tgas) TMin_SD93 = Tgas;

      SD93_Tgas[count] = Tgas;
      SD93_CoolRate[count] = cool_rate;

      count++;
    }

    fclose(fp);
  }
  return 1;
}

double ChemistryContinuumOpacityTables::metal_cooling_rate(double T_gas) const{
  int i_T=0;
  double T_new = SD93_Tgas[0];
  while(T_gas > T_new){
    i_T++;
    T_new = SD93_Tgas[i_T];

    if(i_T == nT_SD93-1) break;
  }

  int i_T_left  = i_T - 1;
  if(i_T_left < 0) i_T_left = 0;
  if(i_T == nT_SD93 - 1) i_T_left = nT_SD93-1;

  double c_T   = (SD93_Tgas[i_T]-T_gas)    / (SD93_Tgas[i_T] -SD93_Tgas[i_T_left]);

  if(i_T == i_T_left)
    c_T = 0.;

  double cool_rate = c_T * SD93_CoolRate[i_T_left] + (1 - c_T) * SD93_CoolRate[i_T];

  if(T_gas > TMax_SD93){
    cool_rate *= pow(T_gas/TMax_SD93, -3);
  }

  if(T_gas < TMin_SD93){
    cool_rate *= exp(-pow(T_gas-TMin_SD93, 2)/200.);
  }

  return cool_rate;
}
double get_cooling_rate_SD93(const ChemistryContinuumOpacityTables *tables, double temperature) {
  return tables->metal_cooling_rate(temperature);
}

#endif



#ifdef CONTINUUM_COOLING_MATSUKOBA



// Immutable fit coefficients and their evaluation; no process-wide initializer.
class ContinuumEmission final {
public:
  static double lambda(double temperature);
  static double cooling(double temperature, double density,
                        const double *abundances, int diagnostics_enabled);
  static int output();
private:
enum ContCooling{
  cH_fb,
  cHM_fb,
#ifdef HELIUM_COOLING
  cHeI_fb,
#endif
  cHM_ff,
  cH_ff,
  cH2H2_CIE,
  cH2He_CIE,
  cEND
};

  static constexpr int N_reaction = cEND - cH_fb;
  inline static constexpr std::array<std::array<double, 5>, 7> A_Cont = [] {
    std::array<std::array<double, 5>, 7> A_Cont{};
  A_Cont[cH_fb][0]    =-22.410, A_Cont[cH_fb][1]    =-0.51833 , A_Cont[cH_fb][2]    = 2.3058e-3, A_Cont[cH_fb][3]    = 0.0000, A_Cont[cH_fb][4]    = 0.0000;
  A_Cont[cHM_fb][0]   =-31.211, A_Cont[cHM_fb][1]   = 1.1678  , A_Cont[cHM_fb][2]   =-3.0473e-2, A_Cont[cHM_fb][3]   = 0.0000, A_Cont[cHM_fb][4]   = 0.0000;
#ifdef HELIUM_COOLING
  A_Cont[cHeI_fb][0]  =-21.4171,A_Cont[cHM_fb][1]   = -0.995727, A_Cont[cHM_fb][2]  = 0.0714639, A_Cont[cHM_fb][3]   = 0.0000, A_Cont[cHM_fb][4]   = 0.0000;
#endif
  A_Cont[cHM_ff][0]   =-38.564, A_Cont[cHM_ff][1]   = 3.4681  , A_Cont[cHM_ff][2]   =-0.17880  , A_Cont[cHM_ff][3]   = 0.0000, A_Cont[cHM_ff][4]   = 0.0000;
  A_Cont[cH_ff][0]    =-28.030, A_Cont[cH_ff][1]    = 0.54235 , A_Cont[cH_ff][2]    =-5.3291e-3, A_Cont[cH_ff][3]    = 0.0000, A_Cont[cH_ff][4]    = 0.0000;
  A_Cont[cH2H2_CIE][0]= 75.232, A_Cont[cH2H2_CIE][1]=-1.4487e2, A_Cont[cH2H2_CIE][2]= 65.126   , A_Cont[cH2H2_CIE][3]=-12.304, A_Cont[cH2H2_CIE][4]= 0.84166;
  A_Cont[cH2He_CIE][0]= 87.458, A_Cont[cH2He_CIE][1]=-1.6412e2, A_Cont[cH2He_CIE][2]= 75.157   , A_Cont[cH2He_CIE][3]=-14.482, A_Cont[cH2He_CIE][4]= 1.0120;

    return A_Cont;
  }();
};

double ContinuumEmission::lambda(double T_gas){
  if(T_gas < 100.) return 0.;

  /************* DEBUGGING... [CHON200417] **************/
  double fac = 1.;
  if(T_gas > 1.e6){
    fac = sqrt(T_gas / 1.e6);
    T_gas = 1.e6;
  }

  double log_T  = log10(T_gas);
  double log_T2 = log_T  * log_T;
  double log_T3 = log_T2 * log_T;
  double log_T4 = log_T3 * log_T;

  double f[N_reaction];

  for(int i=0;i<N_reaction;i++){
    f[i] = A_Cont[i][0] + A_Cont[i][1] * log_T
      + A_Cont[i][2] * log_T2 + A_Cont[i][3] * log_T3 + A_Cont[i][4] * log_T4;
  }

  if(f[cH2H2_CIE] > -34)
    f[cH2H2_CIE] = -34;

  if(f[cH2He_CIE] > -34)
    f[cH2He_CIE] = -34;

  double Lambda = 4 * M_PI * (pow(10., f[cH_fb]) +  pow(10., f[cHM_ff]));
  Lambda *= fac;

  return Lambda;
}

double ContinuumEmission::cooling(double T_gas, double rho_cgs,
                                   const double *abund,
                                   int diagnostics_enabled){
  if(T_gas < 100.) return 0.;

  /************* DEBUGGING... [CHON200417] **************/
  double fac = 1.;
  if(T_gas > 1.e6){
    fac = sqrt(T_gas / 1.e6);
    T_gas = 1.e6;
  }

  double log_T  = log10(T_gas);
  double log_T2 = log_T  * log_T;
  double log_T3 = log_T2 * log_T;
  double log_T4 = log_T3 * log_T;

  double f[N_reaction];

  //printf("T_gas=%g\n", T_gas);
  for(int i=0;i<N_reaction;i++){
    f[i] = A_Cont[i][0] + A_Cont[i][1] * log_T
      + A_Cont[i][2] * log_T2 + A_Cont[i][3] * log_T3 + A_Cont[i][4] * log_T4;
  }

  if(f[cH2H2_CIE] > -34)
    f[cH2H2_CIE] = -34;

  if(f[cH2He_CIE] > -34)
    f[cH2He_CIE] = -34;

  double n_nuclei  = rho_cgs/PROTONMASS;

  double n_elec = n_nuclei * abund[sELEC];
  double n_HI   = n_nuclei * abund[sHI];
  double n_HII  = n_nuclei * abund[sHII];
  double n_H2I  = n_nuclei * abund[sH2I];
#ifdef HELIUM
  double n_HeI  = n_nuclei * abund[sHeI];
#else
  double n_HeI  = n_nuclei * XHe;
  //n_HeI  = 0.;
#endif


  double lambda_tot_thin = 4 * M_PI *
    ( pow(10., f[cH_fb] ) * n_HII * n_elec +
      pow(10., f[cHM_fb]) * n_HI  * n_elec +
#ifdef HELIUM_COOLING
      pow(10., f[cHeI_fb]) * n_HeII* n_elec +
#endif
      pow(10., f[cHM_ff]) * n_HI  * n_elec +
      pow(10., f[cH_ff]) * n_HII * n_elec +
      pow(10., f[cH2H2_CIE]) * n_H2I * n_H2I +
      pow(10., f[cH2He_CIE]) * n_H2I * n_HeI
      );

#ifdef ID_CHEM_DEBUG
  if(diagnostics_enabled)
    printf("cool_rate_cont=%g %g %g %g %g %g\n"
	   , pow(10, f[cH_fb] ) * n_HII * n_elec
	   , pow(10, f[cHM_fb] ) * n_HI * n_elec
	   , pow(10, f[cHM_ff]) * n_HI * n_elec
	   , pow(10, f[cH_ff])  * n_HII * n_elec
	   , pow(10, f[cH2H2_CIE]) * n_H2I * n_H2I
	   , pow(10, f[cH2He_CIE]) * n_HeI * n_H2I
	 );
#endif

  return lambda_tot_thin * fac;
}

int ContinuumEmission::output(){
  int iT;
  const int nT=500;

  FILE *fp;
  if(!(fp = fopen("test_cooling_cont.dat", "w"))){
    printf("Can't open the file.\n");
    return 0;
  }

  double T_gas, log_T, log_T2, log_T3, log_T4;

  for(iT=0;iT<nT;iT++){
    T_gas = 500.*exp(log(2.e5/500.)/nT * iT);
    log_T  = log10(T_gas);
    log_T2 = log_T  * log_T;
    log_T3 = log_T2 * log_T;
    log_T4 = log_T3 * log_T;

    int i;
    double f[N_reaction+1];

    for(i=0;i<N_reaction;i++){
      f[i] = A_Cont[i][0] + A_Cont[i][1] * log_T
	+ A_Cont[i][2] * log_T2 + A_Cont[i][3] * log_T3 + A_Cont[i][4] * log_T4;
    }

    double T3;
    T3 = T_gas/1.e3;
    f[N_reaction] = 1.e-27 * T3 * (1.4924 + 0.07815*T3 + 0.0063*T3*T3) * (1. - 0.1535*sqrt(T3));

    fprintf(fp
	    , "%g %g %g %g %g %g %g %g\n"
	    , T_gas, pow(10,f[0]), pow(10,f[1]), pow(10,f[2]), pow(10,f[3]), pow(10,f[4]), pow(10,f[5]),f[6]/4./M_PI);
  }

  fclose(fp);
  return 1;
}



double continuum_cooling_rate_thin_Lambda(double temperature) {
  return ContinuumEmission::lambda(temperature);
}
double continuum_cooling_rate_thin(double temperature, double density,
                                   const double *abundances, int diagnostics_enabled) {
  return ContinuumEmission::cooling(temperature, density, abundances, diagnostics_enabled);
}
int output_cooling_rate(void) { return ContinuumEmission::output(); }

double ChemistryContinuumOpacityTables::calc_continuum_opacities(double T_gas,
    double rho_cgs, double l_Jeans, struct opacities *opac) const {
  const ChemistryContinuumOpacityTables *tables = this;
  int ifunc;
  double tau_cont, f_reduce, x = 0.;
  double tau_Rossland=0., tau_Planck=0.;
  double kappa_Rossland=0., kappa_Planck=0.;

  if(rho_cgs / PROTONMASS > 1.e2){
    tables->get_opacities_log(T_gas, rho_cgs, &kappa_Rossland, &kappa_Planck);

    tau_Rossland = rho_cgs * kappa_Rossland * l_Jeans;
    tau_Planck   = rho_cgs * kappa_Planck   * l_Jeans;

    x        = tau_Planck * (1. + 3./4.*tau_Rossland);
    f_reduce = 1./ (1. + x);
    tau_cont = sqrt(tau_Rossland*tau_Planck);

    //printf("f_red=%g tau=%g %g %g l_Jeans=%g rho=%g kappa=%g %g\n"
    //	   , f_reduce, tau_Rossland, tau_Planck, tau_cont, l_Jeans, rho_cgs, kappa_Rossland, kappa_Planck);
  }else{
    f_reduce = 1.;
    tau_cont = 0.;
  }

  opac -> tau = tau_cont;
  opac -> tau_Rossland = tau_Rossland;
  opac -> tau_Planck   = tau_Planck;
  opac -> kappa_Rossland = kappa_Rossland;
  opac -> kappa_Planck   = kappa_Planck;
  opac -> x   = x;
  opac -> f_reduce = f_reduce;


  return f_reduce;
}

double calc_continuum_opacities(
    const ChemistryContinuumOpacityTables *tables, double T_gas,
    double rho_cgs, double l_Jeans, struct opacities *opac) { return tables->calc_continuum_opacities(T_gas, rho_cgs, l_Jeans, opac); }



static const char fname_Rossland[]="opacities_Rossland.dat";
static const char fname_Planck[]="opacities_Planck.dat";

int ChemistryContinuumOpacityTables::opacity_upper_index(const double *grid, int size, double value) {
  if(value <= grid[0]) return 0;
  if(value >= grid[size-1]) return size-1;

  int left = 0;
  int right = size - 1;
  while(right - left > 1) {
    int middle = left + (right - left) / 2;
    if(value > grid[middle])
      left = middle;
    else
      right = middle;
  }
  return right;
}

void ChemistryContinuumOpacityTables::lntrim(char *str) {
  char *p;
  p = strchr(str, '\n');
  if(p != NULL) {
    *p = '\0';
  }
}


int ChemistryContinuumOpacityTables::read_opacities(const char *table_root) {
  ChemistryContinuumOpacityTables *tables = this;
  double *Kappa_Tgas = NULL;
  double *Kappa_rho = NULL;
  double *Kappa_Rossland = NULL;
  double *Kappa_Planck = NULL;
  double RhoMin_tau = 0.;
  double RhoMax_tau = 0.;
  double TMin_tau = 0.;
  double TMax_tau = 0.;
  int n_T = 0;
  int n_rho = 0;
  FILE *fp_r = NULL;
  FILE *fp_p = NULL;

  {
    if(!(fp_r = open_chemistry_table(table_root, fname_Rossland))){
      goto fail;
    }

    if(!(fp_p = open_chemistry_table(table_root, fname_Planck))){
      goto fail;
    }

    char buf[200];
    char *tp;

    int count_T, count_rho;
    count_T = 0;
    while(fgets(buf, sizeof(buf), fp_r) != NULL){
      lntrim(buf);

      count_rho=0;

      tp = strtok(buf, " ");
      do{
	count_rho++;
      }while((tp = strtok(NULL, " ")));
      count_T++;
    }

    rewind(fp_p),rewind(fp_r);

    n_T   = count_T-1;
    n_rho = count_rho-1;

    if(!(Kappa_Tgas = static_cast<double *>(malloc(sizeof(double)*n_T)))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    if(!(Kappa_rho = static_cast<double *>(malloc(sizeof(double)*n_rho)))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    if(!(Kappa_Rossland = static_cast<double *>(malloc(sizeof(double)*(n_T)*(n_rho))))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    if(!(Kappa_Planck = static_cast<double *>(malloc(sizeof(double)*(n_T)*(n_rho))))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    count_T = count_rho = 0;
    RhoMax_tau = TMax_tau = -1.e30;
    RhoMin_tau = TMin_tau = 1.e30;

    double tau;

    while(fgets(buf, sizeof(buf), fp_r) != NULL){
      lntrim(buf);
      count_rho=0;

      tp = strtok(buf, " ");
      do{
	if(tp == NULL) break;

	tau = atof(tp);
	if(count_T > 0 && count_rho > 0)
	  Kappa_Rossland[n_rho * (count_T-1) + (count_rho-1)] = pow(10, tau);

	if(count_rho == 0 && count_T > 0){
	  double temp;

	  temp = pow(10, tau);
	  Kappa_Tgas[count_T-1] = temp;

	  if(TMax_tau > pow(10,tau)) TMax_tau = temp;
	  if(TMin_tau < pow(10,tau)) TMin_tau = temp;
	}

	if(count_T == 0){
	  double rho;

	  //rho = pow(10, tau) / PROTONMASS;
	  rho = pow(10, tau);
	  Kappa_rho[count_rho] = rho;

	  if(RhoMax_tau > pow(10,tau)) RhoMax_tau = rho;
	  if(RhoMin_tau < pow(10,tau)) RhoMin_tau = rho;
	}

	count_rho++;
      }while((tp = strtok(NULL, " ")));

      count_T++;
    }

    count_T=0;
    while(fgets(buf, sizeof(buf), fp_p) != NULL){
      lntrim(buf);
      count_rho=0;

      tp = strtok(buf, " ");
      do{
	if(tp == NULL) break;

	tau = atof(tp);
	if(count_T > 0 && count_rho > 0)
	  Kappa_Planck[n_rho * (count_T-1) + (count_rho-1)] = pow(10, tau);

	count_rho++;
      }while((tp = strtok(NULL, " ")));
      count_T++;
    }

    fclose(fp_p),fclose(fp_r);
    fp_p = fp_r = NULL;
  }

  tables->temperature_count = n_T;
  tables->density_count = n_rho;
  tables->temperature = Kappa_Tgas;
  tables->density = Kappa_rho;
  tables->rosseland = Kappa_Rossland;
  tables->planck = Kappa_Planck;
  tables->density_min = RhoMin_tau;
  tables->density_max = RhoMax_tau;
  tables->temperature_min = TMin_tau;
  tables->temperature_max = TMax_tau;
  return 1;

fail:
  if(fp_p != NULL)
    fclose(fp_p);
  if(fp_r != NULL)
    fclose(fp_r);
  free(Kappa_Tgas);
  free(Kappa_rho);
  free(Kappa_Rossland);
  free(Kappa_Planck);
  return 0;
}

void ChemistryContinuumOpacityTables::get_opacities(double T_gas,
    double rho_cgs, double *kappa_Rossland, double *kappa_Planck) const {
  const ChemistryContinuumOpacityTables *tables = this;
  const double *Kappa_Tgas = tables->temperature;
  const double *Kappa_rho = tables->density;
  const double *Kappa_Rossland = tables->rosseland;
  const double *Kappa_Planck = tables->planck;
  const int n_T = tables->temperature_count;
  const int n_rho = tables->density_count;
  const double TMax_tau = tables->temperature_max;
  int i_T, i_rho;

  i_T = opacity_upper_index(Kappa_Tgas, n_T, T_gas);
  i_rho = opacity_upper_index(Kappa_rho, n_rho, rho_cgs);

  int i_T_left, i_rho_left;

  i_T_left   = i_T - 1;
  i_rho_left = i_rho - 1;
  if(i_T_left < 0) i_T_left = 0;
  if(i_rho_left < 0) i_rho_left = 0;

  if(i_T == n_T - 1) i_T_left = n_T-1;
  if(i_rho == n_rho - 1) i_rho_left = n_rho-1;


  double c_T, c_rho;
  double kappa_r, kappa_p;

  c_T   = (Kappa_Tgas[i_T]-T_gas)    / (Kappa_Tgas[i_T] -Kappa_Tgas[i_T_left]);
  c_rho = (Kappa_rho[i_rho]-rho_cgs) / (Kappa_rho[i_rho] -Kappa_rho[i_rho_left]);

  if(i_T == i_T_left)
    c_T = 0.;

  if(i_rho == i_rho_left)
    c_rho = 0.;

  kappa_r =
    c_T       * c_rho       * Kappa_Rossland[i_T_left * n_rho + i_rho_left] +
    c_T       * (1 - c_rho) * Kappa_Rossland[i_T_left * n_rho + i_rho] +
    (1 - c_T) * c_rho       * Kappa_Rossland[i_T      * n_rho + i_rho_left] +
    (1 - c_T) * (1 - c_rho) * Kappa_Rossland[i_T      * n_rho + i_rho];

  kappa_p =
    c_T       * c_rho       * Kappa_Planck[i_T_left * n_rho + i_rho_left] +
    c_T       * (1 - c_rho) * Kappa_Planck[i_T_left * n_rho + i_rho] +
    (1 - c_T) * c_rho       * Kappa_Planck[i_T      * n_rho + i_rho_left] +
    (1 - c_T) * (1 - c_rho) * Kappa_Planck[i_T      * n_rho + i_rho];

  if(T_gas > TMax_tau){
    kappa_p *= pow(T_gas/TMax_tau, -3);
  }

  *kappa_Rossland = kappa_r;
  *kappa_Planck   = kappa_p;
}

void ChemistryContinuumOpacityTables::get_opacities_log(double T_gas,
    double rho_cgs, double *kappa_Rossland, double *kappa_Planck) const {
  const ChemistryContinuumOpacityTables *tables = this;
  const double *Kappa_Tgas = tables->temperature;
  const double *Kappa_rho = tables->density;
  const double *Kappa_Rossland = tables->rosseland;
  const double *Kappa_Planck = tables->planck;
  const int n_T = tables->temperature_count;
  const int n_rho = tables->density_count;
  int i_T, i_rho;

  i_T = opacity_upper_index(Kappa_Tgas, n_T, T_gas);
  i_rho = opacity_upper_index(Kappa_rho, n_rho, rho_cgs);

  int i_T_left, i_rho_left;

  i_T_left   = i_T - 1;
  i_rho_left = i_rho - 1;
  if(i_T_left < 0) i_T_left = 0;
  if(i_rho_left < 0) i_rho_left = 0;

  if(i_T == n_T - 1) i_T_left = n_T-1;
  if(i_rho == n_rho - 1) i_rho_left = n_rho-1;


  double c_T, c_rho;
  double kappa_r, kappa_p;

  c_T   = (log10(Kappa_Tgas[i_T])-log10(T_gas))    / (log10(Kappa_Tgas[i_T])  - log10(Kappa_Tgas[i_T_left]));
  c_rho = (log10(Kappa_rho[i_rho])-log10(rho_cgs)) / (log10(Kappa_rho[i_rho]) - log10(Kappa_rho[i_rho_left]));

  if(i_T == i_T_left)
    c_T = 0.;

  if(i_rho == i_rho_left)
    c_rho = 0.;

  kappa_r = pow(10,
		c_T       * c_rho       * log10(Kappa_Rossland[i_T_left * n_rho + i_rho_left]) +
		c_T       * (1 - c_rho) * log10(Kappa_Rossland[i_T_left * n_rho + i_rho]) +
		(1 - c_T) * c_rho       * log10(Kappa_Rossland[i_T      * n_rho + i_rho_left]) +
		(1 - c_T) * (1 - c_rho) * log10(Kappa_Rossland[i_T      * n_rho + i_rho]));

  kappa_p = pow(10,
		c_T       * c_rho       * log10(Kappa_Planck[i_T_left * n_rho + i_rho_left]) +
		c_T       * (1 - c_rho) * log10(Kappa_Planck[i_T_left * n_rho + i_rho]) +
		(1 - c_T) * c_rho       * log10(Kappa_Planck[i_T      * n_rho + i_rho_left]) +
		(1 - c_T) * (1 - c_rho) * log10(Kappa_Planck[i_T      * n_rho + i_rho]));

  *kappa_Rossland = kappa_r;
  *kappa_Planck   = kappa_p;

}
#endif

double  beta_esc(double tau_L,double tau_C){
  double beta_esc;

  if(tau_L < 0.)
    beta_esc=1.;
  else if(tau_L < 1.e-5)
    beta_esc = exp(-tau_C);
  else
    beta_esc = exp(-tau_C)*(1.-exp(-tau_L))/tau_L;

  return beta_esc;
}


double H_criticaldensity(double xnH, double T_K,double xNc_H, double y_H,double y_e){
  double n_crit, n_crit_2,n_crit_3;
  double ncrit_e_2,ncrit_H_2,ncrit_e_3,ncrit_H_3;
  double A_21;
  double gamma_H_e_21,gamma_H_H_21,gamma_H_e_31,gamma_H_H_31, gamma_H_e_32,gamma_H_H_32;
  double esc_21,esc_31,esc_32;
  double tau_line_21,tau_line_31,tau_line_32;
  double fraction_1,fraction_2,fraction_3;
  double v_D,Q,r_up,r_low;
  double n_e,ln57T,C_2s2p,ratio2s_2p;
  double g_2s=2.,g_2p=6.,g1=2.,g2=8.,g3=18.;
  double A_2s1s=8.23,A_2p1s=6.27e8,A_31=5.575e7,A_32=4.410e7;
  double DE_21=10.198,DE_31=12.088,DE_32=1.890;
  double alpha_21=24.0,alpha_31=22.0,alpha_32=67.0;
  double chi_21=0.28,chi_31=0.37,chi_32=0.30;
  double f_21=0.4162,f_31=7.910e-2,f_32=0.6407;
  double xk_B=1.381e-16,h_Pla=6.626e-27,c_ph=2.998e10;

  n_e        = xnH*y_e;
  ln57T      = log(5.7*T_K);
  C_2s2p     = 6.21e-4*ln57T*(1.+0.78/ln57T)*n_e/sqrt(T_K);
  ratio2s_2p = g_2s/g_2p * C_2s2p/(C_2s2p+A_2s1s);
  A_21       = A_2p1s / (1.+ratio2s_2p);

  //- optical depth for line emission (level 2->1)
  v_D = sqrt(2.*xk_B*T_K/PROTONMASS);
  tau_line_21 = A_21/8./M_PI * pow(h_Pla*c_ph/(DE_21*eV_to_erg),3.)
    * g2/g1 * xNc_H/v_D;
  //- escape probability
  esc_21 = beta_esc(tau_line_21,0.);


  //- collision with electron
  r_up  = 2.;
  r_low = 1.;
  Q     = DE_21*eV_to_erg / (xk_B*T_K);
  gamma_H_e_21 = 1.e-8 * pow(r_low*r_low/(r_up*r_up-r_low*r_low),1.5)
    * r_low*r_low*r_low/r_up/r_up * alpha_21 * sqrt(Q*(Q+1.))/(Q+chi_21);
  ncrit_e_2 = esc_21*A_21 / gamma_H_e_21;
  //- collision with atomic hydrogen
  gamma_H_H_21 = 7.89e-15 * (r_low*r_low/r_up/r_up)
    * pow(1./r_low/r_low-1./r_up/r_up,-2.)
    * f_21 * sqrt(T_K)
    * (1.+1.27e-5    /(1./r_low/r_low-1./r_up/r_up)     *T_K)
    / (1.+4.76e-17*pow(1./r_low/r_low-1./r_up/r_up,-2.0)*T_K);
  ncrit_H_2 = esc_21*A_21 / gamma_H_H_21;

  n_crit_2  = 1. / (y_H/ncrit_H_2 + y_e/ncrit_e_2);

  fraction_2 = xnH/(n_crit_2+xnH) * g2/g1 * exp(-Q);
  fraction_1 = 1. - fraction_2;

  if(n_crit_2 < xnH){
    tau_line_21 = A_21/8./M_PI * pow(h_Pla*c_ph/(DE_21*eV_to_erg), 3.0)
      * (g2/g1*fraction_1-fraction_2) * (xNc_H/v_D);
    esc_21    = beta_esc(tau_line_21,0.);
    ncrit_e_2 = esc_21*A_21 / gamma_H_e_21;
    ncrit_H_2 = esc_21*A_21 / gamma_H_H_21;

    n_crit_2   = 1. / (y_H/ncrit_H_2 + y_e/ncrit_e_2);
    fraction_2 = xnH/(n_crit_2+xnH) * g2/g1 * exp(-Q);
    fraction_1 = 1. - fraction_2;
  }

  tau_line_31 = (A_31/(8.*M_PI)) * pow(PLANCK*c_ph/(DE_31*eV_to_erg),3) * (g3/g1)*fraction_1 * (xNc_H/v_D);
  tau_line_32 = (A_32/(8.*M_PI)) * pow(PLANCK*c_ph/(DE_32*eV_to_erg),3) * (g3/g2)*fraction_2 * (xNc_H/v_D);

  esc_31 = beta_esc(tau_line_31,0.);
  esc_32 = beta_esc(tau_line_32,0.);

  r_up  = 3.;
  r_low = 1.;
  Q     = DE_31*eV_to_erg / (xk_B*T_K);
  gamma_H_e_31 = 1.e-8 * pow(r_low*r_low/(r_up*r_up-r_low*r_low),1.5) * pow(r_low,3)/pow(r_up,2) * alpha_31 * sqrt(Q*(Q+1.))/(Q+chi_31);

  r_low = 2.;
  Q     = DE_32*eV_to_erg / (xk_B*T_K);
  gamma_H_e_32 = 1.e-8 * pow(r_low*r_low/(r_up*r_up-r_low*r_low),1.5) * pow(r_low,3)/pow(r_up,2) * alpha_32 * sqrt(Q*(Q+1.))/(Q+chi_32);

  ncrit_e_3 = (esc_31*A_31+esc_32*A_32) / (gamma_H_e_31+gamma_H_e_32);

  gamma_H_H_31 =
    7.89e-15 * (r_low*r_low/r_up/r_up)
    * pow(1./r_low/r_low-1./r_up/r_up,-2)
    * f_31 * sqrt(T_K)
    * (1. + 1.27e-5 / (1./(r_low*r_low) - 1./(r_up*r_up)) * T_K)
    / (1. + 4.76e-17 * pow(1./r_low/r_low-1./r_up/r_up, -2.0) * T_K);

  r_low = 2.;
  gamma_H_H_32 = 7.89e-15 * (r_low*r_low/r_up/r_up)
    * pow(1./r_low/r_low-1./r_up/r_up,-2.0)
    * f_32 * sqrt(T_K)
    * (1.+1.27e-5/(1./r_low/r_low-1./r_up/r_up)*T_K)
    / (1.+4.76e-17*pow(1./r_low/r_low-1./r_up/r_up,-2.0)*T_K);
  ncrit_H_3 = (esc_31*A_31+esc_32*A_32) / (gamma_H_H_31+gamma_H_H_32);

  n_crit_3  = 1. / (y_H/ncrit_H_3 + y_e/ncrit_e_3);

  fraction_3 = fraction_2 * (g3/g2) * exp(-Q);
  fraction_1 = 1. - (fraction_2+fraction_3);
  // in the case of xnH > xnH_crit_3
  if(n_crit_3 < xnH){
    tau_line_31 = (A_31/(8.*M_PI)) * pow(PLANCK*c_ph/(DE_31*eV_to_erg),3.0)
      * ((g3/g1)*fraction_1-fraction_3) * (xNc_H/v_D);
    tau_line_32 = (A_32/(8.*M_PI)) * pow(PLANCK*c_ph/(DE_32*eV_to_erg),3.0)
      * ((g3/g2)*fraction_2-fraction_3) * (xNc_H/v_D);
    esc_31    = beta_esc(tau_line_31,0.);
    esc_32    = beta_esc(tau_line_32,0.);
    ncrit_e_3 = (esc_31*A_31+esc_32*A_32) / (gamma_H_e_31+gamma_H_e_32);
    ncrit_H_3 = (esc_31*A_31+esc_32*A_32) / (gamma_H_H_31+gamma_H_H_32);

    n_crit_3  = 1. / (y_H/ncrit_H_3 + y_e/ncrit_e_3);
  }

  n_crit = n_crit_3;

  return n_crit;
}

double fesc_Lya(double xnH,double T_K,double y_H,double y_e,double rr,double vel,
                double *cooling_rate_Lya, double *cooling_rate_2ph,
                int diagnostics_enabled){
  double fesc_Lya;
  double p_21,A_21,A_2ph,C_21,C_12;
  double r_2s2p,C_2s2p;
  double N_ref,tau_21,sigma_Ly,T_4,V_ratio;
  double ln57T,n_H,n_e;
  double gamma_H_H_21,gamma_H_e_21,Q_21;
  double g_2=8., g_1=2.;
  double g_2p=6.,g_2s=2.;
  double A_2p1s=6.27e8,A_2s1s=8.23;
  double r_up=2.,r_low=1.;
  double f_21=0.4162,alpha_21=24.0;
  double DE_21=10.198,chi_21=0.28;
  double boltz=1.3807e-16,c_ph=2.9979e10;

  n_H    = y_H*xnH;
  n_e    = y_e*xnH;
  ln57T  = log(5.7*T_K);
  C_2s2p = 6.21e-4*ln57T*(1.+0.78/ln57T)*n_e/sqrt(T_K);
  r_2s2p = g_2s/g_2p * C_2s2p/(C_2s2p+A_2s1s);

  Q_21   = DE_21*eV_to_erg / (boltz*T_K);
  gamma_H_H_21 = 7.89e-15 * (r_low*r_low/r_up/r_up)
    * pow(1./r_low/r_low-1./r_up/r_up, -2.0)
    * f_21 * sqrt(T_K)
    * (1.+1.27e-5 *pow(1./r_low/r_low-1./r_up/r_up,-1.0)*T_K)
    / (1.+4.76e-17*pow(1./r_low/r_low-1./r_up/r_up,-2.0)*T_K*T_K);

  gamma_H_e_21 = 1.e-8 * pow(r_low*r_low/(r_up*r_up-r_low*r_low),1.5)
    * r_low*r_low*r_low/r_up/r_up * alpha_21 *
    sqrt(Q_21*(Q_21+1.))/(Q_21+chi_21);

  A_21  = A_2p1s * 1.     / (1.+r_2s2p);
  A_2ph = A_2s1s * r_2s2p / (1.+r_2s2p);
  C_21  = n_e*gamma_H_e_21 + n_H*gamma_H_H_21;
  C_12  = C_21 * g_2 / g_1 * exp(-Q_21);

  p_21  = A_21 / (A_21+A_2ph+C_21);


  T_4      = T_K/1.e4;
  sigma_Ly = 5.9e-14/sqrt(T_4);
  tau_21   = sigma_Ly*n_H*rr;

  N_ref = 15.*pow(tau_21/pow(10.,5.5), 1.0/3.0);

  V_ratio = vel/c_ph;
  fesc_Lya = pow(p_21, N_ref)/(1.+N_ref)*exp(-3.*N_ref*V_ratio);

  double r_21;
  double cool_A21, cool_2ph;

  r_21 = C_12 / (A_21*fesc_Lya + A_2ph + C_21);

  cool_A21 = -A_21  * DE_21 * eV_to_erg * fesc_Lya * n_H * r_21/(1. + r_21);
  cool_2ph = -A_2ph * DE_21 * eV_to_erg            * n_H * r_21/(1. + r_21);

#ifdef ID_CHEM_DEBUG
  if(diagnostics_enabled)
    printf("[HI cooling] nH=%g HI=%g e=%g C_21=%g r21=%g cool_A21=%g\n"
	   , n_H, y_H, y_e, C_21, r_21, cool_A21);
#endif


  //printf("C_21=%g A_21=%g A_2ph=%g\n", C_21, A_21, A_2ph);
  //printf("tau21=%g N_ref=%g exp=%g p_21=%g r_21=%g r_2s2p=%g C_2s2p=%g\n", tau_21, N_ref, exp(-3.*N_ref*V_ratio), p_21, r_21, r_2s2p, C_2s2p);
  //printf("cool_A21=%g cool_2ph=%g A_21=%g A_2ph=%g\n", cool_A21, cool_2ph, A_21, A_2ph);

  *cooling_rate_Lya = cool_A21;
  *cooling_rate_2ph = cool_2ph;

  return fesc_Lya;
}
#endif
