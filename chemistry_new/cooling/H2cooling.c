#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "allvars.h"
#include "proto.h"
#include "chemistry.h"
#include "chemistry_constants.h"
#include "chemistry_diagnostics.h"
#include "H2cooling.h"

double H2line_LTE(double T);
double H2cool_LTEfit(double T);

double H2cool_G15(double xnH,double T_K,double y_H,double y_H2,double y_e,double y_Hp);
double H2cool_LowD_G15(double xnH,double T_K,double y_H,double y_H2,double y_e,double y_Hp,double y_He);

/*-----------------------------------------------------------------------
  ! LTE H2 cooling fit by Glover(2015) and fesc fit by Omukai-san
  !
  !     xLmbd_H2: H2 cooling rate per unit volume
  !     skeMeMo:I use this for whole region of density
  !             Of course, it is mainly for n > 1.d8
  !-----------------------------------------------------------------------*/
double H2_cooling_G15(double nH, double T_K,double yHI,double yH2,double yHeI,double yHII, double y_e){
  double xLmbd_H2;
  double T_rad = 0.;

  if(T_K > 2.e4)
    T_K = 2.e4;

  //lineがthickになった後に、opt-thin low-density cooling rateが現れないように場合分け
  if(nH < 1.e8){
    //double xLdH2_GP12 = H2_cooling_GP12(nH, T_K, yHI, yH2, yHeI, yHII, y_e);
    double xLdH2_n0  = H2cool_LowD_G15(nH, T_K, yHI, yH2, y_e, yHII, yHeI);

    //if(FlagPrint)
    //      printf("H2_cool_rate: GP12=%g G15=%g\n", xLdH2_GP12, xLdH2_n0);

    double xLdH2_LTE = H2cool_LTEfit(T_K) * nH * yH2;
    xLmbd_H2  = xLdH2_LTE / (1. + xLdH2_LTE/xLdH2_n0);
  }else{
    double xLdH2_LTE = H2cool_LTEfit(T_K) * nH * yH2;
    xLmbd_H2  = xLdH2_LTE;
  }

  return -1 * xLmbd_H2;
}




double H2cool_LTEfit(double T_K){
  double xLmbd_H2;
  double xLmbd,xLmbdr=0.;
  double xlgT3,a;
  
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
    double ar=-20.584225
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
double fesc_fit(double T_K, double xNc_H2){
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


double fesc_unity(double T_K, double xNc_H2){
  return 1.;
}


double H2_cooling_GP12(double xnH,double T_K,double y_H,double y_H2,double y_He,double y_Hp,double y_e)
{
  double cool_rate_G08_H2H, cool_rate_G08_H2H2, cool_rate_G08_H2He, cool_rate_G08_H2p, cool_rate_G08_H2e;
  double cool_rate_LOW_H2, cool_rate_LTE_H2, cool_rate_H2;

  /* At T > 10000K, assume rate same as at T=10000K */
  if(T_K > 1.e4)
    T_K = 1.e4;

  double T3       = T_K/1.e3;
  double logT3    = log10(T3);
  double logT3_p2 = logT3    * logT3;
  double logT3_p3 = logT3_p2 * logT3;
  double logT3_p4 = logT3_p3 * logT3;
  double logT3_p5 = logT3_p4 * logT3;

  double log10T = log10(T_K);

  /* ----- Glover&Abel'08 (Table 8) : H2 cooling rate in low-density limit, assumes 3:1 o:p ratio */

  /* ----- H2-H : based on Wrathmall&Flower'07 */
  if(T_K > 1.e3){
    cool_rate_G08_H2H = pow(10.0, -24.311209
			    + 4.6450521 * logT3
			    - 3.7209846 * logT3_p2
			    + 5.9369081 * logT3_p3
			    - 5.5108047 * logT3_p4
			    + 1.5538288 * logT3_p5);
  }else if(T_K > 1.e2){
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
  cool_rate_G08_H2H2 = pow(10.0, -23.962112
			   + 2.09433740  * logT3
			   - 0.77151436  * logT3_p2
			   + 0.43693353  * logT3_p3
			   - 0.14913216  * logT3_p4
			   - 0.033638326 * logT3_p5);
  


  /* ----- H2-He : based on Flower+'98, Balakrishnan+'99 */
  cool_rate_G08_H2He = pow(10.0, -23.689237
			   + 2.1892372  * logT3
			   - 0.81520438 * logT3_p2
			   + 0.29036281 * logT3_p3
			   - 0.16596184 * logT3_p4
			   + 0.19191375 * logT3_p5);
  
  
  /* ----- H2-H+ : based on Gerlich'90 and Krstic'02; Fit accurate to within 5%. */
  cool_rate_G08_H2p = pow(10.0, -21.716699
			  + 1.3865783   * logT3
			  - 0.37915285  * logT3_p2
			  + 0.11453688  * logT3_p3
			  - 0.23214154  * logT3_p4
			  + 0.058538864 * logT3_p5);  
  
  /* ----- H2-e : based on Draine+'83*/
  if(T_K > 2.e2){
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
  cool_rate_H2 = (cool_rate_G08_H2H  * y_H
		  + cool_rate_G08_H2H2 * y_H2
		  + cool_rate_G08_H2He * y_He
		  + cool_rate_G08_H2p  * y_Hp
		  + cool_rate_G08_H2e  * y_e);

  double xLmbd_H2 = cool_rate_H2 * xnH * xnH * y_H2;
  
  return xLmbd_H2;
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


double HD_cooling_GP12(double nH, double T_gas, double yHI, double yHDI)
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
