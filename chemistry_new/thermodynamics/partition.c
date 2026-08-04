#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "allvars.h"
#include "chemistry.h"
#include "chemistry_constants.h"
#include "partition.h"


/* rossi partition function for H2 */
const double ra0 = 4.55642e-1;
const double ra1 = 6.98477e-1;
const double ra2 = -6.10627e-2;
const double ra3 = -1.26614e-2;
const double ra4 = 2.21442e-3;
const double ra5 = 2.3324e-1;
const double ra6 = 1.66574;


/* irwin partition function for H2 */
const double  b0 =  1.69179;
const double  b1 = -1.7227;
const double  b2 =  7.98033e-1;
const double  b3 = -1.57089e-1;
const double  b4 = -5.35313e-1;
const double  b5 =  1.75818;
const double  b6 = -2.63895;
const double  b7 =  1.35708;

double E_H2_BFM(int iv,int J);

void check_partition_functions_H2(void){
  FILE *fp = fopen("partition_functions_H2.dat", "w");
  if(fp == NULL)
    return;
  
  const int imax=1000;
  int i;
  double T_max, T_min, dln_T;
  double T_gas, p1, p2, p3, p4;

  T_max=1.e7;
  T_min=0.1;

  dln_T = log(T_max/T_min)/(double)imax;

  for(i=0;i<imax;i++){
    T_gas = T_min * exp(dln_T * i);

    p1 = partition_H2_Rossi1985(T_gas);
    p2 = partition_H2_Irwin1987(T_gas);
    p3 = partition_H2_EBM1989(T_gas);
    p4 = partition_H2_EBM1989_fitted_by_Matsukoba(T_gas);

    fprintf(fp, "%g %g %g %g %g\n", T_gas, p1, p2, p3, p4);
  }

  fclose(fp);
}

double partition_H2_Rossi1985(double T_gas)
{
  double z_H2, Z, lnQ;
  
  /* rossi partition function 1985 */    
  Z = T_gas/1000.0;
  
  lnQ = ra0*log(Z) + ra1*Z/2.0 + ra2*Z*Z/6.0 + ra3*Z*Z*Z/12.0 + ra4*Z*Z*Z*Z/20.0 - ra5/Z + ra6;
  
  z_H2=exp(lnQ);
  
  if(T_gas > 6000.0) /* rossi fit diverges for T>1.0e4, my new fit */
    z_H2 = 1.0635e-5 * pow(T_gas, 1.8) * exp(-T_gas/5.0e4);

  if(T_gas < 600.0) /* rossi fit diverges for T>1.0e4, my new fit */
    z_H2 = 1.;

  return z_H2;
}

double partition_H2_Irwin1987(double T_gas)
{
  /* Irwin partition function 1987 */    
  double Z = log10(5040.0/T_gas);  
  double logQ = b0 + b1*Z + b2*Z*Z + b3*Z*Z*Z + b4*pow(Z, 4) + b5*pow(Z, 5) + b6*pow(Z, 6) + b7*pow(Z, 7); 
  double z_H2 = pow(10.0,logQ);

  return z_H2;
}

double partition_H2_EBM1989(double T_K){
  double pf;
  const int iv_max=5,j_max=25;
  double ET[iv_max+1][j_max+1];
  const double ET_diss=5.196e4;
  double z_p,z_o;
  int iv,j;

  for(iv=0;iv<=iv_max;iv++)
    for(j=0;j<=j_max;j++)
      ET[iv][j] = E_H2_BFM(iv,j);

  z_p = z_o = 0.;
  for(iv=0;iv<=iv_max;iv++){
    for(j=0;j<=j_max;j++){
      if(ET[iv][j]-ET[0][0] <= ET_diss){
         if(j%2 == 0) 
	   z_p += (2*j+1) * exp(-(ET[iv][j]-ET[0][0])/T_K);
	 else
	   z_o += (2*j+1) * exp(-(ET[iv][j]-ET[0][1])/T_K);
      }
    }
  }
   

  pf = 0.25*z_p + 0.75*z_o;
  
  return pf;
}

double partition_H2_EBM1989_fitted_by_Matsukoba(double T_K){
  double pf_z_H2;
  double a0,a1,a2,a3,a4,a5;
  double lgT1,lgT2,lgT3,lgT4,lgT5;
  double lg_pf;
    
  a0 = -20.2721571192088;
  a1 =  33.6402868140006;
  a2 = -21.086890724525;
  a3 =  6.27559221791495;
  a4 = -0.871080363948635;
  a5 =  0.0455727749861326;

  lgT1 = log10(T_K);
  lgT2 = lgT1*lgT1;
  lgT3 = lgT2*lgT1;
  lgT4 = lgT3*lgT1;
  lgT5 = lgT4*lgT1;

  lg_pf = a0 + a1*lgT1 + a2*lgT2 + a3*lgT3 + a4*lgT4 + a5*lgT5;
  pf_z_H2 = pow(10.0,lg_pf);

  if(pf_z_H2 < 2.5)  pf_z_H2 = 2.5;
  if(pf_z_H2 > 2.e3) pf_z_H2 = 2.e3;

  return pf_z_H2;
}



double E_H2_BFM(int iv,int J){
  /* ------------------------------------------------- !
  !  Borysow, Frmmhold & Moraldi 1989, ApJ, 336, 495  !
  !    Equation (A11) (in K)                          !
  ! ------------------------------------------------- */


  double E_BFM;

  double vv,rrot,Ev,Ev0,Ev1,Ev2,Ev3,Ev4,Bv,Bv0,Bv1,Bv2,Bv3,Bv4;
  double Dv,Dv0,Dv1,Dv2,Dv3,Dv4,Fv,Fv0,Fv1,Fv2,Fv3,Fv4,Gv,Gv0,Gv1,Gv2,Gv3,Gv4;
  double Hv,Hv0,Hv1,Hv2,Hv3,Hv4,Ov,Ov0,Ov1,Ov2,Ov3,Ov4,Pv,Pv0,Pv1,Pv2,Pv3,Pv4;

  Ev0=0.38496;
  Ev1=-0.04609;
  Ev2=0.00178;
  Ev3=-7.e-5;
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

  vv=iv+0.5;
  Ev=Ev0+Ev1*vv+Ev2*pow(vv,2)+Ev3*pow(vv,3)+Ev4*pow(vv,4);
  Bv=Bv0+Bv1*vv+Bv2*pow(vv,2)+Bv3*pow(vv,3)+Bv4*pow(vv,4);
  Dv=Dv0+Dv1*vv+Dv2*pow(vv,2)+Dv3*pow(vv,3)+Dv4*pow(vv,4);
  Fv=Fv0+Fv1*vv+Fv2*pow(vv,2)+Fv3*pow(vv,3)+Fv4*pow(vv,4);
  Gv=Gv0+Gv1*vv+Gv2*pow(vv,2)+Gv3*pow(vv,3)+Gv4*pow(vv,4);
  Hv=Hv0+Hv1*vv+Hv2*pow(vv,2)+Hv3*pow(vv,3)+Hv4*pow(vv,4);
  Ov=Ov0+Ov1*vv+Ov2*pow(vv,2)+Ov3*pow(vv,3)+Ov4*pow(vv,4);
  Pv=Pv0+Pv1*vv+Pv2*pow(vv,2)+Pv3*pow(vv,3)+Pv4*pow(vv,4);

  rrot = (double)J*(J+1);
  E_BFM=
    1.43879*(-Ev*1.e5
	     +Bv*rrot
	     -Dv*1.e-2*pow(rrot,2)
	     +Fv*1.e-5*pow(rrot,3)
	     -Gv*1.e-8*pow(rrot,4)
	     +Hv*1.e-11*pow(rrot,5)
	     -Ov*1.e-14*pow(rrot,6)
	     +Pv*1.e-17*pow(rrot,7));

  return E_BFM;
}




double partition_H(double T_gas)
{
  double z_H = 2.0;
  
  if(T_gas > 1.0e4){
    z_H=0.0;
    int nlevel;
    for(nlevel=1;nlevel<=10;nlevel++){
      double dn = (double) nlevel;
      double En = (1.0-1.0/(dn*dn)) * 13.598 * eV_to_erg;
      z_H += 2.0 * dn*dn * exp(- En/(BOLTZMANN * T_gas));
    }
  }
  
  return z_H;
}
  

double pf_z_H(double T_K){
  double pf;
  double dE;
  const int level=10;
  double g[level],E[level],f_l[level];
  int i;

  for(i=1;i<=level;i++){
    g[i] = 2.*i*i;
    E[i] = 13.58/((double)i*i)*eV_to_erg;
  }

  for(i=1,pf=0.;i<=level;i++){
    dE = E[1] - E[i];
    f_l[i] = g[i]*exp(-dE/(BOLTZMANN*T_K));
    pf += f_l[i];
  }

  return pf;
}

double linear(double *xa,double *ya,int m,double x);

double pf_z_H2p(double T_K){
  double xTBC[] = 
    {1.00000e0, 1.30000e0, 1.70000e0, 2.00000e0, 3.00000e0,  
     5.00000e0, 7.00000e0, 1.00000e1, 1.50000e1, 2.00000e1,  
     3.00000e1, 5.00000e1, 7.00000e1, 1.00000e2, 1.30000e2,  
     1.70000e2, 2.00000e2, 2.50000e2, 3.00000e2, 5.00000e2,  
     7.00000e2, 1.00000e3, 1.50000e3, 2.00000e3, 3.00000e3,  
     4.00000e3, 5.00000e3, 6.00000e3, 7.00000e3, 8.00000e3,  
     9.00000e3, 1.00000e4};
  double xqh2p[] = 
    { 5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1, 5.00000e-1,  
      5.00000e-1, 5.00024e-1, 5.00921e-1, 5.15632e-1, 5.64397e-1,  
      7.65738e-1, 1.33880e0 , 1.91022e0 , 2.68536e0 , 3.40918e0,  
      4.35256e0 , 5.05717e0 , 6.23172e0 , 7.40726e0 , 1.21361e1,  
      1.70214e1 , 2.50338e1 , 4.11228e1 , 6.13523e1 , 1.16205e2,  
      1.93054e2 , 2.94603e2 , 4.22340e2 , 5.76489e2 , 7.56297e2,  
      9.60417e2 , 1.18728e3};

  double pf_z_H2p = linear(xTBC,xqh2p,32,T_K);

  return pf_z_H2p;
}

/*-----------------------------------------------------------------------
  ! liner interpolation
  !-----------------------------------------------------------------------  */
double linear(double *xa,double *ya,int m,double x){
  double y;
  int ms, i;
  double t, y1, y2;
      
  ms = 0;
  for(i=0;i<m;i++){
    if(x-xa[i] <= 0. || i==m-1){
      ms = i;
      break;
    }
  }

  if(ms==0) ms=1;
  y1 = ya[ms-1];
  y2 = ya[ms];

  t  = (x - xa[ms-1])/(xa[ms] - xa[ms-1]);
  y  = (1.-t)*y1 + t*y2;

  return y;  
}
