#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"
#include "chemistry_tables.h"
#include "saha.h"
#include "chemistry_models.h"
#include "partition.h"

int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


/******************************************/
/* ====================================== */
/*   Saha's value of atomic hydrogen      */
/*       H, H2, Hp, e                     */
/*   (RM wrote: 2018/12/10)               */
/******************************************/
double chem_saha(double xnH,double T_K,double *y_H,double *y_H2,double *y_e,double *y_Hp){
  double x_u,x_l,x_m,f_u,f_l,f_m;
  double g_ion,g_dis;
  double Z_pf,Mass,Ene,T_erg;
  double z_H,z_H2,m_H2,m_Hp,m_Hm;
  double err;
  double err_cnv = 1.e-5;
  const int itrmax = 1000;
  double m_H, m_e=9.1094e-28;
  double z_Hp=1.0,z_e=2.0;
  double h_Pla=6.625e-27,xk_B=1.381e-16, chi_H=2.179e-11,chi_H2=7.177e-12;
  int itr;

  z_H  = partition_H(T_K);
  z_H2 = chemistry_models.h2_partition(T_K);
  m_H = m_Hp = m_Hm = PROTONMASS;
  m_H2 = 2 * m_H;
  T_erg = xk_B*T_K;
  Ene   = 2.0*M_PI*T_erg/(h_Pla*h_Pla);

  /* set the function(g) depending on only temperature from Saha-Boltzmann equation*/
  /*-  H   ->  H+  +   e  */
  Z_pf   = z_Hp*z_e/z_H;
  Mass   = m_Hp*m_e/m_H;
  g_ion  = Z_pf * pow(Mass,1.5) * pow(Ene,1.5) * exp(-chi_H/T_erg);
  /*-  H2  ->  H   +   H  */
  Z_pf   = z_H*z_H/z_H2;
  Mass   = m_H*m_H/m_H2;
  g_dis  = Z_pf * pow(Mass,1.5) * pow(Ene,1.5) * exp(-chi_H2/T_erg);


  /** ------------------------------------------------- **/
  /** ---------          bisection          ----------- **/
  /** ------------------------------------------------- **/
  /** - initial setup   (u:upper limit, l:lower limit)- **/
  x_u = 1.;
  x_l = 0.;
  f_u = 2.*(x_u*x_u)*xnH/g_dis + x_u + sqrt(x_u*g_ion/xnH) - 1.;
  f_l = 2.*(x_l*x_l)*xnH/g_dis + x_l + sqrt(x_l*g_ion/xnH) - 1.;

  if(f_l*f_u >= 0.){
    printf("not found Saha value (chem_saha)\nstopping...\n");
    return NAN;
  }

  for(itr=0;itr<itrmax;itr++){
    x_m = 0.5*(x_u+x_l);

    err = fabs((x_u-x_m)/x_m);
    if(itr <= (itrmax/2))
      err_cnv = 1.e-5;
    else
      err_cnv = 1.e-3;

    if(err <= err_cnv) break;

    f_m = 2.*(x_m*x_m)*xnH/g_dis + x_m + sqrt(x_m*g_ion/xnH) - 1.;

    if(f_u*f_m <= 0.){
      x_l = x_m;
      f_l = f_m;
    }else{
      x_u = x_m;
      f_u = f_m;
    }

    /**** abundance become >1.0 or negative  ***/
    if(x_u > 1. || x_l < 0.){
      printf("not found Saha value (chem_saha)\nstopping...\n");
      return NAN;
    }
    /**** NaN/Inf check ****/
    if(isnan(f_m) || isinf(x_m)){
      printf("not found Saha value (chem_saha)\n");
      printf("f_m or x_m  NaN/Inf\nstopping...\n");
      return NAN;
    }

    if(itr > itrmax/2){
      printf("x_u=%g x_l=%g f_l=%g f_u=%g err=%g\n", x_u, x_l, f_l, f_u, err);
    }

  }

  *y_H  = x_m;
  *y_H2 = (x_m*x_m)*xnH/g_dis;
  *y_Hp = sqrt(x_m*g_ion/xnH);
  *y_e  = (*y_Hp);

  return err;
}


/* =======================================================*/
/*  calculate chemistry using saha equation( H & D )      */
/*     (  H, H2, e, H+, H- )                              */
/*       ske MeMo : I made this bymyself                  */
/*                  this will be used @ n > n_saha        */
/* =======================================================*/

double chem_saha_5sp(const ChemistryRateTables *tables, double xnH,
                     double T_K, double *y){
  double z_H,z_H2,z_e,z_Hp, z_Hm;
  double m_H,m_H2,m_e,m_Hp, m_Hm;

  int n_sp;
  enum SPECIES species;
  n_sp = sEND - sELEC;

  double Z_pf[n_sp], Mass[n_sp], g[n_sp], dE[n_sp];
  double f[n_sp];
  double y_init[n_sp], dy[n_sp], dy_norm[n_sp];
  double AA[n_sp*n_sp];
  double y_H, y_H2,  y_e, y_Hp, y_Hm;
  double Ene, T_erg, err, err_max;

  int i, j;
  int itr, info;
  const int itr_max = 100;
  int inc = 1;
  int ipiv[n_sp];
  const double *chem_e = tables->formation_energy_kj_mol;


  for(i=0;i<n_sp;i++) y_init[i] = y[i];

  //printf("y_init=%g %g %g %g %g\n", y[0],y[1],y[2],y[3],y[4]);

  T_erg = BOLTZMANN * T_K;
  Ene   = 2.*M_PI*T_erg/PLANCK/PLANCK;

  z_H   = partition_H(T_K);
  z_H2  = chemistry_models.h2_partition(T_K);
  z_e   = 2.;
  z_Hp  = 1.;
  z_Hm  = 1.;

  m_H = m_Hp = m_Hm = PROTONMASS;
  m_H2  = 2.*m_H;
  m_e   = ELECTRONMASS;

  /*---------- saha equation--------------------------*/

  /* 0) H + e <-> H+ + 2e */
  dE[0]   = (chem_e[sHII]+2*chem_e[sELEC]-chem_e[sHI]-chem_e[sELEC])*kJ_per_mol_to_eV *eV_to_erg;
  Z_pf[0] = z_Hp*z_e/z_H;
  Mass[0] = m_Hp*m_e/m_H;
  g[0]    = Z_pf[0] * pow(Mass[0],1.5) * pow(Ene,1.5) * exp(-dE[0]/T_erg);

  /* 1) H + H2 <-> H + H + H */
  dE[1]   = (3*chem_e[sHI]-chem_e[sHI]-chem_e[sH2I])*kJ_per_mol_to_eV *eV_to_erg;
  Z_pf[1] = z_H*z_H/z_H2;
  Mass[1] = m_H*m_H/m_H2;
  g[1]    = Z_pf[1] * pow(Mass[1],1.5) * pow(Ene,1.5) * exp(-dE[1]/T_erg);

  /* 2) H- + e <-> H + 2e */
  dE[2]   = (2*chem_e[sELEC] + chem_e[sHI] - chem_e[sHM] - chem_e[sELEC])*kJ_per_mol_to_eV *eV_to_erg;
  Z_pf[2] = z_H*z_e/z_Hm;
  Mass[2] = m_H*m_e/m_Hm;
  g[2]    = Z_pf[2] * pow(Mass[2],1.5) * pow(Ene,1.5) * exp(-dE[2]/T_erg);

  //printf("T_gas=%g g=%g %g %g\n", T_K, g[0], g[1], g[2]);

  /*---------- newton method----------------------------------------- */
  for(itr=0;itr<itr_max;itr++){
    // chemical abundance
    y_H = y[sHI] ; y_H2 = y[sH2I] ; y_e  = y[sELEC] ; y_Hp = y[sHII] ; y_Hm = y[sHM];

    // define function of simultaneous equation (f = 0)
    /* 0) H + e <-> H+ + 2e */
    f[0] = g[0]/xnH*y_H - y_Hp*y_e;
    /* 1) H + H2 <-> H + H + H */
    f[1] = g[1]/xnH*y_H2 - y_H*y_H;
    /* 2) H- + e <-> H + 2e */
    f[2] = g[2]/xnH*y_Hm - y_H*y_e;
    /* H abundance */
    f[3] = 1. - (y_H + y_Hp + y_Hm + 2.*y_H2);
    /* charge conservation */
    f[4] = y_Hp - y_e - y_Hm;

    /* calculation of Jacobian A */
    for(i=0;i<n_sp;i++)
      for(j=0;j<n_sp;j++)
	AA[i*n_sp+j] = 0.;

    /*AA[0*n_sp+sHI] = g[0]/xnH  ; AA[0*n_sp+sELEC] = -y_Hp   ; AA[0*n_sp+sHII] = -y_e;
    AA[1*n_sp+sHI] = -2.*y_H   ; AA[1*n_sp+sH2I]  = g[1]/xnH;
    AA[2*n_sp+sHI] = -y_e  ; AA[2*n_sp+sELEC] = -y_H    ; AA[2*n_sp+sHM]  = g[2]/xnH;

    AA[3*n_sp+sHI] = -1. ; AA[3*n_sp+sH2I] = -2. ; AA[3*n_sp+sHII]  = -1. ; AA[3*n_sp+sHM] = -1.;
    AA[4*n_sp+sHII] = 1. ; AA[4*n_sp+sHM]  = -1. ; AA[4*n_sp+sELEC] = -1. ;*/

    AA[0+n_sp*sHI] = g[0]/xnH        ; AA[0+n_sp*sELEC] = -y_Hp ; AA[0+n_sp*sHII] = -y_e;
    AA[1+n_sp*sHI] = -2.*y_H ; AA[1+n_sp*sH2I]  = g[1]/xnH      ;
    AA[2+n_sp*sHI] = -y_e    ; AA[2+n_sp*sELEC] = -y_H  ; AA[2+n_sp*sHM]  = g[2]/xnH;

    AA[3+n_sp*sHI] = -1. ; AA[3+n_sp*sH2I] = -2. ; AA[3+n_sp*sHII]  = -1. ; AA[3+n_sp*sHM] = -1.;
    AA[4+n_sp*sHII] = 1. ; AA[4+n_sp*sHM]  = -1. ; AA[4+n_sp*sELEC] = -1. ;

    /* solve A*dy = -f with gauss elimination */
    /* use LAPACK solver                      */
    for(i=0;i<n_sp;i++)  dy[i] = -f[i];

    dgesv_(&n_sp, &inc, AA, &n_sp, ipiv, dy, &n_sp, &info);


    /* update y */
    for(i=0;i<n_sp;i++)  y[i] = dy[i] + y[i];

    err_max = 0.;
    for(i=0;i<n_sp;i++){
      if(y[i] > 1.e-15 || y_init[i] > 1.e-15){
	err = fabs(dy[i]/y_init[i]);
	err_max = dmax(err,err_max);
      }

      dy_norm[i] = fabs((y_init[i]-y[i])/y_init[i]);
    }


    if(err_max < 1.e-10)
      break;
  }

  /* ----------------   iterationが収束したかの判定   ------------------- */

  if(itr == itr_max){
    printf("didn't converged...\n");
  }

  return err_max;
}






int heap_sort(unsigned long n, double ra[])
{
  int ifunc=1;
  unsigned long i,ir,j,l;
  double rra;

  if (n < 2) return n;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
  }

  return (ifunc);
}
