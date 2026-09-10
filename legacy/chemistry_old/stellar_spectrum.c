#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <gsl/gsl_math.h>
#include "../../src/standalone/allvars.h"
#include "chemistry_old_api.h"

//#ifdef LW_RADIATION_LOCAL
#define eV_to_erg  1.60184e-12

double Hminus_detachment_John1988(double hv);
double calc_detachment_rate(double Trad);
double calc_G0_rate(double Trad);

static const int nbin_Teff = 200;
static const double Tmin = 3000, Tmax = 1.e6;
static double dln_T;

static double *Table_Teff;
static double *Table_Hminus_pd_rate;
static double *Table_G0;

void make_rate_tables_LW_radiation(void){
  Table_Hminus_pd_rate = malloc(sizeof(double) * nbin_Teff);
  Table_G0 = malloc(sizeof(double) * nbin_Teff);
  Table_Teff = malloc(sizeof(double) * nbin_Teff);

  dln_T = log(Tmax/Tmin)/(double) nbin_Teff;
  
  if(ThisTask == 0){
    char fname[]= "./Hminus_detachment_rate.dat";
    FILE *fp;
    if(!(fp = fopen(fname, "w"))){
      printf("Can't open th file, %s!\n", fname);
      exit(1);
    }
  
    int i;
    for(i=0;i<nbin_Teff;i++){
      double Teff = Tmin * exp(dln_T * i);
      Table_Teff[i] = Teff;
      
      double rate = calc_detachment_rate(Teff);
      Table_Hminus_pd_rate[i] = rate;

      double g0 = calc_G0_rate(Teff);
      Table_G0[i] = g0;    
        
      fprintf(fp, "%g %g %g\n", Teff, rate, g0);
    }

    fclose(fp);
  }

#ifdef HAVE_MPI
  if(NTask > 1){
    MPI_Bcast(Table_Teff           , nbin_Teff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Table_Hminus_pd_rate , nbin_Teff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Table_G0             , nbin_Teff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif
}

void get_Hminus_PD_rate(double Teff, double *Hminus_pd_rate, double *G0){
  double x_in = dmin(dmax(Teff,Table_Teff[0]), Table_Teff[nbin_Teff-1]);
  int ix = (int)(log(x_in/Tmin) / dln_T);

  //printf("ix=%d Teff=%g x_in=%g nbin_Teff=%d\n", ix, Teff, x_in, nbin_Teff);
  
  if((ix < nbin_Teff-1 && x_in > Table_Teff[ix+1]) || x_in < Table_Teff[ix]){
    printf("something wrong in reading tables...\n");
  }
  
  if(ix < 0)
    ix = 0;
  else if(ix > nbin_Teff - 2)
    ix = nbin_Teff - 2;
   
  double x1 = Table_Teff[ix];
  double x2 = Table_Teff[ix+1];

  double tap1 = (x2-x_in)/(x2-x1);
  double tap2 = (x_in-x1)/(x2-x1);

  double rate = tap1*Table_Hminus_pd_rate[ix] + tap2*Table_Hminus_pd_rate[ix+1];
  double g0   = tap1*Table_G0[ix]             + tap2*Table_G0[ix+1];

  *Hminus_pd_rate = rate;
  *G0 = g0;
}


double calc_G0_rate(double Trad){
  const double nu_min=5.12, nu_max= 13.6;
  const int n_bin=500;

  double hc = PLANCK * C_LIGHT;
  double dln_nu = log(nu_max/nu_min) / (double)n_bin;
  double nu_old = nu_min, J_nu_old = 0.;

  int i;
  double sum=0.;
  for(i=0;i<n_bin;i++){
    double nu   = nu_min * exp(i*dln_nu);
    double d_nu = (nu - nu_old) * eV_to_Hz;
    double J_nu = 2 * pow(nu * eV_to_erg, 3) / (hc * hc) / (exp(nu*eV_to_erg/BOLTZMANN/Trad) - 1);
    
    if(i > 0)
      sum += 0.5 * (J_nu_old + J_nu) * d_nu;
    
    nu_old = nu;
    J_nu_old = J_nu;
  }

  double L_per_area = 4*M_PI*sum;

  return L_per_area; //unit in erg/s/cm^2
}

double calc_detachment_rate(double Trad){
  const double nu_min=0.76, nu_max= 13.6;
  const int n_bin=500;

  double hc = PLANCK * C_LIGHT;
  double dln_nu = log(nu_max/nu_min) / (double)n_bin;
  double nu_old = nu_min, rate_old = 0.;

  int i;
  double sum=0.;
  for(i=0;i<n_bin;i++){
    double nu   = nu_min * exp(i*dln_nu);
    double d_nu = (nu - nu_old) * eV_to_Hz;
    double J_nu = 2 * pow(nu * eV_to_erg, 3) / (hc * hc) / (exp(nu*eV_to_erg/BOLTZMANN/Trad) - 1);
    double sigma = Hminus_detachment_John1988(nu);
    double rate = 4 * M_PI / (nu * eV_to_erg) * J_nu * sigma;;
    
    if(i > 0)      
      sum += 0.5 * (rate_old + rate) * d_nu;    
    
    nu_old = nu;
    rate_old = rate;
  }

  double L_per_area = sum;

  return L_per_area; //unit in erg/s/cm^2                                                                                                                                                                                                                     
}

const double C_n[6] = {152.519, 49.534, -118.858, 92.536, -34.194, 4.982};
const double lambda_0 = 1.6419;

double Hminus_detachment_John1988(double hv){
  double lambda = PLANCK * C_LIGHT / (hv * eV_to_erg) * 1.e4;

  if(lambda > lambda_0) return 0.;

  int i;
  double f_lambda=0.;
  for(i=0;i<6;i++)
    f_lambda += C_n[i] * pow(1./lambda - 1./lambda_0, 0.5*i);

  return 1.e-18 * pow(lambda, 3) * pow(1./lambda - 1./lambda_0, 1.5) * f_lambda;
}
//#endif
