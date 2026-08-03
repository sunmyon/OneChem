#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "../allvars.h"
#include "../proto.h"
#include "chemistry_internal.h"

#ifdef CHEMISTRY_SC
#ifdef DUST_COOLING
static int flag_dust_init=0;

void get_dust_opacities(double T_dust, double *kappa_Rossland, double *kappa_Planck);
void read_dust_opacity_tables(void);

double dust_cross_section(double *S){
  double sigma_eff, S_eff;
  const double a_min = 5.e-7; // 5nm
  const double a_mid = 1.e-4; // 1um
  const double a_max = 5.e-4; // 5um
  const double rho_dust = 3.; // 3g/cc
  const int n_bin=10000;

  int i;
  double dln_a, a, f, sigma, S_i;
  double a_old, da;
  double sum1, sum2, norm;

  dln_a = log(a_max/a_min)/n_bin;
  a_old = 0.;

  for(i=0,sum1=sum2=norm=0.;i<n_bin;i++){
    a = a_min * exp(i * dln_a);
    if(a < a_mid)
      f = pow(a/a_min, -3.5);
    else
      f = pow(a_mid/a_min, -3.5) * pow(a/a_mid, -5.5);
    
    sigma = M_PI * a * a;
    S_i   = sigma / (4.*M_PI/3. * a * a * a * rho_dust);

    da    = a - a_old;

    sum1  += sigma * f * da;
    sum2  += S_i   * f * da;
    norm  += f * da;

    a_old = a;
  }
  
  sigma_eff = sum1 / norm;
  S_eff = sum2 / norm;
  //printf("sigma_eff=%g [cm^2] S_eff=%g [cm^2/g] rad=%g [um] sum=%g norm=%g\n", sigma_eff, S_eff, sqrt(sigma_eff/M_PI)*1.e5, sum1, norm);

  *S = S_eff;
  return sigma_eff;
}


double dust_temperature(double rho_cgs, double T_gas, double metallicity, double l_Jeans, double stellar_flux, struct opacities *opac, double *edot_dust){
  double tau_g = opac -> tau;

  if(flag_dust_init == 0)
    read_dust_opacity_tables();
  
  double tcoll_inv = rho_cgs / PROTONMASS * 1.1e-5 / DUST_TO_GAS_MASS_RATIO * sqrt(T_gas / 1000.) * (1. - 0.8 * exp(-75./T_gas));  

  double T_dust;
  double Td_L, Td_S;
  double kappa_P, kappa_R;
  double tau, beta_esc;
  double Lambda_coll, Lambda_dust, resid;


  Td_L = 2000.0;
  Td_S = 0.;
  T_dust = (Td_S+Td_L)/2.0;

  //printf("T_gas=%g n=%g\n", T_gas, rho_cgs/PROTONMASS);

  int itr;
  for(itr=0;itr<30;itr++){
    get_dust_opacities(T_dust, &kappa_R, &kappa_P);
    
    kappa_P /= DUST_TO_GAS_MASS_RATIO; //kappa for unit dust mass (independent of metallicity)
    kappa_R /= DUST_TO_GAS_MASS_RATIO; //kappa for unit dust mass (independent of metallicity)

    tau = tau_g + kappa_P * rho_cgs *  metallicity * DUST_TO_GAS_MASS_RATIO * l_Jeans;
    beta_esc = 1./(1 + tau*tau);

    double f_over_rho_g = 9.67e-4;
    if(T_dust < 159.)
      f_over_rho_g = 5.33e-3;
    else if(T_dust < 275.)
      f_over_rho_g = 4.04e-3;
    else if(T_dust < 425.)
      f_over_rho_g = 3.438e-3;
    else if(T_dust < 680.)
      f_over_rho_g = 1.085e-3;
    else if(T_dust < 1230.)
      f_over_rho_g = 9.67e-4;
    
    Lambda_coll = tcoll_inv * f_over_rho_g * (T_gas - T_dust); //Lambda_coll for unit dust mass
    Lambda_dust = (4 * STEFAN * (pow(T_dust,4) - pow(All.T_CMB,4)) - stellar_flux) * kappa_P * beta_esc;
    //Lambda_cust for unit dust mass    

    resid = Lambda_dust - Lambda_coll;
    
    if(resid > 0.) 
      Td_L = T_dust;
    else
      Td_S = T_dust;
    
    if(itr%2000==1999 || FlagPrint)
      printf("[%02d] T_dust=%g dT=%g T_CMB=%g Trad=%G L_coll=%g L_dust=%g kapa=%g beta_esc=%g\n"
	     ,itr , T_dust, T_gas - T_dust, All.T_CMB, pow(stellar_flux/4./STEFAN, 0.25), Lambda_coll, Lambda_dust, kappa_P, beta_esc);

    T_dust = (Td_S+Td_L)/2.0;
  }

  opac -> kappa_dust_Rossland = kappa_R * DUST_TO_GAS_MASS_RATIO * metallicity;
  opac -> kappa_dust_Planck = kappa_P * DUST_TO_GAS_MASS_RATIO * metallicity;
  
  opac -> tau = tau;
  opac -> f_reduce = 1./(1. + tau*tau);

  *edot_dust = -Lambda_dust * rho_cgs * metallicity * DUST_TO_GAS_MASS_RATIO; 
  
  return T_dust;
}

#ifdef LYMANWERNER
double dust_PEheating(double G0, double nH, double y_e, double T_gas, double metallicity){
  double edot_PE;
  double epsilon, n_e;
  double fac1, beta;

  if(G0 > 1.e5) G0 = 1.e5; //added [CHON240107]

  n_e  = nH * y_e;
  fac1 = G0 * sqrt(T_gas) / n_e;

  epsilon = 4.9e-2 / (1. + 4.e-3 * pow(fac1, 0.73))
    + 3.7e-2 * pow(T_gas/1.e4 , 0.7) / (1. + 2.e-4 * fac1);
  beta    = 0.735 * pow(T_gas, -0.068);

  edot_PE  = 1.e-24 * epsilon * G0 * nH - 4.65e-30 * pow(T_gas, 0.94) * pow(fac1, beta) * nH * n_e;
  edot_PE *= metallicity;

  if(FlagPrint)
    printf("edot_PE=%g G0=%g epsilon=%g nH=%g n_e=%g fac1=%g beta=%g heat=%g cool=%g\n"
	   , edot_PE, G0, epsilon, nH, n_e, fac1, beta
	   , 1.e-24 * epsilon * G0 * nH * metallicity
	   , 4.65e-30 * pow(T_gas, 0.94) * pow(fac1, beta) * nH * n_e * metallicity);

  return edot_PE;
}
#endif

#ifdef H2_FORMATION_ON_DUST
// H2 formation rate on grains (Schneider et al. 2006)
double H2form_on_grains(double rho_cgs, double T_gas, double T_dust, double metallicity)
{
  double T2       = 0.01*T_gas;
  double S_H      = 1.0/( 1.0 + 0.4*sqrt(T2+0.01*T_dust) + 0.2*T2 + 0.08*T2*T2);  
  double k_H2form = 6.0e-17 * sqrt(T_gas/300.) /(1. + exp(7.5e2*(1./75.-1./T_dust))) * S_H * metallicity * rho_cgs / PROTONMASS;

  //printf("k_H2form=%g epsilon=%g S_H=%g rho_d=%g sigma=%g vH_avr=%g\n", k_H2form, epsilon_H2, S_H, rho_dust, sigma_dust, vH_avr);

  return k_H2form;
}
#endif


char fname_dust_Rossland[]="kR.out";
char fname_dust_Planck[]="kP.out";

int n_Tdust;
double *Kappa_Tdust;
double *Kappa_Dust_Rossland;
double *Kappa_Dust_Planck;
double TMax_dust_tau, TMin_dust_tau;

void read_dust_opacity_tables(void){
  if(ThisTask == 0){
    FILE *fp_r, *fp_p;

    if(!(fp_r = open_chemistry_table(fname_dust_Rossland))){
      exit(1);
    }
    
    if(!(fp_p = open_chemistry_table(fname_dust_Planck))){
      exit(1);
    }

    char buf[200];

    int count_r, count_p;
    count_r = count_p = 0;

    int n;
    double T_dust, kappa;
    while(fgets(buf, sizeof(buf), fp_r) != NULL){
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2)
	count_r++;      
    }

    while(fgets(buf, sizeof(buf), fp_p) != NULL){
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2)
	count_p++;      
    }

    if(count_r != count_p){
      printf("The number of raws of the opacity tables %s and %s is different.\n"
	     "Plank=%d Rossland=%d\n"
	     , fname_dust_Rossland, fname_dust_Planck, count_p, count_r);

      exit(1);
    }

    rewind(fp_r), rewind(fp_p);
    n_Tdust   = count_r;

    if(!(Kappa_Tdust = malloc(sizeof(double)*n_Tdust))){
      printf("Can't allocate memory.\n");
      exit(1);
    }

    if(!(Kappa_Dust_Rossland = malloc(sizeof(double)*n_Tdust))){
      printf("Can't allocate memory.\n");
      exit(1);
    }

    if(!(Kappa_Dust_Planck = malloc(sizeof(double)*n_Tdust))){
      printf("Can't allocate memory.\n");
      exit(1);
    }

    count_r = count_p = 0;
    TMax_dust_tau = -1.e30;
    TMin_dust_tau = 1.e30;
        
    while(fgets(buf, sizeof(buf), fp_r) != NULL){    
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2){
	Kappa_Tdust[count_r] = T_dust;
	Kappa_Dust_Rossland[count_r] = kappa;	
	count_r++;

	//printf("[Ross]T_dust=%g kappa=%g\n", T_dust, kappa);
      }

      if(T_dust > TMax_dust_tau) TMax_dust_tau = T_dust;
      if(T_dust < TMin_dust_tau) TMin_dust_tau = T_dust;
    }

    while(fgets(buf, sizeof(buf), fp_p) != NULL){    
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2){
	if(fabs(Kappa_Tdust[count_p]/T_dust-1.) > 1.e-3){
	  printf("Something wrong in the dust opacity table.\n"
		 "Temperature value is different Tdust[%d]=%g while %g.\n"
		 ,count_p, Kappa_Tdust[count_p], T_dust);

	  exit(1);
	}

	Kappa_Dust_Planck[count_p] = kappa;	
	count_p++;

	//printf("[Planck]T_dust=%g kappa=%g\n", T_dust, kappa);
      }
    }

    fclose(fp_p),fclose(fp_r);
  }

#ifdef HAVE_MPI
  if(NTask > 1){
    MPI_Bcast(&n_Tdust, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(ThisTask != 0){
      if(!(Kappa_Tdust         = malloc(sizeof(double)*n_Tdust)) ||
	 !(Kappa_Dust_Rossland = malloc(sizeof(double)*n_Tdust)) ||
	 !(Kappa_Dust_Planck   = malloc(sizeof(double)*n_Tdust))
	 ){
	printf("Can't allocate memory.\n");
	exit(1);
      }
    }

    MPI_Bcast(Kappa_Tdust        , n_Tdust, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Kappa_Dust_Rossland, n_Tdust, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Kappa_Dust_Planck  , n_Tdust, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif
  
  flag_dust_init = 1;
}

void get_dust_opacities(double T_dust, double *kappa_Rossland, double *kappa_Planck){
  int i_T;
  double T_new;

  i_T=0, T_new = Kappa_Tdust[0];
  while(T_dust > T_new){
    i_T++;
    T_new = Kappa_Tdust[i_T];

    if(i_T == n_Tdust-1) break;
  }

  int i_T_left;

  i_T_left   = i_T - 1;
  if(i_T_left < 0) i_T_left = 0;
  if(i_T == n_Tdust - 1) i_T_left = n_Tdust-1;

  double c_T;
  double kappa_r, kappa_p;

  c_T   = (Kappa_Tdust[i_T]-T_dust)    / (Kappa_Tdust[i_T] -Kappa_Tdust[i_T_left]);

  if(i_T == i_T_left)
    c_T = 0.;

  kappa_r = 
    c_T       * Kappa_Dust_Rossland[i_T_left] + 
    (1 - c_T) * Kappa_Dust_Rossland[i_T];

  kappa_p = 
    c_T       * Kappa_Dust_Planck[i_T_left] + 
    (1 - c_T) * Kappa_Dust_Planck[i_T];

  *kappa_Rossland = kappa_r;
  *kappa_Planck   = kappa_p;

  /*printf("[T=%g] kappa_r=%g k=%g %g\n"
	 , T_dust, kappa_r
	 , Kappa_Dust_Rossland[i_T_left], Kappa_Dust_Rossland[i_T]
	 );*/
}

#endif
#endif
