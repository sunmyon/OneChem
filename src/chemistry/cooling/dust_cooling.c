#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chemistry_constants.h"
#include "chemistry_diagnostics.h"
#include "chemistry_tables.h"
#include "dust_cooling.h"
#include "chemistry_profile.h"

#ifdef CHEMISTRY_SC
#ifdef DUST_COOLING
struct ChemistryDustTables {
  int temperature_count;
  double *temperature;
  double *rosseland;
  double *planck;
  double log_temperature_min;
  double inverse_log_temperature_step;
};

void get_dust_opacities(const ChemistryDustTables *tables, double T_dust,
                        double *kappa_Rossland, double *kappa_Planck);
static void get_dust_opacities_with_derivative(
    const ChemistryDustTables *tables, double T_dust,
    double *kappa_Rossland, double *kappa_Planck,
    double *kappa_Planck_derivative);

typedef struct {
  double kappa_P;
  double kappa_R;
  double tau;
  double beta_esc;
  double Lambda_coll;
  double Lambda_dust;
  double residual_derivative;
} DustBalance;

static double dust_collision_factor(double T_dust, double *derivative)
{
  static const double boundaries[] = {159., 275., 425., 680.};
  static const double levels[] = {5.33e-3, 4.04e-3, 3.438e-3, 1.085e-3, 9.67e-4};
  const double relative_width = 1.e-3;

  double factor = levels[0];
  *derivative = 0.;
  for(int i = 0; i < 4; i++){
    double width = relative_width * boundaries[i];
    double x = (T_dust - boundaries[i]) / width;
    double transition;
    double transition_derivative;
    if(x <= -20.){
      transition = 0.;
      transition_derivative = 0.;
    } else if(x >= 20.){
      transition = 1.;
      transition_derivative = 0.;
    } else {
      double tanh_x = tanh(x);
      transition = 0.5 * (1. + tanh_x);
      transition_derivative = 0.5 * (1. - tanh_x * tanh_x) / width;
    }
    double level_change = levels[i + 1] - levels[i];
    factor += level_change * transition;
    *derivative += level_change * transition_derivative;
  }

  return factor;
}

static double evaluate_dust_balance(const ChemistryDustTables *tables,
                                    double T_dust, double T_gas,
                                    double radiation_temperature,
                                    double rho_cgs, double metallicity,
                                    double l_Jeans, double stellar_flux,
                                    double tau_g, double tcoll_inv,
                                    DustBalance *balance)
{
  double kappa_P_derivative;
  get_dust_opacities_with_derivative(tables, T_dust, &balance->kappa_R,
                                     &balance->kappa_P,
                                     &kappa_P_derivative);

  balance->kappa_P /= DUST_TO_GAS_MASS_RATIO;
  balance->kappa_R /= DUST_TO_GAS_MASS_RATIO;
  kappa_P_derivative /= DUST_TO_GAS_MASS_RATIO;

  balance->tau = tau_g + balance->kappa_P * rho_cgs * metallicity * DUST_TO_GAS_MASS_RATIO * l_Jeans;
  balance->beta_esc = 1. / (1. + balance->tau * balance->tau);
  double tau_derivative = kappa_P_derivative * rho_cgs * metallicity * DUST_TO_GAS_MASS_RATIO * l_Jeans;
  double beta_derivative = -2. * balance->tau * tau_derivative * balance->beta_esc * balance->beta_esc;

  double f_over_rho_g_derivative;
  double f_over_rho_g = dust_collision_factor(T_dust, &f_over_rho_g_derivative);

  balance->Lambda_coll = tcoll_inv * f_over_rho_g * (T_gas - T_dust);

  double T_radiation = pow(pow(radiation_temperature, 4) + stellar_flux / (4. * STEFAN), 0.25);
  double T4_minus_Tradiation4 = (T_dust - T_radiation) * (T_dust + T_radiation)
    * (T_dust * T_dust + T_radiation * T_radiation);
  double radiative_flux = 4. * STEFAN * T4_minus_Tradiation4;
  balance->Lambda_dust = radiative_flux * balance->kappa_P * balance->beta_esc;
  balance->residual_derivative = 16. * STEFAN * T_dust * T_dust * T_dust * balance->kappa_P * balance->beta_esc
    + radiative_flux * kappa_P_derivative * balance->beta_esc
    + radiative_flux * balance->kappa_P * beta_derivative
    + tcoll_inv * f_over_rho_g
    - tcoll_inv * f_over_rho_g_derivative * (T_gas - T_dust);

  return balance->Lambda_dust - balance->Lambda_coll;
}
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



double dust_temperature(const ChemistryDustTables *tables, double rho_cgs,
                        double T_gas, double radiation_temperature,
                        double metallicity, double l_Jeans,
                        double stellar_flux, struct opacities *opac,
                        double *edot_dust, int diagnostics_enabled){
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  double tau_g = opac -> tau;

  double tcoll_inv = rho_cgs / PROTONMASS * 1.1e-5 / DUST_TO_GAS_MASS_RATIO * sqrt(T_gas / 1000.) * (1. - 0.8 * exp(-75./T_gas));  

  double T_radiation = pow(pow(radiation_temperature, 4) + stellar_flux / (4. * STEFAN), 0.25);
  double Td_S = fmin(2000., fmax(0., fmin(T_gas, T_radiation)));
  double Td_L = fmin(2000., fmax(T_gas, T_radiation));
  double T_dust = fmin(Td_L, fmax(Td_S, T_radiation));
  double previous_step = Td_L - Td_S;
  DustBalance balance;

  //printf("T_gas=%g n=%g\n", T_gas, rho_cgs/PROTONMASS);

  int itr;
  for(itr=0;itr<30;itr++){
    double resid = evaluate_dust_balance(tables, T_dust, T_gas,
        radiation_temperature, rho_cgs, metallicity, l_Jeans, stellar_flux,
        tau_g, tcoll_inv, &balance);
    double residual_scale = fabs(balance.Lambda_dust) + fabs(balance.Lambda_coll);
    if(fabs(resid) <= 1.e-5 * residual_scale)
      break;
    
    if(resid > 0.){
      Td_L = T_dust;
    } else {
      Td_S = T_dust;
    }

    double temperature_tolerance = fmax(1.e-4, 1.e-5 * fmax(T_dust, 1.));
    if(Td_L - Td_S <= temperature_tolerance){
      break;
    }
    
    if(itr%2000==1999 || diagnostics_enabled)
      printf("[%02d] T_dust=%g dT=%g T_CMB=%g Trad=%G L_coll=%g L_dust=%g kapa=%g beta_esc=%g\n"
	     ,itr , T_dust, T_gas - T_dust, radiation_temperature, pow(stellar_flux/4./STEFAN, 0.25), balance.Lambda_coll, balance.Lambda_dust, balance.kappa_P, balance.beta_esc);

    double T_next = 0.5 * (Td_S + Td_L);
    int use_newton = 0;

    if(itr >= 0){
      double derivative = balance.residual_derivative;
      if(isfinite(derivative) && derivative != 0.){
        double T_newton = T_dust - resid / derivative;
        double newton_step = fabs(T_newton - T_dust);
        if(isfinite(T_newton) && T_newton > Td_S && T_newton < Td_L && newton_step < 0.5 * previous_step){
          T_next = T_newton;
          previous_step = newton_step;
          use_newton = 1;
        }
      }
    }

    if(!use_newton)
      previous_step = 0.5 * (Td_L - Td_S);

    T_dust = T_next;
  }

  evaluate_dust_balance(tables, T_dust, T_gas, radiation_temperature,
      rho_cgs, metallicity, l_Jeans, stellar_flux, tau_g, tcoll_inv,
      &balance);

  opac -> kappa_dust_Rossland = balance.kappa_R * DUST_TO_GAS_MASS_RATIO * metallicity;
  opac -> kappa_dust_Planck = balance.kappa_P * DUST_TO_GAS_MASS_RATIO * metallicity;
  
  opac -> tau = balance.tau;
  opac -> f_reduce = balance.beta_esc;

  double delta_coll = fabs(T_gas - T_dust);
  double delta_rad = fabs(T_dust - T_radiation);
  double Lambda_gas_dust = delta_coll >= delta_rad ? balance.Lambda_coll : balance.Lambda_dust;

  *edot_dust = -Lambda_gas_dust * rho_cgs * metallicity * DUST_TO_GAS_MASS_RATIO;
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_DUST_TEMPERATURE, profile_timer);
#endif
  return T_dust;
}

#ifdef LYMANWERNER
double dust_PEheating(double G0, double nH, double y_e, double T_gas,
                      double metallicity, int diagnostics_enabled){
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

  if(diagnostics_enabled)
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


static const char fname_dust_Rossland[]="kR.out";
static const char fname_dust_Planck[]="kP.out";

int dust_opacity_tables_create(const char *table_root,
                               ChemistryDustTables **output_tables){
  if(output_tables == NULL)
    return 0;
  *output_tables = NULL;

  ChemistryDustTables *tables = calloc(1, sizeof(*tables));
  if(tables == NULL){
    printf("Can't allocate memory.\n");
    return 0;
  }

  FILE *fp_r = NULL;
  FILE *fp_p = NULL;

  {
    if(!(fp_r = open_chemistry_table(table_root, fname_dust_Rossland))){
      goto fail;
    }
    
    if(!(fp_p = open_chemistry_table(table_root, fname_dust_Planck))){
      goto fail;
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

    if(count_r != count_p || count_r < 2){
      printf("The number of raws of the opacity tables %s and %s is different.\n"
	     "Plank=%d Rossland=%d\n"
	     , fname_dust_Rossland, fname_dust_Planck, count_p, count_r);

      goto fail;
    }

    rewind(fp_r), rewind(fp_p);
    tables->temperature_count = count_r;

    if(!(tables->temperature = malloc(sizeof(double)*tables->temperature_count))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    if(!(tables->rosseland = malloc(sizeof(double)*tables->temperature_count))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    if(!(tables->planck = malloc(sizeof(double)*tables->temperature_count))){
      printf("Can't allocate memory.\n");
      goto fail;
    }

    count_r = count_p = 0;
        
    while(fgets(buf, sizeof(buf), fp_r) != NULL){    
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2){
	tables->temperature[count_r] = T_dust;
	tables->rosseland[count_r] = kappa;	
	count_r++;

	//printf("[Ross]T_dust=%g kappa=%g\n", T_dust, kappa);
      }

    }

    while(fgets(buf, sizeof(buf), fp_p) != NULL){    
      n = sscanf(buf, "%le %le", &T_dust, &kappa);
      if(n == 2){
	if(fabs(tables->temperature[count_p]/T_dust-1.) > 1.e-3){
	  printf("Something wrong in the dust opacity table.\n"
		 "Temperature value is different Tdust[%d]=%g while %g.\n"
		 ,count_p, tables->temperature[count_p], T_dust);

	  goto fail;
	}

	tables->planck[count_p] = kappa;	
	count_p++;

	//printf("[Planck]T_dust=%g kappa=%g\n", T_dust, kappa);
      }
    }

    fclose(fp_p),fclose(fp_r);
    fp_p = fp_r = NULL;
  }

  tables->log_temperature_min = log(tables->temperature[0]);
  tables->inverse_log_temperature_step =
      (tables->temperature_count - 1) /
      (log(tables->temperature[tables->temperature_count - 1]) -
       tables->log_temperature_min);
  *output_tables = tables;
  return 1;

fail:
  if(fp_p != NULL)
    fclose(fp_p);
  if(fp_r != NULL)
    fclose(fp_r);
  free(tables->temperature);
  free(tables->rosseland);
  free(tables->planck);
  free(tables);
  return 0;
}

void dust_opacity_tables_destroy(ChemistryDustTables *tables)
{
  if(tables == NULL)
    return;
  free(tables->temperature);
  free(tables->rosseland);
  free(tables->planck);
  free(tables);
}

static void get_dust_opacities_with_derivative(
    const ChemistryDustTables *tables, double T_dust,
    double *kappa_Rossland, double *kappa_Planck,
    double *kappa_Planck_derivative){
  const int n_Tdust = tables->temperature_count;
  const double *Kappa_Tdust = tables->temperature;
  const double *Kappa_Dust_Rossland = tables->rosseland;
  const double *Kappa_Dust_Planck = tables->planck;
  const double dust_log_T_min = tables->log_temperature_min;
  const double dust_inv_dlog_T = tables->inverse_log_temperature_step;
  int i_T, i_T_left;

  if(T_dust <= Kappa_Tdust[0]) {
    i_T = i_T_left = 0;
  } else if(T_dust >= Kappa_Tdust[n_Tdust-1]) {
    i_T = i_T_left = n_Tdust-1;
  } else {
    i_T_left = (int)((log(T_dust) - dust_log_T_min) * dust_inv_dlog_T);
    if(i_T_left < 0) i_T_left = 0;
    if(i_T_left > n_Tdust - 2) i_T_left = n_Tdust - 2;

    while(i_T_left > 0 && T_dust <= Kappa_Tdust[i_T_left])
      i_T_left--;
    while(i_T_left < n_Tdust - 2 && T_dust > Kappa_Tdust[i_T_left + 1])
      i_T_left++;

    i_T = i_T_left + 1;
    if(i_T == n_Tdust-1){
      i_T_left = i_T;
    }
  }

  double c_T = 0.;
  double kappa_r, kappa_p;

  if(i_T != i_T_left)
    c_T = (Kappa_Tdust[i_T]-T_dust) / (Kappa_Tdust[i_T] - Kappa_Tdust[i_T_left]);

  kappa_r = 
    c_T       * Kappa_Dust_Rossland[i_T_left] + 
    (1 - c_T) * Kappa_Dust_Rossland[i_T];

  kappa_p = 
    c_T       * Kappa_Dust_Planck[i_T_left] + 
    (1 - c_T) * Kappa_Dust_Planck[i_T];

  *kappa_Rossland = kappa_r;
  *kappa_Planck   = kappa_p;
  *kappa_Planck_derivative = 0.;
  if(i_T != i_T_left)
    *kappa_Planck_derivative = (Kappa_Dust_Planck[i_T] - Kappa_Dust_Planck[i_T_left])
      / (Kappa_Tdust[i_T] - Kappa_Tdust[i_T_left]);

  /*printf("[T=%g] kappa_r=%g k=%g %g\n"
	 , T_dust, kappa_r
	 , Kappa_Dust_Rossland[i_T_left], Kappa_Dust_Rossland[i_T]
	 );*/
}

void get_dust_opacities(const ChemistryDustTables *tables, double T_dust,
                        double *kappa_Rossland, double *kappa_Planck){
  double derivative;
  get_dust_opacities_with_derivative(tables, T_dust, kappa_Rossland,
                                     kappa_Planck, &derivative);
}

#endif
#endif
