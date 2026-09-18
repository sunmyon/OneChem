#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<limits.h>
#include<math.h>
#include<ctype.h>
#include<time.h>
#include<sys/time.h>
#include<sys/stat.h>

#include "allvars.h"
#include "radiation/radiation_rates.h"
#include "standalone_host.h"
#include "standalone_constants.h"
#include "standalone_math.h"
#ifdef CHEMISTRY_OLD
#include "legacy/chemistry_old/chemistry_old_api.h"
#elif defined(CHEMISTRY_PORTABLE)
#include "chemistry/thermodynamics/thermodynamics.h"
#else
#include "legacy/chemistry_new/chemistry_api.h"
#include "legacy/chemistry_new/thermodynamics/thermodynamics.h"
#endif
#ifndef CHEMISTRY_OLD
#include "chemistry_profile.h"
#include "species.h"
#endif
#ifdef CHEMISTRY_PORTABLE
#include "chemistry_init.h"
#include "chemistry_standalone_adapter.h"
#include "chemistry_standalone_types.h"
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define nH_MAX 1.0e16

int set_units();
int initialize_cells(double z, double theta);
void Init_test(int ithis,double *a_start, double *dt_sec);
static void set_chemistry_table_dir(const char *argv0);
static int survey_first_chemistry_failures;
#ifndef CHEMISTRY_OLD
static int run_single_chemistry_step(const char *path, int diagnostic_mode,
                                     const char *guess_path);
static void print_single_chemistry_step_input(const struct chemdata_in *input,
                                              const struct chemdata_out *output,
                                              double temperature, double dt);
#endif

double gettimeofday_sec(void);
#ifdef CHEMISTRY_PORTABLE
int InitTable(ChemRunContext *run_context);
static void prepare_portable_chemistry_step(void);
#else
int InitTable(void);
#endif
void give_file_name(char *option);
void dump_data(int idump, char *outname, int mode);
void close_fp(void);

double update_phys(int ithis, double dt, double SHR0);
double timestep(double dt_pre);

double JeansMass(double T,double nH,double mu);
double compute_collapse_time(double SHR, double rho_gas);
int compute_Jeans_length(int ithis, double T_gas);

#ifndef CHEMISTRY_OLD
double compute_temp(int ithis)
{
  double abund[TGCHEM_NUM_ABUNDANCES] = {0};
  abund[sELEC] = SphP[ithis].elec;
  abund[sHI] = SphP[ithis].HI;
  abund[sHII] = SphP[ithis].HII;
  abund[sH2I] = SphP[ithis].H2I;
  abund[sHM] = SphP[ithis].HM;
#ifdef HELIUM
  abund[sHeI] = SphP[ithis].HeI;
  abund[sHeII] = SphP[ithis].HeII;
  abund[sHeIII] = SphP[ithis].HeIII;
#endif
#ifdef H2_PLUS
  abund[sH2II] = SphP[ithis].H2II;
#endif
#ifdef HD
  abund[sDI] = SphP[ithis].DI;
  abund[sDII] = SphP[ithis].DII;
  abund[sHDI] = SphP[ithis].HDI;
#ifdef HD_PLUS
  abund[sHDII] = SphP[ithis].HDII;
#endif
#ifdef DMINUS
  abund[sDM] = SphP[ithis].DM;
#endif
#endif
#if defined(CARBON_OXYGEN_CHEMISTRY)
  abund[sCI] = SphP[ithis].CI;
  abund[sCII] = SphP[ithis].CII;
  abund[sCO] = SphP[ithis].CO;
  abund[sCO2] = SphP[ithis].CO2;
  abund[sOI] = SphP[ithis].OI;
  abund[sOH] = SphP[ithis].OH;
  abund[sH2O] = SphP[ithis].H2O;
  abund[sO2] = SphP[ithis].O2;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  abund[sCOII] = SphP[ithis].COII;
  abund[sCH] = SphP[ithis].CH;
  abund[sCH2] = SphP[ithis].CH2;
  abund[sOII] = SphP[ithis].OII;
  abund[sOHII] = SphP[ithis].OHII;
  abund[sH2OII] = SphP[ithis].H2OII;
  abund[sH3OII] = SphP[ithis].H3OII;
  abund[sO2II] = SphP[ithis].O2II;
#endif
#endif

  double utherm_cgs = SphP[ithis].ThermalEnergy * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  return chemistry_temperature(abund, SphP[ithis].Gamma, utherm_cgs);
}

static char *trim_text(char *text)
{
  while(isspace((unsigned char)*text))
    text++;
  char *end = text + strlen(text);
  while(end > text && isspace((unsigned char)end[-1]))
    *--end = '\0';
  return text;
}

static int run_single_chemistry_step(const char *path, int diagnostic_mode,
                                     const char *guess_path)
{
  FILE *input = fopen(path, "r");
  if(!input) {
    fprintf(stderr, "cannot open chemistry-step input '%s': %s\n", path, strerror(errno));
    return 2;
  }

  double nH = -1.0, n_nuclei = -1.0, temperature = -1.0, metallicity = -1.0, dt = -1.0;
  double cmb_temperature = -1.0, utherm_input_cgs = -1.0, gamma_input = -1.0;
  double simulation_time = All.Time, scale_factor = All.cf_atime;
  int comoving_integration = All.ComovingIntegrationOn;
  MyIDType diagnostic_id = 0;
#ifdef LW_RADIATION_LOCAL
  double local_J21 = 0.0, local_rate_Hminus_pd = 0.0;
  double local_rate_H2_diss = 0.0, local_G0 = 0.0;
#ifdef HD
  double local_rate_HD_diss = 0.0;
#endif
#endif
#ifdef BACKGROUND_LW_RADIATION
  double background_J21 = All.J21;
  double background_rate_Hminus_pd = All.rate_Hminus_pd;
  double background_rate_H2_diss = All.rate_H2_diss;
#ifdef HD
  double background_rate_HD_diss = All.rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double background_G0 = All.G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
  double photoionization_rate_HI = 0.0, photoheating_rate_HI = 0.0;
#endif
#ifdef THERMAL_FEEDBACK
  double stellar_radiation_flux = 0.0;
#endif
  double abundance[TGCHEM_NUM_ABUNDANCES] = {0};
  int abundance_set[TGCHEM_NUM_ABUNDANCES] = {0};
  char line[4096];
  int line_number = 0;
  while(fgets(line, sizeof(line), input)) {
    line_number++;
    char *text = trim_text(line);
    if(*text == '\0' || *text == '#')
      continue;
    if(strcmp(text, "CHEMISTRY_RERUN_BEGIN") == 0 ||
       strcmp(text, "CHEMISTRY_RERUN_END") == 0)
      continue;
    const char *rerun_marker = "/*CHEMISTRY_RERUN_V1*/";
    char *marker = strstr(text, rerun_marker);
    if(marker)
      text = marker + strlen(rerun_marker);
    char *assignment = strtok(text, ";");
    while(assignment) {
    text = trim_text(assignment);
    char *equals = strchr(text, '=');
    if(!equals) {
      fprintf(stderr, "%s:%d: expected key=value\n", path, line_number);
      fclose(input);
      return 2;
    }
    *equals = '\0';
    char *key = trim_text(text);
    char *value_text = trim_text(equals + 1);
    char *value_end = NULL;
    errno = 0;
    double value = strtod(value_text, &value_end);
    if(errno || value_end == value_text || *trim_text(value_end) != '\0' || !isfinite(value)) {
      fprintf(stderr, "%s:%d: invalid value '%s'\n", path, line_number, value_text);
      fclose(input);
      return 2;
    }

    if(strcmp(key, "nH_cm3") == 0) nH = value;
    else if(strcmp(key, "n_nuclei_cm3") == 0 || strcmp(key, "n_nuclei") == 0) n_nuclei = value;
    else if(strcmp(key, "T_K") == 0 || strcmp(key, "T_gas") == 0) temperature = value;
    else if(strcmp(key, "T_CMB_K") == 0 || strcmp(key, "cmb_temperature") == 0) cmb_temperature = value;
    else if(strcmp(key, "Z") == 0 || strcmp(key, "metallicity") == 0) metallicity = value;
    else if(strcmp(key, "dt_s") == 0 || strcmp(key, "dt") == 0) dt = value;
    else if(strcmp(key, "Utherm_cgs") == 0 || strcmp(key, "utherm_cgs") == 0) utherm_input_cgs = value;
    else if(strcmp(key, "Gamma") == 0 || strcmp(key, "gamma") == 0) gamma_input = value;
    else if(strcmp(key, "diagnostic_id") == 0) diagnostic_id = (MyIDType)value;
    else if(strcmp(key, "comoving_integration") == 0) comoving_integration = value != 0.0;
    else if(strcmp(key, "simulation_time") == 0) simulation_time = value;
    else if(strcmp(key, "scale_factor") == 0) scale_factor = value;
#ifdef LW_RADIATION_LOCAL
    else if(strcmp(key, "local_J21") == 0) local_J21 = value;
    else if(strcmp(key, "local_rate_Hminus_pd_s-1") == 0) local_rate_Hminus_pd = value;
    else if(strcmp(key, "local_rate_H2_diss_s-1") == 0) local_rate_H2_diss = value;
    else if(strcmp(key, "local_rate_Hminus_pd") == 0) local_rate_Hminus_pd = value;
    else if(strcmp(key, "local_rate_H2_diss") == 0) local_rate_H2_diss = value;
#ifdef HD
    else if(strcmp(key, "local_rate_HD_diss") == 0) local_rate_HD_diss = value;
#endif
    else if(strcmp(key, "local_G0") == 0) local_G0 = value;
#endif
#ifdef BACKGROUND_LW_RADIATION
    else if(strcmp(key, "background_J21") == 0) background_J21 = value;
    else if(strcmp(key, "background_rate_Hminus_pd_s-1") == 0) background_rate_Hminus_pd = value;
    else if(strcmp(key, "background_rate_H2_diss_s-1") == 0) background_rate_H2_diss = value;
    else if(strcmp(key, "background_rate_Hminus_pd") == 0) background_rate_Hminus_pd = value;
    else if(strcmp(key, "background_rate_H2_diss") == 0) background_rate_H2_diss = value;
#ifdef HD
    else if(strcmp(key, "background_rate_HD_diss") == 0) background_rate_HD_diss = value;
#endif
#ifdef DUST_COOLING
    else if(strcmp(key, "background_G0") == 0) background_G0 = value;
#endif
#endif
#ifdef IONIZATION_PHOTON
    else if(strcmp(key, "photoionization_rate_HI_s-1") == 0) photoionization_rate_HI = value;
    else if(strcmp(key, "photoheating_rate_HI_erg_s-1") == 0) photoheating_rate_HI = value;
    else if(strcmp(key, "photoionization_rate_HI") == 0) photoionization_rate_HI = value;
    else if(strcmp(key, "photoheating_rate_HI") == 0) photoheating_rate_HI = value;
#endif
#ifdef THERMAL_FEEDBACK
    else if(strcmp(key, "stellar_radiation_flux") == 0) stellar_radiation_flux = value;
#endif
    else {
      const char *species_name = strncmp(key, "x0_", 3) == 0 ? key + 3 : key;
      char species_buffer[64];
      int stored_abundance_assignment = 0;
      if(strncmp(key, "yi[s", 4) == 0) {
        const char *closing = strchr(key + 4, ']');
        size_t length = closing ? (size_t)(closing - (key + 4)) : 0;
        if(!closing || closing[1] != '\0' || length == 0 || length >= sizeof(species_buffer)) {
          fprintf(stderr, "%s:%d: invalid species assignment '%s'\n", path, line_number, key);
          fclose(input);
          return 2;
        }
        memcpy(species_buffer, key + 4, length);
        species_buffer[length] = '\0';
        species_name = species_buffer;
        stored_abundance_assignment = 1;
      }
      enum SPECIES species;
      if(!chemistry_species_from_name(species_name, &species) || species >= sSTORED_END) {
        fprintf(stderr, "%s:%d: unknown key '%s'\n", path, line_number, key);
        fclose(input);
        return 2;
      }
      abundance[species] = stored_abundance_assignment ? value / XH : value;
      abundance_set[species] = 1;
    }
    assignment = strtok(NULL, ";");
    }
  }
  fclose(input);

  if(n_nuclei <= 0.0 && nH > 0.0)
    n_nuclei = nH / XH;
  if(n_nuclei <= 0.0 || temperature <= 0.0 || metallicity < 0.0 || dt <= 0.0) {
    fprintf(stderr, "chemistry-step input requires nH_cm3 or n_nuclei_cm3, T_K, Z, and dt_s\n");
    return 2;
  }
  if(cmb_temperature > 0.0)
    All.T_CMB = cmb_temperature;
  All.Time = simulation_time;
  All.ComovingIntegrationOn = comoving_integration;
  if(scale_factor > 0.0) {
    All.ScaleFactor = scale_factor;
    All.cf_atime = scale_factor;
    All.cf_a3inv = comoving_integration ?
        1.0 / (scale_factor * scale_factor * scale_factor) : 1.0;
  }
#ifdef BACKGROUND_LW_RADIATION
  All.J21 = background_J21;
  All.rate_Hminus_pd = background_rate_Hminus_pd;
  All.rate_H2_diss = background_rate_H2_diss;
#ifdef HD
  All.rate_HD_diss = background_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  All.G0 = background_G0;
#endif
#endif
#ifdef CHEMISTRY_PORTABLE
  prepare_portable_chemistry_step();
#endif
  for(int species = 0; species < sSTORED_END; species++)
    if(!abundance_set[species]) {
      fprintf(stderr, "chemistry-step input requires x0_%s\n", chemistry_species_definition(species)->name);
      return 2;
    }

  struct chemdata_in sph_in;
  struct chemdata_out sph_out;
  memset(&sph_in, 0, sizeof(sph_in));
  memset(&sph_out, 0, sizeof(sph_out));
  sph_in.ID = diagnostic_id;
  sph_in.rho_cgs = n_nuclei * PROTONMASS;
  sph_in.Density = sph_in.rho_cgs / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv);
#ifdef METALLICITY
  sph_in.Metallicity = metallicity;
#endif
#ifdef LW_RADIATION_LOCAL
  sph_in.J21 = local_J21;
  sph_in.rate_Hminus_pd = local_rate_Hminus_pd;
  sph_in.rate_H2_diss = local_rate_H2_diss;
#ifdef HD
  sph_in.rate_HD_diss = local_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  sph_in.G0 = local_G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
  sph_in.k22 = photoionization_rate_HI;
  sph_in.piHI = photoheating_rate_HI;
#endif
#ifdef THERMAL_FEEDBACK
  sph_in.stellar_rad_flux = stellar_radiation_flux;
#endif
  for(int species = 0; species < sSTORED_END; species++)
    sph_out.Abund[species] = XH * abundance[species];
  sph_out.Gamma = gamma_input > 1.0 ? gamma_input :
      compute_gamma_from_abund(sph_out.Abund, temperature);
  double molecular_weight = chemistry_molecular_weight(sph_out.Abund);
  double utherm_cgs = BOLTZMANN * temperature /
      ((sph_out.Gamma - 1.0) * molecular_weight * PROTONMASS);
  if(utherm_input_cgs > 0.0)
    utherm_cgs = utherm_input_cgs;
  sph_out.Utherm = utherm_cgs / (All.UnitEnergy_in_cgs / All.UnitMass_in_g);
  sph_out.t_chem = sph_out.t_cool = MAX_REAL_NUMBER;

  if(diagnostic_mode == 2) {
    double full_abundance[N_SP];
    if(compute_equilibrium_abundances_d(&sph_in, &sph_out, temperature, dt, full_abundance))
      for(int species = 0; species < N_SP; species++)
        printf("ENERGY_STEP_INPUT_%s=%.17g\n",
               chemistry_species_definition(species)->name, full_abundance[species]);
  }

#ifdef CARBON_OXYGEN_CHEMISTRY
  if(diagnostic_mode == 1) {
    const int grid_side = 21;
    const unsigned long molecular_carbon_count = 3 + (unsigned long)grid_side * grid_side;
    const int carbon_reservoir_modes = 2;
    const int oxygen_reservoir_modes = 4;
    const unsigned long carbon_oxygen_count = molecular_carbon_count *
        carbon_reservoir_modes * oxygen_reservoir_modes;
    const int hydrogen_modes = 4;
#ifdef HD
    const int deuterium_modes = 4;
#else
    const int deuterium_modes = 1;
#endif
    const unsigned long grid_count = carbon_oxygen_count * hydrogen_modes * deuterium_modes;
    const double carbon_total = XH * A_C_sun * metallicity;
    const double oxygen_total = XH * A_O_sun * metallicity;
    double started = gettimeofday_sec();

    for(unsigned long guess_index = 0; guess_index < grid_count; guess_index++) {
      struct chemdata_out trial = sph_out;
      double guess[TGCHEM_NUM_ABUNDANCES];
      memcpy(guess, sph_out.Abund, sizeof(guess));
      unsigned long mode_index = guess_index / carbon_oxygen_count;
      unsigned long carbon_oxygen_index = guess_index % carbon_oxygen_count;
      unsigned long reservoir_index = carbon_oxygen_index / molecular_carbon_count;
      unsigned long molecular_carbon_index = carbon_oxygen_index % molecular_carbon_count;
      int oxygen_reservoir_mode = (int)(reservoir_index % oxygen_reservoir_modes);
      int carbon_reservoir_mode = (int)(reservoir_index / oxygen_reservoir_modes);
      int hydrogen_mode = (int)(mode_index % hydrogen_modes);
      int deuterium_mode = (int)(mode_index / hydrogen_modes);

      if(hydrogen_mode != 1) {
        guess[sHI] = guess[sHII] = guess[sH2I] = guess[sHM] = 0.0;
      }
#ifdef HD
      if(deuterium_mode != 1) {
        guess[sDI] = guess[sDII] = guess[sHDI] = 0.0;
        if(deuterium_mode == 0)
          guess[sHDI] = D_number_fraction;
        else if(deuterium_mode == 2)
          guess[sDI] = D_number_fraction;
        else
          guess[sDII] = D_number_fraction;
      }
#endif
      guess[sCI] = guess[sCII] = guess[sCO] = guess[sCO2] = 0.0;
      guess[sOI] = guess[sOH] = guess[sH2O] = guess[sO2] = 0.0;

      double carbon_left = carbon_total;
      double oxygen_left = oxygen_total;
      if(molecular_carbon_index == 0) {
        guess[sCO] = standalone_min(carbon_left, oxygen_left);
      } else if(molecular_carbon_index == 1) {
        guess[sCO2] = standalone_min(carbon_left, 0.5 * oxygen_left);
      } else if(molecular_carbon_index >= 3) {
        unsigned long grid_index = molecular_carbon_index - 3;
        double co2_fraction = (double)(grid_index / grid_side) / (double)(grid_side - 1);
        double co_fraction = (double)(grid_index % grid_side) / (double)(grid_side - 1);
        guess[sCO2] = standalone_min(carbon_left, 0.5 * oxygen_left) * co2_fraction;
        carbon_left -= guess[sCO2];
        oxygen_left -= 2.0 * guess[sCO2];
        guess[sCO] = standalone_min(carbon_left, oxygen_left) * co_fraction;
      }
      carbon_left = carbon_total - guess[sCO] - guess[sCO2];
      oxygen_left = oxygen_total - guess[sCO] - 2.0 * guess[sCO2];
      if(carbon_reservoir_mode == 0)
        guess[sCII] = standalone_max(carbon_left, 0.0);
      else
        guess[sCI] = standalone_max(carbon_left, 0.0);
      if(oxygen_reservoir_mode == 0)
        guess[sH2O] = standalone_max(oxygen_left, 0.0);
      else if(oxygen_reservoir_mode == 1)
        guess[sO2] = 0.5 * standalone_max(oxygen_left, 0.0);
      else if(oxygen_reservoir_mode == 2)
        guess[sOI] = standalone_max(oxygen_left, 0.0);
      else
        guess[sOH] = standalone_max(oxygen_left, 0.0);

      if(hydrogen_mode != 1) {
        double hydrogen_in_other_species = 0.0;
        for(int species = 0; species < sSTORED_END; species++)
          if(species != sHI && species != sHII && species != sH2I)
            hydrogen_in_other_species +=
                chemistry_species_definition(species)->nuclei[CHEM_ELEMENT_H] *
                guess[species];
        double hydrogen_available = standalone_max(XH - hydrogen_in_other_species, 0.0);
        if(hydrogen_mode == 0)
          guess[sH2I] = 0.5 * hydrogen_available;
        else if(hydrogen_mode == 2)
          guess[sHI] = hydrogen_available;
        else
          guess[sHII] = hydrogen_available;
      }

      int guess_iteration_limit = molecular_carbon_index < 3 ? 1000 : 100;
      int status = compute_chemistry_implicit_guess_d(
          &sph_in, &trial, temperature, dt, guess, guess_iteration_limit);
      int nonnegative_root = 1;
      for(int species = 0; species < sSTORED_END; species++)
        if(trial.Abund[species] < 0.0)
          nonnegative_root = 0;
      if(status == 1 && nonnegative_root) {
        printf("CHEM_GRID_FOUND guess=%lu hydrogen_mode=%d deuterium_mode=%d carbon_reservoir_mode=%d oxygen_reservoir_mode=%d molecular_carbon_guess=%lu elapsed_s=%.9g\n",
               guess_index, hydrogen_mode, deuterium_mode, carbon_reservoir_mode,
               oxygen_reservoir_mode, molecular_carbon_index,
               gettimeofday_sec() - started);
        for(int species = 0; species < sSTORED_END; species++)
          printf("xguess_%s=%.17g\n", chemistry_species_definition(species)->name,
                 guess[species] / XH);
        for(int species = 0; species < sSTORED_END; species++)
          printf("root_%s=%.17g\n", chemistry_species_definition(species)->name,
                 trial.Abund[species] / XH);
        return 0;
      }
      if((guess_index + 1) % 1000 == 0)
        fprintf(stderr, "CHEM_GRID_PROGRESS guesses=%lu elapsed_s=%.9g\n",
                guess_index + 1, gettimeofday_sec() - started);
    }
    fprintf(stderr, "CHEM_GRID_NOT_FOUND guesses=%lu elapsed_s=%.9g\n",
            grid_count, gettimeofday_sec() - started);
    return 1;
  }
#else
  (void)diagnostic_mode;
#endif

  double started = gettimeofday_sec();
  if(diagnostic_mode == 2)
    sph_in.flag_print = 1;
  int status;
#ifdef CHEMISTRY_PORTABLE
  prepare_portable_chemistry_step();
#endif
  if(diagnostic_mode == 2) {
    status = compute_abundances_d(1, &sph_in, &sph_out, dt);
  } else if(diagnostic_mode == 3 || diagnostic_mode == 4 || diagnostic_mode == 5) {
    FILE *guess_input = fopen(guess_path, "r");
    if(!guess_input) {
      fprintf(stderr, "cannot open chemistry guess '%s': %s\n",
              guess_path, strerror(errno));
      return 2;
    }
    double newton_guess[TGCHEM_NUM_ABUNDANCES] = {0};
    int guess_set[TGCHEM_NUM_ABUNDANCES] = {0};
    int guess_line_number = 0;
    while(fgets(line, sizeof(line), guess_input)) {
      guess_line_number++;
      char *text = trim_text(line);
      if(*text == '\0' || *text == '#')
        continue;
      char *equals = strchr(text, '=');
      if(!equals)
        continue;
      *equals = '\0';
      char *key = trim_text(text);
      if(strncmp(key, "guess_", 6) != 0)
        continue;
      char *value_text = trim_text(equals + 1);
      char *value_end = NULL;
      errno = 0;
      double value = strtod(value_text, &value_end);
      enum SPECIES species;
      if(errno || value_end == value_text || *trim_text(value_end) != '\0' ||
         !isfinite(value) ||
         !chemistry_species_from_name(key + 6, &species) || species >= sSTORED_END) {
        fprintf(stderr, "%s:%d: invalid chemistry guess '%s'\n",
                guess_path, guess_line_number, key);
        fclose(guess_input);
        return 2;
      }
      newton_guess[species] = XH * value;
      guess_set[species] = 1;
    }
    fclose(guess_input);
    for(int species = 0; species < sSTORED_END; species++)
      if(!guess_set[species]) {
        fprintf(stderr, "chemistry guess requires guess_%s\n",
                chemistry_species_definition(species)->name);
        return 2;
      }
    if(diagnostic_mode == 3) {
      status = compute_chemistry_implicit_guess_d(
          &sph_in, &sph_out, temperature, dt, newton_guess, 1000);
    } else if(diagnostic_mode == 5) {
      const int grid_intervals = 100;
      double grid_started = gettimeofday_sec();
      unsigned long grid_trials = 0;
      for(int hydrogen_index = 0; hydrogen_index <= grid_intervals;
          hydrogen_index++) {
        double hydrogen_fraction = (double)hydrogen_index / grid_intervals;
        for(int carbon_oxygen_index = 0; carbon_oxygen_index <= grid_intervals;
            carbon_oxygen_index++) {
          double carbon_oxygen_fraction =
              (double)carbon_oxygen_index / grid_intervals;
          double trial_guess[TGCHEM_NUM_ABUNDANCES];
          memcpy(trial_guess, sph_out.Abund, sizeof(trial_guess));
          const int hydrogen_species[] = {sELEC, sHI, sHII, sH2I, sHM};
          for(size_t member = 0;
              member < sizeof(hydrogen_species) / sizeof(hydrogen_species[0]);
              member++) {
            int species = hydrogen_species[member];
            trial_guess[species] =
                (1.0 - hydrogen_fraction) * sph_out.Abund[species] +
                hydrogen_fraction * newton_guess[species];
          }
#ifdef HD
          const int deuterium_species[] = {sDI, sDII, sHDI};
          for(size_t member = 0;
              member < sizeof(deuterium_species) / sizeof(deuterium_species[0]);
              member++) {
            int species = deuterium_species[member];
            trial_guess[species] =
                (1.0 - hydrogen_fraction) * sph_out.Abund[species] +
                hydrogen_fraction * newton_guess[species];
          }
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
          for(int species = sCI; species < sSTORED_END; species++)
            trial_guess[species] =
                (1.0 - carbon_oxygen_fraction) * sph_out.Abund[species] +
                carbon_oxygen_fraction * newton_guess[species];
#endif
          struct chemdata_out trial = sph_out;
          int trial_status = compute_chemistry_implicit_guess_d(
              &sph_in, &trial, temperature, dt, trial_guess, 1000);
          grid_trials++;
          if(trial_status == 1) {
            printf("CHEM_GUESS_GRID2_FOUND trials=%lu H_fraction=%.2f CO_fraction=%.2f elapsed_s=%.9g\n",
                   grid_trials, hydrogen_fraction, carbon_oxygen_fraction,
                   gettimeofday_sec() - grid_started);
            for(int species = 0; species < sSTORED_END; species++)
              printf("guess_%s=%.17g\n",
                     chemistry_species_definition(species)->name,
                     trial_guess[species] / XH);
            for(int species = 0; species < sSTORED_END; species++)
              printf("root_%s=%.17g\n",
                     chemistry_species_definition(species)->name,
                     trial.Abund[species] / XH);
            return 0;
          }
        }
      }
      printf("CHEM_GUESS_GRID2_NOT_FOUND trials=%lu elapsed_s=%.9g\n",
             grid_trials, gettimeofday_sec() - grid_started);
      return 1;
    } else {
      int success_count = 0;
      double first_success = -1.0, last_success = -1.0;
      double scan_started = gettimeofday_sec();
      const int scan_intervals = 1000;
      for(int scan_index = 0; scan_index <= scan_intervals; scan_index++) {
        double fraction = (double)scan_index / scan_intervals;
        double trial_guess[TGCHEM_NUM_ABUNDANCES];
        for(int species = 0; species < sSTORED_END; species++)
          trial_guess[species] =
              (1.0 - fraction) * sph_out.Abund[species] +
              fraction * newton_guess[species];
        struct chemdata_out trial = sph_out;
        int trial_status = compute_chemistry_implicit_guess_d(
            &sph_in, &trial, temperature, dt, trial_guess, 1000);
        if(trial_status == 1) {
          if(first_success < 0.0)
            first_success = fraction;
          last_success = fraction;
          success_count++;
          double scan_dii = 0.0;
#ifdef HD
          scan_dii = trial_guess[sDII] / XH;
#endif
          printf("CHEM_GUESS_SCAN_SUCCESS fraction=%.3f HII=%.17g DII=%.17g",
                 fraction, trial_guess[sHII] / XH, scan_dii);
#ifdef CARBON_OXYGEN_CHEMISTRY
          printf(" CO=%.17g CO2=%.17g",
                 trial_guess[sCO] / XH, trial_guess[sCO2] / XH);
#endif
          printf("\n");
        }
      }
      printf("CHEM_GUESS_SCAN_SUMMARY successes=%d first=%.17g last=%.17g elapsed_s=%.9g\n",
             success_count, first_success, last_success,
             gettimeofday_sec() - scan_started);

      for(int group_mask = 1; group_mask < 16; group_mask++) {
        double trial_guess[TGCHEM_NUM_ABUNDANCES];
        memcpy(trial_guess, sph_out.Abund, sizeof(trial_guess));
        if(group_mask & 1) {
          trial_guess[sELEC] = newton_guess[sELEC];
          trial_guess[sHI] = newton_guess[sHI];
          trial_guess[sHII] = newton_guess[sHII];
          trial_guess[sH2I] = newton_guess[sH2I];
          trial_guess[sHM] = newton_guess[sHM];
        }
#ifdef HD
        if(group_mask & 2) {
          trial_guess[sDI] = newton_guess[sDI];
          trial_guess[sDII] = newton_guess[sDII];
          trial_guess[sHDI] = newton_guess[sHDI];
        }
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
        if(group_mask & 4)
          for(int species = sCI; species < sCOII; species++)
            trial_guess[species] = newton_guess[species];
        if(group_mask & 8)
          for(int species = sCOII; species < sSTORED_END; species++)
            trial_guess[species] = newton_guess[species];
#endif
        struct chemdata_out trial = sph_out;
        int trial_status = compute_chemistry_implicit_guess_d(
            &sph_in, &trial, temperature, dt, trial_guess, 1000);
        printf("CHEM_GUESS_GROUPS H=%d D=%d CO_stored=%d CO_intermediate=%d converged=%d\n",
               !!(group_mask & 1), !!(group_mask & 2), !!(group_mask & 4),
               !!(group_mask & 8),
               trial_status == 1);
      }
      for(int reverted_species = 0; reverted_species < sSTORED_END;
          reverted_species++) {
        double trial_guess[TGCHEM_NUM_ABUNDANCES];
        memcpy(trial_guess, newton_guess, sizeof(trial_guess));
        trial_guess[reverted_species] = sph_out.Abund[reverted_species];
        struct chemdata_out trial = sph_out;
        int trial_status = compute_chemistry_implicit_guess_d(
            &sph_in, &trial, temperature, dt, trial_guess, 1000);
        printf("CHEM_GUESS_REVERT species=%s converged=%d physical=%.17g successful=%.17g\n",
               chemistry_species_definition(reverted_species)->name,
               trial_status == 1,
               sph_out.Abund[reverted_species] / XH,
               newton_guess[reverted_species] / XH);
      }
      {
        double trial_guess[TGCHEM_NUM_ABUNDANCES];
        memcpy(trial_guess, sph_out.Abund, sizeof(trial_guess));
        const int selected_species[] = {
          sHII, sH2I, sHM,
#ifdef CARBON_OXYGEN_CHEMISTRY
          sCI, sCII, sCO, sOH, sOII, sOHII,
#endif
        };
        const int selected_count = sizeof(selected_species) / sizeof(selected_species[0]);
        for(int selected = 0; selected < selected_count; selected++) {
          int species = selected_species[selected];
          trial_guess[species] = newton_guess[species];
        }
        struct chemdata_out trial = sph_out;
        int trial_status = compute_chemistry_implicit_guess_d(
            &sph_in, &trial, temperature, dt, trial_guess, 1000);
        printf("CHEM_GUESS_SELECTED converged=%d count=%d\n",
               trial_status == 1, selected_count);
      }
      return success_count > 0 ? 0 : 1;
    }
  } else {
    status = compute_chemistry_implicit_d(&sph_in, &sph_out, temperature, dt);
  }
  int converged = diagnostic_mode == 2 ? status == 0 : status == 1;
  double elapsed = gettimeofday_sec() - started;
  double final_temperature = chemistry_temperature(sph_out.Abund, sph_out.Gamma,
      sph_out.Utherm * (All.UnitEnergy_in_cgs / All.UnitMass_in_g));
  printf("CHEM_STEP converged=%d status=%d elapsed_s=%.9g n_nuclei_cm3=%.17g T_initial_K=%.17g T_final_K=%.17g Z=%.17g dt_s=%.17g\n",
         converged, status, elapsed, n_nuclei, temperature, final_temperature, metallicity, dt);
  for(int species = 0; species < sSTORED_END; species++)
    printf("%s=%.17g\n", chemistry_species_definition(species)->name,
           sph_out.Abund[species] / XH);
  return converged ? 0 : 1;
}

static void print_single_chemistry_step_input(const struct chemdata_in *input,
                                              const struct chemdata_out *output,
                                              double temperature, double dt)
{
  fprintf(stderr, "CHEM_STEP_INPUT_BEGIN\n");
  fprintf(stderr, "n_nuclei_cm3=%.17g\n", input->rho_cgs / PROTONMASS);
  fprintf(stderr, "T_K=%.17g\n", temperature);
  fprintf(stderr, "T_CMB_K=%.17g\n", All.T_CMB);
#ifdef METALLICITY
  fprintf(stderr, "Z=%.17g\n", input->Metallicity);
#else
  fprintf(stderr, "Z=0\n");
#endif
  fprintf(stderr, "dt_s=%.17g\n", dt);
  for(int species = 0; species < sSTORED_END; species++)
    fprintf(stderr, "x0_%s=%.17g\n", chemistry_species_definition(species)->name,
            output->Abund[species] / XH);
  fprintf(stderr, "CHEM_STEP_INPUT_END\n");
}
#endif

#if defined(IONIZATION_PHOTON) || defined(X_RAY)
int init_ionization_rad(void);
int ionization_rates(void);
int make_ionization_rate_table(void);
int init_ionization_XRay(void);
void calc_phis_coll(double Ncoll_H, double Ncoll_He, double *phi_H, double *phi_He, double *kion1, double *kion2);

int init_ionization_XRay_2nd(void);
void calc_phis_coll_2nd(double Ncoll_H, double Ncoll_He, double xion, double *phiH_2nd, double *phiHe_2nd, double *gammaH_2nd, double *gammaHe_2nd);
#endif

#ifdef BACKGROUND_LW_RADIATION
int set_J21(const RadiationRateContext *radiation_context);
#endif

int main(int argc, char **argv){
  int  ifunc, i;
  int  istep_max = 10000000;
  double  theta;
  char   outputdir[255],outname[255],option[255];
#ifdef BACKGROUND_LW_RADIATION
  RadiationRateContext radiation_context;
#endif
#ifdef CHEMISTRY_PORTABLE
  ChemRunContext chemistry_run = {0};
#endif

#ifdef HIGH_Z
  theta = 1.99 * M_PI;
  StandaloneCollapseRedshift = 30.0;
  All.ExpansionParam = 1.0 / (1.0 + StandaloneCollapseRedshift);  // S12
#else
  StandaloneCollapseRedshift = 0.0;
  All.ExpansionParam = 1.0;
#endif

  double atime = All.ExpansionParam,  a_end = atime + 1.0e-10;
  All.ScaleFactor = All.ExpansionParam;
  
  All.Time     = atime;	
  All.cf_atime = All.ScaleFactor;
  All.cf_a3inv = 1.;

  All.T_CMB = STANDALONE_CMB_TEMPERATURE_0 / All.ScaleFactor;

  ifunc = set_units();
#ifndef CHEMISTRY_OLD
  set_chemistry_table_dir(argv[0]);
#endif
#ifdef CHEMISTRY_PORTABLE
  ifunc = InitTable(&chemistry_run);
#else
  ifunc = InitTable();
#endif
  if(ifunc != 0) {
    fprintf(stderr, "Failed to initialize chemistry (status=%d).\n", ifunc);
    endrun();
  }

#ifndef CHEMISTRY_OLD
  if(argc == 3 && strcmp(argv[1], "--chem-step") == 0)
    return run_single_chemistry_step(argv[2], 0, NULL);
  if(argc == 3 && strcmp(argv[1], "--chem-grid") == 0)
    return run_single_chemistry_step(argv[2], 1, NULL);
  if(argc == 3 && strcmp(argv[1], "--energy-step") == 0)
    return run_single_chemistry_step(argv[2], 2, NULL);
  if(argc == 4 && strcmp(argv[1], "--chem-guess") == 0)
    return run_single_chemistry_step(argv[2], 3, argv[3]);
  if(argc == 4 && strcmp(argv[1], "--chem-guess-scan") == 0)
    return run_single_chemistry_step(argv[2], 4, argv[3]);
  if(argc == 4 && strcmp(argv[1], "--chem-guess-grid2") == 0)
    return run_single_chemistry_step(argv[2], 5, argv[3]);
  if(argc == 2 && strcmp(argv[1], "--survey-first-failures") == 0)
    survey_first_chemistry_failures = 1;
#endif

#ifdef BACKGROUND_LW_RADIATION
  if(!radiation_rate_context_init(&radiation_context)) {
    fprintf(stderr, "Failed to initialize the radiation-rate context.\n");
    return EXIT_FAILURE;
  }
#endif

#ifdef METALLICITY
  All.Metallicity = pow(10., METALLICITY);
#endif

  const char *output_root = getenv("ONEZONE_OUTPUT_ROOT");
  if(!output_root || output_root[0] == '\0')
    output_root = "./output";
  double run_nH_max = nH_MAX;
  const char *maximum_density_text = getenv("ONEZONE_MAX_NH");
  if(maximum_density_text && maximum_density_text[0] != '\0') {
    char *maximum_density_end = NULL;
    double requested_maximum = strtod(maximum_density_text, &maximum_density_end);
    if(maximum_density_end != maximum_density_text &&
       *maximum_density_end == '\0' && requested_maximum > 0.0)
      run_nH_max = standalone_min(run_nH_max, requested_maximum);
  }
  int requested_redshift_index = -1;
  int requested_metallicity_index = -1;
  int requested_J21_index = -1;
  const char *redshift_index_text = getenv("ONEZONE_REDSHIFT_INDEX");
  const char *metallicity_index_text = getenv("ONEZONE_METALLICITY_INDEX");
  const char *J21_index_text = getenv("ONEZONE_J21_INDEX");
  if(redshift_index_text)
    requested_redshift_index = atoi(redshift_index_text);
  if(metallicity_index_text)
    requested_metallicity_index = atoi(metallicity_index_text);
  if(J21_index_text)
    requested_J21_index = atoi(J21_index_text);
  sprintf(outputdir, "%s", output_root);
  mkdir(outputdir, 0755);
  
#ifndef METALLICITY
  sprintf(outputdir,"%s/nomet",output_root);
#else
  sprintf(outputdir,"%s/metals",output_root);
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  sprintf(outputdir,"%s_CO",outputdir);
#ifdef CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
  sprintf(outputdir,"%s_EQ_FULL",outputdir);
#elif defined(CHEM_INTERMEDIATE_EQUILIBRIUM)
  sprintf(outputdir,"%s_EQ_REDUCED",outputdir);
#else
  sprintf(outputdir,"%s_NEQ",outputdir);
#endif
#endif

#ifdef HELIUM
  sprintf(outputdir,"%s_He",outputdir);
#endif

#ifdef HD
  sprintf(outputdir,"%s_HD",outputdir);
#endif

#ifdef LYMANWERNER
  sprintf(outputdir,"%s_LW",outputdir);
#endif

#ifdef CMB_HEATING
#ifdef SURVEY_CMB_REDSHIFT
  sprintf(outputdir,"%s_CMB",outputdir);
#else
  sprintf(outputdir, "%s_CMBz%d", outputdir, (int) StandaloneCollapseRedshift);
#endif
#endif

  mkdir(outputdir, 0755);

  double time1 = gettimeofday_sec();

#ifdef SURVEY_LYMANWERNER
  const int nl = 9;
  const double LW[9] = {
    0.0, 1.e-3, 1.e-2, 1.e-1, 1.0, 1.e1, 1.e2, 1.e3, 1.e4
  };

  for(int il=0;il<nl;il++)
#endif
    {
#ifdef SURVEY_LYMANWERNER
      if(requested_J21_index >= 0 && il != requested_J21_index)
        continue;
#endif
#ifdef SURVEY_CMB_REDSHIFT
      const int nz=6;
      double Redshifts[6]={0., 5., 10., 15., 20., 30.};
      
      for(int iz=0;iz<nz;iz++)
#endif
	{
#ifdef SURVEY_CMB_REDSHIFT
	  if(requested_redshift_index >= 0 && iz != requested_redshift_index)
	    continue;
#endif
#ifdef SURVEY_LYMANWERNER
	  All.LymanWerner = LW[il];
#ifdef BACKGROUND_LW_RADIATION
	  All.J21 = All.LymanWerner;
#endif
#endif

#ifdef SURVEY_CMB_REDSHIFT
	  All.ExpansionParam = All.ScaleFactor = 1./(1 + Redshifts[iz]);
	  All.T_CMB = STANDALONE_CMB_TEMPERATURE_0 / All.ScaleFactor;
	  
	  printf("[iz=%d] atime=%g\n", iz, All.ScaleFactor);
#endif

#ifdef SURVEY_METALLICITY
		  const int nmet=6;
		  double metallicity_survey[6] = {1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.};
		  for(int imet=0;imet<nmet;imet++)
#endif
		    {
#ifdef SURVEY_METALLICITY
		      if(requested_metallicity_index >= 0 && imet != requested_metallicity_index)
		        continue;
#endif
#ifdef SURVEY_METALLICITY
		      All.Metallicity = metallicity_survey[imet];
		      //All.Metallicity = pow(10., -3.25-0.05 * imet);
		      printf("[imet=%d] Metallicity=%g\n", imet, All.Metallicity);
#endif

#ifdef BACKGROUND_LW_RADIATION
	      if(!set_J21(&radiation_context)) {
	        fprintf(stderr, "Failed to prepare the uniform LW background.\n");
	        radiation_rate_context_destroy(&radiation_context);
	        return EXIT_FAILURE;
	      }
#endif
#ifdef CHEMISTRY_PORTABLE
              prepare_portable_chemistry_step();
#endif
	      
	      /****** give initial condition ******/
	      ifunc = initialize_cells(StandaloneCollapseRedshift, theta);
  
	      // INITIAL DEFORMATION PARAMETER ///// 
	      double dt_phys_init = 1.0e-20;
	      double t_phys = 0;
	      SphP[0].t_cool = MAX_REAL_NUMBER;
	      SphP[0].t_chem = MAX_REAL_NUMBER;

	      /***************************************************************************************************/
	      /******************** determine the output directory and file names ********************************/
	      /***************************************************************************************************/

	      sprintf(option,"%s","output");  
	      give_file_name(option);
	      sprintf(outname,"%s/%s.dat",outputdir,option);

	      printf("output file name is %s\n", outname);
	      /***************************************************************************************************/
	      /***************************************************************************************************/
  
#ifdef MODE_TEST
	      Init_test(i=0,&t_phys,&dt_phys_init);
	
	      sprintf(outname,"%s/%s","./","aiu");
#ifndef MODE_REPEAT
	      istep_max = 1;
#else
	      istep_max = 100;
#endif
#endif  

	      double time01 = gettimeofday_sec();
	
	      ////////////////////////////////////////////////////////////////////////////////////////////////////
	      // BEGIN MAIN LOOP /////////////////////////////////////////////////////////////////////////////////
	      ////////////////////////////////////////////////////////////////////////////////////////////////////
	      int idump = 0, istep =  0;

	      double nH_init = XH * SphP[0].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs/PROTONMASS;
	      double SHR = 1.666;
	      double rho_gas  = SphP[0].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs*All.cf_a3inv;
	      SphP[0].t_col = compute_collapse_time(1.66,rho_gas); // in sec (if rho0 is in cgs)    

	      printf("[INIT] nH=%g(rho=%g) temp=%g gamma=%g Eth=%g dt_ent=%g dt=%g l_Jeans=%g\n"
		     , nH_init, SphP[0].Density, SphP[0].temperature, SphP[0].Gamma
		     , SphP[0].ThermalEnergy,  SphP[0].DtEntropy, dt_phys_init, SphP[0].Jeans_Length);
	      fflush(stdout);	

	      double dt_phys_pre = dt_phys_init;
	      while(1){
#ifdef MODE_REPEAT
		struct sph_particle_data sph_old = SphP[0];
#endif		
		double dt_phys = timestep(dt_phys_pre);
		double dt_in_sec = dt_phys * All.UnitTime_in_s/All.HubbleParam;
		double t_end = t_phys + dt_phys;      
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE loop start: istep=%d dt=%g t=%g\n", istep, dt_phys, t_phys);
#endif
#ifndef MODE_TEST
		if(istep > 0)
		  SHR = update_phys(i=0, dt_in_sec, SHR);
#endif

		SphP[0].SHR = SHR;
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE before compute_temp\n");
#endif
		SphP[0].temperature = compute_temp(0);
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE after compute_temp: T=%g\n", SphP[0].temperature);
#endif
		
		double rho_gas  = SphP[0].Density * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs*All.cf_a3inv;		
		double nH = XH*rho_gas/PROTONMASS;


#ifdef CHEMISTRY_OLD
		int flag_implicit, flag_print=0, mode=1;
		ifunc = compute_abundances_d(mode=1, i=0, t_phys, t_end, &flag_implicit, flag_print);
#else
			struct chemdata_in sph_in;
			struct chemdata_out sph_out;
			const int ithis = 0;
			memset(&sph_in, 0, sizeof(sph_in));
			memset(&sph_out, 0, sizeof(sph_out));
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE before particle2in: i=%d dt=%g T=%g Eth=%g elec=%g HI=%g HII=%g H2I=%g HM=%g\n",
			ithis, dt_phys, SphP[ithis].temperature, SphP[ithis].ThermalEnergy, SphP[ithis].elec,
			SphP[ithis].HI, SphP[ithis].HII, SphP[ithis].H2I, SphP[ithis].HM);
#endif
			sph_in.ID = P[ithis].ID;
			sph_in.Density = SphP[ithis].Density;
			sph_in.rho_cgs = SphP[ithis].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
#ifdef METALLICITY
			sph_in.Metallicity = SphP[ithis].Metallicity;
#endif
#if defined(LYMANWERNER) && defined(LW_RADIATION_LOCAL)
			sph_in.rate_Hminus_pd = SphP[ithis].rate_Hminus_pd;
			sph_in.rate_H2_diss = SphP[ithis].rate_H2_diss;
			sph_in.J21 = SphP[ithis].J21;
#ifdef HD
			sph_in.rate_HD_diss = SphP[ithis].rate_HD_diss;
#endif
#ifdef DUST_COOLING
			sph_in.G0 = SphP[ithis].G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
			sph_in.k22 = SphP[ithis].k22;
			sph_in.piHI = SphP[ithis].piHI;
#endif
#ifdef THERMAL_FEEDBACK
			sph_in.stellar_rad_flux = SphP[ithis].stellar_radiation_flux;
#endif
			SphP[ithis].n_nuclei = sph_in.rho_cgs / PROTONMASS;
			SphP[ithis].ThermalEnergy = SphP[ithis].Entropy / (SphP[ithis].Gamma - 1.0) * pow(SphP[ithis].Density * All.cf_a3inv, SphP[ithis].Gamma - 1.0);
			sph_out.Utherm = SphP[ithis].ThermalEnergy;
			sph_out.Abund[sELEC] = SphP[ithis].elec;
			sph_out.Abund[sHI] = SphP[ithis].HI;
			sph_out.Abund[sHII] = SphP[ithis].HII;
			sph_out.Abund[sH2I] = SphP[ithis].H2I;
			sph_out.Abund[sHM] = SphP[ithis].HM;
#ifdef HELIUM
			sph_out.Abund[sHeI] = SphP[ithis].HeI;
			sph_out.Abund[sHeII] = SphP[ithis].HeII;
			sph_out.Abund[sHeIII] = SphP[ithis].HeIII;
#endif
#ifdef H2_PLUS
			sph_out.Abund[sH2II] = SphP[ithis].H2II;
#endif
#ifdef HD
			sph_out.Abund[sDI] = SphP[ithis].DI;
			sph_out.Abund[sDII] = SphP[ithis].DII;
			sph_out.Abund[sHDI] = SphP[ithis].HDI;
#ifdef HD_PLUS
			sph_out.Abund[sHDII] = SphP[ithis].HDII;
#endif
#ifdef DMINUS
			sph_out.Abund[sDM] = SphP[ithis].DM;
#endif
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
			sph_out.Abund[sCI] = SphP[ithis].CI;
			sph_out.Abund[sCII] = SphP[ithis].CII;
			sph_out.Abund[sCO] = SphP[ithis].CO;
			sph_out.Abund[sCO2] = SphP[ithis].CO2;
			sph_out.Abund[sOI] = SphP[ithis].OI;
			sph_out.Abund[sOH] = SphP[ithis].OH;
			sph_out.Abund[sH2O] = SphP[ithis].H2O;
			sph_out.Abund[sO2] = SphP[ithis].O2;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
			sph_out.Abund[sCOII] = SphP[ithis].COII;
			sph_out.Abund[sCH] = SphP[ithis].CH;
			sph_out.Abund[sCH2] = SphP[ithis].CH2;
			sph_out.Abund[sOII] = SphP[ithis].OII;
			sph_out.Abund[sOHII] = SphP[ithis].OHII;
			sph_out.Abund[sH2OII] = SphP[ithis].H2OII;
			sph_out.Abund[sH3OII] = SphP[ithis].H3OII;
			sph_out.Abund[sO2II] = SphP[ithis].O2II;
#endif
#endif
			sph_out.Gamma = SphP[ithis].Gamma;
			sph_out.t_chem = SphP[ithis].t_chem;
			sph_out.t_cool = SphP[ithis].t_cool;
			sph_out.Tdust = SphP[ithis].Tdust;
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE before compute: rho=%g U=%g gamma=%g Abund=%g %g %g %g %g\n",
			sph_in.rho_cgs, sph_out.Utherm, sph_out.Gamma,
			sph_out.Abund[sELEC], sph_out.Abund[sHI], sph_out.Abund[sHII], sph_out.Abund[sH2I], sph_out.Abund[sHM]);
#endif
      
		struct chemdata_out sph_out_initial = sph_out;
		double chemistry_input_temperature = SphP[ithis].temperature;
#ifdef CHEMISTRY_PORTABLE
		prepare_portable_chemistry_step();
#endif
		ifunc = compute_abundances_d(1, &sph_in, &sph_out, dt_in_sec);
		if(ifunc != 0)
		  {
		    fprintf(stderr,
		            "CHEM_SURVEY_FAILURE Z=%.17g T_CMB_K=%.17g istep=%d status=%d\n",
		            All.Metallicity, All.T_CMB, istep, ifunc);
		    print_single_chemistry_step_input(&sph_in, &sph_out_initial,
		                                      chemistry_input_temperature, dt_in_sec);
		    if(survey_first_chemistry_failures)
		      break;
		  }
#ifdef TRACE_1ZONE_CHEM
		fprintf(stderr, "TRACE after compute: ifunc=%d U=%g gamma=%g Abund=%g %g %g %g %g\n",
			ifunc, sph_out.Utherm, sph_out.Gamma,
			sph_out.Abund[sELEC], sph_out.Abund[sHI], sph_out.Abund[sHII], sph_out.Abund[sH2I], sph_out.Abund[sHM]);
#endif
		if(ifunc == 1){
		  sph_in.flag_print = 1;
		  ifunc = compute_abundances_d(1, &sph_in, &sph_out, dt_in_sec);
		}
		
			SphP[ithis].ThermalEnergy = sph_out.Utherm;
			SphP[ithis].elec = sph_out.Abund[sELEC];
			SphP[ithis].HI = sph_out.Abund[sHI];
			SphP[ithis].HII = sph_out.Abund[sHII];
			SphP[ithis].H2I = sph_out.Abund[sH2I];
			SphP[ithis].HM = sph_out.Abund[sHM];
#ifdef HELIUM
			SphP[ithis].HeI = sph_out.Abund[sHeI];
			SphP[ithis].HeII = sph_out.Abund[sHeII];
			SphP[ithis].HeIII = sph_out.Abund[sHeIII];
#endif
#ifdef H2_PLUS
			SphP[ithis].H2II = sph_out.Abund[sH2II];
#endif
#ifdef HD
			SphP[ithis].DI = sph_out.Abund[sDI];
			SphP[ithis].DII = sph_out.Abund[sDII];
			SphP[ithis].HDI = sph_out.Abund[sHDI];
#ifdef HD_PLUS
			SphP[ithis].HDII = sph_out.Abund[sHDII];
#endif
#ifdef DMINUS
			SphP[ithis].DM = sph_out.Abund[sDM];
#endif
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
			SphP[ithis].CI = sph_out.Abund[sCI];
			SphP[ithis].CII = sph_out.Abund[sCII];
			SphP[ithis].CO = sph_out.Abund[sCO];
			SphP[ithis].CO2 = sph_out.Abund[sCO2];
			SphP[ithis].OI = sph_out.Abund[sOI];
			SphP[ithis].OH = sph_out.Abund[sOH];
			SphP[ithis].H2O = sph_out.Abund[sH2O];
			SphP[ithis].O2 = sph_out.Abund[sO2];
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
			SphP[ithis].COII = sph_out.Abund[sCOII];
			SphP[ithis].CH = sph_out.Abund[sCH];
			SphP[ithis].CH2 = sph_out.Abund[sCH2];
			SphP[ithis].OII = sph_out.Abund[sOII];
			SphP[ithis].OHII = sph_out.Abund[sOHII];
			SphP[ithis].H2OII = sph_out.Abund[sH2OII];
			SphP[ithis].H3OII = sph_out.Abund[sH3OII];
			SphP[ithis].O2II = sph_out.Abund[sO2II];
#endif
#endif
			SphP[ithis].Gamma = sph_out.Gamma;
			SphP[ithis].t_chem = sph_out.t_chem;
			SphP[ithis].t_cool = sph_out.t_cool;
			SphP[ithis].Tdust = sph_out.Tdust;
#ifdef CHEMISTRY_STANDALONE_HOST
			SphP[ithis].tau_cont = sph_out.tau_cont;
			SphP[ithis].f_reduce = sph_out.continuum_reduction_factor;
			SphP[ithis].fesc_H2 = sph_out.diagnostics.H2_escape_fraction;
#ifdef HD
			SphP[ithis].fesc_HD = sph_out.diagnostics.HD_escape_fraction;
#endif
			for(int irate = 0; irate < k_MAX; irate++) {
			  SphP[ithis].k[irate] = sph_out.forward_rates[irate];
			  SphP[ithis].k_inv[irate] = sph_out.inverse_rates[irate];
			}
			SphP[ithis].edot_rad = sph_out.diagnostics.edot_rad;
			SphP[ithis].edot_chem = sph_out.diagnostics.edot_chem;
			SphP[ithis].edot_H2 = sph_out.diagnostics.edot_H2;
			SphP[ithis].edot_comp = sph_out.diagnostics.edot_comp;
			SphP[ithis].edot_ceHI = sph_out.diagnostics.edot_ceHI;
			SphP[ithis].edot_radHI = sph_out.diagnostics.edot_radHI;
			SphP[ithis].edot_ciHI = sph_out.diagnostics.edot_ciHI;
			SphP[ithis].edot_H2form_gas = sph_out.diagnostics.edot_H2form_gas;
			SphP[ithis].edot_H2form = sph_out.diagnostics.edot_H2form_gas;
			SphP[ithis].edot_HMform = sph_out.diagnostics.edot_HMform;
			SphP[ithis].edot_H2form_HM = sph_out.diagnostics.edot_H2form_HM;
			SphP[ithis].edot_H2form_dust = sph_out.diagnostics.edot_H2form_dust;
			SphP[ithis].edot_H2form_3body = sph_out.diagnostics.edot_H2form_3body;
			SphP[ithis].n_crit = sph_out.diagnostics.n_crit;
#ifdef HELIUM
			SphP[ithis].edot_ceHeI = sph_out.diagnostics.edot_ceHeI;
			SphP[ithis].edot_ceHeII = sph_out.diagnostics.edot_ceHeII;
			SphP[ithis].edot_ciHeI = sph_out.diagnostics.edot_ciHeI;
			SphP[ithis].edot_ciHeII = sph_out.diagnostics.edot_ciHeII;
			SphP[ithis].edot_ciHeIS = sph_out.diagnostics.edot_ciHeIS;
			SphP[ithis].edot_reHeII1 = sph_out.diagnostics.edot_reHeII1;
			SphP[ithis].edot_reHeII2 = sph_out.diagnostics.edot_reHeII2;
			SphP[ithis].edot_reHeIII = sph_out.diagnostics.edot_reHeIII;
#endif
#ifdef HD
			SphP[ithis].edot_HD = sph_out.diagnostics.edot_HD;
#endif
#ifdef METALLICITY
			SphP[ithis].edot_CI = sph_out.diagnostics.edot_CI;
			SphP[ithis].edot_CII = sph_out.diagnostics.edot_CII;
			SphP[ithis].edot_OI = sph_out.diagnostics.edot_OI;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
			SphP[ithis].edot_CO = sph_out.diagnostics.edot_CO;
			SphP[ithis].edot_OH = sph_out.diagnostics.edot_OH;
			SphP[ithis].edot_H2O = sph_out.diagnostics.edot_H2O;
#endif
#ifdef DUST_COOLING
			SphP[ithis].edot_dust = sph_out.diagnostics.edot_dust;
#ifdef LYMANWERNER
			SphP[ithis].edot_PE = sph_out.diagnostics.edot_PE;
#endif
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
			SphP[ithis].edot_cont = sph_out.diagnostics.edot_cont;
#else
			SphP[ithis].edot_reHII = sph_out.diagnostics.edot_reHI;
			SphP[ithis].edot_brem = sph_out.diagnostics.edot_brem;
			SphP[ithis].edot_CIE = sph_out.diagnostics.edot_CIE;
			SphP[ithis].edot_cont = sph_out.diagnostics.edot_cont;
#endif
#endif
			SphP[ithis].Entropy = SphP[ithis].ThermalEnergy * (SphP[ithis].Gamma - 1.0) / pow(SphP[ithis].Density * All.cf_a3inv, SphP[ithis].Gamma - 1.0);
#endif
		
		SphP[0].temperature = compute_temp(0);
		
		t_phys += dt_phys;
		All.Time = t_phys;
    
		if(istep%100==0) 
		  {
		    printf("[step%d] (%gsec) t=%g, nH=%g(rho=%g) temp=%g tdust=%g gamma=%g SHR=%g Eth=%g dt_ent=%g dt=%g elec=%g H2I=%g\n"
			   , istep
			   , gettimeofday_sec() - time1
			   , t_phys * All.UnitTime_in_s/All.HubbleParam, nH, SphP[0].Density, SphP[0].temperature,SphP[0].Tdust, SphP[0].Gamma
			   , SHR,SphP[0].ThermalEnergy,  SphP[0].DtEntropy, dt_phys, SphP[0].elec
			   , SphP[0].H2I);
		  }

		if(nH >= pow(10.0,idump/10.0-1.0) * nH_init){
		  dump_data(idump++, outname, 0);
		}

		dt_phys_pre = dt_phys;
		istep++;
  
		if(istep == istep_max) break;
		if(nH > run_nH_max || SphP[0].temperature > 2.e6)
		  break;
		
#ifdef MODE_REPEAT
		SphP[0] = sph_old;
#endif
	      }
      
	      close_fp();
	    }
	}
    }

  double time2 = gettimeofday_sec();
  
  printf("time:%g\n",time2 - time1);
#if defined(CHEMISTRY_PROFILE) && !defined(CHEMISTRY_OLD)
  chemistry_profile_report();
#endif
#ifdef CHEMISTRY_PORTABLE
  DestroyChem(&chemistry_run);
#endif
#ifdef BACKGROUND_LW_RADIATION
  radiation_rate_context_destroy(&radiation_context);
#endif
}

void give_file_name(char *option){
#ifdef ENZO_CHEM
  sprintf(option,"%s_%s",option,"EZ");
#elif MARTIN97
  sprintf(option,"%s_%s",option,"MT");
#elif GA08
  sprintf(option,"%s_%s",option,"GA");
#endif
	
#ifdef HD
  sprintf(option,"%s_%s",option,"HD");
#endif

#ifdef HLINE
  sprintf(option,"%s_%s",option,"Hline");
#endif
	
#ifdef METAL
  sprintf(option,"%s_%s%g",option,"ME", METAL);
#endif

#ifdef METALLICITY
  sprintf(option,"%s_%s%g",option,"ME", (double) (log10(All.Metallicity)));

#ifndef DUST_COOLING
  sprintf(option,"%s_%s",option,"nodust");
#endif

#ifndef METAL_COOLING_LW
  sprintf(option,"%s_%s",option,"noline");
#endif
#endif
	
#ifdef HYBRID_IMPLICIT
  sprintf(option,"%s_%s",option,"HYB");
#elif ENERGY_IMPLICIT
  sprintf(option,"%s_%s",option,"EIMP");
#elif CHEMISTRY_IMPLICIT
  sprintf(option,"%s_%s",option,"CIMP");
#else
  sprintf(option,"%s_%s",option,"EXP");
#endif
	
#ifdef LYMANWERNER  
#ifdef TYPE_A
  sprintf(option,"%s_A_%g",option,All.LymanWerner);
#elif TYPE_B
  sprintf(option,"%s_B_%g",option,All.LymanWerner);
#elif TYPE_3E4
  sprintf(option,"%s_Teff_3E4_%g",option,All.LymanWerner);
#elif TYPE_POW
  sprintf(option,"%s_P_%g",option,All.LymanWerner);
#endif	

#ifndef SHIELD_LW
  sprintf(option,"%s_noshield",option);
#endif

#ifdef NO_H2_DESTRUCTION
  sprintf(option,"%s_NOH2Diss",option);
#endif
#endif

#ifdef CMB_HEATING
  sprintf(option,"%s_CMBz%d",option, (int)(1./All.ExpansionParam - 1));
#endif


#ifdef SLOW_COLLAPSE
  sprintf(option, "%s_BETA_%g", option, (double)SLOW_COLLAPSE);
#endif

#ifdef X_RAY
  sprintf(option,"%s_XR_%g",option,X_RAY);
#endif

#ifdef HIGH_TEMP
  sprintf(option,"%s_highT",option);
#endif

#ifdef FULLY_MOLECULAR
  sprintf(option,"%s_Molecular",option);
#endif
}



static FILE *fpout, *fcool, *fchem, *fchem_inv, *fH2, *fHM, *fS, *felec, *fHI;
static int flag_file_open = 0;
static int f_count = 0;

static FILE *open_output_file(const char *path, const char *mode){
  FILE *file = fopen(path, mode);
  if(file == NULL){
    fprintf(stderr, "Failed to open output file '%s': %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

void dump_data(int idump, char *outname, int mode){

  int i;
	
  static double  dens_out, nHout, Tout, Tdust, ab_H2I, ab_HDI, ab_elec, time;
  static double  ab_HI, ab_HII, ab_DI, ab_DII;
  static double  ab_H2II, ab_HDII, ab_HM, ab_DM;
  static double  ab_HeI, ab_HeII, ab_HeIII;
  static double  fesc_line, fesc_cont;
  static double  AdiaIndex_out, MMolWeigh_out;
  static double  M_Jeans;
  static double  cool_ceHI, cool_ceHeI, cool_ceHeII;
  static double  cool_ciHI, cool_ciHeI, cool_ciHeII, cool_ciHeIS;
  static double  cool_ciHI0;
  static double  cool_reHII, cool_reHeII1, cool_reHeII2, cool_reHeIII;
  static double  cool_comp, cool_brem, coolRad;
  static double  cool_cont, cool_2ph;
  static double  cool_CI, cool_CII, cool_OI, cool_CO, cool_OH, cool_H2O, cool_dust, cool_PE;
  static double  heat_radHI, heat_radHeI, heat_radHeII;
  static double  adia, coolH2, coolHD, coolCIE, coolChem;
  static double  Escape, tau_cont;
  static double  gamma;
  double  cool_tot;
  double  coolH2form_HM, coolH2form_H2II, coolH2form_3body, coolH2form_dust, SHR, t_ff_out;
  
  char output_cool[255], output_chem[255], output_chem_inv[255], output_H2I[255], output_HM[255], output_HI[255];
	
  if(mode == 1 && f_count > 0){
    fclose(fpout),fclose(fcool),fclose(fchem),fclose(fchem_inv),fclose(fH2),fclose(fHM),fclose(felec),fclose(fS);
    fclose(fHI);
    return; //if mode is 1, just close the file.
  }
  
  if(idump != f_count) printf("something wrong in dump_data!\n");
  if(f_count == 0){
    fpout = open_output_file(outname,"w");
		
    sprintf(output_cool,"%s_cool",outname);
    fcool = open_output_file(output_cool,"w");

    sprintf(output_chem,"%s_chem",outname);
    fchem = open_output_file(output_chem,"w");

    sprintf(output_chem_inv,"%s_chem_inv",outname);
    fchem_inv = open_output_file(output_chem_inv,"w");

    sprintf(output_H2I,"%s_H2I",outname);
    fH2 = open_output_file(output_H2I,"w");

    sprintf(output_HI,"%s_HI",outname);
    fHI = open_output_file(output_HI,"w");

    sprintf(output_HM,"%s_HM",outname);
    fHM = open_output_file(output_HM,"w");

    sprintf(output_HM,"%s_elec",outname);
    felec = open_output_file(output_HM,"w");

    sprintf(output_HM,"%s_S",outname);
    fS = open_output_file(output_HM,"w");

    flag_file_open = 1;
  }else if(flag_file_open == 0){
    fpout = open_output_file(outname,"a");
		
    sprintf(output_cool,"%s_cool",outname);
    fcool = open_output_file(output_cool,"a");

    sprintf(output_chem,"%s_chem",outname);
    fchem = open_output_file(output_chem,"a");

    sprintf(output_chem_inv,"%s_chem_inv",outname);
    fchem_inv = open_output_file(output_chem_inv,"a");

    sprintf(output_H2I,"%s_H2I",outname);
    fH2 = open_output_file(output_H2I,"a");

    sprintf(output_HI,"%s_HI",outname);
    fHI = open_output_file(output_HI,"a");

    sprintf(output_HM,"%s_HM",outname);
    fHM = open_output_file(output_HM,"a");

    sprintf(output_HM,"%s_elec",outname);
    felec = open_output_file(output_HM,"a");

    sprintf(output_HM,"%s_S",outname);
    fS = open_output_file(output_HM,"a");

    flag_file_open = 1;
  }
	
  dens_out = SphP[0].Density;
  nHout    = XH*SphP[0].Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam/PROTONMASS;
  Tout     = SphP[0].temperature;
  Tdust    = SphP[0].Tdust;
  gamma    = SphP[0].Gamma;
  SHR      = SphP[0].SHR;

  ab_elec  = SphP[0].elec /XH;
  ab_HI    = SphP[0].HI   /XH;
  ab_HII   = SphP[0].HII  /XH;
#ifdef HELIUM
  ab_HeI   = SphP[0].HeI  /XH;
  ab_HeII  = SphP[0].HeII /XH;
  ab_HeIII = SphP[0].HeIII/XH;
#else
  ab_HeI   = 0.;
  ab_HeII  = 0.;
  ab_HeIII = 0.;
#endif
  ab_H2I   = SphP[0].H2I  /XH;
#ifdef H2_PLUS
  ab_H2II  = SphP[0].H2II /XH;
#endif
  ab_HM    = SphP[0].HM   /XH;
#ifdef HD
  ab_DI    = SphP[0].DI   /XH;
  ab_DII   = SphP[0].DII  /XH;	
  ab_HDI   = SphP[0].HDI  /XH;
#ifdef HD_PLUS
  ab_HDII  = SphP[0].HDII /XH;
#endif
#ifdef DMINUS
  ab_DM    = SphP[0].DM   /XH;
#endif
#else
  ab_DI    = 0.;
  ab_DII   = 0.;
  ab_HDI   = 0.;
  ab_HDII  = 0.;
  ab_DM    = 0.;
#endif

  fesc_cont = SphP[0].f_reduce;
  tau_cont = SphP[0].tau_cont;

  double sum_fraction, molecular_weight;

  sum_fraction  =  SphP[0].elec
    + SphP[0].HI + SphP[0].HII + SphP[0].HM + SphP[0].H2I 
#ifdef H2_PLUS
    + SphP[0].H2II
#endif
#ifdef HD
    + SphP[0].DI + SphP[0].DII + SphP[0].HDI 
#ifdef HD_PLUS
    + SphP[0].HDII
#endif
#ifdef DMINUS
    + SphP[0].DM
#endif
#endif
#ifdef HELIUM
    + SphP[0].HeI + SphP[0].HeII + SphP[0].HeIII
#else
    + 0.0625
#endif
    ;
  
  molecular_weight = 1.0/sum_fraction;

  double unit = 1/(All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv);
  /*conversion factor edot/density -> cgs (!edot is already in cgs unit)*/
	
	
  cool_ceHI   = fabs(SphP[0].edot_ceHI)    / SphP[0].Density * unit;
  cool_ceHeI  = fabs(SphP[0].edot_ceHeI)   / SphP[0].Density * unit;
  cool_ceHeII = fabs(SphP[0].edot_ceHeII)  / SphP[0].Density * unit;
  
  cool_ciHI   = fabs(SphP[0].edot_ciHI)    / SphP[0].Density * unit;
  cool_ciHeI  = fabs(SphP[0].edot_ciHeI)   / SphP[0].Density * unit;
  cool_ciHeII = fabs(SphP[0].edot_ciHeII)  / SphP[0].Density * unit;
  cool_ciHeIS = fabs(SphP[0].edot_ciHeIS)  / SphP[0].Density * unit;
	
  cool_reHII   = fabs(SphP[0].edot_reHII)  / SphP[0].Density * unit;
  cool_reHeII1 = fabs(SphP[0].edot_reHeII1)/ SphP[0].Density * unit;
  cool_reHeII2 = fabs(SphP[0].edot_reHeII2)/ SphP[0].Density * unit;
  cool_reHeIII = fabs(SphP[0].edot_reHeIII)/ SphP[0].Density * unit;
	
  cool_comp = fabs(SphP[0].edot_comp)      / SphP[0].Density * unit;
  cool_brem = fabs(SphP[0].edot_brem)      / SphP[0].Density * unit;
	
  adia     = SphP[0].edot_adia             / SphP[0].Density * unit;
  //adia     = SphP[0].edot_grav/SphP[0].Density * unit;
  coolH2   = fabs(SphP[0].edot_H2)         / SphP[0].Density * unit;

  coolHD   = fabs(SphP[0].edot_HD)         / SphP[0].Density * unit;
  coolCIE  = fabs(SphP[0].edot_CIE)        / SphP[0].Density * unit;
  
  coolRad  = fabs(SphP[0].edot_rad)        / SphP[0].Density * unit;
  coolChem = fabs(SphP[0].edot_chem)       / SphP[0].Density * unit;

  cool_2ph  = fabs(SphP[0].edot_2ph)        / SphP[0].Density * unit;
  cool_cont = fabs(SphP[0].edot_cont)       / SphP[0].Density * unit;

  cool_CI   = fabs(SphP[0].edot_CI)         / SphP[0].Density * unit;
  cool_CII  = fabs(SphP[0].edot_CII)        / SphP[0].Density * unit;
  cool_OI   = fabs(SphP[0].edot_OI)         / SphP[0].Density * unit;
#ifdef CARBON_OXYGEN_CHEMISTRY
  cool_CO   = fabs(SphP[0].edot_CO)         / SphP[0].Density * unit;
  cool_OH   = fabs(SphP[0].edot_OH)         / SphP[0].Density * unit;
  cool_H2O  = fabs(SphP[0].edot_H2O)        / SphP[0].Density * unit;
#else
  cool_CO = cool_OH = cool_H2O = 0.;
#endif
  cool_dust = fabs(SphP[0].edot_dust)       / SphP[0].Density * unit;
  cool_PE   = fabs(SphP[0].edot_PE)         / SphP[0].Density * unit;

  heat_radHI = SphP[0].edot_radHI / SphP[0].Density * unit;
  heat_radHeI = SphP[0].edot_radHeI / SphP[0].Density * unit;
  heat_radHeII = SphP[0].edot_radHeII / SphP[0].Density * unit;
	
  coolH2form_HM     = fabs(SphP[0].edot_H2form_HM   ) /SphP[0].Density * unit;
  coolH2form_H2II   = fabs(SphP[0].edot_H2form_H2II ) /SphP[0].Density * unit;
  coolH2form_dust   = fabs(SphP[0].edot_H2form_dust ) /SphP[0].Density * unit;
  coolH2form_3body  = fabs(SphP[0].edot_H2form_3body) /SphP[0].Density * unit;
  
  cool_tot = fabs(SphP[0].edot_rad + SphP[0].edot_chem)/SphP[0].Density * unit;
	
  SHR = SphP[0].SHR;
  t_ff_out = compute_freefall_time(SphP[0].Density) /
      STANDALONE_YEAR_IN_SECONDS;
  
  AdiaIndex_out = SphP[0].Gamma;
  MMolWeigh_out = SphP[0].molecular_weight;
  M_Jeans = JeansMass(Tout,nHout,MMolWeigh_out);
  Escape = SphP[0].fesc_H2;
  
  // WRITE OUTPUT FILES //////////////////////////////////////////////////////////////////////////////
  fprintf(fpout,"%20.10e %20.10e %20.10e"
	  " %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e"
	  " %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e"
	  " %20.10e %20.10e %20.10e %20.10e",
	  SphP[0].n_nuclei        , Tout                    , adia, // 1 2 3
	  gamma                   , SHR                     , // 4 5
	  molecular_weight        , fesc_cont               , tau_cont, // 6 7 8
	  0.                      , 0.                      , // 9 10
	  ab_elec                 , ab_HI                   , ab_HII, // 11 12 13 
	  ab_HeI                  , ab_HeII                 , ab_HeIII, // 14 15 16 
	  ab_H2I                  , ab_H2II                 , ab_HM, // 17 18 19 
	  ab_DI                   , ab_DII                  , ab_HDI, // 20 21 22 
	  ab_HDII                 , ab_DM                   , 0. , // 23 24 25 
	  Escape, SHR, t_ff_out   , Tdust); // 26 27 28 29
#ifdef CARBON_OXYGEN_CHEMISTRY
  fprintf(fpout, " %20.10e %20.10e %20.10e %20.10e"
                 " %20.10e %20.10e %20.10e %20.10e",
          SphP[0].CI / XH, SphP[0].CII / XH, SphP[0].CO / XH, SphP[0].CO2 / XH,
          SphP[0].OI / XH, SphP[0].OH / XH, SphP[0].H2O / XH, SphP[0].O2 / XH); // 30--37
#endif
  fputc('\n', fpout);
  
  fprintf(fcool,"%20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e "
	  "%20.10e "
	  "%20.10e "
	  "%20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e %20.10e\n",
	  SphP[0].n_nuclei, coolRad,  // 1 2
	  cool_ceHI , cool_ceHeI  , cool_ceHeII , // 3 4 5
	  cool_ciHI , cool_ciHeI  , cool_ciHeII , cool_ciHeIS , // 6 7 8 9
	  cool_reHII, cool_reHeII1, cool_reHeII2, cool_reHeIII, // 10 11 12 13
	  cool_comp , cool_brem   , // 14 15
	  coolH2    , coolHD      , coolCIE     , // 16 17 18
	  coolH2form_HM   , // 19
	  coolH2form_H2II , // 20
	  coolH2form_3body, // 21
	  coolChem        , // 22
	  cool_2ph  , cool_cont, //23 24
	  cool_CI, cool_CII  , cool_OI, cool_dust, cool_PE, //25 26 27 28 29
	  heat_radHI, heat_radHeI, heat_radHeII, // 30 31 32
	  adia            ,  // 33 
	  SphP[0].n_crit  , cool_tot, // 34 35
	  cool_CO, cool_OH, cool_H2O); // 36 37 38

  fprintf(fchem,
	  "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e\n",
	  SphP[0].n_nuclei, Tout, SphP[0].k[1], SphP[0].k[2]  , SphP[0].k[3] , SphP[0].k[4], 
	  SphP[0].k[5]  , SphP[0].k[6]   , SphP[0].k[7]  , SphP[0].k[8] , SphP[0].k[9],
	  SphP[0].k[10] , SphP[0].k[11]  , SphP[0].k[12] , SphP[0].k[13], SphP[0].k[14],
	  SphP[0].k[15] , SphP[0].k[16]  , SphP[0].k[17] , SphP[0].k[18], SphP[0].k[19],
	  SphP[0].k[20] , SphP[0].k[21]  , SphP[0].k[22] , SphP[0].k[23], SphP[0].k[24],
	  SphP[0].k[25] , SphP[0].k[26]  , SphP[0].k[27] , SphP[0].k[28], SphP[0].k[29],
	  SphP[0].k[30] , SphP[0].k[31]  , SphP[0].k[32] , SphP[0].k[33], SphP[0].k[34],
	  SphP[0].elec  , SphP[0].HII    , SphP[0].H2I   , SphP[0].HM   , ab_H2II);

  fprintf(fchem_inv,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e\n",
	  SphP[0].n_nuclei , Tout, SphP[0].k_inv[1], SphP[0].k_inv[2]  , SphP[0].k_inv[3] , SphP[0].k_inv[4], 
	  SphP[0].k_inv[5]  , SphP[0].k_inv[6]   , SphP[0].k_inv[7]  , SphP[0].k_inv[8] , SphP[0].k_inv[9],
	  SphP[0].k_inv[10] , SphP[0].k_inv[11]  , SphP[0].k_inv[12] , SphP[0].k_inv[13], SphP[0].k_inv[14],
	  SphP[0].k_inv[15] , SphP[0].k_inv[16]  , SphP[0].k_inv[17] , SphP[0].k_inv[18], SphP[0].k_inv[19],
	  SphP[0].k_inv[20] , SphP[0].k_inv[21]  , SphP[0].k_inv[22] , SphP[0].k_inv[23], SphP[0].k_inv[24],
	  SphP[0].k_inv[25] , SphP[0].k_inv[26]  , SphP[0].k_inv[27] , SphP[0].k_inv[28], SphP[0].k_inv[29],
	  SphP[0].k_inv[30] , SphP[0].k_inv[31]  , SphP[0].k_inv[32] , SphP[0].k_inv[33], SphP[0].k_inv[34]);


  fprintf(fHI,
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e "
	  "%20.10e %20.10e"
	  "\n",
	  SphP[0].n_nuclei,
	  SphP[0].k[1]  * SphP[0].HI   * SphP[0].elec   ,
	  SphP[0].k[2]  * SphP[0].HII  * SphP[0].elec   ,
	  SphP[0].k[7]  * SphP[0].HI   * SphP[0].elec   , 
	  SphP[0].k[8]  * SphP[0].HI   * SphP[0].HM     , 
#ifdef H2_PLUS
	  SphP[0].k[9]  * SphP[0].HI   * SphP[0].HII    ,
	  SphP[0].k[10] * SphP[0].HI   * SphP[0].H2II   ,
#else
	  0., 0.,
#endif
	  SphP[0].k[11] * SphP[0].HI   * SphP[0].H2I    ,
#ifdef H2_PLUS
	  SphP[0].k[12] * SphP[0].HII  * SphP[0].H2I    ,
#else
	  0.,
#endif
	  SphP[0].k[13] * SphP[0].H2I  * SphP[0].elec   ,
#ifdef HELIUM
	  SphP[0].k[14] * SphP[0].HI   * SphP[0].HeI    ,
	  SphP[0].k[15] * SphP[0].HII  * SphP[0].HeI    ,
#else
	  0., 0.,
#endif
	  SphP[0].k[16] * SphP[0].HM   * SphP[0].elec   ,
	  SphP[0].k[17] * SphP[0].HI   * SphP[0].HM     ,
	  SphP[0].k[18] * SphP[0].HII  * SphP[0].HM     ,
#ifdef H2_PLUS
	  SphP[0].k[20] * SphP[0].H2II * SphP[0].elec   ,
	  SphP[0].k[21] * SphP[0].H2II * SphP[0].HM     ,
#else
	  0., 0.,
#endif
	  SphP[0].k[22] * SphP[0].HI   / SphP[0].n_nuclei,
	  SphP[0].k[29] * SphP[0].H2I  / SphP[0].n_nuclei,
#ifdef H2_FORMATION_ON_DUST
	  SphP[0].k[51] * SphP[0].HI   / SphP[0].n_nuclei,
#else
	  0.,
#endif
	  SphP[0].k[55] * SphP[0].HI   * SphP[0].HI   * SphP[0].HI * SphP[0].n_nuclei,
	  SphP[0].k[57] * SphP[0].H2I  * SphP[0].HI   * SphP[0].HI * SphP[0].n_nuclei,
	  SphP[0].k[58] * SphP[0].H2I  * SphP[0].H2I  ,
	  SphP[0].k_inv[1]  * SphP[0].HII  * SphP[0].elec  * SphP[0].elec * SphP[0].n_nuclei ,
	  SphP[0].k_inv[2]  * SphP[0].HI   / SphP[0].n_nuclei ,
	  SphP[0].k_inv[7]  * SphP[0].HM   / SphP[0].n_nuclei , 
	  SphP[0].k_inv[8]  * SphP[0].H2I  * SphP[0].elec     , 
	  SphP[0].k_inv[13] * SphP[0].HI   * SphP[0].HI * SphP[0].elec * SphP[0].n_nuclei
	  );

  fprintf(fH2,
	  "%20.10le %20.10le %20.10le "
	  "%20.10le %20.10le "
	  "%20.10le %20.10le "
	  "%20.10le %20.10le "
	  "%20.10le %20.10le "
	  "%20.10le %20.10le "
	  "%20.10le %20.10le"
	  "\n",
	  SphP[0].n_nuclei,
	  SphP[0].k[8]  * SphP[0].HM   * SphP[0].HI   ,
#ifdef H2_PLUS
	  SphP[0].k[10] * SphP[0].H2II * SphP[0].HI   ,
#else
	  0.,
#endif
	  SphP[0].k[11] * SphP[0].H2I  * SphP[0].HI   , 
#ifdef H2_PLUS
	  SphP[0].k[12] * SphP[0].H2I  * SphP[0].HII  , 
#else
	  0.,
#endif
	  SphP[0].k[13] * SphP[0].H2I  * SphP[0].elec ,
#ifdef H2_PLUS
	  SphP[0].k[21] * SphP[0].H2II * SphP[0].HM   ,
#else
	  0.,
#endif
	  SphP[0].k[29] * SphP[0].H2I  / SphP[0].n_nuclei,
	  SphP[0].k[55] * SphP[0].HI   * SphP[0].HI   * SphP[0].HI * SphP[0].n_nuclei,
	  SphP[0].k[57] * SphP[0].H2I  * SphP[0].HI   * SphP[0].HI * SphP[0].n_nuclei,
	  SphP[0].k[58] * SphP[0].H2I  * SphP[0].H2I  ,
#if CHEM_RATE_COUNT > 60
	  SphP[0].k[60] * SphP[0].H2I  * SphP[0].elec ,
#else
	  0.,
#endif
	  SphP[0].k[51] * SphP[0].HI / SphP[0].n_nuclei,
	  SphP[0].k_inv[8] * SphP[0].H2I * SphP[0].elec   ,
	  SphP[0].k_inv[13] * SphP[0].HI * SphP[0].HI * SphP[0].elec * SphP[0].n_nuclei
	  );
  
  double HM_equ;
  HM_equ = 
    (SphP[0].k[7] * SphP[0].HI * SphP[0].elec 
     + SphP[0].k_inv[8]  * SphP[0].H2I * SphP[0].elec 
     + SphP[0].k_inv[18] * SphP[0].HI  * SphP[0].HI )
    /(SphP[0].k_inv[7] / SphP[0].n_nuclei
      + SphP[0].k[8] *SphP[0].HI 
      + SphP[0].k[18]*SphP[0].HII + SphP[0].k[25] / SphP[0].n_nuclei);

  fprintf(fHM,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e"
	  "\n",
	  SphP[0].n_nuclei , Tout,
	  SphP[0].HM, 
	  HM_equ, 
	  SphP[0].k[7]  * SphP[0].elec * SphP[0].HI  ,
	  SphP[0].k[8]  * SphP[0].HM   * SphP[0].HI  , 
	  SphP[0].k[16] * SphP[0].HM   * SphP[0].elec,
	  SphP[0].k[17] * SphP[0].HM   * SphP[0].HI  , 
	  SphP[0].k[18] * SphP[0].HM   * SphP[0].HII , 
#ifdef H2_PLUS
	  SphP[0].k[19] * SphP[0].HII  * SphP[0].HM  ,
	  SphP[0].k[21] * SphP[0].H2II * SphP[0].HM  , 
#else
	  0., 0.,
#endif
	  SphP[0].k_inv[7]  * SphP[0].HM  / SphP[0].n_nuclei  ,
	  SphP[0].k_inv[8]  * SphP[0].H2I * SphP[0].elec  , 
	  SphP[0].k_inv[18] * SphP[0].HI  * SphP[0].HI  , 
	  SphP[0].k[25] * SphP[0].HM   / SphP[0].n_nuclei
	  );

  fprintf(felec,
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e " 
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e "
	  "%20.10e %20.10e %20.10e\n",
	  SphP[0].n_nuclei, Tout,
	  SphP[0].k[1]  * SphP[0].elec  * SphP[0].HI   ,
	  SphP[0].k[2]  * SphP[0].elec  * SphP[0].HII  ,
#ifdef HELIUM
	  SphP[0].k[3]  * SphP[0].HeI   * SphP[0].elec ,
	  SphP[0].k[4]  * SphP[0].HeII  * SphP[0].elec ,
	  SphP[0].k[5]  * SphP[0].HeII  * SphP[0].elec ,
	  SphP[0].k[6]  * SphP[0].HeIII * SphP[0].elec ,
#else
	  0.,0.,0.,0.,
#endif
	  SphP[0].k[7]  * SphP[0].HI    * SphP[0].elec ,
	  SphP[0].k[8]  * SphP[0].HM    * SphP[0].HI   ,
	  SphP[0].k[16] * SphP[0].HM    * SphP[0].elec ,
	  SphP[0].k[17] * SphP[0].HM    * SphP[0].HI   ,
#ifdef H2_PLUS
	  SphP[0].k[19] * SphP[0].HM    * SphP[0].H2I  , 
	  SphP[0].k[20] * SphP[0].H2II  * SphP[0].elec ,
#else
	  0., 0.,
#endif
	  SphP[0].k[22] * SphP[0].HI    / SphP[0].n_nuclei,
#ifdef HELIUM
	  SphP[0].k[24] * SphP[0].HeI   / SphP[0].n_nuclei,
#else
	  0.,
#endif
	  SphP[0].k[25] * SphP[0].HM    / SphP[0].n_nuclei,
	  SphP[0].k[30] * SphP[0].HI    * SphP[0].HI,
	  SphP[0].k_inv[1]  * SphP[0].HII * SphP[0].elec * SphP[0].elec * SphP[0].n_nuclei,
	  SphP[0].k_inv[2]  * SphP[0].HI  / SphP[0].n_nuclei  ,
	  SphP[0].k_inv[7]  * SphP[0].HM  / SphP[0].n_nuclei  ,
	  SphP[0].k_inv[8]  * SphP[0].H2I * SphP[0].elec  ,
	  SphP[0].k[1]  * SphP[0].elec  * SphP[0].HI   
	  + SphP[0].k[8]  * SphP[0].HM    * SphP[0].HI
	  + SphP[0].k[25] * SphP[0].HM    / SphP[0].n_nuclei
	  + SphP[0].k[30] * SphP[0].HI    * SphP[0].HI
	  + SphP[0].k_inv[2]  * SphP[0].HI  / SphP[0].n_nuclei  
	  + SphP[0].k_inv[7]  * SphP[0].HM  / SphP[0].n_nuclei,
	  SphP[0].k[2]  * SphP[0].elec  * SphP[0].HII  
	  + SphP[0].k[7]  * SphP[0].HI    * SphP[0].elec 
	  + SphP[0].k_inv[1]  * SphP[0].HII * SphP[0].elec * SphP[0].elec * SphP[0].n_nuclei
	  + SphP[0].k_inv[8]  * SphP[0].H2I * SphP[0].elec 
	  );
	
  f_count++;

  if(f_count > 0 && f_count % 1000 == 0){
    fclose(fpout),fclose(fcool),fclose(fchem),fclose(fchem_inv),fclose(fH2),fclose(fHM),fclose(felec),fclose(fS);
    fclose(fHI);

    flag_file_open = 0;
  }
}

void close_fp(void){
  if(flag_file_open){
    fclose(fpout),fclose(fcool),fclose(fchem),fclose(fchem_inv),fclose(fH2),fclose(fHM),fclose(felec),fclose(fS);
    fclose(fHI);
    
    flag_file_open = 0;
  }
  f_count = 0;
}

#ifdef CHEMISTRY_PORTABLE
int InitTable(ChemRunContext *run_context){
#else
int InitTable(void){
#endif
  int ifunc = 0;
	
#ifdef CHEMISTRY_OLD
  ifunc = H2line_tabulate();
#endif
	
#ifdef METAL_COOLING
  ifunc = H2O_cooling_table();
  ifunc = orthoH2O_cooling_table();
  ifunc = paraH2O_cooling_table();
  ifunc = CO_cooling_table();
  ifunc = OH_cooling_table();
#endif
	
#if defined(IONIZATION_PHOTON) || defined(X_RAY)
  init_ionization_rad();
#endif

#ifdef CHEMISTRY_PORTABLE
  const ChemistryInitConfig chemistry_config = {All.ChemistryTableDir};
  ChemistryStatus chemistry_status = InitChem(&chemistry_config, run_context);
  if(chemistry_status != CHEMISTRY_STATUS_OK)
    return (int)chemistry_status;
  run_context->timescale_safety_factor = All.Epsilon;
  const ChemStandaloneUnits units = {
      All.UnitEnergy_in_cgs / All.UnitMass_in_g};
  chemistry_standalone_configure(run_context, &units);
  prepare_portable_chemistry_step();
#else
  ifunc = InitChem(); // tabulate reaction rates
#endif

#ifdef IONIZATION_PHOTON
  make_ionization_rate_table();
#else
#ifdef X_RAY
  init_ionization_XRay();
#ifdef X_RAY_2nd
  init_ionization_XRay_2nd();
#endif
#endif
#endif

  return ifunc;
}

#ifdef CHEMISTRY_PORTABLE
static void prepare_portable_chemistry_step(void)
{
  ChemGlobalStep step = {0};
  step.cmb_temperature = All.T_CMB;
  step.scale_factor = All.cf_atime;
  step.simulation_time = All.Time;
  step.comoving_integration = All.ComovingIntegrationOn;
#ifdef BACKGROUND_LW_RADIATION
  step.background_J21 = All.J21;
  step.background_rate_Hminus_pd = All.rate_Hminus_pd;
  step.background_rate_H2_diss = All.rate_H2_diss;
#ifdef HD
  step.background_rate_HD_diss = All.rate_HD_diss;
#endif
#ifdef DUST_COOLING
  step.background_G0 = All.G0;
#endif
#endif
  chemistry_standalone_prepare_step(&step);
}
#endif


#ifdef BACKGROUND_LW_RADIATION
int set_J21(const RadiationRateContext *radiation_context){
#ifdef TYPE_B
  const double Teff = 1.e4;
#else
  const double Teff = 1.e5;
#endif
  const RadiationBlackbodyInput background = {All.J21, Teff};
  RadiationRates rates;

  if(!radiation_rate_uniform_background(
         radiation_context, &background, &rates))
    return 0;

  All.rate_Hminus_pd = rates.hminus_photodetachment_rate;
#ifdef DUST_COOLING
#ifdef CHEMISTRY_OLD
  All.G0 = rates.g0 / 1.6e-3;
#else
  All.G0 = rates.g0;
#endif
#endif
  All.rate_H2_diss = rates.h2_photodissociation_rate;
#ifdef HD
  All.rate_HD_diss = rates.hd_photodissociation_rate;
#endif
  return 1;
}
#endif


int initialize_cells(double z, double theta){
  double rho_inst, Temp_inst, rho0, Tvir, n_nuclei0;
  int j, ifunc;
  double ab_H, ab_He, ab_elec, ab_HI, ab_HII, ab_HeI, ab_HeII, ab_HeIII, ab_HM, ab_DI, ab_DII, ab_HDI, ab_HDII, ab_DM;
  double ab_H2I, ab_H2II, ab_HD;
#if defined(METAL) || defined(CARBON_OXYGEN_CHEMISTRY)
  double ab_CII, ab_OI;
#endif
  double sum_fraction, egyspec_to_temp;
  
  //HERE Density MUST BE GIVEN IN PHYSICAL CGS UNITS
  for(int i = 0; i < STANDALONE_CELL_COUNT; i++){
    /* Every survey point is an independent collapse calculation. */
    memset(&SphP[i], 0, sizeof(SphP[i]));

    SphP[i].DMdensity= (9.0/2.0) * pow(theta-sin(theta),2.0) / All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * RHOCRITICAL * All.Omega0 
      * pow((1.0+z),3)/pow((1.0-cos(theta)),3);
    if(All.ComovingIntegrationOn)
      SphP[i].DMdensity= 200.0 / All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam * RHOCRITICAL * All.Omega0;
    else
      SphP[i].DMdensity= 200.0 / All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam * RHOCRITICAL * All.Omega0 * pow((1.0+z),3);
		
    StandaloneVirialDarkMatterDensity = SphP[0].DMdensity;
    SphP[i].Density = StandaloneVirialDarkMatterDensity * 0.16; // gram/cc
		
    //SphP[i].Density  = 4.5e-3 * PROTONMASS / XH;
  
    rho0 = SphP[i].Density;

    Tvir = 21.0;
    SphP[i].temperature= Tvir;

    //T0 = Tvir*10;
    n_nuclei0 = 8.9e-3;
    double T0 = 21.;

#ifdef OMUKAI05
    n_nuclei0 = 0.1;
    T0 = 300.;    
#endif

    n_nuclei0 = 0.001; 


#ifdef HIGH_TEMP
    n_nuclei0 = 1.e12;
    T0 = 1.e4;
#endif

#ifdef FULLY_MOLECULAR
    n_nuclei0 = 1.e3;
    T0=20.;
#endif

    //n_nuclei0 = 1.;
    //n_nuclei0 = 2.e6;
    rho0 = PROTONMASS * n_nuclei0;
    SphP[i].Density = rho0 / All.UnitDensity_in_cgs;

    //T0 = 6000.;
    //T0 = 8000.;
    SphP[i].temp_hold  =SphP[i].temperature;

    // PRIMORDIAL ELEMENTS /////////////////////////////////////////////////////////////////////////////
    //ab_HII  = 3.7e-4 * H_number_fraction;
    ab_HII  = 2.e-4;
#ifdef HIGH_TEMP
    ab_HII  = 2.e-1;
#endif

    //ab_HII  = 1.e-4;
    //ab_HII  = 1.0e-2 * H_number_fraction;
    ab_HM   = 1.0e-12 * H_number_fraction;
    //ab_H2I  = 2.0e-6 * H_number_fraction;
    //ab_H2I  = 2.0e-6;
    //ab_H2I  = 1.0e-6;
    ab_H2I  = 1.0e-6;

    //ab_H2I  = 3.0e-9;
    ab_H2II = 3.7e-4 * ab_H2I;

#ifdef FULLY_MOLECULAR
    ab_H2I = (0.5 - 1.e-4) * H_number_fraction;
    ab_HI  = 1.e-4 * H_number_fraction;
    ab_HII = 1.e-6 * H_number_fraction;
#endif

#ifdef HD
    ab_DII  = 1.0e-20;
    ab_DM   = 1.0e-20;
    ab_HDI  = 1.0e-9 * H_number_fraction;
    ab_HDII = 1.0e-20;
    ab_DI   = D_number_fraction - ab_HDI - ab_HDII - ab_DII - ab_DM;
#endif
    ab_HI   = H_number_fraction - 1.0*ab_HII - 2.0*ab_H2I;
		
#if defined(METAL) || defined(CARBON_OXYGEN_CHEMISTRY)
    ab_CII  = XH * A_C_sun * All.Metallicity;
    ab_OI   = XH * A_O_sun * All.Metallicity;
#endif

#ifdef METALLICITY
    SphP[i].Metallicity = All.Metallicity;
#endif
		
    ab_HeII = 1.0e-10;
    ab_HeIII= 1.0e-20;
    ab_HeI  = He_number_fraction - ab_HeII - ab_HeIII;

    SphP[i].HI   = ab_HI;
    SphP[i].HII  = ab_HII;
    SphP[i].H2I  = ab_H2I;
#ifdef H2_PLUS
    SphP[i].H2II = ab_H2II;
#endif
    SphP[i].HM   = ab_HM;

#ifdef HELIUM
    SphP[i].HeI  = ab_HeI;
    SphP[i].HeII = ab_HeII;
    SphP[i].HeIII= ab_HeIII;
#endif

#ifdef HD
    SphP[i].DI   = ab_DI;
    SphP[i].DII  = ab_DII;
    SphP[i].HDI  = ab_HDI;
#ifdef HD_PLUS
    SphP[i].HDII = ab_HDII;
#endif
#ifdef DMINUS
    SphP[i].DM   = ab_DM;
#endif
#endif
		
#ifdef METAL
    SphP[i].CII = ab_CII;
    SphP[i].OI  = ab_OI;
    SphP[i].CI  = SphP[i].OII = SphP[i].O2 = SphP[i].O2II = 1.e-10 * ab_CII;
    SphP[i].CO  = SphP[i].COII = SphP[i].CO2 = 1.e-10 * ab_CII;
    SphP[i].OH  = SphP[i].OHII = SphP[i].H2O = SphP[i].H2OII = 1.e-10 * ab_CII;
    SphP[i].H3OII = SphP[i].CH = SphP[i].CH2 = 1.e-10 * ab_CII;
#endif		
#ifdef CARBON_OXYGEN_CHEMISTRY
    SphP[i].CI = SphP[i].CO = SphP[i].CO2 = 0.0;
    SphP[i].CII = ab_CII;
    SphP[i].OI = ab_OI;
    SphP[i].OH = SphP[i].H2O = SphP[i].O2 = 0.0;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
    SphP[i].COII = SphP[i].CH = SphP[i].CH2 = 0.0;
    SphP[i].OII = SphP[i].OHII = SphP[i].H2OII = SphP[i].H3OII = SphP[i].O2II = 0.0;
#endif
#endif
		
    // compute density
    ab_elec = ab_HII + ab_HeII + 2.0 * ab_HeIII + ab_H2II  - ab_HM 
#ifdef HD
      + ab_DII + ab_HDII - ab_DM
#endif
#if defined(METAL) || defined(CARBON_OXYGEN_CHEMISTRY)
      + ab_CII
#endif			
      ;
		
    // charge neutrality
    SphP[i].elec = ab_elec;

    sum_fraction    =  
      SphP[i].elec + SphP[i].HI  + SphP[i].HII + SphP[i].HM  + SphP[i].H2I
#ifdef H2_PLUS
      + SphP[i].H2II
#endif
#ifdef HELIUM
      + SphP[i].HeI + SphP[i].HeII + SphP[i].HeIII 
#endif
#ifdef HD
      + SphP[i].DI  + SphP[i].DII  + SphP[i].HDI 
#ifdef HD_PLUS      
      + SphP[i].HDII
#endif
#ifdef DMINUS
      + SphP[i].DM
#endif
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
      + SphP[i].CI + SphP[i].CII + SphP[i].CO + SphP[i].CO2
      + SphP[i].OI + SphP[i].OH + SphP[i].H2O + SphP[i].O2
#endif
      ;
			
    SphP[i].molecular_weight = 1.0/sum_fraction;

    SphP[i].Gamma = 1.6666;
#ifdef CHEMISTRY_OLD
    ifunc = compute_gamma(i);
#else
    {
      double abund[TGCHEM_NUM_ABUNDANCES] = {0};
      abund[sELEC] = SphP[i].elec;
      abund[sHI] = SphP[i].HI;
      abund[sHII] = SphP[i].HII;
      abund[sH2I] = SphP[i].H2I;
      abund[sHM] = SphP[i].HM;
#ifdef HELIUM
      abund[sHeI] = SphP[i].HeI;
      abund[sHeII] = SphP[i].HeII;
      abund[sHeIII] = SphP[i].HeIII;
#endif
#ifdef H2_PLUS
      abund[sH2II] = SphP[i].H2II;
#endif
#ifdef HD
      abund[sDI] = SphP[i].DI;
      abund[sDII] = SphP[i].DII;
      abund[sHDI] = SphP[i].HDI;
#ifdef HD_PLUS
      abund[sHDII] = SphP[i].HDII;
#endif
#ifdef DMINUS
      abund[sDM] = SphP[i].DM;
#endif
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
      abund[sCI] = SphP[i].CI;
      abund[sCII] = SphP[i].CII;
      abund[sCO] = SphP[i].CO;
      abund[sCO2] = SphP[i].CO2;
      abund[sOI] = SphP[i].OI;
      abund[sOH] = SphP[i].OH;
      abund[sH2O] = SphP[i].H2O;
      abund[sO2] = SphP[i].O2;
#endif
      SphP[i].Gamma = compute_gamma_from_abund(abund, T0);
    }
#endif

    SphP[i].pressure = All.UnitMass_in_g/All.UnitEnergy_in_cgs * rho0 /PROTONMASS /sum_fraction * BOLTZMANN * T0 ;
    SphP[i].ThermalEnergy = SphP[i].pressure/rho0/(SphP[i].Gamma-1.0);
    SphP[i].Entropy = SphP[i].ThermalEnergy * (SphP[i].Gamma-1.0) / pow(SphP[i].Density, SphP[i].Gamma-1.0);
    

    egyspec_to_temp  = (SphP[i].Gamma-1.0)/BOLTZMANN * PROTONMASS * SphP[i].molecular_weight;
    SphP[i].ThermalEnergy = T0 / egyspec_to_temp / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);
    
    SphP[i].Density = rho0 /  All.UnitDensity_in_cgs / (All.HubbleParam * All.HubbleParam);
    SphP[i].Entropy = SphP[i].ThermalEnergy * (SphP[i].Gamma-1.0) / pow(SphP[i].Density, SphP[i].Gamma-1.0);
    
    P[i].ID = i;
    for(j=0;j<3;j++) P[i].Pos[j] = 0.0;


    double column_dens_H2, column_dens_HI, column_dens_HD;
		
    ifunc = compute_Jeans_length(i, T0);
    column_dens_H2 = SphP[i].Jeans_Length * SphP[i].H2I * n_nuclei0;
    column_dens_HI = SphP[i].Jeans_Length * SphP[i].HI  * n_nuclei0;
#ifdef HD
    column_dens_HD = SphP[i].Jeans_Length * SphP[i].HDI * n_nuclei0;
#endif

#ifdef SHIELD_6RAY				
    for(j=0;j<6;j++){
      SphP[i].ColumnDens_H2[j] = column_dens_H2;
      SphP[i].ColumnDens_HI[j] = column_dens_HI;
#ifdef HD 
     SphP[i].ColumnDens_HD[j] = column_dens_HD;
#endif
    }
#endif    

#ifdef LYMANWERNER
    All.J21 = All.LymanWerner;
    SphP[i].J21 = All.J21;
    
#ifdef DUST_COOLING
    SphP[i].G0 = All.G0;
#endif
#endif

#ifdef X_RAY
    double k22, k24;
    double phi_H, phi_He;
    double column_dens_HeI;

    All.X_J21 = X_RAY * 1e-21;

    column_dens_HeI = SphP[i].Jeans_Length * SphP[i].HeI * n_nuclei0;

    calc_phis_coll(column_dens_HI, column_dens_HeI, &phi_H, &phi_He, &k22, &k24);
    SphP[i].k22 = k22 * All.X_J21;
    SphP[i].k24 = k24 * All.X_J21;
  
    SphP[i].piHI = phi_H * All.X_J21;
    SphP[i].piHeI= phi_He * All.X_J21;
#ifdef X_RAY_2nd
    double phiH_2nd, phiHe_2nd, gammaH_2nd, gammaHe_2nd;
    
    calc_phis_coll_2nd(column_dens_HI, column_dens_HeI, SphP[0].elec, &phiH_2nd, &phiHe_2nd, &gammaH_2nd, &gammaHe_2nd);
    
    SphP[i].k22_2nd = phiH_2nd;
    SphP[i].k24_2nd = phiHe_2nd;
    SphP[i].piHI  = gammaH_2nd * All.X_J21;
    SphP[i].piHeI = gammaHe_2nd * All.X_J21;
#endif
#endif
	
		
  }
  //not used in non-collapse run
  //StandaloneVirialDarkMatterDensity = 8.0 * SphP[0].DMdensity;
  
  return 0;
}

#ifdef MODE_TEST
void Init_test(int ithis,double *a_start,double *dt_sec){
  int i=0, ifunc;
  int flag_column = 0;
  double fac, J21;
  double rho, T_gas, y[14], yi[TGCHEM_NUM_ABUNDANCES];
  double y_HI,y_HII,y_HM,y_HeI,y_HeII,y_HeIII,y_elec,y_H2;
  double x, gamma_H2_minus1_inv, denominator,sum_y;
  double n_nuclei, n_H, dt;
  double egyspec_to_temp, sum_fraction, molecular_weight;
  double rho_cgs, a_start0;
  double gamma, Metallicity = 0.0;
#ifdef CHEMISTRY_PORTABLE
  MyIDType diagnostic_id;
  int comoving_integration;
  double simulation_time, scale_factor, cmb_temperature, metallicity;
  double utherm_cgs;
#ifdef LW_RADIATION_LOCAL
  double local_J21, local_rate_Hminus_pd, local_rate_H2_diss;
#ifdef HD
  double local_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double local_G0;
#endif
#endif
#ifdef BACKGROUND_LW_RADIATION
  double background_J21, background_rate_Hminus_pd, background_rate_H2_diss;
#ifdef HD
  double background_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double background_G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
  double photoionization_rate_HI, photoheating_rate_HI;
#endif
#ifdef THERMAL_FEEDBACK
  double stellar_radiation_flux;
#endif
#endif
  
#ifdef METALLICITY
  Metallicity = pow(10., METALLICITY);
#endif

  a_start0 = *a_start;

#ifndef CHEMISTRY_PORTABLE
  All.ComovingIntegrationOn=0, All.Time=a_start0=304000, dt=8.64578e+07, n_nuclei=5.89166e+07, rho=1.65878e-10, T_gas=207.72, gamma=1.4186;
  yi[0]=4.22719e-09, yi[1]=1.0251e-10, yi[2]=4.22718e-09;
  yi[3]=0.375, yi[4]=1.24678e-19;
  yi[5]=2.25552e-16, yi[6]=1.20735e-14, yi[7]=2.25e-05;
  SphP[ithis].l_core=1.78848e+16;
  SphP[ithis].tau_ion=0, SphP[ithis].k22=154.376, SphP[ithis].piHI=9.44418e-10;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304001, dt=2.70181e+06, n_nuclei=7.08012e+11, rho=1.99339e-06, T_gas=95.2861, gamma=1.4186;
  yi[0]=2.31543e-07, yi[1]=7.5e-13, yi[2]=2.31543e-07;
  yi[3]=0.375, yi[4]=4.70937e-44;
  yi[5]=7.69221e-14, yi[6]=1.05697e-13, yi[7]=2.25e-05;
  SphP[ithis].l_core=5.08141e+16;
  SphP[ithis].tau_ion=4.54305e+06, SphP[ithis].k22=2.39405, SphP[ithis].piHI=1.27227e-10;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304001, dt=2.70181e+06, n_nuclei=7.08012e+11, rho=1.99339e-06, T_gas=95.2861, gamma=1.4186;
  yi[0]=2.31543e-07, yi[1]=7.5e-13, yi[2]=2.31543e-07;
  yi[3]=0.375, yi[4]=4.70937e-44;
  yi[5]=7.69221e-14, yi[6]=1.05697e-13, yi[7]=2.25e-05;
  SphP[ithis].l_core=5.08141e+16;
  SphP[ithis].tau_ion=4.54305e+06, SphP[ithis].k22=2.39405, SphP[ithis].piHI=1.27227e-10;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304000, dt=4.32289e+07, n_nuclei=1.62925e+08, rho=4.58711e-10, T_gas=325.112, gamma=1.4186;
  yi[0]=1.07368e-12, yi[1]=8.07211e-18, yi[2]=1.07368e-12;
  yi[3]=0.375, yi[4]=1.5498e-21;
  yi[5]=3.86509e-22, yi[6]=6.80239e-18, yi[7]=2.25e-05;
  SphP[ithis].l_core=7.36039e+14;
  SphP[ithis].tau_ion=0.0191457, SphP[ithis].k22=2.76558e2, SphP[ithis].piHI=0.000103011e-5;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304005, dt=7.7812e+08, n_nuclei=3.7906e+07, rho=1.06723e-10, T_gas=288.388, gamma=1.4186;
  yi[0]=1.84249e-07, yi[1]=7.5e-22, yi[2]=1.84248e-07;
  yi[3]=0.375, yi[4]=7.5e-56;
  yi[5]=1.34227e-25, yi[6]=7.81063e-13, yi[7]=2.25e-05;
  SphP[ithis].l_core=1.13428e+15;
  SphP[ithis].tau_ion=9.17998, SphP[ithis].k22=3.39036e-07, SphP[ithis].piHI=4.43952e-18;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304080, dt=2638.48, n_nuclei=4.77603e+12, rho=1.34468e-05, T_gas=16020.1, gamma=1.66667;
  yi[0]=0.7425, yi[1]=7.5e-51, yi[2]=38305.5;
  yi[3]=7.5e-51, yi[4]=7.5e-56;
  yi[5]=7.5e-26, yi[6]=0.75, yi[7]=7.5e-31;
  SphP[ithis].l_core=8.38663e+12;
  SphP[ithis].tau_ion=1.9119e+08, SphP[ithis].k22=0, SphP[ithis].piHI=0;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304002, dt=675452, n_nuclei=1.02894e+13, rho=2.89695e-05, T_gas=3381.43, gamma=1.45568;
  yi[0]=7.2272e-40, yi[1]=3.75566e-34, yi[2]=1.27733e-20;
  yi[3]=0.249128, yi[4]=0.251746;
  yi[5]=7.5e-26, yi[6]=8.12874e-24, yi[7]=2.25001e-05;
  SphP[ithis].l_core=6.09392e+12;
  SphP[ithis].tau_ion=1.30246e+09, SphP[ithis].k22=0, SphP[ithis].piHI=0;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304002, dt=675452, n_nuclei=9.05849e+12, rho=2.5504e-05, T_gas=2734.52, gamma=1.36035;
  yi[0]=1.09624e-09, yi[1]=0.00375566, yi[2]=1.09624e-09;
  yi[3]=0.373122, yi[4]=2.76348e-14;
  yi[5]=9.15798e-08, yi[6]=2.4576e-14, yi[7]=2.24084e-05;
  SphP[ithis].l_core=5.98626e+12;
  SphP[ithis].tau_ion=1.30246e+09, SphP[ithis].k22=0, SphP[ithis].piHI=0;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304000, dt=6.91663e+08, n_nuclei=256813, rho=7.2305e-13, T_gas=27.7766, gamma=1.4186;
  yi[0]=6.58471e-12, yi[1]=1.51714e-16, yi[2]=6.58471e-12;
  yi[3]=0.375, yi[4]=2.40255e-26;
  yi[5]=2.41759e-18, yi[6]=1.02781e-23, yi[7]=2.25e-05;
  SphP[ithis].l_core=9.686e+15;
  SphP[ithis].tau_ion=134.729, SphP[ithis].k22=5.62975, SphP[ithis].piHI=1.80881e-10;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304000, dt=6.91663e+08, n_nuclei=256813, rho=7.2305e-13, T_gas=27.7766, gamma=1.4186;
  yi[0]=6.58471e-12, yi[1]=1.51714e-16, yi[2]=6.58471e-12;
  yi[3]=0.375, yi[4]=2.40255e-26;
  yi[5]=2.41759e-18, yi[6]=1.02781e-23, yi[7]=2.25e-05;
  SphP[ithis].l_core=9.686e+15;
  SphP[ithis].tau_ion=134.729, SphP[ithis].k22=5.62975, SphP[ithis].piHI=1.80881e-10;


  All.ComovingIntegrationOn=0, All.Time=a_start0=304000, dt=1.72916e+08, n_nuclei=1.05692e+07, rho=2.97574e-11, T_gas=155.732, gamma=1.4186;
  yi[0]=6.66367e-09, yi[1]=1.76201e-10, yi[2]=6.66366e-09;
  yi[3]=0.375, yi[4]=1.39049e-19;
  yi[5]=4.99784e-17, yi[6]=8.74007e-15, yi[7]=2.25e-05;
  SphP[ithis].l_core=1.82843e+15;
  SphP[ithis].tau_ion=1231.19, SphP[ithis].k22=1.76591e-15, SphP[ithis].piHI=1.73441e-25;

  All.ComovingIntegrationOn=0, All.Time=a_start0=304009, dt=1.38333e+09, n_nuclei=90853.9, rho=2.55797e-13, T_gas=2.34716, gamma=1.41862;
  yi[0]=2.04592e-05, yi[1]=2.03116e-06, yi[2]=2.04588e-05;
  yi[3]=0.374978, yi[4]=3.12239e-18;
  yi[5]=8.68692e-11, yi[6]=4.50413e-10, yi[7]=3.51422e-07;
  SphP[ithis].l_core=3.64763e+15;
  SphP[ithis].tau_ion=1.18097e+08, SphP[ithis].k22=0, SphP[ithis].piHI=0;

  All.ComovingIntegrationOn=0, All.Time=a_start0=991479, dt=1.39326e+06, n_nuclei=8.94671e+12, rho=2.51893e-05, T_gas=2159.57, gamma=1.45345;
  yi[0]=2.96356e-10, yi[1]=0.240009, yi[2]=2.96351e-10;
  yi[3]=0.254995, yi[4]=2.7915e-16;
  yi[5]=6.22067e-06, yi[6]=7.31872e-15, yi[7]=1.62793e-05;
  SphP[ithis].l_core=7.78957e+12;
  SphP[ithis].tau_ion=0, SphP[ithis].k22=0, SphP[ithis].piHI=0;
#else
  All.ComovingIntegrationOn = 0;
  All.Time = a_start0 = 304000.0;
  dt = 1.0e8;
  n_nuclei = 1.0e8;
  rho = n_nuclei * PROTONMASS;
  T_gas = 200.0;
  gamma = 5.0 / 3.0;
  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    yi[i] = 0.0;
  yi[sELEC] = 1.0e-4;
  yi[sHI] = XH - 1.0e-4 - 2.0e-8;
  yi[sHII] = 1.0e-4;
  yi[sH2I] = 1.0e-8;
  yi[sHM] = 1.0e-20;
#ifdef HELIUM
  yi[sHeI] = XHe;
  yi[sHeII] = 0.0;
  yi[sHeIII] = 0.0;
#endif
#ifdef H2_PLUS
  yi[sH2II] = 0.0;
#endif
#ifdef HD
  yi[sDI] = D_number_fraction;
  yi[sDII] = 0.0;
  yi[sHDI] = 0.0;
#ifdef HD_PLUS
  yi[sHDII] = 0.0;
#endif
#ifdef DMINUS
  yi[sDM] = 0.0;
#endif
#endif
#endif

#ifdef CHEMISTRY_PORTABLE
  diagnostic_id = P[ithis].ID;
  comoving_integration = All.ComovingIntegrationOn;
  simulation_time = All.Time;
  scale_factor = All.cf_atime;
  cmb_temperature = All.T_CMB;
  metallicity = Metallicity;
  utherm_cgs = 0.0;
#ifdef LW_RADIATION_LOCAL
  local_J21 = 0.0;
  local_rate_Hminus_pd = 0.0;
  local_rate_H2_diss = 0.0;
#ifdef HD
  local_rate_HD_diss = 0.0;
#endif
#ifdef DUST_COOLING
  local_G0 = 0.0;
#endif
#endif
#ifdef BACKGROUND_LW_RADIATION
  background_J21 = All.J21;
  background_rate_Hminus_pd = All.rate_Hminus_pd;
  background_rate_H2_diss = All.rate_H2_diss;
#ifdef HD
  background_rate_HD_diss = All.rate_HD_diss;
#endif
#ifdef DUST_COOLING
  background_G0 = All.G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
  photoionization_rate_HI = SphP[ithis].k22;
  photoheating_rate_HI = SphP[ithis].piHI;
#endif
#ifdef THERMAL_FEEDBACK
  stellar_radiation_flux = SphP[ithis].stellar_radiation_flux;
#endif

  /* Paste one complete CHEMISTRY_RERUN_V1 log line immediately below. */

  P[ithis].ID = diagnostic_id;
  All.ComovingIntegrationOn = comoving_integration;
  All.Time = a_start0 = simulation_time;
  All.ScaleFactor = All.cf_atime = scale_factor;
  All.cf_a3inv = comoving_integration ?
      1.0 / (scale_factor * scale_factor * scale_factor) : 1.0;
  All.T_CMB = cmb_temperature;
  Metallicity = metallicity;
#ifdef BACKGROUND_LW_RADIATION
  All.J21 = background_J21;
  All.rate_Hminus_pd = background_rate_Hminus_pd;
  All.rate_H2_diss = background_rate_H2_diss;
#ifdef HD
  All.rate_HD_diss = background_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  All.G0 = background_G0;
#endif
#endif
#ifdef LW_RADIATION_LOCAL
  SphP[ithis].J21 = local_J21;
  SphP[ithis].rate_Hminus_pd = local_rate_Hminus_pd;
  SphP[ithis].rate_H2_diss = local_rate_H2_diss;
#ifdef HD
  SphP[ithis].rate_HD_diss = local_rate_HD_diss;
#endif
#ifdef DUST_COOLING
  SphP[ithis].G0 = local_G0;
#endif
#endif
#ifdef IONIZATION_PHOTON
  SphP[ithis].k22 = photoionization_rate_HI;
  SphP[ithis].piHI = photoheating_rate_HI;
#endif
#ifdef THERMAL_FEEDBACK
  SphP[ithis].stellar_radiation_flux = stellar_radiation_flux;
#endif
  (void)utherm_cgs;
#endif

#ifndef CHEMISTRY_PORTABLE
  if(All.ComovingIntegrationOn){
    a3inv=  1/(a_start0*a_start0*a_start0);
  }else{
    a3inv=  1;
  }
#endif
	
  SphP[ithis].n_nuclei = n_nuclei;
  SphP[ithis].Density  = rho;
  SphP[ithis].Gamma    = gamma;  
  
  i=0;
  ///////set abundances //////////////////////
  SphP[ithis].elec = yi[i++];
  SphP[ithis].HI   = yi[i++];
  SphP[ithis].HII  = yi[i++];
#ifdef HELIUM
  SphP[ithis].HeI  = yi[i++];
  SphP[ithis].HeII = yi[i++];
  SphP[ithis].HeIII= yi[i++];
#endif
  SphP[ithis].H2I  = yi[i++];
#ifdef H2_PLUS
  SphP[ithis].H2II = yi[i++];
#endif
  SphP[ithis].HM   = yi[i++];
#ifdef HD
  SphP[ithis].DI   = yi[i++];
  SphP[ithis].DII  = yi[i++];
  SphP[ithis].HDI  = yi[i++];
#ifdef HD_PLUS
  SphP[ithis].HDII = yi[i++];
#endif
#ifdef DMINUS
  SphP[ithis].DM   = yi[i++];
#endif
#endif
#if defined(CHEMISTRY_PORTABLE) && defined(CARBON_OXYGEN_CHEMISTRY)
  SphP[ithis].CI = yi[sCI];
  SphP[ithis].CII = yi[sCII];
  SphP[ithis].CO = yi[sCO];
  SphP[ithis].CO2 = yi[sCO2];
  SphP[ithis].OI = yi[sOI];
  SphP[ithis].OH = yi[sOH];
  SphP[ithis].H2O = yi[sH2O];
  SphP[ithis].O2 = yi[sO2];
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  SphP[ithis].COII = yi[sCOII];
  SphP[ithis].CH = yi[sCH];
  SphP[ithis].CH2 = yi[sCH2];
  SphP[ithis].OII = yi[sOII];
  SphP[ithis].OHII = yi[sOHII];
  SphP[ithis].H2OII = yi[sH2OII];
  SphP[ithis].H3OII = yi[sH3OII];
  SphP[ithis].O2II = yi[sO2II];
#endif
#endif

#ifdef METALLICITY
  SphP[ithis].Metallicity = Metallicity;
#endif

    
  sum_fraction  =  SphP[ithis].elec
    + SphP[ithis].HI + SphP[ithis].HII + SphP[ithis].HM + SphP[ithis].H2I 
#ifdef H2_PLUS
    + SphP[ithis].H2II
#endif
#ifdef HD
    + SphP[ithis].DI + SphP[ithis].DII + SphP[ithis].HDI 
#ifdef HD_PLUS
    + SphP[ithis].HDII
#endif
#ifdef DMINUS
    + SphP[ithis].DM
#endif
#endif
#ifdef HELIUM
    + SphP[ithis].HeI + SphP[ithis].HeII + SphP[ithis].HeIII
#else
    + 0.0625
#endif
    ;
  
  molecular_weight = 1.0/sum_fraction;
  egyspec_to_temp  = (SphP[ithis].Gamma-1.0)/BOLTZMANN * PROTONMASS * molecular_weight;
  
  SphP[ithis].ThermalEnergy = T_gas / egyspec_to_temp / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);
#ifdef CHEMISTRY_PORTABLE
  if(utherm_cgs > 0.0)
    SphP[ithis].ThermalEnergy = utherm_cgs /
        (All.UnitEnergy_in_cgs / All.UnitMass_in_g);
#endif
  
  //////////////set density and entropy////////////////////////////
  rho_cgs = n_nuclei * PROTONMASS;
  SphP[ithis].Density = rho_cgs /  All.UnitDensity_in_cgs / (All.HubbleParam * All.HubbleParam * All.cf_a3inv);
  SphP[ithis].Entropy = SphP[ithis].ThermalEnergy * (SphP[ithis].Gamma-1.0) / pow(SphP[ithis].Density*All.cf_a3inv, SphP[ithis].Gamma-1.0);

  rho_cgs = All.UnitDensity_in_cgs * SphP[ithis].Density * All.HubbleParam * All.HubbleParam * All.cf_a3inv; /* physical gram / cm3 */
  
  /* here, mean molecular weight must be taken into account */
	
  printf("rho=%g\n",SphP[ithis].Density);
	

#ifdef SHIELD_6RAY
  double column_dens_H2, column_dens_HI, column_dens_HD;
  ifunc = compute_Jeans_length(ithis, T_gas);

  if(flag_column == 0){
    column_dens_H2 = SphP[ithis].Jeans_Length * SphP[ithis].H2I * n_nuclei;
    column_dens_HI = SphP[ithis].Jeans_Length * SphP[ithis].HI  * n_nuclei;
    column_dens_HD = SphP[ithis].Jeans_Length * SphP[ithis].HDI * n_nuclei;
    
    for(i=0;i<6;i++){
      SphP[ithis].ColumnDens_H2[i] = column_dens_H2;
      SphP[ithis].ColumnDens_HI[i] = column_dens_HI;
      SphP[ithis].ColumnDens_HD[i] = column_dens_HD;
    }
  }
#endif
	
  if(All.ComovingIntegrationOn){
    double hubble_a;
		
    hubble_a = All.Hubble*sqrt(All.Omega0/(a_start0*a_start0*a_start0) + (1-All.Omega0-All.OmegaLambda)/(a_start0*a_start0) + All.OmegaLambda);
		
    dt  *= (a_start0 * hubble_a);   
    dt /= All.UnitTime_in_s;  
    dt *= All.HubbleParam;     /* in second */
		
    *dt_sec = dt;
  }else{
    *dt_sec = dt / All.UnitTime_in_s * All.HubbleParam;
  }
	
	
  *a_start = a_start0;
}
#endif




double JeansMass(double T, double nH, double mu)
{
  double M_Jeans;

//  M_Jeans = 700.0 * pow(SphP[0].temperature/200.0,1.5) * pow(0.162*nH/1.0e4,-0.5);
   // Dopcke et al.(2011)
  M_Jeans = 1.69e3 * pow(mu/1.23,-1.5) * pow(T/200.0,1.5) * pow(nH/1.0e4,-0.5);
   // Schneider et al.(2012)

  return M_Jeans;
}

double gettimeofday_sec(){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
}

int set_units()    /* ... set some units */
{
  //All.UnitLength_in_cm= 3.085678e21; /* 1.0 kpc */
  //All.UnitLength_in_cm= 3.085678e18; /* 1.0 pc */
  All.UnitLength_in_cm= 1.496e13; /* 1.0 au */

  //All.UnitMass_in_g=    1.989e43;    /* 1.0e10 solar masses */ 
  All.UnitMass_in_g=    1.989e33;    /* 1 solar masses */ 

  All.UnitVelocity_in_cm_per_s=1.e5;  /* 1 km/sec */

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitTime_in_s= All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);
  
  All.G = GRAVITY/pow(All.UnitLength_in_cm,3)*All.UnitMass_in_g*pow(All.UnitTime_in_s,2);
  
  //All.HubbleParam = H_PARAM;
  All.HubbleParam = 1.;

  All.Omega0 = OMEGA0;
  All.OmegaLambda = OMEGALAMBDA;
  All.OmegaBaryon = OMEGABARYON;
  All.Hubble = HUBBLE * All.UnitTime_in_s;
  
  All.ComovingIntegrationOn = 0;
  All.Epsilon = 0.01;
  
  return 1;
}

static void set_chemistry_table_dir(const char *argv0)
{
  char exe_path[PATH_MAX];

  if(realpath(argv0, exe_path) == NULL){
    printf("Can't resolve executable path from argv[0]=%s.\n", argv0);
    endrun();
  }

  char *slash = strrchr(exe_path, '/');
  if(slash == NULL){
    printf("Can't determine executable directory from %s.\n", exe_path);
    endrun();
  }

  *slash = '\0';
#ifdef CHEMISTRY_OLD
  snprintf(All.ChemistryTableDir, sizeof(All.ChemistryTableDir), "%s/legacy/chemistry_old", exe_path);
#elif defined(CHEMISTRY_PORTABLE)
  snprintf(All.ChemistryTableDir, sizeof(All.ChemistryTableDir), "%s/src/chemistry", exe_path);
#else
  snprintf(All.ChemistryTableDir, sizeof(All.ChemistryTableDir), "%s/legacy/chemistry_new", exe_path);
#endif
}


void endrun(void){
  exit(1);
}
 
