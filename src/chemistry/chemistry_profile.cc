#include <stdio.h>
#include <string.h>
#include <time.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"
#include "chemistry_core.h"
#include "species.h"

#ifdef CHEMISTRY_PROFILE

typedef struct {
  double temperature, dt, density, metallicity;
  double energy_old, energy_new, residual;
  double chemical, radiative;
  double chem_H2, chem_Hion, chem_Hminus;
  double rad_H2, rad_HD, rad_HI, rad_CI, rad_CII, rad_OI;
  double rad_CO, rad_OH, rad_H2O, rad_dust, rad_continuum, rad_other;
  double CI, CII, CO, CO2, OI, OH, H2O, O2;
  double COII, CH, CH2, OII, OHII, H2OII, H3OII, O2II;
  int chemistry_status;
} ResidualTrace;

#define TRACE_COUNT 96
static double profile_seconds[CHEM_PROFILE_COUNT];
static unsigned long profile_calls[CHEM_PROFILE_COUNT];
static ResidualTrace traces[TRACE_COUNT];
static unsigned long trace_next;
static int last_failure;
static unsigned long equilibrium_initializations;
static unsigned long limited_equilibrium_initializations;
static double minimum_equilibrium_initialization_fraction = 1.0;
static unsigned long equilibrium_electron_ratio_gt_2;
static unsigned long equilibrium_electron_ratio_gt_10;
static unsigned long equilibrium_electron_ratio_gt_1e3;
static double maximum_equilibrium_electron_ratio = 1.0;
static unsigned long newton_solve_histogram[6];
static unsigned long newton_solves;
static unsigned long newton_solves_with_reset;
static unsigned long newton_solves_with_clipping;
static unsigned long newton_reset_count;
static unsigned long newton_clipped_steps;
static unsigned long newton_max_iterations;
static unsigned long clipped_species_count[N_SP];
static unsigned long newton_limiter_count[N_SP];
static unsigned long high_density_newton_limiter_count[N_SP];

double profile_start(void)
{
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);
  return time.tv_sec + 1.e-9 * time.tv_nsec;
}

void profile_end(ChemistryProfileSection section, double started_at)
{
  profile_seconds[section] += profile_start() - started_at;
  profile_calls[section]++;
}

void chemistry_profile_set_failure(int reason)
{
  last_failure = reason;
}

int chemistry_profile_last_failure(void)
{
  return last_failure;
}

void chemistry_profile_record_equilibrium_initialization(double fraction,
                                                         double electron_before,
                                                         double electron_after)
{
  equilibrium_initializations++;
  if(fraction < 1.0) {
    limited_equilibrium_initializations++;
    if(fraction < minimum_equilibrium_initialization_fraction)
      minimum_equilibrium_initialization_fraction = fraction;
  }
  double ratio = electron_after / dmax(electron_before, 1.e-300);
  if(ratio > 2.0)
    equilibrium_electron_ratio_gt_2++;
  if(ratio > 10.0)
    equilibrium_electron_ratio_gt_10++;
  if(ratio > 1.e3)
    equilibrium_electron_ratio_gt_1e3++;
  if(ratio > maximum_equilibrium_electron_ratio)
    maximum_equilibrium_electron_ratio = ratio;
}

void chemistry_profile_record_newton_solve(int iterations, int reset_count,
                                           int clipped_steps)
{
  const int boundaries[] = {2, 4, 8, 16, 69};
  int bucket = 0;
  while(bucket < 5 && iterations > boundaries[bucket])
    bucket++;
  newton_solve_histogram[bucket]++;
  newton_solves++;
  if(reset_count > 0)
    newton_solves_with_reset++;
  if(clipped_steps > 0)
    newton_solves_with_clipping++;
  newton_reset_count += (unsigned long)reset_count;
  newton_clipped_steps += (unsigned long)clipped_steps;
  if((unsigned long)iterations > newton_max_iterations)
    newton_max_iterations = (unsigned long)iterations;
}

void chemistry_profile_record_clipped_species(int species)
{
  if(species >= 0 && species < N_SP)
    clipped_species_count[species]++;
}

void chemistry_profile_record_newton_limiter(int species, double nH)
{
  if(species < 0 || species >= N_SP)
    return;
  newton_limiter_count[species]++;
  if(nH >= 1.e12)
    high_density_newton_limiter_count[species]++;
}

void chemistry_profile_record_residual(const ChemEnvironment *environment, const ChemState *state,
                                       const ChemWorkspace *work, double temperature, double dt,
                                       double energy_old, double energy_new, double residual,
                                       int chemistry_status)
{
  ResidualTrace *trace = &traces[trace_next++ % TRACE_COUNT];
  memset(trace, 0, sizeof(*trace));
  trace->temperature = temperature;
  trace->dt = dt;
  trace->density = environment->rho_cgs / PROTONMASS;
  trace->metallicity = environment->metallicity;
  trace->energy_old = energy_old;
  trace->energy_new = energy_new;
  trace->residual = residual;
  trace->chemical = work->evaluation.rhoU_dot_chem;
  trace->radiative = work->evaluation.rhoU_dot_rad;
  trace->chem_H2 = work->diagnostics.edot_H2form_gas;
  trace->chem_Hion = work->diagnostics.edot_ciHI;
  trace->chem_Hminus = work->diagnostics.edot_HMform;
  trace->rad_H2 = work->diagnostics.edot_H2;
#ifdef HD
  trace->rad_HD = work->diagnostics.edot_HD;
#endif
  trace->rad_HI = work->diagnostics.edot_ceHI;
#ifdef METALLICITY
  trace->rad_CI = work->diagnostics.edot_CI;
  trace->rad_CII = work->diagnostics.edot_CII;
  trace->rad_OI = work->diagnostics.edot_OI;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  trace->rad_CO = work->diagnostics.edot_CO;
  trace->rad_OH = work->diagnostics.edot_OH;
  trace->rad_H2O = work->diagnostics.edot_H2O;
#endif
#ifdef DUST_COOLING
  trace->rad_dust = work->diagnostics.edot_dust;
#endif
  trace->rad_continuum = work->diagnostics.edot_cont;
  trace->rad_other = trace->radiative - trace->rad_H2 - trace->rad_HD - trace->rad_HI -
                     trace->rad_CI - trace->rad_CII - trace->rad_OI - trace->rad_CO -
                     trace->rad_OH - trace->rad_H2O - trace->rad_dust - trace->rad_continuum;
#ifdef CARBON_OXYGEN_CHEMISTRY
  trace->CI = state->abund[sCI]; trace->CII = state->abund[sCII];
  trace->CO = state->abund[sCO]; trace->CO2 = state->abund[sCO2];
  trace->OI = state->abund[sOI]; trace->OH = state->abund[sOH];
  trace->H2O = state->abund[sH2O]; trace->O2 = state->abund[sO2];
  trace->COII = state->abund[sCOII]; trace->CH = state->abund[sCH];
  trace->CH2 = state->abund[sCH2]; trace->OII = state->abund[sOII];
  trace->OHII = state->abund[sOHII]; trace->H2OII = state->abund[sH2OII];
  trace->H3OII = state->abund[sH3OII]; trace->O2II = state->abund[sO2II];
#endif
  trace->chemistry_status = chemistry_status;
}

void chemistry_profile_report_energy_warning(void)
{
  unsigned long count = trace_next < TRACE_COUNT ? trace_next : TRACE_COUNT;
  unsigned long first = trace_next - count;
  fprintf(stderr, "PROFILE_RESIDUAL_BEGIN count=%lu\n", count);
  fprintf(stderr, "# T status g eold enew chem H2chem Hionchem HMchem rad H2 HD HI CI CII OI CO OH H2O dust cont other n Z dt\n");
  for(unsigned long index = first; index < trace_next; index++) {
    ResidualTrace *t = &traces[index % TRACE_COUNT];
    fprintf(stderr, "RES %.17g %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            t->temperature, t->chemistry_status, t->residual, t->energy_old, t->energy_new,
            t->chemical, t->chem_H2, t->chem_Hion, t->chem_Hminus, t->radiative,
            t->rad_H2, t->rad_HD, t->rad_HI, t->rad_CI, t->rad_CII, t->rad_OI,
            t->rad_CO, t->rad_OH, t->rad_H2O, t->rad_dust, t->rad_continuum,
            t->rad_other, t->density, t->metallicity, t->dt);
    fprintf(stderr, "ABUND %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            t->CI, t->CII, t->CO, t->CO2, t->OI, t->OH, t->H2O, t->O2,
            t->COII, t->CH, t->CH2, t->OII, t->OHII, t->H2OII, t->H3OII, t->O2II);
  }
  fprintf(stderr, "PROFILE_RESIDUAL_END\n");
}

void chemistry_profile_report(void)
{
  static const char *names[CHEM_PROFILE_COUNT] = {
    "chemistry total", "energy solver", "energy residual", "implicit chemistry",
    "reaction rates", "intermediate solve", "intermediate assembly", "intermediate linear",
    "Newton RHS", "Newton Jacobian", "Jacobian reaction assembly", "Jacobian response",
    "Jacobian intermediate", "Jacobian Schur update",
    "Newton linear solve", "Newton step", "chemical cooling", "radiative cooling",
    "molecular cooling", "dust temperature"
  };
  double total = profile_seconds[CHEM_PROFILE_TOTAL];
  fprintf(stderr, "PROFILE_TIMING_BEGIN\n");
  for(int section = 0; section < CHEM_PROFILE_COUNT; section++)
    fprintf(stderr, "PROFILE_TIME %-22s %.9g s %8.3f %% calls=%lu\n", names[section],
            profile_seconds[section], total > 0.0 ? 100.0 * profile_seconds[section] / total : 0.0,
            profile_calls[section]);
  fprintf(stderr,
          "PROFILE_EQUILIBRIUM_INITIALIZATION calls=%lu limited=%lu min_fraction=%.17g\n",
          equilibrium_initializations, limited_equilibrium_initializations,
          minimum_equilibrium_initialization_fraction);
  fprintf(stderr,
          "PROFILE_EQUILIBRIUM_ELECTRON ratio_gt_2=%lu ratio_gt_10=%lu ratio_gt_1e3=%lu max_ratio=%.17g\n",
          equilibrium_electron_ratio_gt_2, equilibrium_electron_ratio_gt_10,
          equilibrium_electron_ratio_gt_1e3,
          maximum_equilibrium_electron_ratio);
  fprintf(stderr,
          "PROFILE_NEWTON_SOLVES calls=%lu with_reset=%lu resets=%lu with_clipping=%lu clipped_steps=%lu max_iterations=%lu\n",
          newton_solves, newton_solves_with_reset, newton_reset_count,
          newton_solves_with_clipping, newton_clipped_steps,
          newton_max_iterations);
  fprintf(stderr,
          "PROFILE_NEWTON_HIST le2=%lu le4=%lu le8=%lu le16=%lu le69=%lu ge70=%lu\n",
          newton_solve_histogram[0], newton_solve_histogram[1],
          newton_solve_histogram[2], newton_solve_histogram[3],
          newton_solve_histogram[4], newton_solve_histogram[5]);
  for(int species = 0; species < N_SP; species++)
    if(clipped_species_count[species] > 0)
      fprintf(stderr, "PROFILE_CLIPPED_SPECIES name=%s count=%lu\n",
              chemistry_species_definition(static_cast<SPECIES>(species))->name,
              clipped_species_count[species]);
  for(int species = 0; species < N_SP; species++)
    if(newton_limiter_count[species] > 0)
      fprintf(stderr,
              "PROFILE_NEWTON_LIMITER name=%s count=%lu high_density_count=%lu\n",
              chemistry_species_definition(static_cast<SPECIES>(species))->name,
              newton_limiter_count[species],
              high_density_newton_limiter_count[species]);
  fprintf(stderr, "PROFILE_TIMING_END\n");
}

#endif
