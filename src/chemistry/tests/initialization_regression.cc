#include "solver/fixed_rate_solver.hpp"
#include "chemistry_model.hpp"
#include "chemistry_core.h"
/* Internal regression: run with the same compile-time species as the core. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "chemistry_constants.h"
#include "chemistry_core.h"
#include "chemistry_types.h"
#include "chemistry_solver.h"
#include "chemistry_tables.h"
#include "forward_rates.h"
#include "species.h"
#include "thermodynamics.h"

static int check_short_step(const ChemRunView *run, double ion_floor)
{
  ChemGlobalStep global{};
  global.cmb_temperature = 2.7255;
  global.scale_factor = 1.;
  ChemEnvironment environment{};
  const ChemRunView &view = *run;
  environment.run = &view;
  environment.global = &global;
  environment.rho_cgs = PROTONMASS / XH;
  ChemState state{};
  state.abund[sHII] = 1.e-4 * XH;
  state.abund[sH2I] = 1.e-3 * XH;
  state.abund[sHM] = 1.e-8 * XH;
  state.abund[sHI] = XH - state.abund[sHII] -
      2. * state.abund[sH2I] - state.abund[sHM];
  state.abund[sELEC] = state.abund[sHII] - state.abund[sHM];
#ifdef HELIUM
  state.abund[sHeII] = ion_floor * XH;
  state.abund[sHeIII] = ion_floor * XH;
  state.abund[sHeI] = XHe - state.abund[sHeII] - state.abund[sHeIII];
  state.abund[sELEC] += state.abund[sHeII] + 2. * state.abund[sHeIII];
#endif
#ifdef H2_PLUS
  state.abund[sH2II] = ion_floor * XH;
  state.abund[sHI] -= 2. * state.abund[sH2II];
  state.abund[sELEC] += state.abund[sH2II];
#endif
#ifdef HD
  state.abund[sDI] = D_number_fraction;
#endif
  const double temperature = 300.;
  const double dt = 1.e-12;
  state.gamma = chemistry_gamma_from_abundances(state.abund, temperature);
  const double mu = chemistry_molecular_weight(state.abund);
  state.utherm_cgs = temperature * BOLTZMANN /
      ((state.gamma - 1.) * PROTONMASS * mu);
  ChemState initial = state;
  ChemWorkspace workspace{};
  workspace.thermodynamics.temperature = temperature;
  workspace.thermodynamics.molecular_weight = mu;
  workspace.evaluation.t_cool_rad = MAX_REAL_NUMBER;
  workspace.evaluation.t_cool_chem = MAX_REAL_NUMBER;
  workspace.evaluation.t_chem = MAX_REAL_NUMBER;
  ChemStepData step;
  chemistry::Driver::prepare_cell(environment, state, temperature, dt, step);
  if(chemistry::FixedRateSolver{}.advance(environment, step, state,
                                   workspace, dt, {}) != CHEMISTRY_STATUS_OK ||
     workspace.status != CHEMISTRY_STATUS_OK || workspace.evaluation.err_chem){
    fprintf(stderr, "short-step solve failed: status=%d error=%g failure=%d\n", workspace.status, workspace.evaluation.err_chem, chemistry_profile_last_failure());
    return 0;
  }
  for(int species = 0; species < N_SP; species++){
    const auto &definition = *chemistry_species_definition(static_cast<SPECIES>(species));
    if((definition.nuclei[CHEM_ELEMENT_C] || definition.nuclei[CHEM_ELEMENT_O]) &&
       state.abund[species] != 0.) {
      fprintf(stderr, "export floor introduced metals into a zero-metallicity state\n");
      return 0;
    }
    if(!isfinite(state.abund[species]) || state.abund[species] < 0. ||
       fabs(state.abund[species] - initial.abund[species]) > 1.e-11 * XH){
      fprintf(stderr, "old-state regression: species=%s initial=%.17g final=%.17g\n",
              chemistry_species_definition(static_cast<SPECIES>(species))->name,
              initial.abund[species], state.abund[species]);
      return 0;
    }
  }
  return 1;
}

int main(int argc, char **argv)
{
  ChemistryConfig config{};
  chemistry_config_defaults(&config);
  config.table_root = argc > 1 ? argv[1] : "src/chemistry";
  auto created = chemistry::Model::create(config);
  if(created.status != CHEMISTRY_STATUS_OK) return 1;
  const auto &run = created.model->run();
  int valid = 1;
#ifdef H2_PLUS
  for(int i = 0; i < N_T; i++){
    const double rate = run.rate_tables->forward[10][i];
    if(!isfinite(rate) || rate < 0.){
      fprintf(stderr, "nonfinite k10 at T=%.17g\n",run.rate_tables->temperature[i]);
      valid = 0;
      break;
    }
  }
#endif
  valid &= check_short_step(&run, 0.);
  valid &= check_short_step(&run, 1.e-110);
  valid &= check_short_step(&run, 1.e-100);
  if(valid) puts("PASS: finite k10 table and physical old-state preservation");
  return valid ? 0 : 1;
}
