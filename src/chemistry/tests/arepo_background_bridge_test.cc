#include <math.h>
#include <string.h>

extern "C" {
#include "allvars.h"
}
#include "chemistry_arepo.h"
#include "chemistry_model.hpp"
#include <assert.h>
#include <stdlib.h>

struct global_data_all_processes All;

static int advance_calls;
struct sph_particle_data particle_state;
struct sph_particle_data *SphP = &particle_state;
struct particle_data particle;
struct particle_data *P = &particle;
static int active_index;
struct mock_timebins TimeBinsHydro = {1, &active_index};
double CPU_Step[5];
int ThisTask, NTask = 1, PTask;
int *Send_count, *Recv_count, *Send_offset, *Recv_offset;
double measure_time(void) { return 0.; }
double second(void) { return 0.; }
double timediff(double a, double b) { return b-a; }
void mpi_printf(const char *format, ...) { (void)format; }
double chemistry_temperature(const double *a, double g, double u) {
  (void)a; (void)g; return u;
}
double chemistry_molecular_weight(const double *a) { (void)a; return 1.; }
extern "C" void cooling_only(void);
#ifdef CHEM_BALANCE
size_t FreeBytes(void) { return 1024*1024; }
void *mymalloc(const char *label, size_t size) { (void)label; return malloc(size); }
void myfree(void *p) { free(p); }
#include "../adapters/arepo/chem_balance.cc"
#endif

void chemistry_config_defaults(ChemistryConfig *config) {
  memset(config, 0, sizeof(*config));
  config->abi_version = CHEMISTRY_ABI_VERSION;
}
namespace chemistry {
Model::Creation Model::create(const ChemistryConfig &config) noexcept {
  assert(config.timescale_safety_factor == All.Epsilon);
  return {Model::Handle(new Model), CHEMISTRY_STATUS_OK};
}
Model::~Model() noexcept = default;
RateTables::RateTables() noexcept = default;
RateTables::~RateTables() noexcept = default;
void RateTables::Deleter::operator()(ChemistryRateTables *) const noexcept {}
void ModelResourceDeleter::operator()(ChemReactionNetwork *) const noexcept {}
void ModelResourceDeleter::operator()(ChemistryCoolingTables *) const noexcept {}
void ModelResourceDeleter::operator()(ChemistryDustTables *) const noexcept {}
void ModelResourceDeleter::operator()(ChemistryContinuumOpacityTables *) const noexcept {}
}
ChemistryStatus chemistry::Driver::advance(const Model &model, const ChemistryStepContext &step_input,
    ChemWorkspace &workspace, const ChemistryCellEnvironment &cell_input, double dt,
    const ChemistryAdvanceOptions &options_input, ChemistryState &state_input,
    ChemistryResult &result_output) const {
  (void)model; (void)workspace;
  const auto *step = &step_input;
  const auto *cell = &cell_input;
  const auto *options = &options_input;
  auto *state = &state_input;
  auto *result = &result_output;

  assert(step && options);
  assert(cell->density_cgs == 160. && dt == 3.);
  assert(state->specific_thermal_energy_cgs == 15.);
  assert(state->abundances[sHI] == 0.76);
  assert(state->abundance_count == TGCHEM_NUM_ABUNDANCES);
  memset(result, 0, sizeof(*result));
  if(advance_calls++ == 0) {
    assert(!cell->diagnostics_enabled);
    return CHEMISTRY_STATUS_NOT_CONVERGED;
  }
  assert(cell->diagnostics_enabled);
  state->specific_thermal_energy_cgs = 21.;
  state->abundances[sHI] = 0.7;
  state->gamma = 1.4;
  result->chemistry_time_s = 123.;
  result->cooling_time_s = 456.;
  return CHEMISTRY_STATUS_OK;
}

ChemistryStatus chemistry::Driver::prepare(const ChemistryGlobalEnvironment &global,
                                           ChemistryStepContext &step) {
  step = {};
  step.global = global;
  return CHEMISTRY_STATUS_OK;
}

static int nearly_equal(double actual, double expected)
{
  return fabs(actual - expected) <= 1.e-13 * fabs(expected);
}

int main(void)
{
  ChemistryStepContext step;

  All.Epsilon = 0.01;
  All.J21 = 1.0;
  All.T_CMB = 30.0;
  All.cf_atime = 0.1;
  All.Time = 42.0;

  if(!chemistry_arepo_host_init() ||
     !chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_H2_diss, 1.39e-12) ||
     !nearly_equal(step.global.background_rate_HD_diss, 1.55e-12))
    return 1;

  if(!chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_HD_diss, 1.55e-12))
    return 1;

  All.J21 = 2.0;
  if(!chemistry_arepo_host_prepare_step(&step) ||
     !nearly_equal(step.global.background_rate_H2_diss, 2.78e-12) ||
     !nearly_equal(step.global.background_rate_HD_diss, 3.10e-12))
    return 1;

  All.UnitDensity_in_cgs = 2.; All.HubbleParam = 2.; All.cf_a3inv = 4.;
  All.UnitEnergy_in_cgs = 6.; All.UnitMass_in_g = 2.;
  All.cf_atime = 0.5; All.cf_time_hubble_a = 0.5;
  All.Timebase_interval = 1.; All.UnitTime_in_s = 3.;
  particle.Mass = 4.; particle.ID = 123; particle.TimeBinHydro = 1;
  particle_state.Density = 5.; particle_state.Utherm = 5.;
  particle_state.Energy = 10.; particle_state.Gamma = 1.6;
  particle_state.Abund[sHI] = 0.76;
#ifdef CHEM_BALANCE
  struct chemcomm_in input = {};
  struct chemcomm_out output = {};
  ChemDataGet = &input;
  ChemDataResult = &output;
  particle2in_chem(&input, 0);
  assert(chemistry_arepo_host_prepare_step(&step));
  chem_evaluate(0, chemistry_arepo_host_context(), &step);
  out2particle_chem(&output, 0);
#else
  cooling_only();
#endif
  assert(advance_calls == 2);
  assert(particle_state.Utherm == 7. && particle_state.Energy == 12.);
  assert(particle_state.Abund[sHI] == 0.7 && particle_state.AbundMass[sHI] == 2.8);
  assert(particle_state.Gamma == 1.4);
  assert(particle_state.t_chem == 123. && particle_state.t_cool == 456.);
  chemistry_arepo_host_destroy();
  assert(chemistry_arepo_host_context() == NULL);
  return 0;
}
