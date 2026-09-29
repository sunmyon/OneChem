#include "thermodynamics/thermodynamics.h"
#include "chemistry_model.hpp"
#include "solver/fixed_rate_solver.hpp"
#include "chemistry_core.h"
#include "reaction_network/reaction_rates.h"
#include "cooling/chemistry_cooling.h"
#include "species.h"
#include "chemistry_constants.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <type_traits>

static_assert(!std::is_copy_constructible<chemistry::FixedRateProblem>::value,
              "borrowed reference addresses must remain stable");
static_assert(std::is_const<std::remove_pointer_t<decltype(ChemNewtonSystem::initial_abundance)>>::value,
              "Newton must not modify the physical old state");

#define CHECK(condition) do { \
  if(!(condition)) { \
    std::fprintf(stderr, "Fixed-rate contract failed: %s\n", #condition); \
    return 1; \
  } \
} while(0)

int main(int argc, char **argv)
{
  ChemistryConfig config{};
  chemistry_config_defaults(&config);
  config.table_root = argc > 1 ? argv[1] : "src/chemistry";
  auto created = chemistry::Model::create(config);
  CHECK(created.status == CHEMISTRY_STATUS_OK);
  ChemGlobalStep global{};
  global.cmb_temperature = 2.7255;
  global.scale_factor = 1.;
  ChemEnvironment environment{};
  environment.run = &created.model->run();
  environment.global = &global;
  environment.rho_cgs = 1.e3 * PROTONMASS / XH;
  environment.metallicity = 1.e-4;
  ChemState initial{};
  for(std::size_t i = 0; i < chemistry_element_conservation_count; ++i) {
    const auto &element = chemistry_element_conservation[i];
    initial.abund[element.reservoirs[0]] =
        element.total_per_hydrogen * XH *
        (element.metallicity_scaled ? environment.metallicity : 1.);
  }
  initial.abund[sHII] = 1.e-4 * XH;
  initial.abund[sH2I] = 1.e-3 * XH;
  initial.abund[sHM] = 1.e-20 * XH;
  initial.abund[sHI] -= initial.abund[sHII] + 2. * initial.abund[sH2I] + initial.abund[sHM];
  initial.abund[sELEC] = chemistry_net_charge(initial.abund);
  const double temperature = 1.e4, dt = 1.e4;
  initial.gamma = chemistry_gamma_from_abundances(initial.abund, temperature);
  initial.utherm_cgs = temperature * BOLTZMANN /
      ((initial.gamma - 1.) * PROTONMASS * chemistry_molecular_weight(initial.abund));
  ChemWorkspace prepared{};
  prepared.thermodynamics.temperature = temperature;
  ChemStepData step{};
  chemistry::Driver::prepare_cell(environment, initial, temperature, dt, step);
  chemistry::ReactionRateEvaluator evaluator(environment, initial, prepared);
  CHECK(evaluator.evaluate(step) == 0);
  CHECK(evaluator.evaluate_timescale() == 0);
  const auto first_rates = prepared.rates;
  const double first_timescale = prepared.evaluation.t_chem;
  const double first_energy_rate = prepared.evaluation.rhoU_dot_chem;
  auto other_environment = environment;
  other_environment.rho_cgs *= 3.;
  ChemWorkspace other_workspace{};
  other_workspace.thermodynamics.temperature = 300.;
  chemistry::ReactionRateEvaluator other(other_environment, initial, other_workspace);
  CHECK(other.evaluate(step) == 0);
  CHECK(other.evaluate_timescale() == 0);
  CHECK(std::memcmp(&first_rates, &prepared.rates, sizeof(first_rates)) == 0);
  CHECK(prepared.evaluation.t_chem == first_timescale);
  CHECK(prepared.evaluation.rhoU_dot_chem == first_energy_rate);
  // Interleaved evaluation must not retain a different cell's inputs.
  CHECK(evaluator.evaluate(step) == 0);
  CHECK(evaluator.evaluate_timescale() == 0);
  CHECK(std::memcmp(&first_rates, &prepared.rates, sizeof(first_rates)) == 0);
  CHECK(prepared.evaluation.t_chem == first_timescale);
  CHECK(prepared.evaluation.rhoU_dot_chem == first_energy_rate);

  // Cooling sessions borrow independent workspaces and share immutable Model data.
  auto cooling_work = prepared;
  auto other_cooling_work = other_workspace;
  chemistry::CoolingEvaluator cooling(environment, initial, cooling_work);
  chemistry::CoolingEvaluator other_cooling(other_environment, initial, other_cooling_work);
  CHECK(cooling.radiative(step) == CHEMISTRY_STATUS_OK);
  CHECK(cooling.chemical(0., initial.abund) == CHEMISTRY_STATUS_OK);
  const double radiative_rate = cooling_work.evaluation.rhoU_dot_rad;
  const double chemical_rate = cooling_work.evaluation.rhoU_dot_chem;
  CHECK(std::isfinite(radiative_rate) && std::isfinite(chemical_rate));
  CHECK(other_cooling.radiative(step) == CHEMISTRY_STATUS_OK);
  CHECK(other_cooling.chemical(0., initial.abund) == CHEMISTRY_STATUS_OK);
  CHECK(cooling_work.evaluation.rhoU_dot_rad == radiative_rate);
  CHECK(cooling_work.evaluation.rhoU_dot_chem == chemical_rate);
  CHECK(cooling.radiative(step) == CHEMISTRY_STATUS_OK);
  CHECK(cooling.chemical(0., initial.abund) == CHEMISTRY_STATUS_OK);
  CHECK(cooling_work.evaluation.rhoU_dot_rad == radiative_rate);
  CHECK(cooling_work.evaluation.rhoU_dot_chem == chemical_rate);
  CHECK(cooling.chemical(dt, initial.abund) == CHEMISTRY_STATUS_OK);
  CHECK(std::isfinite(cooling_work.evaluation.rhoU_dot_chem));
  other_cooling_work.thermodynamics.temperature = 1.e7;
  CHECK(other_cooling.high_temperature(step) == CHEMISTRY_STATUS_OK);
  CHECK(std::isfinite(other_cooling_work.evaluation.rhoU_dot_rad));
  auto invalid_environment = environment;
  invalid_environment.run = nullptr;
  CHECK(chemistry::CoolingEvaluator(invalid_environment, initial, cooling_work).radiative(step)
        == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  const chemistry::FixedRateProblem problem(
      environment, step, initial, prepared.rates, temperature, dt);
  CHECK(problem.valid());
  const auto reference = problem.reference();
  const auto rates = problem.rates();
  const auto physical_initial = initial;
  // A caller can reuse its inputs without changing the frozen problem.
  initial = {};
  step = {};
  global.cmb_temperature = 999.;
  environment.rho_cgs *= 2.;
  prepared.rates = {};
  CHECK(problem.environment().global->cmb_temperature == 2.7255);
  CHECK(problem.environment().rho_cgs != environment.rho_cgs);
  CHECK(std::memcmp(&rates, &problem.rates(), sizeof(rates)) == 0);

  chemistry::NewtonScratch scratch;
  chemistry::FixedRateSolver solver;
  chemistry::FixedRateOptions options;
  options.report_failure = 0;
  ChemState first{}, second{};
  ChemWorkspace work1{}, work2{};
  CHECK(solver.solve(problem, options, first, work1, scratch) == CHEMISTRY_STATUS_OK);
  CHECK(solver.solve(problem, options, second, work2, scratch) == CHEMISTRY_STATUS_OK);
  CHECK(std::memcmp(first.abund, second.abund, sizeof(first.abund)) == 0);
  CHECK(first.gamma == physical_initial.gamma);
  CHECK(first.utherm_cgs == physical_initial.utherm_cgs);
  CHECK(work1.evaluation.rhoU_dot_chem == work2.evaluation.rhoU_dot_chem);
  CHECK(std::memcmp(&reference, &problem.reference(), sizeof(reference)) == 0);

  options.max_iterations = 0;
  options.allow_restarts = 0;
  CHECK(solver.solve(problem, options, second, work2, scratch) == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  CHECK(std::memcmp(second.abund, physical_initial.abund, sizeof(second.abund)) == 0);
  CHECK(std::memcmp(&reference, &problem.reference(), sizeof(reference)) == 0);
  options.max_iterations = 3000;
  options.allow_restarts = 1;
  CHECK(solver.solve(problem, options, second, work2, scratch) == CHEMISTRY_STATUS_OK);
  CHECK(std::memcmp(first.abund, second.abund, sizeof(first.abund)) == 0);
  auto invalid_initial = physical_initial;
  invalid_initial.abund[sHI] = -1.;
  const chemistry::FixedRateProblem invalid(
      problem.environment(), problem.step(), invalid_initial,
      problem.rates(), temperature, dt);
  CHECK(!invalid.valid());
  CHECK(solver.solve(invalid, options, second, work2, scratch) == CHEMISTRY_STATUS_INCOMPATIBLE_STATE);
  CHECK(std::memcmp(second.abund, invalid_initial.abund, sizeof(second.abund)) == 0);
  std::puts("PASS: fixed-rate snapshots, independent sessions and failure recovery");
  return 0;
}
