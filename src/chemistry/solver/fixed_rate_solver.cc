#include "fixed_rate_solver.hpp"
#include "chemistry_newton_internal.h"
#include "chemistry_newton_control.h"
#include "reaction_network/reaction_network.h"
#include "chemistry_constants.h"
#include "solver/linear_algebra.hpp"
#include "chemistry_newton_reset.h"
#include "chemistry_newton_diagnostics.h"
#include "chemistry_core.h"
#include "reaction_rates.h"
#include <cmath>
#include <cstdio>
#include <cstring>


namespace chemistry {

FixedRateProblem::FixedRateProblem(
    const ChemEnvironment &environment, const ChemStepData &step,
    const ChemState &initial, const ChemRates &rates,
    double temperature, double duration)
    : run_(*environment.run), global_(*environment.global),
      environment_(environment), step_(step), initial_(initial), rates_(rates),
      temperature_(temperature), duration_(duration),
      nuclei_density_(environment.rho_cgs / PROTONMASS),
      hydrogen_density_(XH * nuclei_density_) {
  environment_.run = &run_;
  environment_.global = &global_;
  valid_ = seed_.prepare(
      &environment_, &initial_, &rates_, hydrogen_density_, reference_);
}

void FixedRateProblem::evaluate_rhs(const double *abundance, double *rhs,
                                    double *production, double *loss) const {
  (rates_.network)->rhs_gross( CHEM_REACTION_IMPLICIT,
                              hydrogen_density_, abundance, rates_.forward,
                              rates_.inverse, rhs, production, loss);
}

double FixedRateProblem::evaluate_residual(
    const double *abundance, const double *rhs, const double *production,
    const double *loss, const ChemConservationSystem &conservation) const {
  return chemistry_newton_evaluate_residual(
      duration_, abundance, reference_.initial_abundance, reference_.upper_bound,
      rhs, production, loss, &rates_, &conservation);
}

bool FixedRateProblem::build_jacobian(
    const ChemNewtonSystem &system, const double *rhs, int rebuild,
    int &force_numerical, double jacobian[N_SP][N_SP]) const {
  return system.build_jacobian(
      hydrogen_density_, rhs, &rates_, rebuild,
      &force_numerical, jacobian);
}

} // namespace chemistry

// Newton iteration and recovery for the fixed-rate problem.

#define CHEMISTRY_NEWTON_TRACE_ABUNDANCE 1.e-12

namespace chemistry {
typedef enum {
  NEWTON_RECOVERY_NONE,
  NEWTON_RECOVERY_INVALID_BASIN,
  NEWTON_RECOVERY_INVALID_STEP
} NewtonRecoveryReason;

typedef enum {
  NEWTON_BASIN_RESTART,
  NEWTON_BASIN_CONTINUE,
  NEWTON_BASIN_FAILED
} NewtonBasinOutcome;

typedef enum {
  NEWTON_STEP_CONTINUE,
  NEWTON_STEP_RECOVER_BASIN
} NewtonStepOutcome;



enum {
  NEWTON_REQUIRED_CONFIRMATIONS = 2,
  NEWTON_CLIPPING_HOLD_ITERATIONS = 2,
  NEWTON_BOUNDARY_STEPS_BEFORE_RESELECTION = 2,
  NEWTON_BASIN_ITERATION_LIMIT = 70
};

typedef struct {
  int iteration;
  int basin_start;

  int convergence_count;
  int clipping_hold;
  int boundary_streak;
  int chart_reselection_requested;

  int lu_valid;
  int force_numerical_jacobian;

  ChemNewtonResetIndex reset;
  int deterministic_reset_count;
  int random_restart_count;
  int clipped_step_count;
  int max_iterations;
  int allow_restarts;

  uint64_t random_state;
  double step_error;
  NewtonRecoveryReason recovery;
} NewtonControl;

typedef struct {
  double residual_error;
  ChemNewtonStepResult result;
  double correction[N_SP];
  int limiting_species;
} NewtonIterationData;

static void reset_for_newton_basin(NewtonControl *newton){
  newton->basin_start = newton->iteration;
  newton->convergence_count = 0;
  newton->clipping_hold = 0;
  newton->boundary_streak = 0;
  newton->chart_reselection_requested = 0;
  newton->lu_valid = 0;
  newton->force_numerical_jacobian = 0;
  newton->step_error = 100.0;
  newton->recovery = NEWTON_RECOVERY_NONE;
}

static void reset_after_chart_change(NewtonControl *newton){
  newton->convergence_count = 0;
  newton->boundary_streak = 0;
  newton->chart_reselection_requested = 0;
  newton->lu_valid = 0;
}

static void reset_for_newton_recovery(NewtonControl *newton){
  newton->convergence_count = 0;
  newton->step_error = 100.0;
}

static int chemistry_step_has_relevant_clip(
    const double *abundance, const double *newton_step,
    const int *clipped_species, double convergence_abundance_threshold,
    unsigned char *relevant_clipped_species){
  int relevant_step_was_clipped = 0;

  memset(relevant_clipped_species, 0,
         N_SP * sizeof(*relevant_clipped_species));

  for(int species = 0; species < N_SP; species++){
    const double abundance_before_step =
        abundance[species] - newton_step[species];
    int relevant_clip =
        clipped_species[species] &&
        dmax(fabs(abundance_before_step), fabs(abundance[species])) >=
            convergence_abundance_threshold;

#if defined(CHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION) &&                  \
    !defined(CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE)
    if(chemistry_is_full_implicit_initialized_intermediate(static_cast<SPECIES>(species)))
      relevant_clip = 0;
#endif

    if(relevant_clip){
      relevant_step_was_clipped = 1;
      relevant_clipped_species[species] = 1;
    }
  }

  return relevant_step_was_clipped;
}

class NewtonSession {
public:
  NewtonSession(const FixedRateProblem &problem, const FixedRateOptions &options)
      : problem_(problem), options_(options), system_(problem.seed()) {
    control_.iteration = 1;
    control_.basin_start = 1;
    control_.reset = CHEM_NEWTON_RESET_NONE;
    control_.deterministic_reset_count = CHEM_NEWTON_RESET_LAST;
    control_.max_iterations = options.max_iterations;
    control_.allow_restarts = options.allow_restarts;
    control_.random_state = chemistry_newton_random_seed(
        problem.hydrogen_density(), problem.temperature(), problem.duration());
    control_.step_error = 100.0;
    control_.recovery = NEWTON_RECOVERY_NONE;
  }
  NewtonSession(const NewtonSession &) = delete;
  NewtonSession &operator=(const NewtonSession &) = delete;

private:
  friend class FixedRateSolver;
  int run(ChemState *state, ChemWorkspace *work, NewtonScratch &scratch);
  bool prepare_basin();
  bool reselect_chart() {
    if(!chemistry_newton_reselect_conservation_chart(
        problem_.environment().metallicity, &system_))
      return false;
    reset_after_chart_change(&control_);
    return true;
  }
  NewtonStepOutcome resolve_invalid_step();
  NewtonStepOutcome finish_step(NewtonIterationData *iteration_data);
  NewtonBasinOutcome recover_basin();

  const FixedRateProblem &problem_;
  const FixedRateOptions options_;
  ChemNewtonSystem system_;
  NewtonControl control_{};
};

bool NewtonSession::prepare_basin(){
  const double nH = problem_.hydrogen_density();
  const double metallicity = problem_.environment().metallicity;
  const ChemRates *rates = &problem_.rates();
  auto *system = &system_;
  auto *newton = &control_;
  reset_for_newton_basin(newton);

  if(newton->reset != CHEM_NEWTON_RESET_NONE){
    chemistry_newton_apply_deterministic_reset(
        newton->reset, metallicity, system->initial_abundance,
        system->abundance);

    reselect_chart();

    chemistry_correct_full_implicit_intermediate_guess(
        nH, system->abundance, rates, &system->conservation);

    for(int species = 0; species < N_SP; species++){
      system->total_step[species] =
          system->abundance[species] - system->initial_abundance[species];
    }
  }

  int basin_is_valid = chemistry_newton_normalize_current_state(
      system->abundance, system->upper_bound);

  if(basin_is_valid){
    basin_is_valid = (&system->conservation)->update_algebraic(nH, system->abundance, rates);
  }
  if(basin_is_valid)
    basin_is_valid = chemistry_abundances_are_nonnegative(system->abundance);

  if(!basin_is_valid){
    newton->recovery = NEWTON_RECOVERY_INVALID_BASIN;
    return 0;
  }

  for(int species = 0; species < N_SP; species++){
    if((&system->conservation)->is_algebraic(static_cast<SPECIES>(species)))
      system->total_step[species] = 0.0;
  }

  return 1;
}

NewtonStepOutcome NewtonSession::resolve_invalid_step(){
  auto *newton = &control_;
  if(reselect_chart()){
    newton->iteration++;

    return NEWTON_STEP_CONTINUE;
  }

  newton->recovery = NEWTON_RECOVERY_INVALID_STEP;
  return NEWTON_STEP_RECOVER_BASIN;
}

NewtonStepOutcome NewtonSession::finish_step(NewtonIterationData *iteration_data){
  const double nH = problem_.hydrogen_density();
  auto *system = &system_;
  auto *newton = &control_;
  const int step_was_clipped = iteration_data->result.step_was_clipped;

  if(step_was_clipped)
    newton->lu_valid = 0;

  newton->step_error = chemistry_newton_material_step_norm(
      system->abundance, iteration_data->correction, &system->conservation,
      &iteration_data->limiting_species);

  for(int species = 0; species < N_SP; species++){
    system->total_step[species] =
        system->abundance[species] - system->initial_abundance[species];

    if((&system->conservation)->is_algebraic(static_cast<SPECIES>(species)))
      system->total_step[species] = 0.0;
  }

  unsigned char relevant_clipped_species[N_SP];
  int relevant_step_was_clipped = chemistry_step_has_relevant_clip(
      system->abundance, iteration_data->correction,
      iteration_data->result.clipped_species,
      CHEMISTRY_NEWTON_TRACE_ABUNDANCE, relevant_clipped_species);

  if(relevant_step_was_clipped){
    newton->clipping_hold = NEWTON_CLIPPING_HOLD_ITERATIONS;
    newton->clipped_step_count++;
  }

  newton->boundary_streak =
      step_was_clipped ? newton->boundary_streak + 1 : 0;

  if(newton->boundary_streak >= NEWTON_BOUNDARY_STEPS_BEFORE_RESELECTION)
    newton->chart_reselection_requested = 1;

  if(iteration_data->residual_error > CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT ||
      relevant_step_was_clipped){
    newton->convergence_count = 0;
  }else if(newton->clipping_hold > 0){
    newton->clipping_hold--;
    newton->convergence_count = 0;
  }else if(newton->step_error <= 1.0){
    newton->convergence_count++;
  }else{
    newton->convergence_count = 0;
  }

  /* Profiling only; these calls do not change the Newton state. */
  if(newton->step_error > 0.0)
    chemistry_profile_record_newton_limiter(
        iteration_data->limiting_species, nH);

  for(int species = 0; species < N_SP; species++){
    if(relevant_clipped_species[species])
      chemistry_profile_record_clipped_species(species);
  }

  if(newton->allow_restarts && newton->convergence_count == 0 &&
     newton->iteration - newton->basin_start + 1 >=
         NEWTON_BASIN_ITERATION_LIMIT &&
     newton->reset <= newton->deterministic_reset_count){
    newton->recovery = NEWTON_RECOVERY_NONE;
    return NEWTON_STEP_RECOVER_BASIN;
  }

  newton->iteration++;
  return NEWTON_STEP_CONTINUE;
}

NewtonBasinOutcome NewtonSession::recover_basin(){
  const double nH = problem_.hydrogen_density();
  const double metallicity = problem_.environment().metallicity;
  const ChemRates *rates = &problem_.rates();
  auto *system = &system_;
  auto *newton = &control_;
  reset_for_newton_recovery(newton);

  if(newton->allow_restarts &&
     newton->reset < newton->deterministic_reset_count){
    newton->reset = static_cast<ChemNewtonResetIndex>(newton->reset + 1);
    return NEWTON_BASIN_RESTART;
  }

  if(newton->recovery == NEWTON_RECOVERY_INVALID_STEP)
    chemistry_profile_set_failure(6);

  if(newton->allow_restarts && newton->iteration < newton->max_iterations){
    int found_random_guess = chemistry_find_random_newton_guess(
        nH, metallicity, system->initial_abundance, system->abundance,
        rates, &system->conservation, &newton->random_state,
        &newton->iteration, newton->max_iterations,
        &newton->random_restart_count);

    if(found_random_guess){
      newton->reset = static_cast<ChemNewtonResetIndex>(newton->reset + 1);
      newton->force_numerical_jacobian = 0;
      return NEWTON_BASIN_RESTART;
    }
  }

  if(newton->recovery != NEWTON_RECOVERY_NONE){
    if(newton->recovery == NEWTON_RECOVERY_INVALID_BASIN)
      chemistry_profile_set_failure(3);

    return NEWTON_BASIN_FAILED;
  }

  newton->iteration++;
  return NEWTON_BASIN_CONTINUE;
}

int NewtonSession::run(ChemState *state, ChemWorkspace *work,
                       NewtonScratch &scratch){
  const auto *environment = &problem_.environment();
  const auto *step = &problem_.step();
  const auto &state_old = problem_.initial_state();
  const double dt = problem_.duration();
  const double n_nuclei = problem_.nuclei_density();
  const double nH = problem_.hydrogen_density();
  const double T_gas = problem_.temperature();
  const int max_iterations = options_.max_iterations;
  const int allow_restarts = options_.allow_restarts;
  const int report_failure = options_.report_failure;
  auto &system = system_;
  auto &newton = control_;
  *state = state_old;
  work->rates = problem_.rates();
  work->thermodynamics.temperature = T_gas;
  work->evaluation.err_chem = 1.0;
  if(!problem_.valid() || !system.initialize(
      nH, dt, options_.initial_guess, &problem_.rates())){
    *state = state_old;
    return 0;
  }

  /* Abundance state and physical bounds. */
  double *abundance = system.abundance;
  const double *initial_abundance = system.initial_abundance;
  const double *cooling_initial_abundance = system.cooling_initial_abundance;

  const double *initialization_floor = system.initialization_floor;
  const double *lower_bound = system.newton_lower_bound;
  const double *upper_bound = system.upper_bound;

  /* Reduced-system layout. */
  double *total_step = system.total_step;
  enum SPECIES *solver_species = system.solver_species;
  ChemConservationSystem *conservation = &system.conservation;

  auto &lapack_pivots = scratch.lapack_pivots;
  auto &rhs_jacobian = scratch.rhs_jacobian;
  auto &reaction_rhs = scratch.reaction_rhs;
  auto &species_production = scratch.species_production;
  auto &species_loss = scratch.species_loss;
  auto &newton_matrix = scratch.newton_matrix;
  auto &linear_rhs = scratch.linear_rhs;
  NewtonIterationData iteration_data{};

  ChemNewtonDiagnostics diagnostics;
  NewtonBasinOutcome basin_outcome;

restart_newton_basin:
  basin_outcome = NEWTON_BASIN_CONTINUE;

  if(!prepare_basin()){
    basin_outcome = recover_basin();
  }

  if(basin_outcome == NEWTON_BASIN_CONTINUE){
    /* Diagnostics only: print and snapshot the accepted basin entry. */
    print_newton_basin_entry(work, T_gas, abundance);

#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
    memcpy(diagnostics.basin_entry, abundance, sizeof(diagnostics.basin_entry));
#endif
  }

  /* Main Newton loop.
   * Each pass evaluates the residual, solves the reduced linear system,
   * applies a bounded correction, and updates convergence state. */
  while(basin_outcome == NEWTON_BASIN_CONTINUE &&
        newton.iteration <= max_iterations &&
        newton.convergence_count < NEWTON_REQUIRED_CONFIRMATIONS){
    /* Re-select a conservation reservoir after repeated boundary contact. */
    if(newton.chart_reselection_requested){
      reselect_chart();
    }

    /* Reaction RHS and residual. */
    double rhs_timer = profile_start();

    problem_.evaluate_rhs(abundance, reaction_rhs, species_production, species_loss);

    profile_end(CHEM_PROFILE_NEWTON_RHS, rhs_timer);

    iteration_data.residual_error = problem_.evaluate_residual(abundance, reaction_rhs, species_production,
                                   species_loss, *conservation);

    /* Reduced Jacobian and linear solve. */
    int rebuild_newton_matrix = 1;

#ifdef CHEMISTRY_NEWTON_LU_REUSE
    /* Reuse is limited to the mandatory second convergence confirmation after
     * an already tiny correction. This guard preserves the thermochemical
     * branch selected by a freshly rebuilt Jacobian. */
    rebuild_newton_matrix =
        !newton.lu_valid || newton.convergence_count != 1 ||
        newton.step_error >= 1.0e-12;
#endif

    if(!problem_.build_jacobian(system, reaction_rhs, rebuild_newton_matrix,
                              newton.force_numerical_jacobian, rhs_jacobian)){
      chemistry_profile_set_failure(5);
      *state = state_old;

      return 0;
    }

    /* Diagnostics only. */
    print_newton_iteration(work, newton.iteration, T_gas, dt,
                           newton.step_error,
                           iteration_data.limiting_species, abundance,
                           initial_abundance, reaction_rhs);
    print_newton_jacobian(work, newton.iteration, rhs_jacobian);

    /* Backward Euler Newton system:
     *   (I - dt dR/dy) delta_y = dt R(y) - (y - y_initial). */
    if(rebuild_newton_matrix)
      numerics::assemble_backward_euler_matrix(
          system.solver_species_count, dt,
          [&](int row, int column) {
            return rhs_jacobian[solver_species[column]][solver_species[row]];
          }, newton_matrix);

    numerics::assemble_backward_euler_rhs(
        system.solver_species_count, dt,
        [&](int row) { return reaction_rhs[solver_species[row]]; },
        [&](int row) { return total_step[solver_species[row]]; }, linear_rhs);

    double linear_solve_timer = profile_start();
    int linear_solver_status = 1;
#ifdef CHEMISTRY_NEWTON_LU_REUSE
    if(rebuild_newton_matrix){
      linear_solver_status = numerics::LapackDense::factorize(
          system.solver_species_count, newton_matrix, lapack_pivots);
      newton.lu_valid = linear_solver_status == 0;
    }
    if(newton.lu_valid)
      linear_solver_status = numerics::LapackDense::solve_factored(
          system.solver_species_count, newton_matrix, lapack_pivots, linear_rhs);
#else
    linear_solver_status = numerics::LapackDense::solve(
        system.solver_species_count, newton_matrix, lapack_pivots, linear_rhs);
#endif
    profile_end(CHEM_PROFILE_NEWTON_LINEAR_SOLVE, linear_solve_timer);

    memset(iteration_data.correction, 0,
           sizeof(iteration_data.correction));

    for(int row = 0; row < system.solver_species_count; row++)
      iteration_data.correction[solver_species[row]] = linear_rhs[row];

    /* Apply the bounded correction. */
    newton.step_error = 0.0;

    double newton_step_timer = profile_start();

    int step_is_valid =
        linear_solver_status == 0 &&
        chemistry_newton_apply_step(
            nH, abundance, iteration_data.correction, lower_bound, upper_bound,
            &problem_.rates(), conservation, &iteration_data.result);

    profile_end(CHEM_PROFILE_NEWTON_STEP, newton_step_timer);

    NewtonStepOutcome step_outcome;

    if(step_is_valid){
      step_outcome = finish_step(&iteration_data);
    }else{
      step_outcome = resolve_invalid_step();
    }

    if(step_outcome == NEWTON_STEP_RECOVER_BASIN){
      basin_outcome = recover_basin();
    }
  }

  if(basin_outcome == NEWTON_BASIN_RESTART)
    goto restart_newton_basin;

  if(basin_outcome == NEWTON_BASIN_FAILED){
    *state = state_old;
    return 0;
  }

  /* Phase 3: certify the final iterate and report a failed solve. */
  if(!(conservation)->update_algebraic(nH, abundance, &problem_.rates())){
    chemistry_profile_set_failure(7);
    *state = state_old;

    return 0;
  }

  int converged =
      newton.convergence_count >= NEWTON_REQUIRED_CONFIRMATIONS;

  if(converged){
    chemistry_profile_record_newton_solve(
        newton.iteration - 1, newton.reset, newton.clipped_step_count);
  }

  print_random_restart_success(
      &diagnostics, converged, newton.iteration, newton.reset,
      newton.random_restart_count);

  if(!converged){
    work->evaluation.err_chem = 1.0;

    /* Always emit a reproducible replay line on chemistry non-convergence,
     * even when report_failure is false (the energy-loop trial entry point
     * currently sets report_failure=0 to keep its verbose species dump
     * silent for expected trial-temperature probes, but a rejected
     * temperature should still be reproducible/inspectable on demand). */
    printf("CHEM_NOT_CONVERGED id=%d T_gas=%g dt=%g n_nuclei=%g err=%g"
           " iteration=%d reset=%d random_restarts=%d allow_restarts=%d"
           " report_failure=%d\n",
           environment->ID, T_gas, dt, n_nuclei, newton.step_error,
           newton.iteration, newton.reset, newton.random_restart_count,
           allow_restarts, report_failure);
    chemistry_print_rerun_line(environment, n_nuclei, T_gas, dt, &state_old);

    if(allow_restarts && report_failure){
      chemistry_newton_diagnostics_report_failure(
          environment, work, &state_old, dt, n_nuclei, nH, T_gas,
          newton.step_error, iteration_data.limiting_species,
          newton.iteration, newton.reset, newton.random_restart_count,
          allow_restarts, abundance, initial_abundance,
          iteration_data.correction, CHEMISTRY_NEWTON_TRACE_ABUNDANCE,
          conservation);
    }

    if(allow_restarts)
      chemistry_profile_set_failure(8);

    *state = state_old;

    return 0;
  }

  system.update_timescale(work, dt);

  /* Phase 4: export the converged state and evaluate chemical cooling. */
  work->evaluation.err_chem = 0.0;

  return system.export_solution(
      environment, state, work, &state_old, dt, nH);
}

ChemistryStatus FixedRateSolver::solve(const FixedRateProblem &problem,
                           const FixedRateOptions &options,
                           ChemState &state, ChemWorkspace &work,
                           NewtonScratch &scratch) const {
  NewtonSession session(problem, options);
  state = problem.initial_state();
  if(!problem.valid()) return CHEMISTRY_STATUS_INCOMPATIBLE_STATE;
  if(options.max_iterations <= 0) return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  if(session.run(&state, &work, scratch)) return CHEMISTRY_STATUS_OK;
  return work.status != CHEMISTRY_STATUS_OK ? work.status : CHEMISTRY_STATUS_NOT_CONVERGED;
}

} // namespace chemistry

ChemistryStatus chemistry::FixedRateSolver::advance(const ChemEnvironment &environment,
    const ChemStepData &step, ChemState &state, ChemWorkspace &work, double dt,
    const FixedRateOptions &options) const {
  chemistry_profile_set_failure(0);
  work.evaluation.err_chem = 1.0;
  ReactionRateEvaluator(environment, state, work).evaluate(step);
  const FixedRateProblem problem(environment, step, state, work.rates,
                                 work.thermodynamics.temperature, dt);
  NewtonScratch scratch;
  return solve(problem, options, state, work, scratch);
}
