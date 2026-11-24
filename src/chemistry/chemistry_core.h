#pragma once

#include "chemistry_build_config.h"
#include "chemistry_constants.h"
#include "species.h"

#include "chemistry_types.h"
#include "include/chemistry/chemistry.h"

#ifdef CHEMISTRY_PROFILE

double profile_start(void);
void profile_end(ChemistryProfileSection section, double started_at);
void chemistry_profile_set_failure(int reason);
int chemistry_profile_last_failure(void);
void chemistry_profile_record_residual(const ChemEnvironment *environment,
                                       const ChemState *state,
                                       const ChemWorkspace *work,
                                       double temperature, double dt,
                                       double energy_old, double energy_new,
                                       double residual, int chemistry_status);
void chemistry_profile_report_energy_warning(void);
void chemistry_profile_record_equilibrium_initialization(double fraction,
                                                         double electron_before,
                                                         double electron_after);
void chemistry_profile_record_newton_solve(int iterations, int reset_count,
                                           int clipped_steps);
void chemistry_profile_record_clipped_species(int species);
void chemistry_profile_record_newton_limiter(int species, double nH);
void chemistry_profile_report(void);

#else

static inline double profile_start(void) { return 0.0; }

static inline void profile_end(ChemistryProfileSection section,
                               double started_at) {
  (void)section;
  (void)started_at;
}

static inline void chemistry_profile_set_failure(int reason) { (void)reason; }

static inline int chemistry_profile_last_failure(void) { return 0; }

static inline void chemistry_profile_record_residual(
    const ChemEnvironment *environment, const ChemState *state,
    const ChemWorkspace *work, double temperature, double dt, double energy_old,
    double energy_new, double residual, int chemistry_status) {
  (void)environment;
  (void)state;
  (void)work;
  (void)temperature;
  (void)dt;
  (void)energy_old;
  (void)energy_new;
  (void)residual;
  (void)chemistry_status;
}

static inline void chemistry_profile_report_energy_warning(void) {}

static inline void chemistry_profile_record_equilibrium_initialization(
    double fraction, double electron_before, double electron_after) {
  (void)fraction;
  (void)electron_before;
  (void)electron_after;
}

static inline void chemistry_profile_record_newton_solve(int iterations,
                                                         int reset_count,
                                                         int clipped_steps) {
  (void)iterations;
  (void)reset_count;
  (void)clipped_steps;
}

static inline void chemistry_profile_record_clipped_species(int species) {
  (void)species;
}

static inline void chemistry_profile_record_newton_limiter(int species,
                                                           double nH) {
  (void)species;
  (void)nH;
}

static inline void chemistry_profile_report(void) {}

#endif

#ifdef __cplusplus
namespace chemistry {
class Model;
// Independent orchestration: the caller owns Model, state and workspace.
class Driver final {
public:
  ChemistryStatus advance(const ChemEnvironment &environment, double duration_s,
                          ChemState &state, ChemWorkspace &workspace, int mode = 0) const;
  ChemistryStatus advance(const Model &model, const ChemistryStepContext &step,
                          ChemWorkspace &workspace, const ChemistryCellEnvironment &cell,
                          double duration_s, const ChemistryAdvanceOptions &options,
                          ChemistryState &state, ChemistryResult &result) const;
  static ChemistryStatus prepare(const ChemistryGlobalEnvironment &global,
                                 ChemistryStepContext &step);
  static void prepare_cell(const ChemEnvironment &environment, const ChemState &state,
                           double temperature, double duration_s, ChemStepData &step);
private:
  static void set_high_temperature_ionization_state(ChemState &state);
  static void evaluate_high_temperature_cooling(const ChemEnvironment &environment,
      const ChemStepData &step, const ChemState &state, double temperature, ChemWorkspace &workspace);
  static int evolve_high_temperature(const ChemEnvironment &environment, const ChemStepData &step,
      ChemState &state, double &remaining_time_s, ChemWorkspace &workspace);
  static ChemSolverPlan select_solvers(double n_nuclei, double time_cooling,
      double time_chemistry, double time_remaining);
  static ChemistryStatus evaluate_timesteps(const ChemEnvironment &environment,
      const ChemStepData &step, const ChemState &state, ChemWorkspace &workspace);
  static void import_global_environment(const ChemistryStepContext &step, ChemGlobalStep &global);
  static void import_cell_environment(const Model &model, const ChemGlobalStep &global,
      const ChemistryCellEnvironment &cell, ChemEnvironment &environment);
};
} // namespace chemistry
#endif
