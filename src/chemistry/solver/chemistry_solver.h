#pragma once

#include "chemistry_types.h"

int chemistry_update_intermediate_equilibrium_species(unsigned char mode, double nH,
                                                       double *abundance, const ChemRates *rates);

#if defined(__cplusplus) && defined(ENERGY_IMPLICIT)
namespace chemistry {
class EnergySolver final {
public:
  ChemistryStatus solve(const ChemEnvironment &environment, const ChemStepData &step,
            ChemState &state, ChemWorkspace &workspace,
            ChemistrySolver chemistry_solver, double dt) const;
private:
  struct EnergyInput {
    const ChemEnvironment &environment;
    const ChemStepData &step;
    ChemistrySolver solver;
    double duration;
  };
  enum class ResidualMode { trial, commit, commit_with_report };
  struct EnergyContinuationPoint {
    double temperature;
    double energy;
    double residual;
    ChemState state;
    ChemWorkspace workspace;
  };
  static void restore_workspace(ChemWorkspace *destination,
                              const ChemWorkspace *source);
  static double evaluate_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g, ResidualMode mode, const double *newton_guess, int max_chemistry_iterations);
  static double trial_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g);
  static double residual_with_guess(const EnergyInput &problem, ChemState *state,
                                ChemWorkspace *work, double T_gas_new, double *e_new,
                                double *g, const double *newton_guess,
                                int max_chemistry_iterations);
  static double commit_residual(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_gas_new, double *e_new, double *g);
  static double commit_residual_with_report(const EnergyInput &problem, ChemState *state,
                              ChemWorkspace *work,
                              double T_gas_new, double *e_new, double *g);
  static int evaluate_energy_continuation_point(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace, double temperature,
    const double *newton_guess, EnergyContinuationPoint *point);
  static int refine_energy_continuation_bracket(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace,
    EnergyContinuationPoint lower, EnergyContinuationPoint upper,
    EnergyContinuationPoint *root);
  static int solve_continuation(const EnergyInput &problem, const ChemState *initial_state, const ChemWorkspace *initial_workspace, double temperature_initial,
    double temperature_min, double temperature_max, ChemState *state,
    ChemWorkspace *workspace, double *temperature, double *energy,
    double *residual);
  static double solve_bisection(const EnergyInput &problem, ChemState *state, ChemWorkspace *work, double T_min, double g_min, double T_max, double g_max, double *T_new0, double *e_new0, int max_iter);
};
} // namespace chemistry
#endif

#ifdef __cplusplus
struct ChemConservationSystem;
namespace chemistry {
class ExplicitSolver final {
public:
  ChemistryStatus solve(const ChemEnvironment &environment, const ChemStepData &step,
            ChemState &state, ChemWorkspace &workspace, double dt) const;
private:
  static void update_species(enum SPECIES species, double nH, double dt, const double *initial_abundance, double *abundance, const ChemRates *rates);
  static int update_algebraic_species(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation);
  static int update_simultaneously(
    double nH, double metallicity, double dt, double *abundance,
    const ChemRates *rates);
};
#ifdef SAHA
class SahaSolver final {
public:
  ChemistryStatus solve(const ChemEnvironment &environment, ChemState &state,
            ChemWorkspace &workspace, double dt) const;
};
#endif
} // namespace chemistry
#endif
