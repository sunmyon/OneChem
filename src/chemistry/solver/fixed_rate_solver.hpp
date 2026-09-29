#pragma once


#include "solver/chemistry_newton_internal.h"

namespace chemistry {

// Snapshot of one physical fixed-T, fixed-dt problem. Network/table pointers
// are borrowed from Model; cell/step values, initial state and rates are copied.
// A failed preparation is represented by valid()==false and cannot be solved.
class FixedRateProblem final {
public:
  FixedRateProblem(const ChemEnvironment &environment, const ChemStepData &step,
                   const ChemState &initial, const ChemRates &rates,
                   double temperature, double duration);
  FixedRateProblem(const FixedRateProblem &) = delete;
  FixedRateProblem &operator=(const FixedRateProblem &) = delete;
  FixedRateProblem(FixedRateProblem &&) = delete;
  FixedRateProblem &operator=(FixedRateProblem &&) = delete;

  bool valid() const noexcept { return valid_; }
  const ChemEnvironment &environment() const noexcept { return environment_; }
  const ChemStepData &step() const noexcept { return step_; }
  const ChemState &initial_state() const noexcept { return initial_; }
  const ChemRates &rates() const noexcept { return rates_; }
  const ChemNewtonReference &reference() const noexcept { return reference_; }
  const ChemNewtonSystem &seed() const noexcept { return seed_; }
  double temperature() const noexcept { return temperature_; }
  double duration() const noexcept { return duration_; }
  double nuclei_density() const noexcept { return nuclei_density_; }
  double hydrogen_density() const noexcept { return hydrogen_density_; }

  void evaluate_rhs(const double *abundance, double *rhs,
                    double *production, double *loss) const;
  double evaluate_residual(const double *abundance, const double *rhs,
                           const double *production, const double *loss,
                           const ChemConservationSystem &conservation) const;
  bool build_jacobian(const ChemNewtonSystem &system, const double *rhs,
                      int rebuild, int &force_numerical,
                      double jacobian[N_SP][N_SP]) const;

private:
  const ChemRunView run_;
  const ChemGlobalStep global_;
  ChemEnvironment environment_;
  const ChemStepData step_;
  const ChemState initial_;
  const ChemRates rates_;
  const double temperature_, duration_, nuclei_density_, hydrogen_density_;
  ChemNewtonReference reference_{};
  ChemNewtonSystem seed_{};
  bool valid_;
};

} // namespace chemistry


namespace chemistry {

#ifdef MAX_ITER_CHEM
inline constexpr int default_chemistry_iterations = MAX_ITER_CHEM;
#else
inline constexpr int default_chemistry_iterations = 3000;
#endif
struct FixedRateOptions {
  const double *initial_guess = nullptr; // optional stored units, per baryon
  int max_iterations = default_chemistry_iterations;
  int allow_restarts = 1;
  int report_failure = 1;
};

// Scratch belongs to one simultaneous solve. It contains no model ownership
// or accepted physical state and is not copied with candidate states.
struct NewtonScratch {
  int lapack_pivots[N_SP];
  double rhs_jacobian[N_SP][N_SP];
  double reaction_rhs[N_SP];
  double species_production[N_SP];
  double species_loss[N_SP];
  double newton_matrix[N_SP * N_SP];
  double linear_rhs[N_SP];
};

// Concrete chemistry solver during extraction. Generic numerical machinery
// will be separated after the existing physical contracts are explicit.
class FixedRateSolver {
public:
  ChemistryStatus advance(const ChemEnvironment &environment, const ChemStepData &step,
              ChemState &state, ChemWorkspace &work, double dt,
              const FixedRateOptions &options) const;
  ChemistryStatus solve(const FixedRateProblem &problem, const FixedRateOptions &options,
            ChemState &state, ChemWorkspace &work, NewtonScratch &scratch) const;
};

} // namespace chemistry
