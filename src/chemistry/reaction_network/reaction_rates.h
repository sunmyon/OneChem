#pragma once

#include "chemistry_types.h"




#ifdef __cplusplus
namespace chemistry {

// Short-lived, non-owning evaluation session. Inputs and workspace must outlive
// it; the caller prepares workspace.thermodynamics.temperature before evaluate().
// Model data stays immutable, while all cell-specific outputs go to workspace.
class ReactionRateEvaluator final {
public:
  ReactionRateEvaluator(const ChemEnvironment &environment,
                        const ChemState &state, ChemWorkspace &workspace) noexcept
      : environment_(environment), state_(state), workspace_(workspace) {}
  ReactionRateEvaluator(const ReactionRateEvaluator &) = delete;
  ReactionRateEvaluator &operator=(const ReactionRateEvaluator &) = delete;

  int evaluate(const ChemStepData &step);
  // Uses the coefficients already present in workspace.rates.
  int evaluate_timescale();

private:
  void evaluate_chemical_energy_timescale();
  static double species_timescale(double abundance, double production,
                                 double destruction);
  const ChemEnvironment &environment_;
  const ChemState &state_;
  ChemWorkspace &workspace_;
};

} // namespace chemistry
#endif
