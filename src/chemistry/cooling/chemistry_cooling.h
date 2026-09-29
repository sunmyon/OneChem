#pragma once

#include "chemistry_types.h"



#ifdef __cplusplus
namespace chemistry {
// Borrows one cell's inputs and workspace. The caller prepares temperature,
// reaction rates and continuum state before evaluating the corresponding terms.
class CoolingEvaluator final {
public:
  CoolingEvaluator(const ChemEnvironment &environment, const ChemState &state,
                   ChemWorkspace &workspace) noexcept
      : environment_(environment), state_(state), workspace_(workspace) {}
  CoolingEvaluator(const CoolingEvaluator &) = delete;
  CoolingEvaluator &operator=(const CoolingEvaluator &) = delete;
  ChemistryStatus radiative(const ChemStepData &step);
  ChemistryStatus high_temperature(const ChemStepData &step);
  ChemistryStatus chemical(double dt, const double *old_abundances);
private:
  const ChemEnvironment &environment_;
  const ChemState &state_;
  ChemWorkspace &workspace_;
};
} // namespace chemistry
#endif
