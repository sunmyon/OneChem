#pragma once

#include "chemistry_types.h"
#include "radiation/radiation_rates.h"

#ifdef __cplusplus
extern "C" {
#endif

struct chemdata_in;
struct chemdata_out;

/* Host-unit conversion belongs to this adapter, never to chemistry core. */
typedef struct {
  double specific_energy_cgs;
} ChemStandaloneUnits;

/* Called once per run, after the standalone host has established its units. */
void chemistry_standalone_configure(const ChemRunView *run,
                                    const ChemStandaloneUnits *units);
int chemistry_standalone_set_timescale_safety_factor(double value);

/* Called once whenever CMB, scale factor, or uniform backgrounds change. */
void chemistry_standalone_prepare_step(const ChemGlobalStep *step);

/* Convenience boundary for the standalone host.  The radiation model is
 * host-owned and outlives compute_abundances_d; only its derived rates are
 * copied into the adapter's shared chemistry-step state. */
int chemistry_standalone_prepare_step_with_background(
    const ChemGlobalStep *step,
    const RadiationRateContext *radiation_context,
    const RadiationBlackbodyInput *background);

/* Compatibility entry points used only by the standalone driver. AREPO and
 * other hosts should call their own adapter API instead. */
int compute_abundances_d(int mode, const struct chemdata_in *input,
                         struct chemdata_out *output, double duration_s);
int compute_chemistry_implicit_d(const struct chemdata_in *input,
                                 struct chemdata_out *output,
                                 double temperature, double duration_s);
int compute_chemistry_implicit_guess_d(const struct chemdata_in *input,
                                       struct chemdata_out *output,
                                       double temperature, double duration_s,
                                       const double *newton_guess,
                                       int max_iterations);
int compute_equilibrium_abundances_d(const struct chemdata_in *input,
                                     const struct chemdata_out *output,
                                     double temperature, double duration_s,
                                     double *full_abundance);

#ifdef __cplusplus
}
#endif
