#pragma once

#include "chemistry_types.h"
#include "radiation/radiation_rates.h"

struct chemdata_in;
struct chemdata_out;

/* Host-unit conversion belongs to this adapter, never to chemistry core. */
typedef struct {
  double specific_energy_cgs;
} ChemStandaloneUnits;

// One adapter per independent run. Model data is borrowed and must outlive it.
class ChemistryStandaloneAdapter final {
public:
  void configure(const ChemRunView *run, const ChemStandaloneUnits *units);
  bool set_timescale_safety_factor(double value);
  void prepare_step(const ChemGlobalStep *step);
  bool prepare_step_with_background(const ChemGlobalStep *step,
                                   const RadiationRateContext *radiation_context,
                                   const RadiationBlackbodyInput *background);

  // Host convention: advance returns 0 on success; the other solves return 1.
  // An unprepared adapter fails without changing output values.
  int advance(int mode, const chemdata_in *input, chemdata_out *output,
              double duration_s);
  int solve_fixed(const chemdata_in *input, chemdata_out *output,
                  double temperature, double duration_s);
  int solve_fixed_guess(const chemdata_in *input, chemdata_out *output,
                        double temperature, double duration_s,
                        const double *newton_guess, int max_iterations);
  int equilibrium_abundances(const chemdata_in *input, const chemdata_out *output,
                            double temperature, double duration_s,
                            double *full_abundance);

private:
  bool adapter_is_ready() const;
  void import_environment(const chemdata_in *, ChemEnvironment *);
  void import_state(const chemdata_out *, ChemState *);
  void export_state(const ChemState *, chemdata_out *);
  static void export_workspace(const ChemWorkspace *, chemdata_out *);
  ChemRunView standalone_run{};
  ChemGlobalStep standalone_step{};
  ChemStandaloneUnits standalone_units{};
  bool standalone_run_is_configured = false;
  bool standalone_step_is_prepared = false;
};
