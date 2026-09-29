#pragma once

#include "chemistry_core.h"
#include "chemistry_types.h"

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 6
#else
#define CHEM_MAX_CONSERVATION_CONSTRAINTS 4
#endif

class ChemConservationSystem {
public:
  bool build(double metallicity, const double *abundance);
  bool is_algebraic(enum SPECIES species) const;
  bool project(double *abundance) const;
  bool update_algebraic(double nH, double *abundance, const ChemRates *rates) const;
  bool reduce_jacobian(double nH, const double *abundance, const ChemRates *rates, const enum SPECIES *solver_species, int n_solver_species, double drfdy[N_SP][N_SP]) const;
  int constraint_count() const noexcept { return count; }
  bool dependent(int species) const noexcept { return is_dependent[species]; }
  SPECIES reservoir(int row) const noexcept { return dependent_species[row]; }
  double constraint_total(int row) const noexcept { return total[row]; }
  double weight(int row, int species) const noexcept { return coefficient[row][species]; }
  double inverse_entry(int index) const noexcept { return dependent_inverse[index]; }
private:
  bool invert();
  void add_constraint(SPECIES species, double value);
  int count{};
  double coefficient[CHEM_MAX_CONSERVATION_CONSTRAINTS][N_SP]{};
  double total[CHEM_MAX_CONSERVATION_CONSTRAINTS]{};
  enum SPECIES dependent_species[CHEM_MAX_CONSERVATION_CONSTRAINTS]{};
  double dependent_inverse[
      CHEM_MAX_CONSERVATION_CONSTRAINTS *
      CHEM_MAX_CONSERVATION_CONSTRAINTS]{};
  unsigned char is_dependent[N_SP]{};
};

/* Frozen physical initial state and bounds for a fixed-rate problem. */
typedef struct {
  double initial_abundance[N_SP];
  double cooling_initial_abundance[N_SP];
  double initialization_floor[N_SP];
  double newton_lower_bound[N_SP];
  double upper_bound[N_SP];
} ChemNewtonReference;

/* Mutable iterate/chart plus borrowed, read-only problem data. */
struct ChemNewtonSystem {
  double abundance[N_SP];
  const double *initial_abundance;
  const double *cooling_initial_abundance;
  const double *initialization_floor;
  const double *newton_lower_bound;
  const double *upper_bound;
  double total_step[N_SP];
  enum SPECIES solver_species[N_SP];
  int solver_species_count;
  ChemConservationSystem conservation;
  bool prepare(const ChemEnvironment *environment, const ChemState *state,
               const ChemRates *rates, double nH, ChemNewtonReference &reference);
  bool initialize(double nH, double dt, const double *guess, const ChemRates *rates);
  bool build_jacobian(double nH, const double *rhs, const ChemRates *rates,
                      int rebuild, int *force_numerical, double jacobian[N_SP][N_SP]) const;
  void update_timescale(ChemWorkspace *work, double dt) const;
  bool export_solution(const ChemEnvironment *environment,
                       ChemState *state, ChemWorkspace *work, const ChemState *input,
                       double dt, double nH);
};

int chemistry_species_is_intermediate_equilibrium(enum SPECIES species);
int chemistry_update_intermediate_equilibrium_species(
    unsigned char mode, double nH, double *abundance, const ChemRates *rates);

#ifdef CHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
int chemistry_is_full_implicit_initialized_intermediate(
    enum SPECIES species);
#endif

/* Shared Newton preparation and evaluation helpers. */

void chemistry_print_rerun_line(const ChemEnvironment *environment,
                                       double n_nuclei, double T_gas,
                                       double dt, const ChemState *state);

int chemistry_abundances_are_nonnegative(const double *abundance);

void chemistry_correct_full_implicit_intermediate_guess(
    double nH, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation);

#ifdef __cplusplus
namespace chemistry { using ConservationSystem = ::ChemConservationSystem; }
#endif
