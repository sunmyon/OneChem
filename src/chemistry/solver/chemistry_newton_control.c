#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "chemistry_constants.h"
#include "chemistry_newton_control.h"
#include "reaction_network.h"
#include "species.h"

#define CHEM_NEWTON_MATERIAL_ABUNDANCE 1.e-10
#define CHEM_NEWTON_GROSS_TOLERANCE 1.e-10
#define CHEM_NEWTON_STEP_ABSOLUTE_TOLERANCE 3.e-11
#define CHEM_NEWTON_STEP_RELATIVE_TOLERANCE 1.e-8
#define CHEM_NEWTON_TRACE_ABUNDANCE 1.e-12

static double chemistry_newton_roundoff_factor(size_t operations)
{
  double product = (double)operations * DBL_EPSILON;
  return product < 1.0 ? product / (1.0 - product) : INFINITY;
}

double chemistry_newton_evaluate_residual(
    double dt, const double *abundance, const double *initial_abundance,
    const double *upper_bound, const double *rate,
    const double *production, const double *loss, const ChemRates *rates,
    const ChemConservationSystem *conservation)
{
  double certified_error = 0.0;
  double species_residual[N_SP];

  for(int species = 0; species < N_SP; species++) {
    double change = abundance[species] - initial_abundance[species];
    double residual = fma(rate[species], dt, -change);
    species_residual[species] = residual;
  }

  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_intermediate_equilibrium(species) ||
       conservation->is_dependent[species] ||
       dmax(fabs(initial_abundance[species]), fabs(abundance[species])) <
           CHEM_NEWTON_MATERIAL_ABUNDANCE)
      continue;

    double change = abundance[species] - initial_abundance[species];
    double gross_scale = dt * (production[species] + loss[species]) +
                         fabs(change);
    size_t reaction_terms = chemistry_network_species_reaction_count(
        rates->network, species, CHEM_REACTION_IMPLICIT);
    double rhs_roundoff = chemistry_newton_roundoff_factor(
        reaction_terms + 1) * dt * (production[species] + loss[species]);
    double residual_roundoff = chemistry_newton_roundoff_factor(3) *
        (fabs(rate[species] * dt) + fabs(change));
    double variable_scale = dmax(
        upper_bound[species],
        dmax(fabs(initial_abundance[species]), fabs(abundance[species])));
    double representation_error =
        chemistry_newton_roundoff_factor(4) * variable_scale;
    double roundoff_allowance = rhs_roundoff + residual_roundoff +
                                representation_error;
    double component = fabs(species_residual[species]) /
        (CHEM_NEWTON_GROSS_TOLERANCE * dmax(gross_scale, DBL_MIN) +
         roundoff_allowance);
    if(component > certified_error)
      certified_error = component;
  }

  double constraint_residual[CHEM_MAX_CONSERVATION_CONSTRAINTS] = {0.0};
  double constraint_absolute[CHEM_MAX_CONSERVATION_CONSTRAINTS] = {0.0};
  double constraint_roundoff[CHEM_MAX_CONSERVATION_CONSTRAINTS] = {0.0};
  size_t constraint_terms[CHEM_MAX_CONSERVATION_CONSTRAINTS] = {0};

  for(int species = 0; species < N_SP; species++) {
    if(conservation->is_dependent[species])
      continue;
    double change = abundance[species] - initial_abundance[species];
    size_t reaction_terms = chemistry_network_species_reaction_count(
        rates->network, species, CHEM_REACTION_IMPLICIT);
    double independent_roundoff =
        chemistry_newton_roundoff_factor(reaction_terms + 1) * dt *
            (production[species] + loss[species]) +
        chemistry_newton_roundoff_factor(3) *
            (fabs(rate[species] * dt) + fabs(change)) +
        chemistry_newton_roundoff_factor(4) *
            dmax(upper_bound[species],
                 dmax(fabs(initial_abundance[species]),
                      fabs(abundance[species])));
    for(int constraint = 0; constraint < conservation->count; constraint++) {
      double coefficient = conservation->coefficient[constraint][species];
      if(coefficient == 0.0)
        continue;
      double term = coefficient * species_residual[species];
      constraint_residual[constraint] += term;
      constraint_absolute[constraint] += fabs(term);
      constraint_roundoff[constraint] +=
          fabs(coefficient) * independent_roundoff;
      constraint_terms[constraint]++;
    }
  }
  for(int constraint = 0; constraint < conservation->count; constraint++)
    constraint_roundoff[constraint] +=
        chemistry_newton_roundoff_factor(constraint_terms[constraint] + 1) *
        constraint_absolute[constraint];

  for(int dependent = 0; dependent < conservation->count; dependent++) {
    enum SPECIES species = conservation->dependent_species[dependent];
    if(dmax(fabs(initial_abundance[species]), fabs(abundance[species])) <
       CHEM_NEWTON_MATERIAL_ABUNDANCE)
      continue;
    double implied_residual = 0.0;
    double implied_absolute = 0.0;
    double implied_roundoff = 0.0;
    for(int constraint = 0; constraint < conservation->count; constraint++) {
      double inverse = conservation->dependent_inverse[
          constraint * conservation->count + dependent];
      double term = inverse * constraint_residual[constraint];
      implied_residual -= term;
      implied_absolute += fabs(term);
      implied_roundoff += fabs(inverse) * constraint_roundoff[constraint];
    }
    implied_roundoff += chemistry_newton_roundoff_factor(
        (size_t)conservation->count + 1) * implied_absolute;
    double change = abundance[species] - initial_abundance[species];
    double gross_scale = dt * (production[species] + loss[species]) +
                         fabs(change);
    double component = fabs(implied_residual) /
        (CHEM_NEWTON_GROSS_TOLERANCE * dmax(gross_scale, DBL_MIN) +
         implied_roundoff + DBL_MIN);
    if(component > certified_error)
      certified_error = component;
  }

  return certified_error;
}

double chemistry_newton_material_step_norm(
    const double *abundance, const double *step,
    const ChemConservationSystem *conservation, int *limiting_species)
{
  double norm = 0.0;
  *limiting_species = -1;
  for(int species = 0; species < N_SP; species++) {
    if(conservation->is_dependent[species])
      continue;
    double before = abundance[species] - step[species];
    if(dmax(fabs(before), fabs(abundance[species])) <
       CHEM_NEWTON_MATERIAL_ABUNDANCE)
      continue;
    double scale = CHEM_NEWTON_STEP_ABSOLUTE_TOLERANCE +
        CHEM_NEWTON_STEP_RELATIVE_TOLERANCE *
            dmax(fabs(before), fabs(abundance[species]));
    double component = fabs(step[species]) / scale;
    if(component > norm) {
      norm = component;
      *limiting_species = species;
    }
  }
  return norm;
}

int chemistry_newton_normalize_current_state(
    double *abundance, const double *upper_bound)
{
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_intermediate_equilibrium(species))
      continue;
    if(!isfinite(abundance[species]))
      return 0;

    double tolerance = dmax(
        64.0 * DBL_EPSILON * dmax(1.0, fabs(upper_bound[species])),
        DBL_MIN);
    /* lower_bound is a Newton scaling/floor parameter, not the physical
     * boundary.  Exact zero and positive trace abundances remain valid. */
    if(abundance[species] < 0.0)
      return 0;
    if(abundance[species] > upper_bound[species]) {
      if(abundance[species] > upper_bound[species] + tolerance)
        return 0;
      abundance[species] = upper_bound[species];
    }
  }

  return 1;
}

static double chemistry_newton_dependent_derivative(
    int dependent, enum SPECIES variable_species,
    const ChemConservationSystem *conservation)
{
  double derivative = 0.0;
  for(int constraint = 0; constraint < conservation->count; constraint++)
    derivative -= conservation->dependent_inverse[
        constraint * conservation->count + dependent] *
        conservation->coefficient[constraint][variable_species];
  return derivative;
}

static int chemistry_newton_clip_independent_components(
    const double *initial, double *step,
    const double *lower_bound, const double *upper_bound,
    const ChemConservationSystem *conservation,
    ChemNewtonStepResult *result)
{
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_algebraic(species, conservation))
      continue;
    if(!isfinite(step[species]))
      return 0;

    double proposed = initial[species] + step[species];
    if(proposed < lower_bound[species]) {
      double distance = initial[species] - lower_bound[species];
      double outward = -step[species];
      double retained_distance = distance > 0.0 ?
          distance / (1.0 + outward / distance) : 0.0;
      step[species] =
          lower_bound[species] + retained_distance - initial[species];
    } else if(proposed > upper_bound[species]) {
      double distance = upper_bound[species] - initial[species];
      double outward = step[species];
      double retained_distance = distance > 0.0 ?
          distance / (1.0 + outward / distance) : 0.0;
      step[species] =
          upper_bound[species] - retained_distance - initial[species];
    } else {
      continue;
    }

    result->step_was_clipped = 1;
    result->clipped_species[species] = 1;
  }
  return 1;
}

static void chemistry_newton_form_trial(
    const double *initial, const double *step,
    const ChemConservationSystem *conservation, double *trial)
{
  for(int species = 0; species < N_SP; species++)
    trial[species] = chemistry_species_is_algebraic(species, conservation) ?
        initial[species] : initial[species] + step[species];
}

static int chemistry_newton_trial_is_admissible(
    const double *trial, const double *upper_bound)
{
  for(int species = 0; species < N_SP; species++) {
    if(!isfinite(trial[species]) ||
       trial[species] < 0.0 ||
       (!chemistry_species_is_intermediate_equilibrium(species) &&
        trial[species] > upper_bound[species])) {
      return 0;
    }
  }
  return 1;
}

/* Keep independent corrections component-wise.  A dependent reservoir can
 * reject the combined trial even when unrelated corrections are admissible;
 * never contract all independent species by one common factor. */
static int chemistry_newton_repair_coupled_step(
    double nH, const double *initial, double *step, double *trial,
    const double *upper_bound, const ChemRates *rates,
    const ChemConservationSystem *conservation, ChemNewtonStepResult *result)
{
  for(int pass = 0; pass < N_SP; pass++) {
    int best_species = -1;
    double best_violation = INFINITY;
    int best_is_valid = 0;

    for(int species = 0; species < N_SP; species++) {
      if(chemistry_species_is_algebraic(species, conservation) ||
         step[species] == 0.0)
        continue;

      double candidate[N_SP];
      for(int other = 0; other < N_SP; other++)
        candidate[other] = chemistry_species_is_algebraic(other, conservation) ?
            initial[other] : initial[other] +
            (other == species ? 0.0 : step[other]);
      int algebraic_ok = chemistry_update_implicit_algebraic_species(
          nH, candidate, rates, conservation);
      int admissible = algebraic_ok && chemistry_newton_trial_is_admissible(
          candidate, upper_bound);

      /* Reconstruct conservation reservoirs for a comparable violation score
       * even if the coupled intermediate iteration rejected the candidate. */
      if(!admissible)
        chemistry_update_conservation_dependent_species(candidate, conservation);
      double violation = 0.0;
      for(int other = 0; other < N_SP; other++) {
        if(!isfinite(candidate[other])) {
          violation = INFINITY;
          break;
        }
        double scale = fmax(upper_bound[other], 1.e-30);
        if(candidate[other] < 0.0)
          violation += -candidate[other] / scale;
        else if(!chemistry_species_is_intermediate_equilibrium(other) &&
                candidate[other] > upper_bound[other])
          violation += (candidate[other] - upper_bound[other]) / scale;
      }
      if(!admissible && violation == 0.0)
        violation = INFINITY;
      if(best_species < 0 || (admissible && !best_is_valid) ||
         (admissible == best_is_valid && violation < best_violation)) {
        best_species = species;
        best_violation = violation;
        best_is_valid = admissible;
      }
    }

    if(best_species < 0 || !isfinite(best_violation))
      return 0;

    double original_step = step[best_species];
    if(best_is_valid) {
      double feasible = 0.0, infeasible = 1.0;
      for(int iteration = 0; iteration < 48; iteration++) {
        double fraction = 0.5 * (feasible + infeasible);
        step[best_species] = fraction * original_step;
        chemistry_newton_form_trial(initial, step, conservation, trial);
        int valid = chemistry_update_implicit_algebraic_species(
            nH, trial, rates, conservation);
        if(valid)
          valid = chemistry_newton_trial_is_admissible(
              trial, upper_bound);
        if(valid)
          feasible = fraction;
        else
          infeasible = fraction;
      }
      step[best_species] = 0.99 * feasible * original_step;
    } else {
      step[best_species] = 0.0;
    }
    result->step_was_clipped = 1;
    result->clipped_species[best_species] = 1;
    chemistry_newton_form_trial(initial, step, conservation, trial);
    if(chemistry_update_implicit_algebraic_species(
           nH, trial, rates, conservation) &&
       chemistry_newton_trial_is_admissible(
           trial, upper_bound))
      return 1;
  }
  return 0;
}

static int chemistry_newton_project_negative_charge_carrier(
    double nH, const double *initial, double *step, double *trial,
    const double *lower_bound, const ChemRates *rates,
    const ChemConservationSystem *conservation,
    ChemNewtonStepResult *result)
{
  int charge_dependent = -1;
  for(int dependent = 0; dependent < conservation->count; dependent++) {
    enum SPECIES species = conservation->dependent_species[dependent];
    if(conservation->total[dependent] == 0.0 &&
       chemistry_species_definition(species)->charge < 0) {
      charge_dependent = dependent;
      break;
    }
  }
  if(charge_dependent < 0)
    return 0;

  enum SPECIES dependent_species =
      conservation->dependent_species[charge_dependent];
  if(!isfinite(trial[dependent_species]) ||
     trial[dependent_species] >= lower_bound[dependent_species])
    return 0;

  int projected_species = -1;
  double projected_derivative = 0.0;
  double largest_reducible_charge = -1.0;
  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_algebraic(species, conservation) ||
       chemistry_species_definition(species)->charge >= 0)
      continue;
    double derivative = chemistry_newton_dependent_derivative(
        charge_dependent, species, conservation);
    if(!(derivative < 0.0) || trial[species] <= lower_bound[species])
      continue;
    double reducible_charge =
        fabs((double)chemistry_species_definition(species)->charge) *
        (trial[species] - lower_bound[species]);
    if(reducible_charge > largest_reducible_charge) {
      largest_reducible_charge = reducible_charge;
      projected_species = species;
      projected_derivative = derivative;
    }
  }
  if(projected_species < 0)
    return 0;

  double target = lower_bound[dependent_species];
  double correction =
      (target - trial[dependent_species]) / projected_derivative;
  double projected = dmax(lower_bound[projected_species],
                          trial[projected_species] + correction);
  if(projected >= trial[projected_species])
    return 0;
  if(projected > lower_bound[projected_species])
    projected = nextafter(projected, -INFINITY);
  step[projected_species] = projected - initial[projected_species];
  chemistry_newton_form_trial(initial, step, conservation, trial);
  if(!chemistry_update_implicit_algebraic_species(
         nH, trial, rates, conservation))
    return 0;
  result->step_was_clipped = 1;
  result->clipped_species[projected_species] = 1;
  return 1;
}

int chemistry_newton_apply_step(
    double nH, double *abundance, double *step,
    const double *lower_bound, const double *upper_bound,
    const ChemRates *rates, const ChemConservationSystem *conservation,
    ChemNewtonStepResult *result)
{
  double initial[N_SP], trial[N_SP];
  memcpy(initial, abundance, sizeof(initial));
  memset(result, 0, sizeof(*result));

  if(!chemistry_newton_clip_independent_components(
         initial, step, lower_bound, upper_bound, conservation,
         result))
    return 0;

  chemistry_newton_form_trial(initial, step, conservation, trial);
  int valid = chemistry_update_implicit_algebraic_species(
      nH, trial, rates, conservation);
  if(!valid)
    valid = chemistry_newton_project_negative_charge_carrier(
        nH, initial, step, trial, lower_bound, rates, conservation, result);
  if(valid)
    valid = chemistry_newton_trial_is_admissible(
        trial, upper_bound);

  if(!valid)
    valid = chemistry_newton_repair_coupled_step(
        nH, initial, step, trial, upper_bound, rates, conservation, result);

  if(!valid)
    return 0;
  for(int species = 0; species < N_SP; species++) {
    step[species] = trial[species] - initial[species];
    abundance[species] = trial[species];
  }
  return 1;
}

int chemistry_newton_build_explicit_guess(
    double nH, double dt, const double *initial_abundance,
    double *abundance, const double *lower_bound, const double *upper_bound,
    const ChemRates *rates, const ChemConservationSystem *conservation)
{
  double rhs[N_SP], step[N_SP] = {0.0};
  chemistry_network_rhs(
      rates->network, CHEM_REACTION_IMPLICIT, nH, initial_abundance,
      rates->forward, rates->inverse, rhs);

  for(int species = 0; species < N_SP; species++) {
    if(chemistry_species_is_algebraic(species, conservation))
      continue;
    step[species] = dt * rhs[species];
  }

  memcpy(abundance, initial_abundance, N_SP * sizeof(*abundance));
  ChemNewtonStepResult result;
  if(!chemistry_newton_apply_step(nH, abundance, step, lower_bound,
                                  upper_bound, rates, conservation, &result)) {
    /* A predictor is optional; retain the valid physical seed on failure. */
    memcpy(abundance, initial_abundance, N_SP * sizeof(*abundance));
  }
  return 1;
}

int chemistry_newton_reselect_conservation_chart(
    double metallicity, ChemNewtonSystem *system)
{
  /* Historical note: before compensated summation was applied to the
   * reaction RHS, solver development required a stateful hysteretic recovery
   * mechanism to resolve repeated boundary contact and Newton stagnation.  It
   * retained, for each solve, the attempted C/O pivot pairs, a blacklist of
   * failed pivots, scores for alternative conservation charts, and rollback
   * state for unsuccessful chart changes.  After compensated summation was
   * introduced, the 100k and 1M off-trajectory validation sets both converged
   * completely with the simpler stateless recovery used here, suggesting
   * that cancellation-driven stagnation had motivated much of the historical
   * machinery.  The history is intentionally documented because those tests
   * establish convergence only: they do not prove identical root selection
   * or continuity on every time-dependent trajectory.  Reassess whether the
   * historical hysteretic mechanism is needed if either property regresses.
   */
  ChemConservationSystem replacement;
  if(!chemistry_build_coupled_conservation_system(
         &replacement, metallicity, system->abundance) ||
     replacement.count != system->conservation.count)
    return 0;

  int changed = 0;
  for(int constraint = 0;
      constraint < system->conservation.count;
      constraint++)
    changed |= replacement.dependent_species[constraint] !=
               system->conservation.dependent_species[constraint];
  if(!changed)
    return 0;

  double trial[N_SP];
  memcpy(trial, system->abundance, sizeof(trial));
  if(!chemistry_update_conservation_dependent_species(
         trial, &replacement))
    return 0;
  for(int species = 0; species < N_SP; species++)
    if(!isfinite(trial[species]) || trial[species] < 0.0)
      return 0;

  memcpy(system->abundance, trial, sizeof(trial));
  system->conservation = replacement;
  system->solver_species_count = 0;
  for(int species = 0; species < N_SP; species++) {
    if(!chemistry_species_is_algebraic(species, &system->conservation))
      system->solver_species[system->solver_species_count++] = species;
    system->total_step[species] =
        chemistry_species_is_algebraic(species, &system->conservation) ? 0.0 :
        system->abundance[species] - system->initial_abundance[species];
  }
  return 1;
}
