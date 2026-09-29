#pragma once

#include "chemistry_newton_internal.h"


typedef struct {
  int step_was_clipped;
  int clipped_species[N_SP];
} ChemNewtonStepResult;

double chemistry_newton_evaluate_residual(
    double dt, const double *abundance, const double *initial_abundance,
    const double *upper_bound, const double *rate,
    const double *production, const double *loss, const ChemRates *rates,
    const ChemConservationSystem *conservation);

double chemistry_newton_material_step_norm(
    const double *abundance, const double *step,
    const ChemConservationSystem *conservation, int *limiting_species);

int chemistry_newton_normalize_current_state(
    double *abundance, const double *upper_bound);

int chemistry_newton_apply_step(
    double nH, double *abundance, double *step,
    const double *lower_bound, const double *upper_bound,
    const ChemRates *rates, const ChemConservationSystem *conservation,
    ChemNewtonStepResult *result);

int chemistry_newton_build_explicit_guess(
    double nH, double dt, const double *initial_abundance,
    double *abundance, const double *lower_bound, const double *upper_bound,
    const ChemRates *rates, const ChemConservationSystem *conservation);

int chemistry_newton_reselect_conservation_chart(
    double metallicity, ChemNewtonSystem *system);
