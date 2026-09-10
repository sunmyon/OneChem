#pragma once

#include "chemistry_newton_internal.h"
#include "chemistry_newton_reset.h"
#include "chemistry_types.h"

/* Diagnostic state must remain cheap in production builds.  The abundance
 * snapshot exists only when random-restart diagnostics are explicitly
 * compiled; otherwise this object occupies one byte and stores no species
 * arrays. */
typedef struct {
#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  double basin_entry[N_SP];
#else
  unsigned char unused;
#endif
} ChemNewtonDiagnostics;

void print_newton_basin_entry(const ChemWorkspace *work, double temperature,
                              const double *abundance);

void print_newton_iteration(
    const ChemWorkspace *work, int iteration, double temperature, double dt,
    double newton_step_error, int limiting_species, const double *abundance,
    const double *initial_abundance, const double *reaction_rhs);

void print_newton_jacobian(const ChemWorkspace *work, int iteration,
                           double rhs_jacobian[N_SP][N_SP]);

void print_random_restart_success(
    const ChemNewtonDiagnostics *diagnostics, int converged, int iteration,
    ChemNewtonResetIndex reset, int random_restart_count);

void chemistry_newton_diagnostics_report_failure(
    const ChemEnvironment *environment, const ChemWorkspace *work,
    const ChemState *input_state, double dt, double n_nuclei, double nH,
    double temperature, double newton_step_error, int limiting_species,
    int iteration, ChemNewtonResetIndex reset, int random_restart_count,
    int allow_restarts, const double *abundance,
    const double *initial_abundance, const double *newton_step,
    double convergence_abundance_threshold,
    const ChemConservationSystem *conservation);
