#include "chemistry_newton_diagnostics.h"

#include <math.h>
#include <stdio.h>

#include "chemistry_constants.h"
#include "reaction_network.h"
#include "species.h"

void print_newton_basin_entry(const ChemWorkspace *work, double temperature,
                              const double *abundance) {
  if (!work->diagnostics_enabled)
    return;

  printf("T_gas=%g y=%g %g %g %g %g\n", temperature, abundance[sELEC],
         abundance[sHI], abundance[sHII], abundance[sH2I], abundance[sHM]);
}

void print_newton_iteration(
    const ChemWorkspace *work, int iteration, double temperature, double dt,
    double newton_step_error, int limiting_species, const double *abundance,
    const double *initial_abundance, const double *reaction_rhs) {
  if (!work->diagnostics_enabled)
    return;

  printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g reaction_rhs=%g %g %g "
         "%g %g dt=%g\n",
         iteration, temperature, newton_step_error, limiting_species,
         abundance[sELEC], abundance[sHI], abundance[sHII], abundance[sH2I],
         abundance[sHM], reaction_rhs[sELEC], reaction_rhs[sHI],
         reaction_rhs[sHII], reaction_rhs[sH2I], reaction_rhs[sHM], dt);

  printf("[%d] T_gas=%g y_H2=%g total_step=%g resid=%g\n", iteration,
         temperature, abundance[sH2I],
         abundance[sH2I] - initial_abundance[sH2I], reaction_rhs[sH2I] * dt);

#ifdef ID_CHEM_DEBUG
#ifdef HD
  printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g %g %g %g dt=%g\n",
         iteration, temperature, newton_step_error, limiting_species,
         abundance[sELEC], abundance[sHI], abundance[sHII], abundance[sH2I],
         abundance[sHM], abundance[sDI], abundance[sDII], abundance[sHDI], dt);
#else
  printf("[%d] T_gas=%g err=%g flag=%d y=%g %g %g %g %g\n", iteration,
         temperature, newton_step_error, limiting_species, abundance[sELEC],
         abundance[sHI], abundance[sHII], abundance[sH2I], abundance[sHM]);
#endif

  printf("[%d] T_gas=%g y_H2=%g total_step=%g resid=%g\n", iteration,
         temperature, abundance[sH2I],
         abundance[sH2I] - initial_abundance[sH2I], reaction_rhs[sH2I] * dt);
#endif
}

void print_newton_jacobian(
    const ChemWorkspace *work, int iteration,
    double rhs_jacobian[N_SP][N_SP]) {
#ifdef ID_CHEM_DEBUG
  if (!work->diagnostics_enabled)
    return;

  for (int species = 0; species < N_SP; species++) {
    printf("[%d] rhs_jacobian[%d]=%g %g %g %g %g %g %g %g\n", iteration,
           species, rhs_jacobian[species][0], rhs_jacobian[species][1],
           rhs_jacobian[species][2], rhs_jacobian[species][3],
           rhs_jacobian[species][4], rhs_jacobian[species][5],
           rhs_jacobian[species][6], rhs_jacobian[species][7]);
  }
#else
  (void)work;
  (void)iteration;
  (void)rhs_jacobian;
#endif
}

void print_random_restart_success(
    const ChemNewtonDiagnostics *diagnostics, int converged, int iteration,
    ChemNewtonResetIndex reset, int random_restart_count) {
#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  if (!converged)
    return;

  printf("CHEM_RANDOM_SUCCESS iterations=%d reset=%d random_restarts=%d\n",
         iteration, reset, random_restart_count);

  for (int species = 0; species < N_SP; species++) {
    printf("guess_%s=%.17g\n", chemistry_species_definition(static_cast<SPECIES>(species))->name,
           diagnostics->basin_entry[species]);
  }
#else
  (void)diagnostics;
  (void)converged;
  (void)iteration;
  (void)reset;
  (void)random_restart_count;
#endif
}

static void print_current_abundances(const double *abundance) {
  printf("y[%d]=%g, y[%d]=%g, y[%d]=%g, ", sELEC, abundance[sELEC], sHI,
         abundance[sHI], sHII, abundance[sHII]);

#ifdef HELIUM
  printf("y[%d]=%g, y[%d]=%g, y[%d]=%g, ", sHeI, abundance[sHeI], sHeII,
         abundance[sHeII], sHeIII, abundance[sHeIII]);
#endif

  printf("y[%d]=%g, ", sH2I, abundance[sH2I]);

#ifdef H2_PLUS
  printf("y[%d]=%g, ", sH2II, abundance[sH2II]);
#endif

  printf("y[%d]=%g\n", sHM, abundance[sHM]);

#ifdef HD
  printf("y[%d]=%g, y[%d]=%g y[%d]=%g", sDI, abundance[sDI], sDII,
         abundance[sDII], sHDI, abundance[sHDI]);

#ifdef HD_PLUS
  printf(" y[%d]=%g", sHDII, abundance[sHDII]);
#endif

#ifdef HD_MINUS
  printf(" y[%d]=%g", sDM, abundance[sDM]);
#endif

  printf("\n");
#endif
}

static void print_initial_abundances(const double *initial_abundance) {
  printf("abundance[%d]=%g, abundance[%d]=%g, abundance[%d]=%g, ", sELEC,
         initial_abundance[sELEC], sHI, initial_abundance[sHI], sHII,
         initial_abundance[sHII]);

#ifdef HELIUM
  printf("abundance[%d]=%g, abundance[%d]=%g, abundance[%d]=%g, ", sHeI,
         initial_abundance[sHeI], sHeII, initial_abundance[sHeII], sHeIII,
         initial_abundance[sHeIII]);
#endif

  printf("abundance[%d]=%g, ", sH2I, initial_abundance[sH2I]);

#ifdef H2_PLUS
  printf("abundance[%d]=%g, ", sH2II, initial_abundance[sH2II]);
#endif

  printf("abundance[%d]=%g\n", sHM, initial_abundance[sHM]);

#ifdef HD
  printf("abundance[%d]=%g, abundance[%d]=%g abundance[%d]=%g", sDI,
         initial_abundance[sDI], sDII, initial_abundance[sDII], sHDI,
         initial_abundance[sHDI]);

#ifdef HD_PLUS
  printf(" abundance[%d]=%g", sHDII, initial_abundance[sHDII]);
#endif

#ifdef HD_MINUS
  printf(" abundance[%d]=%g", sDM, initial_abundance[sDM]);
#endif

  printf("\n");
#endif
}

static void print_failure_input(const ChemEnvironment *environment,
                                const ChemState *input_state, double dt,
                                double n_nuclei, double temperature) {
  printf("CHEM_FAILURE_INPUT_BEGIN\n");
  printf("n_nuclei_cm3=%.17g\n", n_nuclei);
  printf("T_K=%.17g\n", temperature);
  printf("T_CMB_K=%.17g\n", environment->global->cmb_temperature);
  printf("Z=%.17g\n", environment->metallicity);
  printf("dt_s=%.17g\n", dt);

  for (int species = 0; species < N_SP; species++) {
    printf("x0_%s=%.17g\n", chemistry_species_definition(static_cast<SPECIES>(species))->name,
           input_state->abund[species] / XH);
  }

  printf("CHEM_FAILURE_INPUT_END\n");
}

static void print_failure_species(
    const ChemWorkspace *work, double nH, double dt, const double *abundance,
    const double *initial_abundance, const double *newton_step,
    double convergence_abundance_threshold,
    const ChemConservationSystem *conservation) {
  double final_rhs[N_SP];

  chemistry_network_rhs(work->rates.network, CHEM_REACTION_IMPLICIT, nH,
                        abundance, work->rates.forward, work->rates.inverse,
                        final_rhs);

  printf("CHEM_FAILURE_SPECIES_BEGIN\n");

  for (int species = 0; species < N_SP; species++) {
    double implicit_residual =
        final_rhs[species] * dt -
        (abundance[species] - initial_abundance[species]);

    double relative_step =
        fabs(newton_step[species]) /
        dmax(fabs(abundance[species]), convergence_abundance_threshold);

    printf("CHEM_FAILURE_SPECIES name=%s abundance=%.17g initial=%.17g"
           " newton_step=%.17g relative_step=%.17g residual=%.17g"
           " algebraic=%d\n",
           chemistry_species_definition(static_cast<SPECIES>(species))->name, abundance[species],
           initial_abundance[species], newton_step[species], relative_step,
           implicit_residual,
           (conservation)->is_algebraic(static_cast<SPECIES>(species)));
  }

  printf("CHEM_FAILURE_SPECIES_END\n");
}

void chemistry_newton_diagnostics_report_failure(
    const ChemEnvironment *environment, const ChemWorkspace *work,
    const ChemState *input_state, double dt, double n_nuclei, double nH,
    double temperature, double newton_step_error, int limiting_species,
    int iteration, ChemNewtonResetIndex reset, int random_restart_count,
    int allow_restarts, const double *abundance,
    const double *initial_abundance, const double *newton_step,
    double convergence_abundance_threshold,
    const ChemConservationSystem *conservation) {
  printf("iteration is not converged err=%le(in %d),iteration=%d,"
         " id=%d, reset=%d, random_restarts=%d, allow_restarts=%d.\n dt=%g, "
         "rho=%g, T_gas=%g, n_nuclei=%g;\n",
         newton_step_error, limiting_species, iteration, environment->ID, reset,
         random_restart_count, allow_restarts, dt, environment->rho_cgs,
         temperature, n_nuclei);

  print_current_abundances(abundance);
  print_initial_abundances(initial_abundance);
  print_failure_input(environment, input_state, dt, n_nuclei, temperature);
  print_failure_species(work, nH, dt, abundance, initial_abundance, newton_step,
                        convergence_abundance_threshold, conservation);
}
