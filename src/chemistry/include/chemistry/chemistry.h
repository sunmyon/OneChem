#ifndef CHEMISTRY_PUBLIC_API_H
#define CHEMISTRY_PUBLIC_API_H

#include <stddef.h>

#include "chemistry_config.h"
#include "chemistry_composition.h"
#include "chemistry_solver.h"
#include "chemistry_state.h"
#include "chemistry_status.h"

typedef struct ChemistryContext ChemistryContext;
typedef struct ChemistryWorkspace ChemistryWorkspace;

/*
 * Explicit, transactional run initialization.  This is the only public call
 * allowed to load tables or allocate run-lifetime chemistry data.  On failure
 * *context is NULL. A successful context is read-only during evaluations and
 * may be shared by callers that use separate workspaces.
 */
ChemistryStatus chemistry_context_create(const ChemistryConfig *config,
                                         ChemistryContext **context);
void chemistry_context_destroy(ChemistryContext *context);

/* Reconcile stored elemental budgets before chemistry: preserve ratios in
 * He, D and C families, then scale oxygen without carbon after accounting for
 * oxygen in carbon-bearing species. Species membership and charge use the
 * shared metadata. Reconstruct electrons from all stored ions. H-only species
 * are unchanged; changes in H bound to other families are not compensated.
 * Missing positive-target families or incompatible budgets/charge are errors.
 * Failure leaves the caller's abundance array unchanged. */
ChemistryStatus chemistry_normalize_abundances(double *abundance, size_t count,
                                               double metallicity);

size_t chemistry_species_count(const ChemistryContext *context);
const char *chemistry_species_name(const ChemistryContext *context,
                                   size_t species_index);

size_t chemistry_workspace_size(const ChemistryContext *context);
ChemistryStatus chemistry_workspace_init(const ChemistryContext *context,
                                          void *storage,
                                          size_t storage_size,
                                          ChemistryWorkspace **workspace);
void chemistry_workspace_reset(ChemistryWorkspace *workspace);

ChemistryStatus chemistry_step_prepare(
    const ChemistryContext *context,
    const ChemistryGlobalEnvironment *global_environment,
    ChemistryStepContext *step);

ChemistryStatus chemistry_dust_opacities(
    const ChemistryContext *context, double temperature,
    double *kappa_rosseland, double *kappa_planck);

/* Performs no file I/O, hidden initialization, or run-lifetime allocation. */
ChemistryStatus chemistry_advance(
    const ChemistryContext *context,
    const ChemistryStepContext *step,
    ChemistryWorkspace *workspace,
    const ChemistryCellEnvironment *cell,
    double duration_s,
    const ChemistryAdvanceOptions *options,
    ChemistryState *state,
    ChemistryResult *result);

#endif
