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
