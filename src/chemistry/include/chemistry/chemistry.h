#ifndef CHEMISTRY_PUBLIC_API_H
#define CHEMISTRY_PUBLIC_API_H

#include <stddef.h>

#include "chemistry_config.h"
#include "chemistry_composition.h"
#include "chemistry_solver.h"
#include "chemistry_state.h"
#include "chemistry_status.h"

// Shared value contracts; C++ callers use Model, Driver and ChemWorkspace.

/* Reconcile stored elemental budgets before chemistry: preserve ratios in
 * He, D and C families, then scale oxygen without carbon after accounting for
 * oxygen in carbon-bearing species. Species membership and charge use the
 * shared metadata. Reconstruct electrons from all stored ions. H-only species
 * are unchanged; changes in H bound to other families are not compensated.
 * Missing positive-target families or incompatible budgets/charge are errors.
 * Failure leaves the caller's abundance array unchanged. */
#ifdef __cplusplus
ChemistryStatus chemistry_normalize_abundances(double *abundance, size_t count,
                                               double metallicity);
#endif

#endif
