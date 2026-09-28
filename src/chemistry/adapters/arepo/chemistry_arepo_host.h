#ifndef CHEMISTRY_AREPO_HOST_H
#define CHEMISTRY_AREPO_HOST_H

#include "chemistry_arepo_adapter.h"

/* Thin AREPO lifecycle hooks. Each MPI process owns its own runtime. */
int chemistry_arepo_host_init(void);
void chemistry_arepo_host_destroy(void);
const ChemArepoContext *chemistry_arepo_host_context(void);
int chemistry_arepo_host_prepare_step(ChemistryStepContext *step);

int chemistry_arepo_blackbody_rates(double effective_temperature,
                                    double *photodetachment_rate, double *g0);
int chemistry_arepo_dust_opacities(double temperature,
                                   double *kappa_rosseland,
                                   double *kappa_planck);

#endif
