#pragma once

#include "../../include/chemistry/chemistry.h"

#ifdef __cplusplus
namespace chemistry { class Model; }
const chemistry::Model *chemistry_arepo_host_context(void);
extern "C" {
#endif
/* AREPO owns one Model per MPI rank. Numerical state stays in the core. */
int chemistry_arepo_host_init(void);
void chemistry_arepo_host_destroy(void);

int chemistry_arepo_host_prepare_step(ChemistryStepContext *step);
int chemistry_arepo_blackbody_rates(double effective_temperature,
                                    double *photodetachment_rate, double *g0);
int chemistry_arepo_dust_opacities(double temperature,
                                   double *kappa_rosseland, double *kappa_planck);

#ifdef __cplusplus
}
#endif
