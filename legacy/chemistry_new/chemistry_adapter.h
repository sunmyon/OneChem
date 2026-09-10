#ifndef CHEMISTRY_ADAPTER_H
#define CHEMISTRY_ADAPTER_H

#ifdef CHEMISTRY_PORTABLE
#include "chemistry_standalone_types.h"
#endif

#ifndef CHEMISTRY_OLD
void particle_to_abundances_1zone(int i, double *abund);
void abundances_to_particle_1zone(int i, const double *abund);
void particle2in_chemdata_1zone(int i, struct chemdata_in *in, struct chemdata_out *out);
void out2particle_chemdata_1zone(int i, const struct chemdata_out *out);
int compute_gamma(int ithis);
double compute_temp(int ithis);
#endif

#endif
