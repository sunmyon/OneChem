#ifndef CHEMISTRY_NEW_API_H
#define CHEMISTRY_NEW_API_H

#include "chemistry_types.h"

/* Transitional host-facing API for chemistry_new. New integrations should
 * use chemistry_portable's public API and a host adapter. */
int InitChem(void);
int compute_abundances_d(int mode, const chemdata_in *input,
                         chemdata_out *output, double duration_s);
int compute_chemistry_implicit_d(const chemdata_in *input,
                                 chemdata_out *output,
                                 double temperature, double duration_s);
int compute_chemistry_implicit_guess_d(const chemdata_in *input,
                                       chemdata_out *output,
                                       double temperature, double duration_s,
                                       const double *newton_guess,
                                       int max_iterations);
int compute_equilibrium_abundances_d(const chemdata_in *input,
                                     const chemdata_out *output,
                                     double temperature, double duration_s,
                                     double *full_abundance);

#endif
