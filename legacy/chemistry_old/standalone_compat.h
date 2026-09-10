#ifndef CHEMISTRY_OLD_STANDALONE_COMPAT_H
#define CHEMISTRY_OLD_STANDALONE_COMPAT_H

/*
 * Everything the historical chemistry implementation expects from the
 * standalone host is quarantined here so the directory can be detached.
 */
#ifdef HD
#define N_SPECIES 14
#else
#define N_SPECIES 9
#endif
#define TGCHEM_NUM_ABUNDANCES N_SPECIES

#ifndef METAL
#define k_MAX 60
#else
#define k_MAX 200
#endif

extern const int ThisTask;
extern const int NTask;

#include "chemistry.h"

#endif
