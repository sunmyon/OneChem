#ifndef CHEMISTRY_MODELS_H
#define CHEMISTRY_MODELS_H

#include "H2cooling.h"
#include "partition.h"

typedef struct {
  H2CoolingFunction h2_cooling;
  PartitionFunction h2_partition;
#ifdef H2_OPACITY
  LineEscapeFunction h2_line_escape;
#ifdef HD
  LineEscapeFunction hd_line_escape;
#endif
#endif
} ChemistryModels;

extern const ChemistryModels chemistry_models;

#endif
