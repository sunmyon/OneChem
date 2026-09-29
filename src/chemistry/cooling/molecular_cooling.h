#pragma once

#include "chemistry_build_config.h"


typedef struct ChemistryCoolingTables ChemistryCoolingTables;

#if CHEMISTRY_MOLECULAR_LINE_COOLING


void evaluate_molecular_line_cooling(const ChemistryCoolingTables *tables,
                                     double temperature,
                                     double radiation_temperature,
                                     double length, double n_nuclei,
                                     const double *abundance,
                                     double *cooling_CO, double *cooling_OH,
                                     double *cooling_H2O);
#endif
