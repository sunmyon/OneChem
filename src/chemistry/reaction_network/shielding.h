#pragma once

#include "chemistry_build_config.h"

#ifdef CHEMISTRY_LW_SELF_SHIELDING
double SelfShield(double column_density, double temperature, double doppler_b);
double EachShield(double column_density, double alpha, double delta, double eta);
#endif
