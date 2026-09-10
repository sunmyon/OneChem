#pragma once

#ifdef SHIELD_LW
double SelfShield(double column_density, double temperature, double doppler_b);
double EachShield(double column_density, double alpha, double delta, double eta);
#endif
