#pragma once

#ifdef MOLECULAR_LINE_COOLING
int molecular_line_cooling_init(void);
void evaluate_molecular_line_cooling(double temperature, double radiation_temperature, double length, double n_nuclei, const double *abundance, double *cooling_CO, double *cooling_OH, double *cooling_H2O);
#endif
