#pragma once

typedef struct ChemistryCoolingTables ChemistryCoolingTables;

#ifdef MOLECULAR_LINE_COOLING
int molecular_line_cooling_create(const char *table_root,
                                  ChemistryCoolingTables **tables);
void molecular_line_cooling_destroy(ChemistryCoolingTables *tables);
void evaluate_molecular_line_cooling(const ChemistryCoolingTables *tables,
                                     double temperature,
                                     double radiation_temperature,
                                     double length, double n_nuclei,
                                     const double *abundance,
                                     double *cooling_CO, double *cooling_OH,
                                     double *cooling_H2O);
#endif
