#ifndef CHEMISTRY_OLD_API_H
#define CHEMISTRY_OLD_API_H

int InitChem(void);
int H2line_tabulate(void);
int compute_abundances_d(int mode, int particle, double start,
                         double end, int *implicit_flag,
                         int diagnostics_enabled);
int compute_gamma(int particle);
double compute_temp(int particle);
double compute_freefall_time(double density_cgs);

#endif
