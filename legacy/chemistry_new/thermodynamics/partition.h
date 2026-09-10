#ifndef PARTITION_H
#define PARTITION_H

typedef double (*PartitionFunction)(double temperature);

double partition_H2_Rossi1985(double temperature);
double partition_H2_Irwin1987(double temperature);
double partition_H2_EBM1989(double temperature);
double partition_H2_EBM1989_fitted_by_Matsukoba(double temperature);
double partition_H(double temperature);
double pf_z_H2p(double temperature);

void check_partition_functions_H2(void);

#endif
