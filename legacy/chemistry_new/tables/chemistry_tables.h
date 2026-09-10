#pragma once

#include <stdio.h>
#include "chemistry_types.h"

#define N_T 10000
#define T_start 0.1
#define T_end 1.0e6

extern double T_table[N_T];
extern double k_table[CHEM_RATE_COUNT][N_T];
extern double k_table_inverse[CHEM_RATE_COUNT][N_T];
extern double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T], H2_CIE[N_T];
extern double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];
extern double dln_T;
extern double chem_e[15];

FILE *open_chemistry_table(const char *filename);
