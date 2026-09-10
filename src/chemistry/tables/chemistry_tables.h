#pragma once

#include <stdio.h>
#include "chemistry_core.h"

#define N_T 10000
#define T_start 0.1
#define T_end 1.0e6

typedef struct ChemistryRateTables {
  double dln_temperature;
  double formation_energy_kj_mol[15];

  double temperature[N_T];
  double forward[CHEM_RATE_COUNT][N_T];
  double inverse[CHEM_RATE_COUNT][N_T];

  double collisional_ionization_HI[N_T];
  double collisional_ionization_HeI[N_T];
  double collisional_ionization_HeII[N_T];
  double collisional_ionization_HeI_singlet[N_T];
  double recombination_HII[N_T];
  double bremsstrahlung[N_T];
  double H2_CIE[N_T];

  double collisional_excitation_HI[N_T];
  double collisional_excitation_HeI[N_T];
  double collisional_excitation_HeII[N_T];
  double recombination_HeII_radiative[N_T];
  double recombination_HeII_dielectronic[N_T];
  double recombination_HeIII[N_T];
} ChemistryRateTables;

FILE *open_chemistry_table(const char *table_root, const char *filename);
