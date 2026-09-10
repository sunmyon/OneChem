#pragma once

typedef struct ChemistryRateTables ChemistryRateTables;

double chem_saha_5sp(const ChemistryRateTables *tables, double nH,
                     double temperature, double *abundances);
