#pragma once

#include "include/chemistry/chemistry_status.h"

typedef struct ChemRunContext ChemRunContext;
typedef struct ChemistryRateTables ChemistryRateTables;

typedef struct {
  /* Directory containing the tables/ subdirectory. */
  const char *table_root;
} ChemistryInitConfig;

/* Transitional eager initializer; replaced by chemistry_context_create. */
ChemistryStatus InitChem(const ChemistryInitConfig *config,
                         ChemRunContext *run_context);
void DestroyChem(ChemRunContext *run_context);
int chem_tabulate(ChemistryRateTables *tables);
