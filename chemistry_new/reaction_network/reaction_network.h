#pragma once

#include <stddef.h>
#include "chemistry.h"

extern int N_SP, N_RE;

typedef enum {
  CHEM_RATE_FORWARD,
  CHEM_RATE_INVERSE
} ChemRateSource;

typedef enum {
  CHEM_REACTION_IMPLICIT = 1 << 0,
  CHEM_REACTION_EXPLICIT = 1 << 1
} ChemReactionMode;

typedef struct {
  enum SPECIES species;
  unsigned char reactant_order;
  signed char net_change;
} ChemReactionTerm;

#define CHEM_MAX_REACTION_TERMS 4

typedef struct {
  ChemRateSource rate_source;
  unsigned short rate_index;
  unsigned char modes;
  unsigned char num_terms;
  ChemReactionTerm terms[CHEM_MAX_REACTION_TERMS];
} ChemReaction;

size_t chemistry_network_reaction_count(void);
void chemistry_network_rhs(double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *rhs);
void chemistry_network_jacobian_column(int differentiated_species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *column);
void chemistry_network_production_destruction(double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction);
