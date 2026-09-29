#pragma once

#include <stddef.h>
#include "chemistry_core.h"


typedef struct ChemReactionNetwork ChemReactionNetwork;

typedef enum {
  CHEM_RATE_FORWARD,
  CHEM_RATE_INVERSE,
  CHEM_RATE_PHOTO
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



size_t chemistry_network_reaction_count(const ChemReactionNetwork *network);
size_t chemistry_network_species_reaction_count(
    const ChemReactionNetwork *network, enum SPECIES species,
    unsigned char mode);
void chemistry_network_rhs(const ChemReactionNetwork *network, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *rhs);
void chemistry_network_rhs_gross(const ChemReactionNetwork *network,
    unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates, double *rhs,
    double *production, double *loss);
void chemistry_network_jacobian_column(const ChemReactionNetwork *network, unsigned char mode, int differentiated_species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *column);
void chemistry_network_jacobian(const ChemReactionNetwork *network, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double jacobian[N_SP][N_SP]);
void chemistry_network_intermediate_system(const ChemReactionNetwork *network, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, const enum SPECIES *intermediate_species, int n_intermediate, double *rhs, double *matrix);
void chemistry_network_production_destruction(const ChemReactionNetwork *network, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction);
void chemistry_network_species_production_destruction(const ChemReactionNetwork *network, enum SPECIES species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction);
void chemistry_network_species_gross_rhs(const ChemReactionNetwork *network,
    enum SPECIES species, unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates,
    double *production, double *loss);
int chemistry_network_intermediate_equilibrium_abundance(const ChemReactionNetwork *network, enum SPECIES species, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *equilibrium_abundance, double *relaxation_time);
