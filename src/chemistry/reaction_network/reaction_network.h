#pragma once

#include <stddef.h>
#include <memory>
#include <vector>
#include "chemistry_model.hpp"
#include "chemistry_core.h"

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

typedef struct {
  size_t num_reactant_reactions;
  std::unique_ptr<size_t[]> reactant_reactions;
  size_t num_affected_reactions;
  std::unique_ptr<size_t[]> affected_reactions;
} ChemSpeciesReactionIndex;

struct ChemReactionNetwork {
public:
  ChemReactionNetwork() = default;
  ChemReactionNetwork(const ChemReactionNetwork &) = delete;
  ChemReactionNetwork &operator=(const ChemReactionNetwork &) = delete;
  ~ChemReactionNetwork();

  static chemistry::ResourceOwner<ChemReactionNetwork> create(const char *table_root) noexcept;

  size_t size() const;

  size_t species_reaction_count(enum SPECIES species,
    unsigned char mode) const;

  int intermediate_equilibrium_abundance(enum SPECIES species, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *equilibrium_abundance, double *relaxation_time) const;

  void rhs(unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates, double *rhs) const;

  void rhs_gross(unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates, double *rhs,
    double *production, double *loss) const;

  void jacobian_column(unsigned char mode, int differentiated_species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *column) const;

  void jacobian(unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double jacobian[N_SP][N_SP]) const;

  void intermediate_system(unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, const enum SPECIES *intermediate_species, int n_intermediate, double *rhs, double *matrix) const;

  void production_destruction(double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction) const;

  void species_production_destruction(enum SPECIES species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction) const;

  void species_gross_rhs(enum SPECIES species, unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates,
    double *production, double *loss) const;

private:
  static ChemReactionTerm *find_or_add_term(ChemReaction *reaction, enum SPECIES species);

  int build_species_reaction_index();

  static char * trim(char *text);

  static int compiled_option_enabled(const char *name);

  static int requirements_enabled(char *requirements);

  static int parse_reaction_side(char *side, int reactant_side, ChemReaction *reaction);

  static int parse_rate_source(const char *text, ChemRateSource *source);

  static int parse_modes(char *text, unsigned char *modes);

  int append_reaction(const ChemReaction *reaction);

  static int reaction_conserves_species_properties(const ChemReaction *reaction);

  static double integer_power(double value, unsigned char exponent);

  static double reaction_coefficient(const ChemReaction *reaction, const double *forward_rates, const double *inverse_rates);

  void rhs_accumulate(unsigned char mode, double nH, const double *abund,
    const double *forward_rates, const double *inverse_rates, double *rhs,
    double *production, double *loss) const;

  std::vector<ChemReaction> reactions;
  ChemSpeciesReactionIndex species_reaction_index[sEND]{};

};
