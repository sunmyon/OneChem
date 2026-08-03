#include <stddef.h>

#include "../allvars.h"
#include "reaction_network.h"

#define TERM(species_, order_, change_) {species_, order_, change_}
#define REACTION(source_, rate_, modes_, count_, ...) {source_, rate_, modes_, count_, {__VA_ARGS__}}
#define FORWARD(rate_, count_, ...) REACTION(CHEM_RATE_FORWARD, rate_, CHEM_REACTION_IMPLICIT | CHEM_REACTION_EXPLICIT, count_, __VA_ARGS__)
#define INVERSE(rate_, count_, ...) REACTION(CHEM_RATE_INVERSE, rate_, CHEM_REACTION_IMPLICIT | CHEM_REACTION_EXPLICIT, count_, __VA_ARGS__)

static const ChemReaction reactions[] = {
  FORWARD(1, 3, TERM(sHI, 1, -1), TERM(sELEC, 1, 1), TERM(sHII, 0, 1)),
  FORWARD(2, 3, TERM(sHII, 1, -1), TERM(sELEC, 1, -1), TERM(sHI, 0, 1)),
#ifdef HELIUM
  FORWARD(3, 3, TERM(sHeI, 1, -1), TERM(sELEC, 1, 1), TERM(sHeII, 0, 1)),
  FORWARD(4, 3, TERM(sHeII, 1, -1), TERM(sELEC, 1, -1), TERM(sHeI, 0, 1)),
  FORWARD(5, 3, TERM(sHeII, 1, -1), TERM(sELEC, 1, 1), TERM(sHeIII, 0, 1)),
  FORWARD(6, 3, TERM(sHeIII, 1, -1), TERM(sELEC, 1, -1), TERM(sHeII, 0, 1)),
#endif
  FORWARD(7, 3, TERM(sHI, 1, -1), TERM(sELEC, 1, -1), TERM(sHM, 0, 1)),
  FORWARD(8, 4, TERM(sHM, 1, -1), TERM(sHI, 1, -1), TERM(sH2I, 0, 1), TERM(sELEC, 0, 1)),
#ifdef H2_PLUS
  FORWARD(9, 3, TERM(sHI, 1, -1), TERM(sHII, 1, -1), TERM(sH2II, 0, 1)),
  FORWARD(10, 4, TERM(sHI, 1, -1), TERM(sH2II, 1, -1), TERM(sH2I, 0, 1), TERM(sHII, 0, 1)),
#endif
  FORWARD(11, 2, TERM(sH2I, 1, -1), TERM(sHI, 1, 2)),
#ifdef H2_PLUS
  FORWARD(12, 4, TERM(sH2I, 1, -1), TERM(sHII, 1, -1), TERM(sH2II, 0, 1), TERM(sHI, 0, 1)),
#endif
  FORWARD(13, 3, TERM(sH2I, 1, -1), TERM(sELEC, 1, 0), TERM(sHI, 0, 2)),
#ifdef HELIUM
  FORWARD(14, 4, TERM(sHeII, 1, -1), TERM(sHI, 1, -1), TERM(sHeI, 0, 1), TERM(sHII, 0, 1)),
  FORWARD(15, 4, TERM(sHeI, 1, -1), TERM(sHII, 1, -1), TERM(sHeII, 0, 1), TERM(sHI, 0, 1)),
#endif
  FORWARD(16, 3, TERM(sHM, 1, -1), TERM(sELEC, 1, 1), TERM(sHI, 0, 1)),
  FORWARD(17, 3, TERM(sHM, 1, -1), TERM(sHI, 1, 1), TERM(sELEC, 0, 1)),
  FORWARD(18, 3, TERM(sHM, 1, -1), TERM(sHII, 1, -1), TERM(sHI, 0, 2)),
#ifdef H2_PLUS
  FORWARD(19, 4, TERM(sHM, 1, -1), TERM(sHII, 1, -1), TERM(sH2II, 0, 1), TERM(sELEC, 0, 1)),
  FORWARD(20, 3, TERM(sH2II, 1, -1), TERM(sELEC, 1, -1), TERM(sHI, 0, 2)),
  FORWARD(21, 4, TERM(sH2II, 1, -1), TERM(sHM, 1, -1), TERM(sHI, 0, 1), TERM(sH2I, 0, 1)),
#endif
#ifdef IONIZATION_PHOTON
  FORWARD(22, 3, TERM(sHI, 1, -1), TERM(sHII, 0, 1), TERM(sELEC, 0, 1)),
#endif
#ifdef LYMANWERNER
  FORWARD(25, 3, TERM(sHM, 1, -1), TERM(sHI, 0, 1), TERM(sELEC, 0, 1)),
  FORWARD(29, 2, TERM(sH2I, 1, -1), TERM(sHI, 0, 2)),
#endif
  FORWARD(30, 3, TERM(sHI, 2, -1), TERM(sELEC, 0, 1), TERM(sHII, 0, 1)),
#ifdef HD
  FORWARD(32, 4, TERM(sDI, 1, -1), TERM(sHII, 1, -1), TERM(sDII, 0, 1), TERM(sHI, 0, 1)),
  FORWARD(33, 4, TERM(sDII, 1, -1), TERM(sHI, 1, -1), TERM(sDI, 0, 1), TERM(sHII, 0, 1)),
  FORWARD(35, 4, TERM(sDI, 1, -1), TERM(sH2I, 1, -1), TERM(sHI, 0, 1), TERM(sHDI, 0, 1)),
  FORWARD(37, 4, TERM(sDII, 1, -1), TERM(sH2I, 1, -1), TERM(sHII, 0, 1), TERM(sHDI, 0, 1)),
  FORWARD(38, 4, TERM(sHDI, 1, -1), TERM(sHI, 1, -1), TERM(sH2I, 0, 1), TERM(sDI, 0, 1)),
  FORWARD(39, 4, TERM(sHDI, 1, -1), TERM(sHII, 1, -1), TERM(sH2I, 0, 1), TERM(sDII, 0, 1)),
#ifdef HD_PLUS
  FORWARD(36, 4, TERM(sHDII, 1, -1), TERM(sHI, 1, -1), TERM(sHII, 0, 1), TERM(sHDI, 0, 1)),
  FORWARD(40, 3, TERM(sDI, 1, -1), TERM(sHII, 1, -1), TERM(sHDII, 0, 1)),
  FORWARD(41, 3, TERM(sDII, 1, -1), TERM(sHI, 1, -1), TERM(sHDII, 0, 1)),
  FORWARD(42, 4, TERM(sHDII, 1, -1), TERM(sELEC, 1, -1), TERM(sHI, 0, 1), TERM(sDI, 0, 1)),
#endif
#ifdef DMINUS
  FORWARD(43, 3, TERM(sDI, 1, -1), TERM(sELEC, 1, -1), TERM(sDM, 0, 1)),
  FORWARD(44, 3, TERM(sDII, 1, -1), TERM(sDM, 1, -1), TERM(sDI, 0, 2)),
  FORWARD(45, 4, TERM(sHII, 1, -1), TERM(sDM, 1, -1), TERM(sHI, 0, 1), TERM(sDI, 0, 1)),
  FORWARD(46, 4, TERM(sHM, 1, -1), TERM(sDI, 1, -1), TERM(sHI, 0, 1), TERM(sDM, 0, 1)),
  FORWARD(47, 4, TERM(sDM, 1, -1), TERM(sHI, 1, -1), TERM(sHM, 0, 1), TERM(sDI, 0, 1)),
  FORWARD(48, 4, TERM(sDM, 1, -1), TERM(sHI, 1, -1), TERM(sHDI, 0, 1), TERM(sELEC, 0, 1)),
#endif
#ifdef LYMANWERNER
  FORWARD(53, 3, TERM(sHDI, 1, -1), TERM(sHI, 0, 1), TERM(sDI, 0, 1)),
#endif
#endif
#ifdef H2_FORMATION_ON_DUST
  FORWARD(51, 2, TERM(sHI, 1, -2), TERM(sH2I, 0, 1)),
#endif
#ifdef THREEBODY
  FORWARD(55, 2, TERM(sHI, 3, -2), TERM(sH2I, 0, 1)),
  FORWARD(57, 2, TERM(sHI, 2, -2), TERM(sH2I, 1, 1)),
  FORWARD(58, 2, TERM(sH2I, 2, -1), TERM(sHI, 0, 2)),
#endif
#ifdef INVERSE_REACTIONS
  INVERSE(1, 3, TERM(sHII, 1, -1), TERM(sELEC, 2, -1), TERM(sHI, 0, 1)),
  INVERSE(2, 3, TERM(sHI, 1, -1), TERM(sHII, 0, 1), TERM(sELEC, 0, 1)),
  INVERSE(7, 3, TERM(sHM, 1, -1), TERM(sHI, 0, 1), TERM(sELEC, 0, 1)),
  INVERSE(8, 4, TERM(sH2I, 1, -1), TERM(sELEC, 1, -1), TERM(sHI, 0, 1), TERM(sHM, 0, 1)),
  INVERSE(13, 3, TERM(sHI, 2, -2), TERM(sELEC, 1, 0), TERM(sH2I, 0, 1)),
  INVERSE(16, 3, TERM(sHI, 1, -1), TERM(sELEC, 2, -1), TERM(sHM, 0, 1)),
  INVERSE(17, 3, TERM(sHI, 2, -1), TERM(sELEC, 1, -1), TERM(sHM, 0, 1)),
  INVERSE(18, 3, TERM(sHI, 2, -2), TERM(sHII, 0, 1), TERM(sHM, 0, 1)),
  INVERSE(30, 3, TERM(sHI, 1, 1), TERM(sHII, 1, -1), TERM(sELEC, 1, -1)),
#endif
};

static double integer_power(double value, unsigned char exponent)
{
  double result = 1.0;
  while(exponent-- > 0)
    result *= value;
  return result;
}

static double reaction_coefficient(const ChemReaction *reaction, const double *forward_rates, const double *inverse_rates)
{
  return reaction->rate_source == CHEM_RATE_FORWARD ?
      forward_rates[reaction->rate_index] : inverse_rates[reaction->rate_index];
}

size_t chemistry_network_reaction_count(void)
{
  return sizeof(reactions) / sizeof(reactions[0]);
}

void chemistry_network_rhs(double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *rhs)
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    rhs[species] = 0.0;
  }

  for(size_t reaction_index = 0; reaction_index < chemistry_network_reaction_count(); reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & CHEM_REACTION_IMPLICIT))
      continue;
    double rate = reaction_coefficient(reaction, forward_rates, inverse_rates);

    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      rate *= integer_power(number_density[term->species], term->reactant_order);
    }

    rate /= nH;
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      rhs[term->species] += term->net_change * rate;
    }
  }
}

void chemistry_network_jacobian_column(int differentiated_species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *column)
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    column[species] = 0.0;
  }

  for(size_t reaction_index = 0; reaction_index < chemistry_network_reaction_count(); reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & CHEM_REACTION_IMPLICIT))
      continue;
    unsigned char differentiated_order = 0;

    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      if(term->species == differentiated_species)
        differentiated_order = term->reactant_order;
    }
    if(differentiated_order == 0)
      continue;

    double derivative = reaction_coefficient(reaction, forward_rates, inverse_rates) * differentiated_order;
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      unsigned char exponent = term->reactant_order;
      if(term->species == differentiated_species)
        exponent--;
      derivative *= integer_power(number_density[term->species], exponent);
    }

    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      column[term->species] += term->net_change * derivative;
    }
  }
}

void chemistry_network_production_destruction(double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction)
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    production[species] = 0.0;
    destruction[species] = 0.0;
  }

  for(size_t reaction_index = 0; reaction_index < chemistry_network_reaction_count(); reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & CHEM_REACTION_EXPLICIT))
      continue;

    double rate = reaction_coefficient(reaction, forward_rates, inverse_rates);
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      rate *= integer_power(number_density[term->species], term->reactant_order);
    }
    rate /= nH;

    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      if(term->net_change > 0) {
        production[term->species] += term->net_change * rate;
      } else if(term->net_change < 0 && abund[term->species] > 0.0) {
        destruction[term->species] += -term->net_change * rate / abund[term->species];
      }
    }
  }
}
