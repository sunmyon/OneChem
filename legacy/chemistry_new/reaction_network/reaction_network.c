#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "chemistry_tables.h"
#include "reaction_network.h"
#include "species.h"

static ChemReaction *reactions;
static size_t reaction_count;

typedef struct {
  size_t num_reactant_reactions;
  size_t *reactant_reactions;
  size_t num_affected_reactions;
  size_t *affected_reactions;
} ChemSpeciesReactionIndex;

static ChemSpeciesReactionIndex species_reaction_index[sEND];

static void free_species_reaction_index(void)
{
  for(int species = 0; species < sEND; species++) {
    free(species_reaction_index[species].reactant_reactions);
    free(species_reaction_index[species].affected_reactions);
    species_reaction_index[species] = (ChemSpeciesReactionIndex){0};
  }
}

static int build_species_reaction_index(void)
{
  size_t reactant_counts[sEND];
  size_t affected_counts[sEND];
  memset(reactant_counts, 0, sizeof(reactant_counts));
  memset(affected_counts, 0, sizeof(affected_counts));

  for(size_t reaction_index = 0; reaction_index < reaction_count; reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      if(term->reactant_order > 0)
        reactant_counts[term->species]++;
      if(term->net_change != 0)
        affected_counts[term->species]++;
    }
  }

  for(int species = 0; species < sEND; species++) {
    ChemSpeciesReactionIndex *index = &species_reaction_index[species];
    if(reactant_counts[species] > 0) {
      index->reactant_reactions = malloc(reactant_counts[species] * sizeof(*index->reactant_reactions));
      if(index->reactant_reactions == NULL)
        goto failure;
    }
    if(affected_counts[species] > 0) {
      index->affected_reactions = malloc(affected_counts[species] * sizeof(*index->affected_reactions));
      if(index->affected_reactions == NULL)
        goto failure;
    }
  }

  for(size_t reaction_index = 0; reaction_index < reaction_count; reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      ChemSpeciesReactionIndex *index = &species_reaction_index[term->species];
      if(term->reactant_order > 0)
        index->reactant_reactions[index->num_reactant_reactions++] = reaction_index;
      if(term->net_change != 0)
        index->affected_reactions[index->num_affected_reactions++] = reaction_index;
    }
  }
  return 1;

failure:
  free_species_reaction_index();
  return 0;
}

static char *trim(char *text)
{
  while(isspace((unsigned char)*text))
    text++;

  char *end = text + strlen(text);
  while(end > text && isspace((unsigned char)end[-1]))
    end--;
  *end = '\0';
  return text;
}

static int compiled_option_enabled(const char *name)
{
  if(strcmp(name, "always") == 0)
    return 1;
#ifdef HELIUM
  if(strcmp(name, "HELIUM") == 0) return 1;
#endif
#ifdef H2_PLUS
  if(strcmp(name, "H2_PLUS") == 0) return 1;
#endif
#ifdef HD
  if(strcmp(name, "HD") == 0) return 1;
#endif
#if defined(HD) && defined(HD_PLUS)
  if(strcmp(name, "HD_PLUS") == 0) return 1;
#endif
#if defined(HD) && defined(DMINUS)
  if(strcmp(name, "DMINUS") == 0) return 1;
#endif
#ifdef IONIZATION_PHOTON
  if(strcmp(name, "IONIZATION_PHOTON") == 0) return 1;
#endif
#ifdef LYMANWERNER
  if(strcmp(name, "LYMANWERNER") == 0) return 1;
#endif
#ifdef H2_FORMATION_ON_DUST
  if(strcmp(name, "H2_FORMATION_ON_DUST") == 0) return 1;
#endif
#ifdef THREEBODY
  if(strcmp(name, "THREEBODY") == 0) return 1;
#endif
#ifdef INVERSE_REACTIONS
  if(strcmp(name, "INVERSE_REACTIONS") == 0) return 1;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  if(strcmp(name, "CARBON_OXYGEN_CHEMISTRY") == 0) return 1;
#endif

  static const char *known_options[] = {
    "HELIUM", "H2_PLUS", "HD", "HD_PLUS", "DMINUS", "IONIZATION_PHOTON",
    "LYMANWERNER", "H2_FORMATION_ON_DUST", "THREEBODY", "INVERSE_REACTIONS",
    "CARBON_OXYGEN_CHEMISTRY"
  };
  for(size_t i = 0; i < sizeof(known_options) / sizeof(known_options[0]); i++)
    if(strcmp(name, known_options[i]) == 0)
      return 0;
  return -1;
}

static int requirements_enabled(char *requirements)
{
  int all_enabled = 1;
  char *saveptr = NULL;
  for(char *name = strtok_r(requirements, ",", &saveptr); name != NULL; name = strtok_r(NULL, ",", &saveptr)) {
    int enabled = compiled_option_enabled(trim(name));
    if(enabled < 0)
      return -1;
    if(!enabled)
      all_enabled = 0;
  }
  return all_enabled;
}

static ChemReactionTerm *find_or_add_term(ChemReaction *reaction, enum SPECIES species)
{
  for(int i = 0; i < reaction->num_terms; i++)
    if(reaction->terms[i].species == species)
      return &reaction->terms[i];

  if(reaction->num_terms == CHEM_MAX_REACTION_TERMS)
    return NULL;

  ChemReactionTerm *term = &reaction->terms[reaction->num_terms++];
  term->species = species;
  term->reactant_order = 0;
  term->net_change = 0;
  return term;
}

static int parse_reaction_side(char *side, int reactant_side, ChemReaction *reaction)
{
  char *saveptr = NULL;
  for(char *token = strtok_r(side, "+", &saveptr); token != NULL; token = strtok_r(NULL, "+", &saveptr)) {
    token = trim(token);
    long coefficient = 1;
    char *species_name = token;

    if(isdigit((unsigned char)*token)) {
      char *end;
      errno = 0;
      coefficient = strtol(token, &end, 10);
      while(isspace((unsigned char)*end)) end++;
      if(errno != 0 || coefficient <= 0 || coefficient > SCHAR_MAX || *end != '*')
        return 0;
      species_name = trim(end + 1);
    }

    long reactant_order = coefficient;
    char *order_separator = strchr(species_name, '@');
    if(order_separator != NULL) {
      if(!reactant_side)
        return 0;
      *order_separator = '\0';
      char *order_end;
      errno = 0;
      reactant_order = strtol(trim(order_separator + 1), &order_end, 10);
      if(errno != 0 || *trim(order_end) != '\0' || reactant_order <= 0 || reactant_order > UCHAR_MAX)
        return 0;
      species_name = trim(species_name);
    }

    enum SPECIES species;
    if(!chemistry_species_from_name(species_name, &species))
      return 0;

    ChemReactionTerm *term = find_or_add_term(reaction, species);
    if(term == NULL)
      return 0;

    if(reactant_side) {
      if(term->reactant_order + reactant_order > UCHAR_MAX || term->net_change - coefficient < SCHAR_MIN)
        return 0;
      term->reactant_order += reactant_order;
      term->net_change -= coefficient;
    } else {
      if(term->net_change + coefficient > SCHAR_MAX)
        return 0;
      term->net_change += coefficient;
    }
  }
  return 1;
}

static int parse_rate_source(const char *text, ChemRateSource *source)
{
  if(strcmp(text, "forward") == 0) *source = CHEM_RATE_FORWARD;
  else if(strcmp(text, "inverse") == 0) *source = CHEM_RATE_INVERSE;
  else if(strcmp(text, "photo") == 0) *source = CHEM_RATE_PHOTO;
  else return 0;
  return 1;
}

static int parse_modes(char *text, unsigned char *modes)
{
  *modes = 0;
  char *saveptr = NULL;
  for(char *mode = strtok_r(text, ",", &saveptr); mode != NULL; mode = strtok_r(NULL, ",", &saveptr)) {
    mode = trim(mode);
    if(strcmp(mode, "implicit") == 0) *modes |= CHEM_REACTION_IMPLICIT;
    else if(strcmp(mode, "explicit") == 0) *modes |= CHEM_REACTION_EXPLICIT;
    else return 0;
  }
  return *modes != 0;
}

static int append_reaction(const ChemReaction *reaction)
{
  ChemReaction *new_reactions = realloc(reactions, (reaction_count + 1) * sizeof(*reactions));
  if(new_reactions == NULL)
    return 0;
  reactions = new_reactions;
  reactions[reaction_count++] = *reaction;
  return 1;
}

static int reaction_conserves_species_properties(const ChemReaction *reaction)
{
  int element_change[CHEM_ELEMENT_COUNT] = {0};
  int charge_change = 0;

  for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
    const ChemReactionTerm *term = &reaction->terms[term_index];
    const ChemSpeciesDefinition *definition = chemistry_species_definition(term->species);
    if(definition == NULL)
      return 0;
    charge_change += term->net_change * definition->charge;
    for(int element = 0; element < CHEM_ELEMENT_COUNT; element++)
      element_change[element] += term->net_change * definition->nuclei[element];
  }

  if(charge_change != 0)
    return 0;
  for(int element = 0; element < CHEM_ELEMENT_COUNT; element++)
    if(element_change[element] != 0)
      return 0;
  return 1;
}

int chemistry_network_init(void)
{
  if(reactions != NULL)
    return 1;

  FILE *table = open_chemistry_table("reaction_network.tsv");
  if(table == NULL)
    return 0;

  char line[2048];
  int line_number = 0;
  while(fgets(line, sizeof(line), table) != NULL) {
    line_number++;
    char *record = trim(line);
    if(*record == '\0' || *record == '#')
      continue;

    char *fields[6];
    for(int field = 0; field < 6; field++) {
      fields[field] = record;
      char *tab = strchr(record, '\t');
      if(field < 5) {
        if(tab == NULL) goto malformed;
        *tab = '\0';
        record = tab + 1;
      } else if(tab != NULL) {
        goto malformed;
      }
      fields[field] = trim(fields[field]);
    }

    if(strcmp(fields[0], "rate_source") == 0)
      continue;

    int enabled = requirements_enabled(fields[3]);
    if(enabled < 0) {
      fprintf(stderr, "Unknown requirement in reaction_network.tsv:%d.\n", line_number);
      goto failure;
    }
    if(!enabled)
      continue;

    ChemReaction reaction = {0};
    if(!parse_rate_source(fields[0], &reaction.rate_source))
      goto malformed;

    char *rate_end;
    errno = 0;
    long rate_index = strtol(fields[1], &rate_end, 10);
    if(errno != 0 || *trim(rate_end) != '\0' || rate_index < 0 || rate_index >= CHEM_RATE_COUNT || rate_index > USHRT_MAX)
      goto malformed;
    reaction.rate_index = rate_index;

    if(!parse_modes(fields[2], &reaction.modes) ||
       !parse_reaction_side(fields[4], 1, &reaction) ||
       !parse_reaction_side(fields[5], 0, &reaction) ||
       reaction.num_terms == 0)
      goto malformed;
    if(!reaction_conserves_species_properties(&reaction)) {
      fprintf(stderr, "Reaction does not conserve elements or charge in reaction_network.tsv:%d.\n", line_number);
      goto failure;
    }
    if(!append_reaction(&reaction))
      goto failure;
    continue;

malformed:
    fprintf(stderr, "Malformed reaction in reaction_network.tsv:%d.\n", line_number);
failure:
    fclose(table);
    free(reactions);
    reactions = NULL;
    reaction_count = 0;
    return 0;
  }

  fclose(table);
  if(reaction_count == 0 || !build_species_reaction_index()) {
    free_species_reaction_index();
    free(reactions);
    reactions = NULL;
    reaction_count = 0;
    return 0;
  }
  return 1;
}

static double integer_power(double value, unsigned char exponent)
{
  double result = 1.0;
  while(exponent-- > 0)
    result *= value;
  return result;
}

static double reaction_coefficient(const ChemReaction *reaction, const double *forward_rates, const double *inverse_rates)
{
  return reaction->rate_source == CHEM_RATE_INVERSE ?
      inverse_rates[reaction->rate_index] : forward_rates[reaction->rate_index];
}

size_t chemistry_network_reaction_count(void)
{
  return reaction_count;
}

int chemistry_network_intermediate_equilibrium_abundance(enum SPECIES species, unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *equilibrium_abundance, double *relaxation_time)
{
  double production = 0.0;
  double destruction = 0.0;

  const ChemSpeciesReactionIndex *index = &species_reaction_index[species];
  for(size_t affected_index = 0; affected_index < index->num_affected_reactions; affected_index++) {
    size_t reaction_index = index->affected_reactions[affected_index];
    const ChemReaction *reaction = &reactions[reaction_index];
    const ChemReactionTerm *equilibrium_term = NULL;

    if(!(reaction->modes & mode))
      continue;

    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      if(reaction->terms[term_index].species == species) {
        equilibrium_term = &reaction->terms[term_index];
        break;
      }
    }

    if(equilibrium_term == NULL || equilibrium_term->net_change == 0)
      continue;

    if((equilibrium_term->net_change > 0 && equilibrium_term->reactant_order != 0) ||
       (equilibrium_term->net_change < 0 && equilibrium_term->reactant_order != 1))
      return 0;

    double rate = reaction_coefficient(reaction, forward_rates, inverse_rates);
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      if(term->species == species)
        continue;
      rate *= integer_power(abund[term->species] * nH, term->reactant_order);
    }

    if(equilibrium_term->net_change > 0) {
      production += equilibrium_term->net_change * rate / nH;
    } else {
      destruction += -equilibrium_term->net_change * rate;
    }
  }

  if(!(destruction > 0.0) || !isfinite(production) || !isfinite(destruction))
    return 0;

  *equilibrium_abundance = production / destruction;
  if(relaxation_time != NULL)
    *relaxation_time = 1.0 / destruction;
  return isfinite(*equilibrium_abundance) && *equilibrium_abundance >= 0.0;
}

void chemistry_network_rhs(unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *rhs)
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    rhs[species] = 0.0;
  }

  for(size_t reaction_index = 0; reaction_index < chemistry_network_reaction_count(); reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & mode))
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

void chemistry_network_jacobian_column(unsigned char mode, int differentiated_species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *column)
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    column[species] = 0.0;
  }

  const ChemSpeciesReactionIndex *index = &species_reaction_index[differentiated_species];
  for(size_t reactant_index = 0; reactant_index < index->num_reactant_reactions; reactant_index++) {
    size_t reaction_index = index->reactant_reactions[reactant_index];
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & mode))
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

void chemistry_network_jacobian(unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double jacobian[N_SP][N_SP])
{
  double number_density[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    memset(jacobian[species], 0, N_SP * sizeof(jacobian[species][0]));
  }

  for(size_t reaction_index = 0; reaction_index < chemistry_network_reaction_count(); reaction_index++) {
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & mode))
      continue;

    double coefficient = reaction_coefficient(reaction, forward_rates, inverse_rates);
    for(int differentiated_index = 0; differentiated_index < reaction->num_terms; differentiated_index++) {
      const ChemReactionTerm *differentiated = &reaction->terms[differentiated_index];
      if(differentiated->reactant_order == 0)
        continue;

      double derivative = coefficient * differentiated->reactant_order;
      for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
        const ChemReactionTerm *term = &reaction->terms[term_index];
        unsigned char exponent = term->reactant_order;
        if(term->species == differentiated->species)
          exponent--;
        derivative *= integer_power(number_density[term->species], exponent);
      }

      for(int affected_index = 0; affected_index < reaction->num_terms; affected_index++) {
        const ChemReactionTerm *affected = &reaction->terms[affected_index];
        jacobian[differentiated->species][affected->species] += affected->net_change * derivative;
      }
    }
  }
}

void chemistry_network_intermediate_system(unsigned char mode, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, const enum SPECIES *intermediate_species, int n_intermediate, double *rhs, double *matrix)
{
  double number_density[sEND];
  int intermediate_position[sEND];
  for(int species = 0; species < sEND; species++) {
    number_density[species] = abund[species] * nH;
    intermediate_position[species] = -1;
  }
  for(int intermediate = 0; intermediate < n_intermediate; intermediate++) {
    intermediate_position[intermediate_species[intermediate]] = intermediate;
    rhs[intermediate] = 0.0;
  }
  memset(matrix, 0, n_intermediate * n_intermediate * sizeof(*matrix));

  for(int row = 0; row < n_intermediate; row++) {
    enum SPECIES row_species = intermediate_species[row];
    const ChemSpeciesReactionIndex *index = &species_reaction_index[row_species];
    for(size_t affected_index = 0; affected_index < index->num_affected_reactions; affected_index++) {
      const ChemReaction *reaction = &reactions[index->affected_reactions[affected_index]];
      if(!(reaction->modes & mode))
        continue;

      const ChemReactionTerm *row_term = NULL;
      double coefficient = reaction_coefficient(reaction, forward_rates, inverse_rates);
      double rate = coefficient;
      for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
        const ChemReactionTerm *term = &reaction->terms[term_index];
        rate *= integer_power(number_density[term->species], term->reactant_order);
        if(term->species == row_species)
          row_term = term;
      }
      if(row_term == NULL || row_term->net_change == 0)
        continue;
      rhs[row] += row_term->net_change * rate / nH;

      for(int differentiated_index = 0; differentiated_index < reaction->num_terms; differentiated_index++) {
        const ChemReactionTerm *differentiated = &reaction->terms[differentiated_index];
        int column = intermediate_position[differentiated->species];
        if(differentiated->reactant_order == 0 || column < 0)
          continue;

        double derivative = coefficient * differentiated->reactant_order;
        for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
          const ChemReactionTerm *term = &reaction->terms[term_index];
          unsigned char exponent = term->reactant_order;
          if(term->species == differentiated->species)
            exponent--;
          derivative *= integer_power(number_density[term->species], exponent);
        }
        matrix[column * n_intermediate + row] += row_term->net_change * derivative;
      }
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

void chemistry_network_species_production_destruction(enum SPECIES species, double nH, const double *abund, const double *forward_rates, const double *inverse_rates, double *production, double *destruction)
{
  *production = 0.0;
  *destruction = 0.0;

  const ChemSpeciesReactionIndex *index = &species_reaction_index[species];
  for(size_t affected_index = 0; affected_index < index->num_affected_reactions; affected_index++) {
    size_t reaction_index = index->affected_reactions[affected_index];
    const ChemReaction *reaction = &reactions[reaction_index];
    if(!(reaction->modes & CHEM_REACTION_EXPLICIT))
      continue;

    const ChemReactionTerm *affected_term = NULL;
    double rate = reaction_coefficient(reaction, forward_rates, inverse_rates);
    for(int term_index = 0; term_index < reaction->num_terms; term_index++) {
      const ChemReactionTerm *term = &reaction->terms[term_index];
      rate *= integer_power(abund[term->species] * nH, term->reactant_order);
      if(term->species == species)
        affected_term = term;
    }
    rate /= nH;

    if(affected_term->net_change > 0) {
      *production += affected_term->net_change * rate;
    } else if(affected_term->net_change < 0 && abund[species] > 0.0) {
      *destruction += -affected_term->net_change * rate / abund[species];
    }
  }
}
