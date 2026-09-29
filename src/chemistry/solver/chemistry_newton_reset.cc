#include "chemistry_newton_reset.h"
#include "chemistry_constants.h"
#include "species.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
static void chemistry_newton_apply_carbon_oxygen_corner_reset(
    ChemNewtonResetIndex reset, double metallicity, double *abundance);
#endif

static uint64_t chemistry_random_u64(uint64_t *state)
{
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);

  return value ^ (value >> 31);
}

static double chemistry_random_unit(uint64_t *state)
{
  return (chemistry_random_u64(state) >> 11) * 0x1.0p-53;
}

static double chemistry_random_weight(uint64_t *state)
{
  return 1.e-12 + chemistry_random_unit(state);
}

static void chemistry_random_partition(uint64_t *random_state, int count,
                                       const double *multiplicity, double total,
                                       double *values)
{
  double weights[4];
  double normalization = 0.0;

  for(int index = 0; index < count; index++){
    weights[index] = chemistry_random_weight(random_state);
    normalization += multiplicity[index] * weights[index];
  }

  for(int index = 0; index < count; index++)
    values[index] = total * weights[index] / normalization;
}

static double chemistry_random_fraction(uint64_t *random_state)
{
  double first = chemistry_random_weight(random_state);
  double second = chemistry_random_weight(random_state);

  return first / (first + second);
}

static double chemistry_random_log_value(uint64_t *random_state,
                                         double minimum, double maximum)
{
  return exp(log(minimum) +
             chemistry_random_unit(random_state) * log(maximum / minimum));
}

static int chemistry_randomize_newton_guess(
    double nH, double metallicity, double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation, uint64_t *random_state)
{
  const double unit_multiplicity[4] = {1.0, 1.0, 1.0, 1.0};
  const double oxygen_multiplicity[4] = {1.0, 1.0, 1.0, 2.0};
  double multiplicity[CHEM_MAX_RANDOM_PARTITION_SPECIES];
  double values[CHEM_MAX_RANDOM_PARTITION_SPECIES];

  for(size_t operation_index = 0;
      operation_index < chemistry_random_initialization_count;
      operation_index++){
    const ChemRandomInitialization *operation =
        &chemistry_random_initialization[operation_index];

    if(operation->kind == CHEM_RANDOM_INITIALIZE_LOG_SPECIES){
      abundance[operation->species] = chemistry_random_log_value(
          random_state, operation->minimum, operation->maximum);
      continue;
    }

    const ChemElementConservationDefinition *element = NULL;

    for(size_t element_index = 0;
        element_index < chemistry_element_conservation_count; element_index++){
      if(chemistry_element_conservation[element_index].element ==
         operation->element){
        element = &chemistry_element_conservation[element_index];
        break;
      }
    }

    if(element == NULL)
      continue;

    int count = element->random_partition_species_count;
    double total = element->total_per_hydrogen;

    if(element->metallicity_scaled)
      total *= metallicity;

    for(int index = 0; index < count; index++){
      enum SPECIES species = element->random_partition_species[index];
      multiplicity[index] =
          chemistry_species_definition(static_cast<SPECIES>(species))->nuclei[element->element];
    }

    chemistry_random_partition(random_state, count, multiplicity, total,
                               values);

    for(int index = 0; index < count; index++)
      abundance[element->random_partition_species[index]] = values[index];
  }

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  double carbon_left = A_C_sun * metallicity;
  double oxygen_left = A_O_sun * metallicity;

  abundance[sCO2] = dmin(carbon_left, 0.5 * oxygen_left) *
                    chemistry_random_fraction(random_state);
  carbon_left -= abundance[sCO2];
  oxygen_left -= 2.0 * abundance[sCO2];

  abundance[sCO] = dmin(carbon_left, oxygen_left) *
                   chemistry_random_fraction(random_state);
  carbon_left -= abundance[sCO];
  oxygen_left -= abundance[sCO];

  chemistry_random_partition(random_state, 2, unit_multiplicity, carbon_left,
                             values);
  abundance[sCI] = values[0];
  abundance[sCII] = values[1];

  chemistry_random_partition(random_state, 4, oxygen_multiplicity, oxygen_left,
                             values);
  abundance[sOI] = values[0];
  abundance[sOH] = values[1];
  abundance[sH2O] = values[2];
  abundance[sO2] = values[3];
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
  abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
  abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
  abundance[sH3OII] = abundance[sO2II] = 0.0;
#endif
#else
  (void)metallicity;
#endif

  return (conservation)->update_algebraic(nH, abundance, rates);
}

uint64_t chemistry_newton_random_seed(double nH, double temperature, double dt)
{
  uint64_t seed = UINT64_C(0x6a09e667f3bcc909);
  uint64_t value_bits;

  memcpy(&value_bits, &nH, sizeof(value_bits));
  seed ^= value_bits;
  memcpy(&value_bits, &temperature, sizeof(value_bits));
  seed ^= value_bits + UINT64_C(0x9e3779b97f4a7c15);
  memcpy(&value_bits, &dt, sizeof(value_bits));
  seed ^= value_bits + UINT64_C(0x3c6ef372fe94f82b);

  return seed;
}

int chemistry_find_random_newton_guess(
    double nH, double metallicity, const double *initial_abundance,
    double *abundance, const ChemRates *rates,
    const ChemConservationSystem *conservation, uint64_t *random_state,
    int *iteration, int max_iterations, int *random_attempt_count)
{
  while(*iteration < max_iterations){
    memcpy(abundance, initial_abundance, N_SP * sizeof(*abundance));
    (*iteration)++;
    (*random_attempt_count)++;

    if(chemistry_randomize_newton_guess(nH, metallicity, abundance, rates,
                                        conservation, random_state))
      return 1;
  }

  return 0;
}

void chemistry_newton_apply_deterministic_reset(ChemNewtonResetIndex reset,
                                                double metallicity,
                                                const double *initial_abundance,
                                                double *abundance) {
  if(reset <= CHEM_NEWTON_RESET_NONE || reset > CHEM_NEWTON_RESET_LAST)
    return;

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  memcpy(abundance, initial_abundance, N_SP * sizeof(*abundance));
#endif

  switch(reset){
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  case CHEM_NEWTON_RESET_INITIAL_H2: {
    abundance[sHII] = 0.0;
    abundance[sH2I] = 0.5 * (1.0 - abundance[sHI]);
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDII] = 0.0;
    abundance[sHDI] = D_number_fraction / XH - abundance[sDI];
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T3127: {
    /* Root basin found for the saved Z=1, T=3127 K non-equilibrium
     * intermediate failure; C/O abundances scale with metallicity. */
    abundance[sELEC] = 0.25179906689485865;
    abundance[sHI] = 0.21111795614880233;
    abundance[sHII] = 0.25179091941005283;
    abundance[sH2I] = 0.268540173142787;
    abundance[sHM] = 2.8871723966669333e-7;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 1.6078502270673519e-5;
    abundance[sDII] = 7.9679348332395664e-6;
    abundance[sHDI] = 5.5588260539816532e-6;
#endif
    abundance[sCI] = 2.39504403722442e-5 * metallicity;
    abundance[sCII] = 4.682672122055309e-7 * metallicity;
    abundance[sCO] = 1.424602206249225e-5 * metallicity;
    abundance[sCO2] = 2.77335270353058e-4 * metallicity;
    abundance[sOI] = 7.6656784138147242e-8 * metallicity;
    abundance[sOH] = 1.6190277220737118e-6 * metallicity;
    abundance[sH2O] = 1.6557922775014028e-6 * metallicity;
    abundance[sO2] = 1.3659802238392435e-6 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P01_T32465: {
    /* High-temperature basin found for the saved z=5, Z=0.01,
     * T=32465 K non-equilibrium intermediate failure. */
    abundance[sELEC] = 0.99990537140681268;
    abundance[sHI] = 1.2607213620080047e-4;
    abundance[sHII] = 0.99987311647556509;
    abundance[sH2I] = 1.75285745333539e-15;
    abundance[sHM] = 2.0672061115436203e-12;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 3.7373480562468117e-9;
    abundance[sDII] = 2.9601525809838385e-5;
    abundance[sHDI] = 1.10545931196041e-19;
#endif
    abundance[sCI] = 9.8954371995761087e-5 * metallicity;
    abundance[sCII] = 3.0644731129628816e-14 * metallicity;
    abundance[sCO] = 1.2445469056970178e-4 * metallicity;
    abundance[sCO2] = 9.2590917520743629e-5 * metallicity;
    abundance[sOI] = 2.2723925332810629e-8 * metallicity;
    abundance[sOH] = 5.0314506352684828e-18 * metallicity;
    abundance[sH2O] = 6.3241931343626324e-19 * metallicity;
    abundance[sO2] = 6.1939284267970712e-18 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 1.255623341210446e-29 * metallicity;
    abundance[sCH] = 1.9883148389096883e-11 * metallicity;
    abundance[sCH2] = 3.7372262679145738e-19 * metallicity;
    abundance[sOII] = 1.8420215400720396e-4 * metallicity;
    abundance[sOHII] = 8.1138596456252812e-5 * metallicity;
    abundance[sH2OII] = 1.3106030161376541e-19 * metallicity;
    abundance[sH3OII] = 2.9572631077328755e-25 * metallicity;
    abundance[sO2II] = 1.6821257889333632e-18 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T34463: {
    /* High-temperature root found for the saved z=20, Z=0.1,
     * T=34463 K non-equilibrium intermediate failure. */
    abundance[sELEC] = 1.00003264289242;
    abundance[sHI] = 1.03013406908217e-5;
    abundance[sHII] = 0.99998204201121632;
    abundance[sH2I] = 1.2172160777026613e-17;
    abundance[sHM] = 1.1701898781268036e-14;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 3.0535703651096367e-10;
    abundance[sDII] = 2.9604957800857676e-5;
    abundance[sHDI] = 5.5231265654500275e-19;
#endif
    abundance[sCI] = 1.2597966882885965e-4 * metallicity;
    abundance[sCII] = 8.3018830260684831e-14 * metallicity;
    abundance[sCO] = 1.5001076262865832e-5 * metallicity;
    abundance[sCO2] = 1.750191666936573e-4 * metallicity;
    abundance[sOI] = 1.3560747480040021e-9 * metallicity;
    abundance[sOH] = 2.1353875395991124e-13 * metallicity;
    abundance[sH2O] = 6.3809201639573046e-18 * metallicity;
    abundance[sO2] = 1.0436809651017052e-17 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 1.6945907606749465e-14 * metallicity;
    abundance[sCH] = 8.0572677676953679e-11 * metallicity;
    abundance[sCH2] = 7.5419748150528964e-12 * metallicity;
    abundance[sOII] = 1.3339284910343401e-4 * metallicity;
    abundance[sOHII] = 7.6566384940889594e-5 * metallicity;
    abundance[sH2OII] = 2.3006907515112232e-16 * metallicity;
    abundance[sH3OII] = 1.405326460540793e-19 * metallicity;
    abundance[sO2II] = 2.9627906658264367e-18 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T30811: {
    /* Root found for the saved z=20, Z=0.1, T=30811 K failure. */
    abundance[sELEC] = 1.0000288624085911;
    abundance[sHI] = 1.3988119506476293e-5;
    abundance[sHII] = 0.99997788016567746;
    abundance[sH2I] = 1.1813936041141837e-14;
    abundance[sHM] = 2.9075900000133324e-14;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 4.1470368754092607e-10;
    abundance[sDII] = 2.9604848454206449e-5;
    abundance[sHDI] = 7.4587283296514796e-19;
#endif
    abundance[sCI] = 1.2979521907697177e-4 * metallicity;
    abundance[sCII] = 7.7535774052324431e-14 * metallicity;
    abundance[sCO] = 1.1185265893047512e-5 * metallicity;
    abundance[sCO2] = 1.7501948788490025e-4 * metallicity;
    abundance[sOI] = 1.8135284293301643e-9 * metallicity;
    abundance[sOH] = 1.3117562321818905e-17 * metallicity;
    abundance[sH2O] = 2.5589089726834677e-18 * metallicity;
    abundance[sO2] = 1.0430683407407197e-17 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 5.0009305674336627e-19 * metallicity;
    abundance[sCH] = 2.7067407023932438e-11 * metallicity;
    abundance[sCH2] = 1.37168559988535e-16 * metallicity;
    abundance[sOII] = 1.3245682424280528e-4 * metallicity;
    abundance[sOHII] = 8.1317120565874268e-5 * metallicity;
    abundance[sH2OII] = 5.2601891224465922e-19 * metallicity;
    abundance[sH3OII] = 4.4300118965131421e-24 * metallicity;
    abundance[sO2II] = 2.7377339204103489e-18 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T1261: {
    /* Root found for the saved z=0, Z=0.1, T=1261 K failure. */
    abundance[sELEC] = 0.46299502150158972;
    abundance[sHI] = 0.23374610242570271;
    abundance[sHII] = 0.46296887826087041;
    abundance[sH2I] = 0.15163314890128671;
    abundance[sHM] = 5.9271967626337258e-7;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 1.683184927734549e-5;
    abundance[sDII] = 1.2773409639027497e-5;
    abundance[sHDI] = 4.2415217524219631e-12;
#endif
    abundance[sCI] = 9.8936286841851588e-5 * metallicity;
    abundance[sCII] = 4.966451336116868e-5 * metallicity;
    abundance[sCO] = 8.674984114183423e-5 * metallicity;
    abundance[sCO2] = 8.064859901241815e-5 * metallicity;
    abundance[sOI] = 3.5921356531547067e-5 * metallicity;
    abundance[sOH] = 4.5881503536558188e-5 * metallicity;
    abundance[sH2O] = 3.1222407874484122e-5 * metallicity;
    abundance[sO2] = 6.1983349342868853e-5 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 0.0;
    abundance[sCH] = 6.6189123140422269e-10 * metallicity;
    abundance[sCH2] = 9.7751495971518061e-11 * metallicity;
    abundance[sOII] = 1.7050635464271139e-5 * metallicity;
    abundance[sOHII] = 7.2860046485103603e-5 * metallicity;
    abundance[sH2OII] = 5.029057419330885e-8 * metallicity;
    abundance[sH3OII] = 2.168143434270531e-11 * metallicity;
    abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T92675: {
    /* Root found for the saved z=20, Z=1, T=92675 K failure. */
    abundance[sELEC] = 1.0000969021535941;
    abundance[sHI] = 1.9564190488324027e-6;
    abundance[sHII] = 0.99991421835454997;
    abundance[sH2I] = 4.4505255218988363e-16;
    abundance[sHM] = 3.4963947708440967e-17;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 5.7952040282222372e-11;
    abundance[sDII] = 2.960520520585443e-5;
    abundance[sHDI] = 2.7857180798267994e-20;
#endif
    abundance[sCI] = 5.5686779410032298e-5 * metallicity;
    abundance[sCII] = 1.3672326857913847e-13 * metallicity;
    abundance[sCO] = 9.870514390705527e-5 * metallicity;
    abundance[sCO2] = 1.6160806521633251e-4 * metallicity;
    abundance[sOI] = 1.3195881492927567e-10 * metallicity;
    abundance[sOH] = 2.4340338772299575e-17 * metallicity;
    abundance[sH2O] = 1.0298898956331353e-18 * metallicity;
    abundance[sO2] = 8.6677703979452596e-18 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 1.7290999112878472e-24 * metallicity;
    abundance[sCH] = 1.1329856501115037e-11 * metallicity;
    abundance[sCH2] = 1.2358743697202434e-19 * metallicity;
    abundance[sOII] = 6.9253378631097865e-5 * metallicity;
    abundance[sOHII] = 8.3825215070314075e-5 * metallicity;
    abundance[sH2OII] = 3.4874029612235333e-19 * metallicity;
    abundance[sH3OII] = 2.8939600484293469e-26 * metallicity;
    abundance[sO2II] = 4.9171363067507986e-18 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T1501: {
    /* Root found for the saved z=0, Z=1, T=1501 K failure. */
    abundance[sELEC] = 9.5276870514691089e-9;
    abundance[sHI] = 0.011505425276457534;
    abundance[sHII] = 6.8250328029227238e-10;
    abundance[sH2I] = 0.4942257502637345;
    abundance[sHM] = 5.2955689998561437e-15;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.8704494388345964e-7;
    abundance[sDII] = 1.4227412396213185e-14;
    abundance[sHDI] = 2.9318218199783867e-5;
#endif
    abundance[sCI] = 3.8783656922327507e-8 * metallicity;
    abundance[sCII] = 1.3763227740127845e-9 * metallicity;
    abundance[sCO] = 6.3451249345320581e-5 * metallicity;
    abundance[sCO2] = 2.5239240133078359e-4 * metallicity;
    abundance[sOI] = 4.2364060284967586e-11 * metallicity;
    abundance[sOH] = 3.2276692703190392e-9 * metallicity;
    abundance[sH2O] = 6.7532087215730819e-6 * metallicity;
    abundance[sO2] = 1.9285299022076825e-13 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 3.964022484995573e-16 * metallicity;
    abundance[sCH] = 9.134910977838556e-9 * metallicity;
    abundance[sCH2] = 1.0705443282522379e-7 * metallicity;
    abundance[sOII] = 3.581165004835212e-18 * metallicity;
    abundance[sOHII] = 2.7651575745559313e-17 * metallicity;
    abundance[sH2OII] = 7.7561751953795306e-14 * metallicity;
    abundance[sH3OII] = 7.4687736386544986e-9 * metallicity;
    abundance[sO2II] = 4.3727921331630416e-16 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T1515: {
    /* Root found for the saved z=20, Z=1, T=1515 K failure. */
    abundance[sELEC] = 0.30236727032410249;
    abundance[sHI] = 0.5795951115982052;
    abundance[sHII] = 0.30233947754074625;
    abundance[sH2I] = 0.059027888873290288;
    abundance[sHM] = 9.7375063629401053e-8;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 1.9322771279217555e-5;
    abundance[sDII] = 1.0282489050219449e-5;
    abundance[sHDI] = 2.8284577354283037e-12;
#endif
    abundance[sCI] = 1.7626131693844425e-5 * metallicity;
    abundance[sCII] = 1.1724828034265817e-5 * metallicity;
    abundance[sCO] = 1.583620530423108e-5 * metallicity;
    abundance[sCO2] = 2.7081155884091612e-4 * metallicity;
    abundance[sOI] = 4.751595719675113e-6 * metallicity;
    abundance[sOH] = 2.2991851135850131e-6 * metallicity;
    abundance[sH2O] = 1.3739858274135889e-6 * metallicity;
    abundance[sO2] = 1.6165345089195457e-6 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 0.0;
    abundance[sCH] = 1.2334072756372903e-9 * metallicity;
    abundance[sCH2] = 4.2719466865926816e-11 * metallicity;
    abundance[sOII] = 1.3990193853007009e-6 * metallicity;
    abundance[sOHII] = 4.480383830739206e-6 * metallicity;
    abundance[sH2OII] = 3.4372276239549154e-9 * metallicity;
    abundance[sH3OII] = 8.9175999873494229e-13 * metallicity;
    abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T1308: {
    /* Molecular root found for the saved z=20, Z=1, T=1308 K failure. */
    abundance[sELEC] = 3.3816627299155693e-8;
    abundance[sHI] = 0.090029236423097111;
    abundance[sHII] = 5.0973636985557156e-9;
    abundance[sH2I] = 0.45495950168833826;
    abundance[sHM] = 1.5124297197929668e-14;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.3271032762120112e-6;
    abundance[sDII] = 1.0732861297366062e-13;
    abundance[sHDI] = 2.7278159774354116e-5;
#endif
    abundance[sCI] = 3.3904397142241463e-8 * metallicity;
    abundance[sCII] = 2.2304322241103668e-9 * metallicity;
    abundance[sCO] = 6.9165214121594593e-5 * metallicity;
    abundance[sCO2] = 2.4679800839385235e-4 * metallicity;
    abundance[sOI] = 2.7774254938040671e-9 * metallicity;
    abundance[sOH] = 2.251429314960598e-8 * metallicity;
    abundance[sH2O] = 1.2186987979377911e-5 * metallicity;
    abundance[sO2] = 3.2604709313020105e-13 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 5.7265524463316557e-16 * metallicity;
    abundance[sCH] = 2.9756172098530313e-10 * metallicity;
    abundance[sCH2] = 3.4509289303746159e-10 * metallicity;
    abundance[sOII] = 5.2325901030760667e-17 * metallicity;
    abundance[sOHII] = 8.0564714123643722e-16 * metallicity;
    abundance[sH2OII] = 1.1360788922976149e-12 * metallicity;
    abundance[sH3OII] = 2.6487600249823491e-8 * metallicity;
    abundance[sO2II] = 1.412829757993675e-15 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T1625: {
    /* Molecular root found for the saved z=15, Z=0.1, T=1625 K
     * energy-trial failure; C/O abundances scale with metallicity. */
    abundance[sELEC] = 4.5382337802913465e-6;
    abundance[sHI] = 0.13212707000017435;
    abundance[sHII] = 4.5367772036671654e-6;
    abundance[sH2I] = 0.43392009380269708;
    abundance[sHM] = 6.6903775210040103e-11;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 3.2407443278998604e-6;
    abundance[sDII] = 9.6176186589400081e-16;
    abundance[sHDI] = 2.6364518829033115e-5;
#endif
    abundance[sCI] = 2.9418491714483299e-8 * metallicity;
    abundance[sCII] = 1.4083976547487872e-9 * metallicity;
    abundance[sCO] = 7.269600596065384e-5 * metallicity;
    abundance[sCO2] = 2.4318552331419869e-4 * metallicity;
    abundance[sOI] = 1.1383006223979571e-9 * metallicity;
    abundance[sOH] = 6.3461958554081927e-8 * metallicity;
    abundance[sH2O] = 9.0685148668328659e-6 * metallicity;
    abundance[sO2] = 3.393002943944483e-6 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 3.0580833826225916e-16 * metallicity;
    abundance[sCH] = 6.8974084235700438e-9 * metallicity;
    abundance[sCH2] = 8.0746427048903346e-8 * metallicity;
    abundance[sOII] = 9.6314986665154312e-15 * metallicity;
    abundance[sOHII] = 2.8133263461932552e-12 * metallicity;
    abundance[sH2OII] = 4.0537034318655057e-11 * metallicity;
    abundance[sH3OII] = 1.3783036094629452e-8 * metallicity;
    abundance[sO2II] = 3.2894633155236885e-16 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T13492: {
    /* High-ionization grid root found for the saved z=0, Z=0.1,
     * T=13492 K energy trial; C/O abundances scale with metallicity. */
    abundance[sELEC] = 0.94198017890404606;
    abundance[sHI] = 0.05805441285023856;
    abundance[sHII] = 0.94192033322456736;
    abundance[sH2I] = 4.3520211046943833e-10;
    abundance[sHM] = 1.3618147972046997e-8;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 1.7239264614029546e-6;
    abundance[sDII] = 2.7881336626535426e-5;
    abundance[sHDI] = 6.9956361354992132e-14;
#endif
    abundance[sCI] = 1.4553278290871734e-4 * metallicity;
    abundance[sCII] = 5.1832200636743317e-14 * metallicity;
    abundance[sCO] = 8.9818522821354362e-5 * metallicity;
    abundance[sCO2] = 8.0648662842740691e-5 * metallicity;
    abundance[sOI] = 4.1045415185868388e-6 * metallicity;
    abundance[sOH] = 2.2327059181335104e-14 * metallicity;
    abundance[sH2O] = 2.5110321146178189e-16 * metallicity;
    abundance[sO2] = 4.5117405652592979e-21 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 3.5097788223583212e-19 * metallicity;
    abundance[sCH] = 3.1373455903840561e-11 * metallicity;
    abundance[sCH2] = 1.8991432439961332e-15 * metallicity;
    abundance[sOII] = 6.7385275627295196e-5 * metallicity;
    abundance[sOHII] = 2.5239433431743299e-4 * metallicity;
    abundance[sH2OII] = 1.272010818433527e-15 * metallicity;
    abundance[sH3OII] = 5.5310536451555839e-20 * metallicity;
    abundance[sO2II] = 2.999343130095177e-15 * metallicity;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T2817: {
    /* Partly ionized root found for the saved z=0, Z=0.1,
     * T=2817 K energy trial; C/O abundances scale with metallicity. */
    abundance[sELEC] = 0.12300980603699925;
    abundance[sHI] = 0.85730009464603152;
    abundance[sHII] = 0.12300585195318735;
    abundance[sH2I] = 0.0098465571490409364;
    abundance[sHM] = 3.0749809613547881e-8;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.5412102424257962e-5;
    abundance[sDII] = 3.5886705848171222e-6;
    abundance[sHDI] = 6.0449014881965254e-7;
#endif
    abundance[sCI] = 1.36278132667118e-11 * metallicity;
    abundance[sCII] = 4.1463444817672722e-7 * metallicity;
    abundance[sCO] = 6.5280680148379464e-5 * metallicity;
    abundance[sCO2] = 2.5030467177313828e-4 * metallicity;
    abundance[sOI] = 5.3695176281078972e-6 * metallicity;
    abundance[sOH] = 7.7137316555990897e-8 * metallicity;
    abundance[sH2O] = 1.0546855146720842e-7 * metallicity;
    abundance[sO2] = 5.4284453534823787e-9 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 2.4922831628171824e-15 * metallicity;
    abundance[sCH] = 0.0;
    abundance[sCH2] = 0.0;
    abundance[sOII] = 7.9884557638556961e-7 * metallicity;
    abundance[sOHII] = 2.7457484316224926e-6 * metallicity;
    abundance[sH2OII] = 2.4011672237160019e-9 * metallicity;
    abundance[sH3OII] = 7.4078187054886648e-13 * metallicity;
    abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z1_T3137: {
    /* Partly ionized root found for the saved z=0, Z=1,
     * T=3137 K energy trial. */
    abundance[sELEC] = 0.24368227658044167;
    abundance[sHI] = 0.74734107188133525;
    abundance[sHII] = 0.24364822196622987;
    abundance[sH2I] = 0.0044765617199466444;
    abundance[sHM] = 2.3168149884005082e-8;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 1.4037980521944389e-5;
    abundance[sDII] = 5.8307747828085568e-6;
    abundance[sHDI] = 9.7365078531417911e-6;
#endif
    abundance[sCI] = 5.5117392151274453e-5 * metallicity;
    abundance[sCII] = 1.9187888638307844e-5 * metallicity;
    abundance[sCO] = 7.9794004867361318e-5 * metallicity;
    abundance[sCO2] = 1.6187818205856807e-4 * metallicity;
    abundance[sOI] = 8.0972686877595746e-5 * metallicity;
    abundance[sOH] = 4.7800103509553789e-5 * metallicity;
    abundance[sH2O] = 1.6282207832983737e-10 * metallicity;
    abundance[sO2] = 1.680877943286188e-5 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 0.0;
    abundance[sCH] = 2.2457184144262032e-8 * metallicity;
    abundance[sCH2] = 7.5100344065396327e-11 * metallicity;
    abundance[sOII] = 9.0591189405508754e-6 * metallicity;
    abundance[sOHII] = 0.0;
    abundance[sH2OII] = 0.0;
    abundance[sH3OII] = 0.0;
    abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_Z0P1_T2055: {
    /* Molecular root found for the saved z=0, Z=0.1,
     * T=2055 K energy trial; C/O abundances scale with metallicity. */
    abundance[sELEC] = 1.7292161219077955e-4;
    abundance[sHI] = 0.15260158073989594;
    abundance[sHII] = 1.7291719422397354e-4;
    abundance[sH2I] = 0.42359976420960388;
    abundance[sHM] = 7.1611168879531515e-10;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 3.7373343442930741e-6;
    abundance[sDII] = 4.0629493336566334e-9;
    abundance[sHDI] = 2.5863865864268007e-5;
#endif
    abundance[sCI] = 1.6482734171834373e-8 * metallicity;
    abundance[sCII] = 4.7025756586749794e-11 * metallicity;
    abundance[sCO] = 5.7563937441663791e-5 * metallicity;
    abundance[sCO2] = 2.5841534867819711e-4 * metallicity;
    abundance[sOI] = 2.0715754337117436e-8 * metallicity;
    abundance[sOH] = 9.3338730701021626e-8 * metallicity;
    abundance[sH2O] = 4.806464506895454e-7 * metallicity;
    abundance[sO2] = 4.9502143698811777e-19 * metallicity;
#ifndef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = 2.9752211192564415e-18 * metallicity;
    abundance[sCH] = 2.0528614763825189e-9 * metallicity;
    abundance[sCH2] = 2.1312587313702209e-9 * metallicity;
    abundance[sOII] = 2.159307716875683e-12 * metallicity;
    abundance[sOHII] = 8.2736084026850442e-11 * metallicity;
    abundance[sH2OII] = 2.1208936678149553e-9 * metallicity;
    abundance[sH3OII] = 8.4584764328809408e-9 * metallicity;
    abundance[sO2II] = 3.5896005208731092e-16 * metallicity;
#endif
    break;
  }
#else
  case CHEM_NEWTON_RESET_INITIAL_H2: {
    abundance[sHI] = 0.99999;
    abundance[sELEC] = 1.e-5;
    abundance[sHII] = 1.e-5;
    abundance[sH2I] = 1.e-10;
    abundance[sHM] = 1.e-15;
    break;
  }

  case CHEM_NEWTON_RESET_FULLY_IONIZED_H: {
    abundance[sHII] = 0.99999;
    abundance[sELEC] = 0.99999;
    abundance[sHI] = 1.e-5;
    abundance[sH2I] = 1.e-10;
    abundance[sHM] = 1.e-15;
    break;
  }

  case CHEM_NEWTON_RESET_MOLECULAR_H: {
    abundance[sHII] = 1.e-10;
    abundance[sELEC] = 1.e-10;
    abundance[sHI] = 1.e-5;
    abundance[sH2I] = 0.499995;
    abundance[sHM] = 1.e-15;
    break;
  }
#endif

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  case CHEM_NEWTON_RESET_CO_RICH:
  case CHEM_NEWTON_RESET_CO2_RICH:
  {
    double carbon_total = A_C_sun * metallicity;
    double oxygen_total = A_O_sun * metallicity;
    abundance[sCI] = abundance[sCII] = abundance[sCO] = abundance[sCO2] = 0.0;
    abundance[sOI] = abundance[sOH] = abundance[sH2O] = abundance[sO2] = 0.0;
    if(reset == CHEM_NEWTON_RESET_CO_RICH) {
      abundance[sCO] = dmin(carbon_total, oxygen_total);
      abundance[sCII] = carbon_total - abundance[sCO];
      abundance[sO2] = 0.5 * (oxygen_total - abundance[sCO]);
    } else if(reset == CHEM_NEWTON_RESET_CO2_RICH) {
      abundance[sCO2] = dmin(carbon_total, 0.5 * oxygen_total);
      abundance[sCII] = carbon_total - abundance[sCO2];
      abundance[sO2] = 0.5 * (oxygen_total - 2.0 * abundance[sCO2]);
    }
    break;
  }
#endif

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  default:
    chemistry_newton_apply_carbon_oxygen_corner_reset(reset, metallicity,
                                                      abundance);
    break;
#else
  default:
    break;
#endif
  }
}

#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
static void chemistry_newton_apply_carbon_oxygen_corner_reset(
    ChemNewtonResetIndex reset, double metallicity, double *abundance) {
  if(reset < CHEM_NEWTON_RESET_MOLECULAR_RESERVOIR_CORNER ||
      reset > CHEM_NEWTON_RESET_LAST)
    return;

  switch(reset){
  case CHEM_NEWTON_RESET_MOLECULAR_RESERVOIR_CORNER: {
    /* Network-independent molecular reservoir corner.  This is generated
     * from the current element totals, not copied from a saved sample. */
    const double reset_trace_fraction = 1.e-10;
    abundance[sELEC] = 1.e-12;
    abundance[sHI] = 1.25e-7;
    abundance[sHII] = 1.e-13;
    abundance[sH2I] = 0.0;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDI] = 7.2e-5 * reset_deuterium_total;
    abundance[sDII] = 1.e-12 * reset_deuterium_total;
    abundance[sHDI] = reset_deuterium_total - abundance[sDI] - abundance[sDII];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    abundance[sCI] = abundance[sCII] =
        reset_trace_fraction * reset_carbon_total;
    abundance[sCO2] = 0.01582 * reset_carbon_total;
    abundance[sCO] =
        reset_carbon_total - abundance[sCI] - abundance[sCII] - abundance[sCO2];
    abundance[sOH] = abundance[sH2O] =
        reset_trace_fraction * reset_oxygen_total;
    abundance[sO2] = 5.e-4 * reset_oxygen_total;
    abundance[sOI] = reset_oxygen_total - abundance[sCO] -
                     2.0 * abundance[sCO2] - abundance[sOH] - abundance[sH2O] -
                     2.0 * abundance[sO2];
    /* Complete the H budget only after every H-bearing D/O species has
     * been assigned.  The previous formula omitted HD/OH/H2O and made HI
     * negative when an inherited HI chart reconstructed the reset. */
    double reset_fixed_hydrogen = abundance[sHI] + abundance[sHII] +
                                  abundance[sHM] + abundance[sOH] +
                                  2.0 * abundance[sH2O];
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    reset_fixed_hydrogen += abundance[sHDI];
#endif
    abundance[sH2I] = 0.5 * (1.0 - reset_fixed_hydrogen);
#ifdef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_CI_H2O_CORNER:
  case CHEM_NEWTON_RESET_CO_OI_CORNER:
  {
    /* Two conservation-valid boundary charts selected by direct Newton
     * replay.  Both are element-total normalized and therefore work at
     * arbitrary metallicity: CI/H2O and CO/OI reservoirs. */
    const double reset_epsilon_co = 1.e-4;
    const double reset_epsilon_h = 1.e-6;
    const double reset_epsilon_d = 1.e-8;
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    double reset_carbon_scale =
        reset_carbon_total / (1.0 + 3.0 * reset_epsilon_co);
    abundance[sCI] = abundance[sCII] = abundance[sCO] = abundance[sCO2] =
        reset_epsilon_co * reset_carbon_scale;
    if(reset == CHEM_NEWTON_RESET_CI_H2O_CORNER)
      abundance[sCI] = reset_carbon_scale;
    else
      abundance[sCO] = reset_carbon_scale;
    double reset_oxygen_free =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    double reset_oxygen_scale =
        reset_oxygen_free / (1.0 + 4.0 * reset_epsilon_co);
    abundance[sOI] = abundance[sOH] = abundance[sH2O] = abundance[sO2] =
        reset_epsilon_co * reset_oxygen_scale;
    if(reset == CHEM_NEWTON_RESET_CI_H2O_CORNER)
      abundance[sH2O] = reset_oxygen_scale;
    else
      abundance[sOI] = reset_oxygen_scale;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    double reset_deuterium_scale =
        reset_deuterium_total / (1.0 + 2.0 * reset_epsilon_d);
    abundance[sDI] = reset_deuterium_scale;
    abundance[sDII] = abundance[sHDI] = reset_epsilon_d * reset_deuterium_scale;
#endif
    double reset_fixed_hydrogen =
        abundance[sHM] + abundance[sOH] + 2.0 * abundance[sH2O];
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    reset_fixed_hydrogen += abundance[sHDI];
#endif
    double reset_hydrogen_scale =
        (1.0 - reset_fixed_hydrogen) / (2.0 + 2.0 * reset_epsilon_h);
    abundance[sHI] = abundance[sHII] = reset_epsilon_h * reset_hydrogen_scale;
    abundance[sH2I] = reset_hydrogen_scale;
    abundance[sELEC] = dmax(abundance[sHII] + abundance[sCII]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                                + abundance[sDII]
#endif
                                - abundance[sHM],
                            1.e-80);
#ifdef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
#endif
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_SAMPLE_88: {
    /* Near-root state found by direct C/O-chart replay for off-trajectory
     * sample 88 at Z=0.01, nH=7.872536470626472e15 cm^-3,
     * T=1361.4753088816972 K and dt=2.1183195127816079e9 s.  Store the
     * C and free-O distributions as element-total fractions so this is a
     * physical reset state rather than a Z-specific abundance vector. */
    abundance[sHI] = 9.740873540219346e-5;
    abundance[sHII] = 1.097372861653367e-16;
    abundance[sH2I] = 0.49993281348707092;
    abundance[sHM] = 1.3578795895763941e-18;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.4678081475222134e-9;
    abundance[sDII] = 3.8769071436259128e-21;
    abundance[sHDI] = 2.9602795349747216e-5;
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    abundance[sCI] = 2.7236981931876617e-5 * reset_carbon_total;
    abundance[sCII] = 9.2726272782868038e-15 * reset_carbon_total;
    abundance[sCO] = 0.057427191664080272 * reset_carbon_total;
    abundance[sCO2] = 0.51325729125156616 * reset_carbon_total;
    abundance[sCOII] = 5.2377565302451654e-22 * reset_carbon_total;
    abundance[sCH] = 0.00034537570760752605 * reset_carbon_total;
    abundance[sCH2] = 0.42894290439480498 * reset_carbon_total;

    double reset_oxygen_total = A_O_sun * metallicity;
    double reset_free_oxygen = reset_oxygen_total - abundance[sCO] -
                               2.0 * abundance[sCO2] - abundance[sCOII];
    abundance[sOI] = 2.7154698158130587e-10 * reset_free_oxygen;
    abundance[sOH] = 2.3053443829496498e-6 * reset_free_oxygen;
    abundance[sH2O] = 0.99999686142535715 * reset_free_oxygen;
    abundance[sO2] = 4.510173344787477e-16 * reset_free_oxygen;
    abundance[sOII] = 2.168339254323813e-25 * reset_free_oxygen;
    abundance[sOHII] = 1.0656894832838836e-21 * reset_free_oxygen;
    abundance[sH2OII] = 1.8247884206553372e-15 * reset_free_oxygen;
    abundance[sH3OII] = 8.3295871001195412e-7 * reset_free_oxygen;
    abundance[sO2II] = 9.2579730581472613e-22 * reset_free_oxygen;
    abundance[sELEC] =
        dmax(abundance[sHII] + abundance[sCII] + abundance[sCOII] +
                 abundance[sOII] + abundance[sOHII] + abundance[sH2OII] +
                 abundance[sH3OII] + abundance[sO2II]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                 + abundance[sDII]
#endif
                 - abundance[sHM],
             1.e-80);
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_SAMPLE_3624: {
    /* Root reached from off-trajectory sample 3624 at Z=0.01,
     * nH=5.1162223511819658e12 cm^-3, T=561.97213743466682 K and
     * dt=2.6882432374101944e10 s.  C/O abundances are stored per unit
     * metallicity; the conservation reconstruction below absorbs their
     * metallicity-dependent H and charge budgets. */
    abundance[sELEC] = 1.7377056872432293e-11;
    abundance[sHI] = 2.4943441957980505e-5;
    abundance[sHII] = 5.5963065945879032e-15;
    abundance[sH2I] = 0.49996462101189643;
    abundance[sHM] = 2.9032775962086918e-18;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 8.2229956853841889e-10;
    abundance[sDII] = 6.9111703489780351e-20;
    abundance[sHDI] = 2.9604440858326133e-5;
#endif
    abundance[sCI] = 7.2671917575693442e-5 * metallicity;
    abundance[sCII] = 2.030799749037776e-14 * metallicity;
    abundance[sCO] = 3.8344117052100503e-6 * metallicity;
    abundance[sCO2] = 1.276808899894878e-10 * metallicity;
    abundance[sOI] = 7.6342134532373436e-16 * metallicity;
    abundance[sOH] = 2.4176987342949408e-13 * metallicity;
    abundance[sH2O] = 5.7116359556324051e-4 * metallicity;
    abundance[sO2] = 6.0208367797959158e-16 * metallicity;
    abundance[sCOII] = 2.1437035682246683e-24 * metallicity;
    abundance[sCH] = 3.1015025340336876e-7 * metallicity;
    abundance[sCH2] = 2.3918339276449513e-4 * metallicity;
    abundance[sOII] = 7.3364211916183326e-30 * metallicity;
    abundance[sOHII] = 5.6474386764171152e-27 * metallicity;
    abundance[sH2OII] = 5.3148833656578055e-17 * metallicity;
    abundance[sH3OII] = 1.7371259788499905e-9 * metallicity;
    abundance[sO2II] = 4.0432125570381598e-21 * metallicity;
    break;
  }
  case CHEM_NEWTON_RESET_SAVED_SAMPLE_16123: {
    /* Root reached from off-trajectory sample 16123 at Z=0.01,
     * nH=1.7367829880121902e6 cm^-3, T=7.3848425884177065 K and
     * dt=4.3867186761266769e14 s. */
    abundance[sELEC] = 6.7702006277014125e-9;
    abundance[sHI] = 0.0060127875042547572;
    abundance[sHII] = 2.8307320932470535e-14;
    abundance[sH2I] = 0.49699165106017362;
    abundance[sHM] = 1.9316995099120814e-17;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.9602232224597413e-5;
    abundance[sDII] = 2.0920073158092711e-21;
    abundance[sHDI] = 3.0309332973224368e-9;
#endif
    abundance[sCI] = 7.6341982761826718e-5 * metallicity;
    abundance[sCII] = 6.7701634129701266e-7 * metallicity;
    abundance[sCO] = 9.838709462709993e-5 * metallicity;
    abundance[sCO2] = 1.4059239013482965e-4 * metallicity;
    abundance[sOI] = 5.8250646224045399e-8 * metallicity;
    abundance[sOH] = 7.0532967886543955e-9 * metallicity;
    abundance[sH2O] = 1.9536281893676159e-4 * metallicity;
    abundance[sO2] = 6.6539702017001291e-13 * metallicity;
    abundance[sCOII] = 8.1535580224251875e-15 * metallicity;
    abundance[sCH] = 1.2824191541349561e-9 * metallicity;
    abundance[sCH2] = 2.3370763899244913e-10 * metallicity;
    abundance[sOII] = 3.2898571868510469e-21 * metallicity;
    abundance[sOHII] = 7.0790840851173877e-21 * metallicity;
    abundance[sH2OII] = 9.2515341675430762e-17 * metallicity;
    abundance[sH3OII] = 8.8442643975149613e-13 * metallicity;
    abundance[sO2II] = 2.134163989138783e-21 * metallicity;
    break;
  }
  case CHEM_NEWTON_RESET_GRID_SAMPLE_24548: {
    /* Direct-Newton grid candidate 81312 for off-trajectory sample 24548
     * at Z=0.01.  It converges under the ordinary maximum-reservoir chart
     * in ten iterations; C/O values are normalized per unit metallicity. */
    abundance[sELEC] = 2.7906630211449072e-6;
    abundance[sHI] = 0.018204671462273718;
    abundance[sHII] = 2.1112122281401323e-10;
    abundance[sH2I] = 0.49089189183000276;
    abundance[sHM] = 4.7800845857078197e-19;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 2.6789152029287045e-5;
    abundance[sDII] = 2.76662035790657e-6;
    abundance[sHDI] = 4.9490770701122862e-8;
#endif
    abundance[sCI] = 3.1344771474477368e-4 * metallicity;
    abundance[sCII] = 2.3831542016001145e-6 * metallicity;
    abundance[sCO] = 9.7140680819436931e-8 * metallicity;
    abundance[sCO2] = 7.1990372806713356e-8 * metallicity;
    abundance[sOI] = 1.3531878007790244e-11 * metallicity;
    abundance[sOH] = 3.925647551083025e-11 * metallicity;
    abundance[sH2O] = 5.7475877181767181e-4 * metallicity;
    abundance[sO2] = 2.6983770881600788e-11 * metallicity;
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    break;
  }
  case CHEM_NEWTON_RESET_H2_DI_CO_OI_CORNER: {
    /* Conservation-valid H2 + DI + CO + OI corner found by direct Newton
     * replay for off-trajectory sample 16475. */
    const double reset_trace_fraction = 1.e-16;
    abundance[sELEC] = 5.0000000327602629e-13;
    abundance[sHI] = 4.9999999999949997e-13;
    abundance[sHII] = 4.9999999999949997e-13;
    abundance[sH2I] = 0.49999999999949996;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDI] = reset_deuterium_total / (1.0 + 2.0 * reset_trace_fraction);
    abundance[sDII] = abundance[sHDI] = reset_trace_fraction * abundance[sDI];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    abundance[sCO] = reset_carbon_total / (1.0 + 3.0 * reset_trace_fraction);
    abundance[sCI] = abundance[sCII] = abundance[sCO2] =
        reset_trace_fraction * abundance[sCO];
    double reset_free_oxygen =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    abundance[sOI] = reset_free_oxygen / (1.0 + 4.0 * reset_trace_fraction);
    abundance[sOH] = abundance[sH2O] = abundance[sO2] =
        reset_trace_fraction * abundance[sOI];
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    break;
  }
  case CHEM_NEWTON_RESET_H2_DI_CI_OI_CORNER: {
    /* Conservation-valid H2 + DI + CI + OI corner found by direct Newton
     * replay for off-trajectory sample 27727. */
    const double reset_trace_fraction = 1.e-16;
    abundance[sELEC] = 4.9999999500032774e-9;
    abundance[sHI] = 4.9999999500000009e-9;
    abundance[sHII] = 4.9999999500000009e-9;
    abundance[sH2I] = 0.49999999500000009;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDI] = reset_deuterium_total / (1.0 + 2.0 * reset_trace_fraction);
    abundance[sDII] = abundance[sHDI] = reset_trace_fraction * abundance[sDI];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    abundance[sCI] = reset_carbon_total / (1.0 + 3.0 * reset_trace_fraction);
    abundance[sCII] = abundance[sCO] = abundance[sCO2] =
        reset_trace_fraction * abundance[sCI];
    double reset_free_oxygen =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    abundance[sOI] = reset_free_oxygen / (1.0 + 4.0 * reset_trace_fraction);
    abundance[sOH] = abundance[sH2O] = abundance[sO2] =
        reset_trace_fraction * abundance[sOI];
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    break;
  }
  case CHEM_NEWTON_RESET_H2_CO_OI_BROAD_D_CORNER: {
    /* H2 + CO + OI corner with a broader D trace.  Direct replay reaches
     * off-trajectory sample 95470 in seven Newton iterations. */
    const double reset_co_trace_fraction = 1.e-16;
    const double reset_d_trace_fraction = 1.e-4;
    abundance[sELEC] = 5.0000000327602629e-13;
    abundance[sHI] = 4.9999999999949997e-13;
    abundance[sHII] = 4.9999999999949997e-13;
    abundance[sH2I] = 0.49999999999949996;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDI] =
        reset_deuterium_total / (1.0 + 2.0 * reset_d_trace_fraction);
    abundance[sDII] = abundance[sHDI] = reset_d_trace_fraction * abundance[sDI];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    abundance[sCO] = reset_carbon_total / (1.0 + 3.0 * reset_co_trace_fraction);
    abundance[sCI] = abundance[sCII] = abundance[sCO2] =
        reset_co_trace_fraction * abundance[sCO];
    double reset_free_oxygen =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    abundance[sOI] = reset_free_oxygen / (1.0 + 4.0 * reset_co_trace_fraction);
    abundance[sOH] = abundance[sH2O] = abundance[sO2] =
        reset_co_trace_fraction * abundance[sOI];
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    break;
  }
  case CHEM_NEWTON_RESET_COLD_CH2_H2O_CORNER: {
    /* Cold, dense CH2/H2O reservoir branch reached from off-trajectory
     * sample 1149 at Z=0.01, nH=9.956139068448704e15 cm^-3,
     * T=328.0969143137591 K and dt=1.693861624390862e8 s.  Store C and
     * free-O distributions as element-total fractions so the candidate is
     * conservation-valid at arbitrary metallicity. */
    abundance[sHI] = 2.0146333439264597e-6;
    abundance[sHII] = 4.0474475656545944e-16;
    abundance[sH2I] = 0.49997548763421412;
    abundance[sHM] = 5.8817281868462176e-17;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 9.1429852689894538e-11;
    abundance[sDII] = 1.7837452166529701e-20;
    abundance[sHDI] = 2.9605171728042033e-5;
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    abundance[sCI] = 0.018750348837055668 * reset_carbon_total;
    abundance[sCII] = 0.00021187284415783157 * reset_carbon_total;
    abundance[sCO] = 0.023110055861679645 * reset_carbon_total;
    abundance[sCO2] = 8.8836870627141572e-6 * reset_carbon_total;
    abundance[sCOII] = 2.4837219534712631e-13 * reset_carbon_total;
    abundance[sCH] = 0.00093439785259926919 * reset_carbon_total;
    abundance[sCH2] = 0.95698444091719648 * reset_carbon_total;

    double reset_oxygen_total = A_O_sun * metallicity;
    double reset_free_oxygen = reset_oxygen_total - abundance[sCO] -
                               2.0 * abundance[sCO2] - abundance[sCOII];
    abundance[sOI] = 2.6424397940109715e-10 * reset_free_oxygen;
    abundance[sOH] = 4.0151998418255464e-10 * reset_free_oxygen;
    abundance[sH2O] = 0.99999999051252841 * reset_free_oxygen;
    abundance[sO2] = 1.9537449426693896e-9 * reset_free_oxygen;
    abundance[sOII] = 9.5418100473327514e-19 * reset_free_oxygen;
    abundance[sOHII] = 7.0003522283248861e-16 * reset_free_oxygen;
    abundance[sH2OII] = 7.5817045187696974e-15 * reset_free_oxygen;
    abundance[sH3OII] = 4.9142059577934431e-9 * reset_free_oxygen;
    abundance[sO2II] = 1.7666211423415686e-15 * reset_free_oxygen;
    abundance[sELEC] =
        dmax(abundance[sHII] + abundance[sCII] + abundance[sCOII] +
                 abundance[sOII] + abundance[sOHII] + abundance[sH2OII] +
                 abundance[sH3OII] + abundance[sO2II]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                 + abundance[sDII]
#endif
                 - abundance[sHM],
             1.e-80);
    break;
  }
  case CHEM_NEWTON_RESET_COLD_CO_O_O2_CORNER: {
    /* Cold molecular CO/O/O2 basin reached by off-trajectory sample 20 at
     * Z=0.01, nH=1.1937948475726013e8 cm^-3, T=46.052149811917154 K
     * and dt=1.8470203016110297e13 s.  Keep this measured root at the end
     * of the deterministic dictionary; C and free O are normalized to
     * their element totals. */
    abundance[sHI] = 0.0010517724897236018;
    abundance[sHII] = 4.0342846736614977e-12;
    abundance[sH2I] = 0.49946226387603804;
    abundance[sHM] = 9.6731775568793928e-19;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    abundance[sDI] = 6.173468146750116e-6;
    abundance[sDII] = 1.0135491570650765e-17;
    abundance[sHDI] = 2.343179501113449e-5;
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    abundance[sCI] = 2.849186761876138e-5 * reset_carbon_total;
    abundance[sCII] = 2.1800652982451151e-11 * reset_carbon_total;
    abundance[sCO] = 0.94388873302540843 * reset_carbon_total;
    abundance[sCO2] = 0.056082771866082584 * reset_carbon_total;
    abundance[sCOII] = 1.4429740979493491e-18 * reset_carbon_total;
    abundance[sCH] = 2.7179738363885452e-9 * reset_carbon_total;
    abundance[sCH2] = 5.0111572053791895e-10 * reset_carbon_total;

    double reset_oxygen_total = A_O_sun * metallicity;
    double reset_free_oxygen = reset_oxygen_total - abundance[sCO] -
                               2.0 * abundance[sCO2] - abundance[sCOII];
    abundance[sOI] = 0.561680492026429 * reset_free_oxygen;
    abundance[sOH] = 2.8102534546831752e-5 * reset_free_oxygen;
    abundance[sH2O] = 0.055502032664592256 * reset_free_oxygen;
    abundance[sO2] = 0.19136545331429355 * reset_free_oxygen;
    abundance[sOII] = 1.2120178454368017e-14 * reset_free_oxygen;
    abundance[sOHII] = 2.0877440075656534e-14 * reset_free_oxygen;
    abundance[sH2OII] = 3.7687841210819245e-12 * reset_free_oxygen;
    abundance[sH3OII] = 7.3390222443628987e-6 * reset_free_oxygen;
    abundance[sO2II] = 2.556355989923708e-5 * reset_free_oxygen;
    abundance[sELEC] =
        dmax(abundance[sHII] + abundance[sCII] + abundance[sCOII] +
                 abundance[sOII] + abundance[sOHII] + abundance[sH2OII] +
                 abundance[sH3OII] + abundance[sO2II]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                 + abundance[sDII]
#endif
                 - abundance[sHM],
             1.e-80);
    break;
  }
  case CHEM_NEWTON_RESET_MINIMUM_COVER_NEUTRAL_H_DII: {
    /* Minimum-set-cover CI/H2O + neutral-H/DII corner.  With the same
     * within-basin recovery policy in teacher and direct replay, it covers
     * samples 5168, 6447, 38946 and 95809 in 14--20 iterations. */
    const double reset_d_trace = 1.e-4;
    const double reset_h_trace = 1.e-16;
    const double reset_co_trace = 1.e-16;
    abundance[sHM] = 1.e-8;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDII] = reset_deuterium_total / (1.0 + 2.0 * reset_d_trace);
    abundance[sDI] = abundance[sHDI] = reset_d_trace * abundance[sDII];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    double reset_oxygen_total = A_O_sun * metallicity;
    abundance[sCI] = reset_carbon_total / (1.0 + 3.0 * reset_co_trace);
    abundance[sCII] = abundance[sCO] = abundance[sCO2] =
        reset_co_trace * abundance[sCI];
    double reset_free_oxygen =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    abundance[sH2O] = reset_free_oxygen / (1.0 + 4.0 * reset_co_trace);
    abundance[sOI] = abundance[sOH] = abundance[sO2] =
        reset_co_trace * abundance[sH2O];
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    double reset_available_hydrogen = 1.0 - abundance[sHM]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                                      - abundance[sHDI]
#endif
                                      - abundance[sOH] - 2.0 * abundance[sH2O];
    abundance[sHI] = reset_available_hydrogen / (1.0 + 3.0 * reset_h_trace);
    abundance[sHII] = abundance[sH2I] = reset_h_trace * abundance[sHI];
    abundance[sELEC] = dmax(abundance[sHII] + abundance[sCII]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                                + abundance[sDII]
#endif
                                - abundance[sHM],
                            1.e-80);
    break;
  }
  case CHEM_NEWTON_RESET_MINIMUM_COVER_SECOND: {
    /* Second minimum-set-cover corner; covers sample 66332 in 12 direct
     * Newton iterations. */
    const double reset_trace_fraction = 1.e-8;
    const double reset_h_trace = 1.e-8;
    abundance[sHM] = 1.e-30;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
    double reset_deuterium_total = D_number_fraction / XH;
    abundance[sDI] = reset_deuterium_total / (1.0 + 2.0e-16);
    abundance[sDII] = abundance[sHDI] = 1.e-16 * abundance[sDI];
#endif
    double reset_carbon_total = A_C_sun * metallicity;
    abundance[sCI] = reset_carbon_total / (1.0 + 3.0 * reset_trace_fraction);
    abundance[sCII] = abundance[sCO] = abundance[sCO2] =
        reset_trace_fraction * abundance[sCI];
    double reset_oxygen_total = A_O_sun * metallicity;
    double reset_free_oxygen =
        reset_oxygen_total - abundance[sCO] - 2.0 * abundance[sCO2];
    abundance[sH2O] = reset_free_oxygen / (1.0 + 4.0 * reset_trace_fraction);
    abundance[sOI] = abundance[sOH] = abundance[sO2] =
        reset_trace_fraction * abundance[sH2O];
    abundance[sCOII] = abundance[sCH] = abundance[sCH2] = 0.0;
    abundance[sOII] = abundance[sOHII] = abundance[sH2OII] = 0.0;
    abundance[sH3OII] = abundance[sO2II] = 0.0;
    double reset_available_hydrogen = 1.0 - abundance[sHM]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                                      - abundance[sHDI]
#endif
                                      - abundance[sOH] - 2.0 * abundance[sH2O];
    abundance[sH2I] = reset_available_hydrogen / (2.0 + 2.0 * reset_h_trace);
    abundance[sHI] = abundance[sHII] = reset_h_trace * abundance[sH2I];
    abundance[sELEC] = dmax(abundance[sHII] + abundance[sCII]
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
                                + abundance[sDII]
#endif
                                - abundance[sHM],
                            1.e-80);
    break;
  }
  default:
    break;
  }
}
#endif
