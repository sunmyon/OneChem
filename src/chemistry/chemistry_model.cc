#include "chemistry_model.hpp"
#include "chemistry_constants.h"
#include "H2cooling.h"
#include "partition.h"
#include "cooling_func.h"
#include "dust_cooling.h"
#include "reaction_network.h"
#include "thermodynamics.h"
#include "forward_rates.h"
#include "molecular_cooling.h"
#include "chemistry_tables.h"

#include <cmath>
#include <new>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

namespace {
void print_configuration(const ChemistryConfig &config) noexcept {
  const auto entry = [](const char *name, const char *value) {
    printf("  %-30s: %s\n", name, value);
  };
  const auto toggle = [&](const char *name, bool enabled) {
    entry(name, enabled ? "on" : "off");
  };
  const auto number = [](const char *name, double value) {
    printf("  %-30s: %.8g\n", name, value);
  };
  printf("Chemistry configuration\n[model inputs]\n");
  entry("table directory", config.table_root);
  number("timestep safety factor", config.timescale_safety_factor);
#ifdef METALLICITY
  entry("metallicity", "host input (solar units)");
#else
  entry("metal-dependent processes", "off");
#endif
  printf("[species]\n");
  entry("hydrogen", "H, H+, H-, H2, electrons");
#ifdef CHEMISTRY_SPECIES_HELIUM
  toggle("helium", true);
#else
  toggle("helium", false);
#endif
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  toggle("deuterium", true);
#else
  toggle("deuterium", false);
#endif
#ifdef CHEMISTRY_SPECIES_H2_PLUS
  toggle("H2+", true);
#else
  toggle("H2+", false);
#endif
#ifdef CHEMISTRY_SPECIES_HD_PLUS
  toggle("HD+", true);
#else
  toggle("HD+", false);
#endif
#ifdef CHEMISTRY_SPECIES_D_MINUS
  toggle("D-", true);
#else
  toggle("D-", false);
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  toggle("carbon / oxygen", true);
#else
  toggle("carbon / oxygen", false);
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION
  toggle("C/O high ionization", true);
#else
  toggle("C/O high ionization", false);
#endif
  printf("[reactions]\n");
  entry("recombination", CHEMISTRY_CASE_B ? "CASE B" : "CASE A");
  entry("three-body H2 formation", !CHEMISTRY_THREE_BODY ? "off" : CHEMISTRY_THREE_BODY_ABEL02 ? "ABEL02" : CHEMISTRY_THREE_BODY_FLOW07 ? "FLOW07" : "FORR13");
  entry("H2 + H -> 3H", CHEMISTRY_MARTIN97 ? "MARTIN97" : "legacy density-dependent fit");
  entry("reaction-rate variant", CHEMISTRY_ENZO ? "ENZO" : CHEMISTRY_HIRANO ? "HIRANO / O05" : CHEMISTRY_O05 ? "O05" : CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA ? "Matsukoba" : "standard");
  toggle("H2 formation on dust", CHEMISTRY_H2_FORMATION_ON_DUST);
#ifdef CHEMISTRY_INVERSE_REACTIONS
  toggle("inverse reactions", true);
#else
  toggle("inverse reactions", false);
#endif
  printf("[cooling / opacity]\n");
  toggle("dust cooling", CHEMISTRY_DUST_COOLING);
  toggle("molecular line cooling", CHEMISTRY_MOLECULAR_LINE_COOLING);
  toggle("metal line cooling (SD93)", CHEMISTRY_METAL_LINE_COOLING_SD93);
  toggle("helium continuum cooling", CHEMISTRY_HELIUM_COOLING);
  toggle("equilibrium helium cooling", CHEMISTRY_HELIUM_EQ_COOLING);
  toggle("optically thin H lines", CHEMISTRY_HLINE_COOLING_THIN);
#ifdef CHEMISTRY_CMB_HEATING
  toggle("CMB heating", true);
#else
  toggle("CMB heating", false);
#endif
  entry("H2 cooling", CHEMISTRY_H2_COOLING_GP12 ? "GP12" : "G15");
  entry("continuum cooling", CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA ? "Matsukoba" : "legacy prescription");
  toggle("continuum opacity", CHEMISTRY_CONT_OPACITY);
  entry("H2 line opacity", CHEMISTRY_H2_OPACITY ? "escape-fraction fit" : "off");
  printf("[radiation / shielding]\n");
#if defined(CHEMISTRY_LW_RADIATION)
  toggle("LW radiation", true);
#else
  toggle("LW radiation", false);
#endif
#if defined(CHEMISTRY_LW_RADIATION) && defined(CHEMISTRY_BACKGROUND_LW_RADIATION)
  toggle("background LW radiation", true);
#else
  toggle("background LW radiation", false);
#endif
#if defined(CHEMISTRY_LW_RADIATION) && defined(LW_RADIATION_LOCAL)
  toggle("local LW radiation", true);
#else
  toggle("local LW radiation", false);
#endif
#if defined(CHEMISTRY_LW_RADIATION) && defined(CHEMISTRY_LW_SELF_SHIELDING)
  toggle("LW self-shielding", true);
#else
  toggle("LW self-shielding", false);
#endif
#if defined(IONIZING_PHOTON)
  toggle("photoionization", true);
#else
  toggle("photoionization", false);
#endif
#if defined(IONIZING_PHOTON) && defined(RYTR_CARBON_OXYGEN_IONIZATION)
  toggle("C/O photoionization coupling", true);
#else
  toggle("C/O photoionization coupling", false);
#endif
#if defined(CHEMISTRY_LW_RADIATION) && defined(RYTR_NCO)
  toggle("CO photodissociation coupling", true);
#else
  toggle("CO photodissociation coupling", false);
#endif
#if defined(CHEMISTRY_LW_RADIATION) && defined(RYTR_NHD)
  toggle("HD photodissociation coupling", true);
#else
  toggle("HD photodissociation coupling", false);
#endif
#if defined(RYTR_TAU_DUST)
  toggle("transport dust attenuation", true);
#else
  toggle("transport dust attenuation", false);
#endif
#if defined(THERMAL_FEEDBACK)
  toggle("stellar thermal feedback", true);
#else
  toggle("stellar thermal feedback", false);
#endif
  printf("[thermodynamics]\n");
  entry("H2 partition function", CHEMISTRY_H2_PARTITION_EBM1989 ? "EBM1989" : CHEMISTRY_H2_PARTITION_I87 ? "I87" : CHEMISTRY_H2_PARTITION_R85 ? "R85" : "EBM1989 (Matsukoba fit)");
  toggle("H2 rotational heat capacity", CHEMISTRY_H2_ROTATIONAL_GAMMA);
  printf("[solver policy]\n");
#ifdef CHEMISTRY_HYBRID_IMPLICIT
  entry("chemistry solver", "hybrid implicit");
  entry("energy solver", "hybrid implicit");
#else
#if defined(CHEMISTRY_IMPLICIT) || defined(CHEMISTRY_ENERGY_IMPLICIT)
  entry("chemistry solver", "implicit");
#else
  entry("chemistry solver", "explicit");
#endif
#ifdef CHEMISTRY_ENERGY_IMPLICIT
  entry("energy solver", "implicit");
#else
  entry("energy solver", "explicit");
#endif
#endif
#ifdef CHEMISTRY_SAHA
  toggle("high-density Saha switch", true);
#else
  toggle("high-density Saha switch", false);
#endif
#ifdef CHEMISTRY_EXPLICIT_SEQUENTIAL
  entry("explicit update", "sequential");
#else
  entry("explicit update", "simultaneous");
#endif
#if defined(CHEMISTRY_DIFF_ANALYTIC) && !defined(CHEMISTRY_DIFF_NUMERICAL)
  entry("Newton Jacobian", "analytic (numerical fallback)");
#else
  entry("Newton Jacobian", "numerical");
#endif
#ifdef CHEMISTRY_INTERMEDIATE_EQUILIBRIUM
  entry("Newton system", "reduced with equilibrium species");
#elif defined(CHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION)
  entry("Newton system", "full with equilibrium initial guess");
#else
  entry("Newton system", "full");
#endif
#ifdef CHEMISTRY_NEWTON_LU_REUSE
  toggle("Newton LU reuse", true);
#else
  toggle("Newton LU reuse", false);
#endif
#ifdef CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE
  toggle("equilibrium convergence", true);
#else
  toggle("equilibrium convergence", false);
#endif
#ifdef CHEMISTRY_ENERGY_BRACKET_FALLBACK
  toggle("energy bracket fallback", true);
#else
  toggle("energy bracket fallback", false);
#endif
#ifdef CHEMISTRY_NORMALIZE_ABUNDANCES
  toggle("input abundance normalization", true);
#else
  toggle("input abundance normalization", false);
#endif
  toggle("update core length", CHEMISTRY_UPDATE_CORE_LENGTH);
  number("Newton residual gate", CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT);
#if defined(CHEMISTRY_SPECIES_CARBON_OXYGEN) && !defined(CHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION)
  toggle("C/O conservation constraints", true);
#else
  toggle("C/O conservation constraints", false);
#endif
#ifdef CHEMISTRY_SAHA
  number("Saha density threshold [cm^-3]", N_EQCHEM);
#endif
#ifdef CHEMISTRY_JEANS_TEMPERATURE_MAX
  number("Jeans temperature limit [K]", CHEMISTRY_JEANS_TEMPERATURE_MAX);
#else
  entry("Jeans temperature limit [K]", "none");
#endif
#ifdef CHEMISTRY_FIXED_TEMPERATURE
  toggle("fixed temperature", true);
#else
  toggle("fixed temperature", false);
#endif
  printf("[diagnostics]\n");
  toggle("detailed diagnostics", CHEMISTRY_DEBUG);
#ifdef CHEMISTRY_PROFILE
  toggle("profiling", true);
#else
  toggle("profiling", false);
#endif
#ifdef CHEMISTRY_RANDOM_RESTART_DIAGNOSTIC
  toggle("random-restart diagnostics", true);
#else
  toggle("random-restart diagnostics", false);
#endif
#ifdef CHEMISTRY_SAHA_DIAGNOSTICS
  toggle("Saha diagnostics", true);
#else
  toggle("Saha diagnostics", false);
#endif
#ifdef CHEMISTRY_WRITE_DIAGNOSTIC_TABLES
  toggle("diagnostic tables", true);
#else
  toggle("diagnostic tables", false);
#endif
  printf("\n");
}
} // namespace

namespace chemistry {
LogGridView::LogGridView(const double *nodes, int count, double log_min,
                         double log_max, double log_step) noexcept
    : nodes_(nodes), count_(count), log_min_(log_min), log_max_(log_max),
      log_step_(log_step) {}

LogGridView::Sample LogGridView::locate(double value) const noexcept {
  double logarithm = log(value);
  if(logarithm < log_min_) logarithm = log_min_;
  if(logarithm > log_max_) logarithm = log_max_;
  int lower = (int)((logarithm - log_min_) / log_step_);
  if(lower < 0) lower = 0;
  if(lower > count_ - 2) lower = count_ - 2;
  return {lower, (value - nodes_[lower]) / (nodes_[lower + 1] - nodes_[lower])};
}

RateTableView::RateTableView(const ChemistryRateTables &tables,
                             double temperature) noexcept
    : tables_(tables), sample_(LogGridView(tables.temperature, N_T,
          log(T_start), log(T_end), tables.dln_temperature).locate(temperature)) {}

RateTables::RateTables() noexcept = default;
RateTables::~RateTables() noexcept = default;
ChemistryStatus RateTables::build() noexcept {
  std::unique_ptr<ChemistryRateTables> candidate(
      new(std::nothrow) ChemistryRateTables{});
  if(!candidate) {
    fprintf(stderr, "Failed to allocate chemistry rate tables.\n");
    return CHEMISTRY_STATUS_OUT_OF_MEMORY;
  }
  if(tabulate(candidate.get()) <= 0) {
    fprintf(stderr, "Failed to tabulate chemistry rates.\n");
    return CHEMISTRY_STATUS_INTERNAL_ERROR;
  }
  data_ = std::move(candidate);
  return CHEMISTRY_STATUS_OK;
}

#if !CHEMISTRY_MOLECULAR_LINE_COOLING
void ModelResourceDeleter::operator()(ChemistryCoolingTables *) const noexcept {}
#endif
#if !CHEMISTRY_DUST_COOLING
void ModelResourceDeleter::operator()(ChemistryDustTables *) const noexcept {}
#endif
#if !CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA && !CHEMISTRY_METAL_LINE_COOLING_SD93
void ModelResourceDeleter::operator()(ChemistryContinuumOpacityTables *) const noexcept {}
#endif

ChemistryStatus Model::initialize(const char *table_root) noexcept {
  // Keep the historical construction order; partial failure unwinds via RAII.
  network_ = load_network(table_root);
  if(!network_) return CHEMISTRY_STATUS_IO_ERROR;
#if CHEMISTRY_MOLECULAR_LINE_COOLING
  molecular_cooling_ = load_molecular_cooling(table_root);
  if(!molecular_cooling_) return CHEMISTRY_STATUS_IO_ERROR;
#endif
  const auto status = rates_.build();
  if(status != CHEMISTRY_STATUS_OK) return status;
#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA || CHEMISTRY_METAL_LINE_COOLING_SD93
  continuum_opacity_ = load_continuum_opacity(table_root);
  if(!continuum_opacity_) return CHEMISTRY_STATUS_IO_ERROR;
#endif
#if CHEMISTRY_DUST_COOLING
  dust_opacity_ = load_dust_opacity(table_root);
  if(!dust_opacity_) return CHEMISTRY_STATUS_IO_ERROR;
#endif
  view_.network = network_.get();
  view_.rate_tables = rates_.data();
  view_.cooling_tables = molecular_cooling_.get();
  view_.dust_tables = dust_opacity_.get();
  view_.continuum_opacity_tables = continuum_opacity_.get();
  return CHEMISTRY_STATUS_OK;
}

Model::Creation Model::create(const ChemistryConfig &config) noexcept {
  if(config.abi_version != CHEMISTRY_ABI_VERSION ||
     config.table_root == nullptr || config.table_root[0] == '\0' ||
     !(config.timescale_safety_factor > 0.) ||
     !std::isfinite(config.timescale_safety_factor))
    return {nullptr, CHEMISTRY_STATUS_INVALID_ARGUMENT};
  std::unique_ptr<Model> model(new(std::nothrow) Model);
  if(!model) return {nullptr, CHEMISTRY_STATUS_OUT_OF_MEMORY};
  const auto status = model->initialize(config.table_root);
  if(status != CHEMISTRY_STATUS_OK) return {nullptr, status};
  model->view_.timescale_safety_factor = config.timescale_safety_factor;
#ifdef CHEMISTRY_FIXED_TEMPERATURE
#ifndef CHEMISTRY_FIXED_TEMPERATURE_K
#define CHEMISTRY_FIXED_TEMPERATURE_K 10000.
#endif
  model->view_.fixed_temperature = CHEMISTRY_FIXED_TEMPERATURE_K;
  // The high-temperature equilibrium shortcut evolves energy; this fixed-T
  // option deliberately supports the low-temperature reaction network only.
  if(!(model->view_.fixed_temperature > 0) || model->view_.fixed_temperature > 1.e5 ||
     !std::isfinite(model->view_.fixed_temperature))
    return {nullptr, CHEMISTRY_STATUS_INVALID_ARGUMENT};
#endif
  if(config.print_configuration) print_configuration(config);
  return {std::move(model), CHEMISTRY_STATUS_OK};
}

Model::~Model() noexcept = default;


} // namespace chemistry


FILE *open_chemistry_table(const char *table_root, const char *filename)
{
  char path[2048];

  if(table_root == NULL || table_root[0] == '\0' || filename == NULL ||
     filename[0] == '\0') {
    fprintf(stderr, "Chemistry table root and filename must be provided.\n");
    return NULL;
  }

  const int written = snprintf(path, sizeof(path), "%s/tables/%s",
                               table_root, filename);
  if(written < 0 || (size_t)written >= sizeof(path)) {
    fprintf(stderr, "Chemistry table path is too long: root=%s file=%s\n",
            table_root, filename);
    return NULL;
  }

  FILE *fp = fopen(path, "r");
  if(fp == NULL)
    fprintf(stderr, "Can't open chemistry table %s.\n", path);

  return fp;
}

  /*--- [List of Chemistry] ---*/
  /*
    ------- 1:	HI    + e   -> HII   + 2e    hydrogen collisional ionization
    ------- 2:	HII   + e   -> HI    + p     hydrogen radiative recombination
    ------- 3:	HeI   + e   -> HeII  + 2e    He collisional ionization
    ------- 4:	HeII  + e   -> HeI   + p     He recombination
    ------- 5:	HeII  + e   -> HeIII + 2e    HeII collisional ionization
    ------- 6:	HeIII + e   -> HeII  + p     HeII recombination
    ------- 7:	HI    + e   -> HM    + p     photo-attachment
    ------- 8:	HM    + HI  -> H2I*  + e     H2 formation path (H-)
    ------- 9:  HI    + HII -> H2II  + p     H2+ channel
    ------- 10:	H2II  + HI  -> H2I*  + HII   H2+ channel
    ------- 11:	H2I   + HI  -> 3HI
    ------- 12:	H2I   + HII -> H2II  + HI    charge transfer
    ------- 13:	H2I   + e   -> 2HI   + e     collisional dissociation
    ------- 14: HeII  + HI  -> HeI   + HII + p  charge transfer
    ------- 15: HeI   + HII -> HeII  + HI
    ------- 16:	HM    + e   -> HI    + 2e
    ------- 17:	HM    + HI  -> 2H    + e
    ------- 18:	HM    + HII -> 2HI
    ------- 19:	HM    + HII -> H2II  + e
    ------- 20:	H2II  + e   -> 2HI
    ------- 21:	H2II  + HM  -> HI    + H2I

    ------- 22:	HI    + p   -> HII   + e     photo-ionization H
    ------- 23:	HeII  + p   -> HeIII + e     HeII ionization
    ------- 24:	HeI   + p   -> HeII  + e     photo-ionization He
    ------- 25:	HM    + p   -> HI    + e     photo-detachment
    ------- 26:	H2II  + p   -> HI    + HII
    ------- 27:	H2I   + p   -> H2II  + e
    ------- 28:	H2II  + p   -> 2HII  + e
    ------- 29:	H2I   + p   -> 2HI           photo-dissociation

    ------- 60,61: H2    + e   -> HM    + HI
    ------- 62,63: HI    + HI  -> H2II  + e

    HD (Nakamura-Umemura complete model)

    ------- 31:  DII  + e    -> DI   + p      radiative recombination
    ------- 32:  DI   + HII  -> DII  + HI     charge transfer
    ------- 33:  DII  + HI   -> DI   + HII    reverse reaction of 32
    ------- 34:  DI   + HI   -> HDI  + p      * neutral formation path
    ------- 35:  DI   + H2I  -> HI   + HDI    * neutral path by H2
    ------- 36:  HDII + HI   -> HII  + HDI    charge transfer
    ------- 37:  DII  + H2I  -> HII  + HDI    * D+ formation path
    ------- 38:  HDI  + HI   -> H2I  + DI
    ------- 39:  HDI  + HII  -> H2I  + DII
    ------- 40:  DI   + HII  -> HDII + p
    ------- 41:  DII  + HI   -> HDII + p
    ------- 42:  HDII + e    -> HI   + DI     collisional dissociation
    ------- 43:  DI   + e    -> DM   + p      D- formation/photo-attachment
    ------- 44:  DII  + DM   -> 2DI           mutual neutrization
    ------- 45:  HII  + DM   -> DI   + HI     mutual neutrization
    ------- 46:  HM   + DI   -> HI   + DM
    ------- 47:  DM   + HI   -> DI   + HM
    ------- 48:  DM   + HI   -> HDI  + e      * D- formation path

    ------- 50:  DI   + p    -> DII  + e      photo-ionization D
    ------- 51:  HDII + p    -> HI   + DII    photo-dissociation HD+
    ------- 52:  HDII + p    -> HII  + DI     photo-dissociation HD+
    ------- 53:  HDI  + p    -> HI   + DI     photo-dissociation HD

    3 body reactions

    ------- 55:  HI + HI + HI -> H2I + HI
    ------- 56:  H2I + HI     -> HI + HI + HI  same as reaction 11
    ------- 57:  HI + HI + H2I-> H2I + H2I
    ------- 58:  H2I + H2I    -> HI + HI + H2I  reverse reaction
  */

  /* probably this part should be parallelized in the near future */
int chemistry::RateTables::tabulate(ChemistryRateTables *tables)
{
  enum SPECIES species;

  if(tables == NULL)
    return 0;

  double *T_table = tables->temperature;
  double (*k_table)[N_T] = tables->forward;
  double (*k_table_inverse)[N_T] = tables->inverse;
  double *chem_e = tables->formation_energy_kj_mol;
  double *ceHIa = tables->collisional_excitation_HI;
  double *ceHeIa = tables->collisional_excitation_HeI;
  double *ceHeIIa = tables->collisional_excitation_HeII;
  double *ciHIa = tables->collisional_ionization_HI;
  double *ciHeIa = tables->collisional_ionization_HeI;
  double *ciHeIIa = tables->collisional_ionization_HeII;
  double *ciHeISa = tables->collisional_ionization_HeI_singlet;
  double *reHIIa = tables->recombination_HII;
  double *reHeII1a = tables->recombination_HeII_radiative;
  double *reHeII2a = tables->recombination_HeII_dielectronic;
  double *reHeIIIa = tables->recombination_HeIII;
  double *brema = tables->bremsstrahlung;
  double *H2_CIE = tables->H2_CIE;

  chem_e[sELEC] = 0.;
  chem_e[sHI]   = 216.034;
  chem_e[sHII]  = 1528.084;
  chem_e[sH2I]  = 0.;
#ifdef CHEMISTRY_SPECIES_H2_PLUS
  chem_e[sH2II] = 1488.365;
#endif
  chem_e[sHM]   = 143.246;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  chem_e[sDI]   = 219.804;
  chem_e[sDII]  = 1532.211;
  chem_e[sHDI]  = 0.328;
#endif

#ifdef CHEMISTRY_WRITE_DIAGNOSTIC_TABLES
  check_partition_functions_H2();
#endif

  /* ----- Compute temperature grid. */
  tables->dln_temperature = (log(T_end) - log(T_start))/N_T;
  for(int i = 0;i< N_T;i++){
    double ln_T = log(T_start) + i*tables->dln_temperature;
    T_table[i] = exp(ln_T);
  }

  /* ----- Initialize arrays. */
  for(int i=0;i<N_T;i++){
    for(int j=0;j<CHEM_RATE_COUNT;j++) k_table[j][i] = 0;
    for(int j=0;j<CHEM_RATE_COUNT;j++) k_table_inverse[j][i] = 0;

    ceHIa[i]    = 0.;
    ceHeIa[i]   = 0.;
    ceHeIIa[i]  = 0.;
    ciHIa[i]    = 0.;
    ciHeIa[i]   = 0.;
    ciHeISa[i]  = 0.;
    ciHeIIa[i]  = 0.;
#if !CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA
    reHIIa[i]   = 0.;
#ifdef CHEMISTRY_SPECIES_HELIUM
    reHeII1a[i] = 0.;
    reHeII2a[i] = 0.;
    reHeIIIa[i] = 0.;
#endif
    brema[i]    = 0.;
    H2_CIE[i]   = 0.;
#endif
  }

  /* ----- Fill in tables over the range T_start to T_end. */
  for(int i=0 ; i<N_T ; i++){
    double z_H2I = partition_H2(T_table[i]);
    double z_HI  = partition_H(T_table[i]);
#ifdef CHEMISTRY_SPECIES_H2_PLUS
    double z_H2II= pf_z_H2p(T_table[i]);
#endif
    double z_HM  = 1.;
    double z_elec = 2.;
    double z_HII = 1.;

    /* Pure forward reaction rates are selected in forward_rate_table. */
    for(int j=0; j<CHEM_RATE_COUNT; j++) {
      ForwardRateFunction evaluate = forward_rate_table[j];
      if(evaluate != NULL)
        k_table[j][i] = evaluate(T_table[i]);
    }

#ifdef CHEMISTRY_INVERSE_REACTIONS
    double K_T;
    double Ene, dE, Z_pf, Mass;
    double m_HI, m_HII, m_HM, m_elec, m_H2I, m_H2II;

    m_HI = m_HII = m_HM = PROTONMASS;
    m_H2I = m_H2II = 2. * m_HI;
    m_elec = ELECTRONMASS;

    Ene = 2 * M_PI * BOLTZMANN * T_table[i] / PLANCK / PLANCK;

    /**** HI + e- -> Hp + 2e-****/
    if(T_table[i] > 1.e3){
      dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHII] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_elec/ z_HII / z_elec / z_elec;
      Mass = m_HI * m_elec/ m_HII / m_elec / m_elec;

      K_T = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[1][i] = k_table[1][i] * K_T;
    }

    /**** Hp + e- -> HI  ****/
    dE   = (chem_e[sHII] + chem_e[sELEC] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HII * z_elec / z_HI;
    Mass = m_HII * m_elec / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[2][i] = k_table[2][i] * K_T; //note this value is optically thick limit.

    /**** HI + e- -> HM  ****/
    dE   = (chem_e[sHI] + chem_e[sELEC] - chem_e[sHM]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_elec / z_HM;
    Mass = m_HI * m_elec / m_HM;
    K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[7][i] = k_table[7][i] * K_T; //note this value is optically thick limit.

    /**** HI + HM -> H2I + e-  ****/
    dE   = (chem_e[sHI] + chem_e[sHM] - chem_e[sH2I] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HI * z_HM / z_H2I / z_elec;
    Mass = m_HI * m_HM / m_H2I / m_elec;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[8][i] = k_table[8][i] * K_T;

#ifdef CHEMISTRY_SPECIES_H2_PLUS
    /**** H2I + Hp -> H2II + HI  ****/
    dE   = (chem_e[sH2I] + chem_e[sHII] - chem_e[sH2II] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_H2I * z_HII / z_HI / z_H2II;
    Mass = m_H2I * m_HII / m_HI / m_H2II;
    /* Combining in log space avoids underflow * overflow at low T.
     * The extreme table tail can exceed double range even in log space;
     * saturate it with headroom for interpolation, not a NaN/inf table entry.
     * This numerical ceiling is not a physical extrapolation of the fit. */
    const double log_k10 = chemistry_log_forward_rate_12(T_table[i]) +
        log(Z_pf) + 1.5 * log(Mass) - dE/T_table[i];
    k_table[10][i] = exp(fmin(log_k10, log(DBL_MAX) - 2.0));
#endif

    if(T_table[i]>200.){
      /**** H2I + e- -> 2HI + e-  ****/
      dE   = (chem_e[sH2I] + chem_e[sELEC] - 2*chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_elec / z_HI / z_HI / z_elec;
      Mass = m_H2I * m_elec / m_HI / m_HI / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[13][i] = k_table[13][i] * K_T;
    }

    if(T_table[i]>200.){
      /**** HM + e- -> HI + 2e-  ****/
      dE   = (chem_e[sHM] + chem_e[sELEC] - chem_e[sHI] - 2*chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HM * z_elec / z_HI / z_elec / z_elec;
      Mass = m_HM * m_elec / m_HI / m_elec / m_elec;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
#ifdef AIU
      k_table_inverse[16][i] = k_table[16][i] * K_T;
#endif
    }

    /**** HM + Hp -> HI + HI  ****/
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    k_table_inverse[18][i] = k_table[18][i] * K_T;

#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA || CHEMISTRY_METAL_LINE_COOLING_SD93
    dE   = (chem_e[sHM] + chem_e[sHII] - 2*chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
    Z_pf = z_HM * z_HII / z_HI / z_HI;
    Mass = m_HM * m_HII / m_HI / m_HI;
    K_T  = Z_pf * pow(Mass, 1.5) * exp(-dE/T_table[i]);
    //k_table_HI_inverse[1][i] = k_table_HI[1][i] * K_T;
#endif

    if(T_table[i] > 1.e3){
      /***** HI + HI -> Hp + HI + e- *****/
      dE   = (2*chem_e[sHI] - chem_e[sHII] - chem_e[sHI] - chem_e[sELEC]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI / z_HII / z_elec / z_HI;
      Mass = m_HI * m_HI / m_HII / m_elec / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table_inverse[30][i] = k_table[30][i] * K_T;
    }

    if(T_table[i] > 500.){
      /****  HI + HI + HI -> HI + H2I  ****/
      dE   = (3*chem_e[sHI] - chem_e[sH2I] - chem_e[sHI]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_HI * z_HI * z_HI / z_H2I / z_HI;
      Mass = m_HI * m_HI * m_HI / m_H2I / m_HI;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, 1.5) * exp(-dE/T_table[i]);
      k_table_inverse[55][i] = k_table[55][i] * K_T;

      /****  H2I + H2I -> H2I + HI + HI  ****/
      dE   = (2*chem_e[sH2I] - 2*chem_e[sHI] - chem_e[sH2I]) * kJ_per_mol_to_eV *eV_to_erg/BOLTZMANN;
      Z_pf = z_H2I * z_H2I / z_HI / z_HI / z_H2I;
      Mass = m_H2I * m_H2I / m_HI / m_HI / m_H2I;
      K_T  = Z_pf * pow(Mass, 1.5) * pow(Ene, -1.5) * exp(-dE/T_table[i]);
      k_table[57][i] = k_table[58][i] * K_T;
    }
#endif
  }

  /* coefficients for cooling and heating */
  for(int i=0; i<N_T; i++){
#if CHEMISTRY_CASE_B
    double lambda_HI   = 2 * 157807.0/T_table[i];
#if defined(CHEMISTRY_SPECIES_HELIUM) || CHEMISTRY_HELIUM_EQ_COOLING
    double lambda_HeI  = 2 * 285335.0/T_table[i];
    double lambda_HeII = 2 * 631515.0/T_table[i];
#endif
#endif

    /* ------- Collisional excitations (Black 1981; Cen 1992) */
    ceHIa[i]   = 7.5e-19*exp( -118348./T_table[i] )
      /(1.+sqrt(T_table[i]/1.e5));

    ceHeIa[i]   = 9.1e-27*exp( -13179./T_table[i] )
      *pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5));

    ceHeIIa[i]  = 5.54e-17*exp( -473638./T_table[i] )
      *pow(T_table[i],-0.397)/(1.+sqrt(T_table[i]/1.e5));

    /* ------- Collisional ionizations (Tom's polynomial fits) */
    ciHIa[i]    = 2.18e-11*k_table[1][i];
    ciHeIa[i]   = 3.94e-11*k_table[3][i];
    ciHeIIa[i]  = 8.72e-11*k_table[5][i];

    ciHeISa[i]  = 5.01e-27*pow(T_table[i],-0.1687)/(1.+sqrt(T_table[i]/1.e5))
      *exp( -55338./T_table[i] );

    /* ------- Recombinations */
#if CHEMISTRY_CASE_B
    reHIIa[i]   = 3.435e-30 * T_table[i] * pow(lambda_HI, 1.97)/pow(1.0+pow(lambda_HI/2.25, 0.376), 3.72);

#if defined(CHEMISTRY_SPECIES_HELIUM) || CHEMISTRY_HELIUM_EQ_COOLING
    reHeII1a[i] = BOLTZMANN * T_table[i] * 1.26e-14 * pow(lambda_HeI, 0.75);

    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */

    reHeIIIa[i] = 8 * 3.435e-30 * T_table[i] * pow(lambda_HeII, 1.97)/pow(1.0+pow(lambda_HeII/2.25, 0.376), 3.72);
#endif
#else
    reHIIa[i]   = 6.50e-27*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/1.e6,0.7));

#if defined(CHEMISTRY_SPECIES_HELIUM) || CHEMISTRY_HELIUM_EQ_COOLING
    reHeII1a[i] = 1.55e-26*pow(T_table[i],0.3647);

    reHeII2a[i] = 1.24e-13*pow(T_table[i],-1.5)*exp( -470000/T_table[i] )
      *(1.+0.3*exp( -94000./T_table[i] ));     /* He dielectric recombination */

    reHeIIIa[i] = 3.48e-26*sqrt(T_table[i])*pow(T_table[i]/1000.0,-0.2)
      / (1.0 + pow(T_table[i]/4.e6,0.7));
#endif

    /* ------- Bremsstrahlung (Black 1981)(Spitzer & Hart 1979) */
    brema[i]    = 1.43e-27*sqrt(T_table[i])
      *(1.1+0.34*exp(-pow(5.5-log10(T_table[i]),2.0)/3.0));

    const double log10_T = log10(T_table[i]);
    H2_CIE[i]   = pow(10.0, - 116.6
		      + 96.3   *log10_T
		      - 47.15  *log10_T*log10_T
		      + 10.744 *log10_T*log10_T*log10_T
		      -  0.9155*log10_T*log10_T*log10_T*log10_T);
#endif
  }

  return N_T;
}


void chemistry_config_defaults(ChemistryConfig *config)
{
  if(config == NULL)
    return;

  memset(config, 0, sizeof(*config));
  config->abi_version = CHEMISTRY_ABI_VERSION;
  config->timescale_safety_factor = 0.1;
  config->print_configuration = 1;
}

const char *chemistry_status_string(ChemistryStatus status)
{
  switch(status) {
    case CHEMISTRY_STATUS_OK:
      return "ok";
    case CHEMISTRY_STATUS_INVALID_ARGUMENT:
      return "invalid argument";
    case CHEMISTRY_STATUS_INCOMPATIBLE_STATE:
      return "incompatible state";
    case CHEMISTRY_STATUS_IO_ERROR:
      return "I/O error";
    case CHEMISTRY_STATUS_OUT_OF_MEMORY:
      return "out of memory";
    case CHEMISTRY_STATUS_NOT_CONVERGED:
      return "not converged";
    case CHEMISTRY_STATUS_INTERNAL_ERROR:
      return "internal error";
  }
  return "unknown chemistry status";
}
