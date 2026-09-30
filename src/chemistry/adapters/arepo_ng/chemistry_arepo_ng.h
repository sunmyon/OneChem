#pragma once

#include "species.h"
#include "include/chemistry/chemistry.h"
#include <mpi.h>
#include <cstdint>
#include <memory>
#include <vector>

namespace chemistry {
constexpr std::uint64_t arepo_ng_network_signature() {
#define CHEMISTRY_NG_NAME(symbol, ...) #symbol ";"
  constexpr char names[] = CHEM_SPECIES_LIST(CHEMISTRY_NG_NAME);
#undef CHEMISTRY_NG_NAME
  std::uint64_t hash=14695981039346656037ULL;
  for(char c:names) { hash^=static_cast<unsigned char>(c); hash*=1099511628211ULL; }
  return hash;
}
// Value storage only: safe in NG particle containers, ghost exchange and restart.
struct ArepoNgCell {
  double abundance[TGCHEM_NUM_ABUNDANCES]; // n_i / (rho / proton mass)
  double species_mass[TGCHEM_NUM_ABUNDANCES]; // gas mass times abundance
  double metallicity; // solar units
  double metal_mass; // gas mass times metallicity
  double gamma;
  std::uint64_t network_signature;
  void initialize(double mass, double metallicity_solar);
  double number_density(int species, double density_cgs) const; // cm^-3, chemistry proton-mass convention
  bool update_thermodynamics(double internal_energy_cgs);
  double energy_at_temperature(double temperature) const;
  void store_conserved(double mass);
  bool recover_composition(double mass);
  bool recover_conserved(double mass, double internal_energy_cgs);
};

// Wire format deliberately contains no pointers, host indices or MPI handles.
struct ArepoNgWork {
  ChemistryCellEnvironment environment;
  double abundance[TGCHEM_NUM_ABUNDANCES];
  double internal_energy_cgs;
  double gamma;
  double duration_s;
  std::uint64_t id;
  ChemistryResult result;
  int status;
};

class ArepoNgAdapter final {
public:
  ArepoNgAdapter();
  ~ArepoNgAdapter();
  ArepoNgAdapter(const ArepoNgAdapter &) = delete;
  ArepoNgAdapter &operator=(const ArepoNgAdapter &) = delete;
  ChemistryStatus initialize(const ChemistryConfig &, MPI_Comm);
  // Collective on the simulation communicator, including ranks with zero work.
  // On failure, every rank retains its original work vector. MPI remains host-owned.
  ChemistryStatus advance(const ChemistryGlobalEnvironment &, std::vector<ArepoNgWork> &, bool balance);
private:
  struct Implementation;
  std::unique_ptr<Implementation> implementation_;
};
} // namespace chemistry

class BaseModule;
struct sph_particle_data;
// Optional host radiation preparation, called before MPI work redistribution.
using ArepoNgEnvironmentProvider = void (*)(const sph_particle_data *, ChemistryCellEnvironment &);
BaseModule *create_arepo_ng_chemistry_module(ArepoNgEnvironmentProvider provider = nullptr);
void chemistry_ng_initialize_particle(sph_particle_data *, bool initialize_energy);
