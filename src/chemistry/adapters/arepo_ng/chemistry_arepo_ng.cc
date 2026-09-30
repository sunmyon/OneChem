#include "chemistry_arepo_ng.h"
#include "chemistry_model.hpp"
#include "thermodynamics/thermodynamics.h"
#include <algorithm>
#include <cmath>
#include <climits>
#include <numeric>
#include <limits>
#include <type_traits>

namespace chemistry {
static_assert(std::is_trivially_copyable<ArepoNgCell>::value, "NG stores cells by value");
static_assert(std::is_trivially_copyable<ArepoNgWork>::value, "MPI wire packet must be a value");

double ArepoNgCell::number_density(int species, double density_cgs) const {
  if(species<0 || species>=TGCHEM_NUM_ABUNDANCES || !std::isfinite(density_cgs) || density_cgs<0)
    return std::numeric_limits<double>::quiet_NaN();
  return abundance[species] * (density_cgs/PROTONMASS);
}

void ArepoNgCell::store_conserved(double mass) {
  for(int i=0; i<TGCHEM_NUM_ABUNDANCES; ++i) species_mass[i] = mass * abundance[i];
  metal_mass = mass * metallicity;
}
void ArepoNgCell::initialize(double mass, double z) {
  *this = {};
  metallicity = z;
  abundance[sHI] = CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON;
#ifdef CHEMISTRY_SPECIES_HELIUM
  abundance[sHeI] = CHEMISTRY_HELIUM_NUCLEI_PER_BARYON;
#endif
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  abundance[sDI] = CHEMISTRY_DEUTERIUM_NUCLEI_PER_BARYON;
#endif
#ifdef CHEMISTRY_SPECIES_CARBON_OXYGEN
  abundance[sCI] = CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON * A_C_sun * z;
  abundance[sOI] = CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON * A_O_sun * z;
#endif
  gamma = 5.0/3.0;
  network_signature=arepo_ng_network_signature();
  store_conserved(mass);
}
bool ArepoNgCell::update_thermodynamics(double u) {
  if(!(u>0) || !std::isfinite(u)) return false;
  double sum=0;
  for(double a: abundance) { if(a<0 || !std::isfinite(a)) return false; sum+=a; }
  if(!(sum>0) || !(abundance[sHI]+abundance[sHII]+abundance[sHM]+2*abundance[sH2I]>0)) return false;
  // Solve T = (gamma(T)-1) mu m_p u/k_B by bracketing; do not rely on stale gamma.
  const double coefficient = PROTONMASS/BOLTZMANN * chemistry_molecular_weight(abundance) * u;
  double lo=0.1*coefficient, hi=coefficient;
  for(int it=0; it<80; ++it) {
    const double t=0.5*(lo+hi);
    const double g=chemistry_gamma_from_abundances(abundance,t);
    if(t > (g-1)*coefficient) hi=t; else lo=t;
  }
  gamma=chemistry_gamma_from_abundances(abundance,0.5*(lo+hi));
  return std::isfinite(gamma) && gamma>1 && gamma<2;
}
double ArepoNgCell::energy_at_temperature(double temperature) const {
  if(temperature<=0) return 0;
  const double g=chemistry_gamma_from_abundances(abundance,temperature);
  return BOLTZMANN*temperature/((g-1)*PROTONMASS*chemistry_molecular_weight(abundance));
}
bool ArepoNgCell::recover_composition(double mass) {
  if(network_signature!=arepo_ng_network_signature()) return false;
  if(!(mass>0) || !std::isfinite(mass)) return false;
  for(int i=0; i<TGCHEM_NUM_ABUNDANCES; ++i) abundance[i]=species_mass[i]/mass;
  metallicity=metal_mass/mass;
  for(double a:abundance) if(!std::isfinite(a) || a<0) return false;
  return std::isfinite(metallicity) && metallicity>=0;
}
bool ArepoNgCell::recover_conserved(double mass, double u) {
  return recover_composition(mass) && update_thermodynamics(u);
}

struct ArepoNgAdapter::Implementation {
  Model::Handle model;
  ChemWorkspace workspace{};
  MPI_Comm comm=MPI_COMM_NULL;
  int rank=0, size=1;
};
ArepoNgAdapter::ArepoNgAdapter(): implementation_(new Implementation) {}
ArepoNgAdapter::~ArepoNgAdapter() = default;
ChemistryStatus ArepoNgAdapter::initialize(const ChemistryConfig &config, MPI_Comm comm) {
  auto &p=*implementation_;
  p.comm=comm;
  MPI_Comm_rank(comm,&p.rank); MPI_Comm_size(comm,&p.size);
  ChemistryConfig local=config; local.print_configuration=config.print_configuration && p.rank==0;
  auto created=Model::create(local);
  int status=created.status, global=0;
  MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,comm);
  if(global) return static_cast<ChemistryStatus>(global);
  p.model=std::move(created.model);
  return CHEMISTRY_STATUS_OK;
}
ChemistryStatus ArepoNgAdapter::advance(const ChemistryGlobalEnvironment &global,
                                      std::vector<ArepoNgWork> &work, bool balance) {
  auto &p=*implementation_;
  ChemistryStepContext step{};
  int status=p.model ? Driver::prepare(global,step) : CHEMISTRY_STATUS_INVALID_ARGUMENT;
  int failed=0;
  MPI_Allreduce(&status,&failed,1,MPI_INT,MPI_MAX,p.comm);
  if(failed) return static_cast<ChemistryStatus>(failed);
  // Limit individual MPI byte counts/displacements independently of local cell count.
  const size_t batch=std::min<size_t>(4096,INT_MAX/sizeof(ArepoNgWork)/p.size);
  if(!batch) return CHEMISTRY_STATUS_INVALID_ARGUMENT;
  auto updated=work;
  unsigned long long local_rounds=(work.size()+batch-1)/batch, rounds=0;
  MPI_Allreduce(&local_rounds,&rounds,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,p.comm);
  auto solve=[&](ArepoNgWork &w) {
    ChemistryState s{w.abundance,TGCHEM_NUM_ABUNDANCES,w.internal_energy_cgs,w.gamma};
    const ChemistryAdvanceOptions options{CHEMISTRY_MODE_DEFAULT,CHEMISTRY_SOLVER_AUTO};
    w.environment.diagnostic_id=w.id;
    w.status=Driver{}.advance(*p.model,step,p.workspace,w.environment,w.duration_s,options,s,w.result);
    w.internal_energy_cgs=s.specific_thermal_energy_cgs; w.gamma=s.gamma;
  };
  for(unsigned long long round=0; round<rounds; ++round) {
    size_t begin=std::min<size_t>(round*batch,work.size());
    int n=static_cast<int>(std::min(batch,work.size()-begin));
    if(!balance) {
      for(int i=0;i<n;++i) { solve(updated[begin+i]); status=std::max(status,updated[begin+i].status); }
    } else {
      std::vector<int> counts(p.size),send(p.size),recv(p.size),sd(p.size),rd(p.size);
      MPI_Allgather(&n,1,MPI_INT,counts.data(),1,MPI_INT,p.comm);
      long long prefix=std::accumulate(counts.begin(),counts.begin()+p.rank,0LL);
      for(int i=0;i<n;++i) ++send[(prefix+i)%p.size];
      MPI_Alltoall(send.data(),1,MPI_INT,recv.data(),1,MPI_INT,p.comm);
      for(int r=1;r<p.size;++r) { sd[r]=sd[r-1]+send[r-1]; rd[r]=rd[r-1]+recv[r-1]; }
      std::vector<ArepoNgWork> outgoing(n),incoming(rd.back()+recv.back());
      auto cursor=sd;
      std::vector<int> positions(n);
      for(int i=0;i<n;++i) { positions[i]=cursor[(prefix+i)%p.size]++; outgoing[positions[i]]=work[begin+i]; }
      for(int r=0;r<p.size;++r) { send[r]*=sizeof(ArepoNgWork); recv[r]*=sizeof(ArepoNgWork); sd[r]*=sizeof(ArepoNgWork); rd[r]*=sizeof(ArepoNgWork); }
      MPI_Alltoallv(outgoing.data(),send.data(),sd.data(),MPI_BYTE,incoming.data(),recv.data(),rd.data(),MPI_BYTE,p.comm);
      for(auto &w:incoming) solve(w);
      MPI_Alltoallv(incoming.data(),recv.data(),rd.data(),MPI_BYTE,outgoing.data(),send.data(),sd.data(),MPI_BYTE,p.comm);
      for(int i=0;i<n;++i) {
        const auto &w=outgoing[positions[i]];
        updated[begin+i]=w;
        status=std::max(status,w.status);
        if(w.id!=work[begin+i].id) status=CHEMISTRY_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  MPI_Allreduce(&status,&failed,1,MPI_INT,MPI_MAX,p.comm);
  if(failed) return static_cast<ChemistryStatus>(failed);
  work.swap(updated);
  return CHEMISTRY_STATUS_OK;
}
} // namespace chemistry
