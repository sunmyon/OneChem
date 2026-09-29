// NG lifecycle/units only. Numerical chemistry and MPI scheduling live in the adapter.
#include "arepoconfig.h"
#include "modules/basemodule.h"
#include "chemistry_arepo_ng.h"
#include <algorithm>
#include <cmath>
#include <cstring>

#if !defined(RIEMANN_HLLD) || defined(ENABLE_ENTROPY_TREATMENT) || defined(COLLAPSE_EOS) || defined(ISOTHERM_EQS) || defined(COOLING) || defined(STARFORMATION)
#error "CHEMISTRY requires HLLD and the thermal-energy equation, without another cooling/EOS module"
#endif

namespace {
double initial_metallicity=0;
double initial_ionized_fraction=1.e-4;
double initial_molecular_fraction=1.e-6;
double unit_energy() { return All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s; }

class ChemistryModule final: public BaseModule {
  chemistry::ArepoNgAdapter adapter;
  char table_root[200]{};
  double cmb_temperature=2.728;
  int distribute=0;
public:
  ChemistryModule() { moduleName="Chemistry"; }
  void onRegisterParameters(tParameters *params) override {
    params->add_param("ChemistryTableRoot",table_root,PARAM_STRING,PARAM_FIXED);
    params->add_param("ChemistryInitialMetallicity",&initial_metallicity,PARAM_DOUBLE,PARAM_FIXED);
    params->add_param("ChemistryInitialIonizedFraction",&initial_ionized_fraction,PARAM_DOUBLE,PARAM_FIXED);
    params->add_param("ChemistryInitialMolecularFraction",&initial_molecular_fraction,PARAM_DOUBLE,PARAM_FIXED);
    params->add_param("ChemistryCMBTemperature",&cmb_temperature,PARAM_DOUBLE,PARAM_FIXED);
    params->add_param("ChemistryBalanceWork",&distribute,PARAM_INT,PARAM_FIXED);
  }
  void onBegrunBeforeInit(sim *simulation) override {
    // Module parameter storage is outside All's automatically broadcast block.
    MPI_Bcast(table_root,sizeof(table_root),MPI_CHAR,0,simulation->Communicator);
    MPI_Bcast(&initial_metallicity,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&initial_ionized_fraction,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&initial_molecular_fraction,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&cmb_temperature,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&distribute,1,MPI_INT,0,simulation->Communicator);
    if(All.ComovingIntegrationOn)
      Terminate("CHEMISTRY: cosmological drift/entropy factors require a variable-gamma formulation; use non-comoving integration");
    if(!std::isfinite(initial_metallicity) || initial_metallicity<0 || !(cmb_temperature>0) || !std::isfinite(cmb_temperature) || (distribute!=0 && distribute!=1))
      Terminate("CHEMISTRY: invalid initial metallicity, CMB temperature or work balance setting");
    if(!std::isfinite(initial_ionized_fraction) || !std::isfinite(initial_molecular_fraction) || initial_ionized_fraction<0 || initial_molecular_fraction<0 || initial_ionized_fraction+2*initial_molecular_fraction>1)
      Terminate("CHEMISTRY: initial n(HII)/nH and n(H2)/nH fractions exceed the hydrogen budget");
    ChemistryConfig config{}; chemistry_config_defaults(&config); config.table_root=table_root;
    const auto status=adapter.initialize(config,simulation->Communicator);
    if(status!=CHEMISTRY_STATUS_OK) Terminate("CHEMISTRY: model initialization failed (%d)",int(status));
  }
  void onBegrunAfterInit(simparticles *particles) override {
    for(gas_iterator p(particles);p();p++) {
      if(!p->Chemistry.recover_conserved(p->getMass(),p->Utherm*unit_energy()))
        Terminate("CHEMISTRY: invalid saved particle composition");
      p->Pressure=p->get_pressure(); p->Csnd=p->get_sound_speed();
    }
  }
  void onSnapshotAddFields(snap_io *snap) override {
    snap->init_field<unsigned long long,sph_particle_data>(
      "ChemistryNetworkSignature","CHNS",GAS_ONLY,READ_IF_PRESENT,1,
      [](sph_particle_data *p,unsigned long long *b) { *b=p->Chemistry.network_signature; },
      [](sph_particle_data *p,unsigned long long *b) { p->Chemistry.network_signature=*b; });
    snap->init_field<FileFloatOrDouble,sph_particle_data>(
      "ChemicalAbundances","CHAB",GAS_ONLY,READ_IF_PRESENT,TGCHEM_NUM_ABUNDANCES,
      [](sph_particle_data *p,FileFloatOrDouble *b) { for(int i=0;i<TGCHEM_NUM_ABUNDANCES;++i) b[i]=p->Chemistry.abundance[i]; },
      [](sph_particle_data *p,FileFloatOrDouble *b) { for(int i=0;i<TGCHEM_NUM_ABUNDANCES;++i) p->Chemistry.abundance[i]=b[i]; });
    snap->init_field<FileFloatOrDouble,sph_particle_data>(
      "ChemistryMetallicity","CHMT",GAS_ONLY,READ_IF_PRESENT,1,
      [](sph_particle_data *p,FileFloatOrDouble *b) { *b=p->Chemistry.metallicity; },
      [](sph_particle_data *p,FileFloatOrDouble *b) { p->Chemistry.metallicity=*b; });
    snap->init_field<FileFloatOrDouble,sph_particle_data>(
      "ChemistryGamma","CHGM",GAS_ONLY,READ_IF_PRESENT,1,
      [](sph_particle_data *p,FileFloatOrDouble *b) { *b=p->Chemistry.gamma; },
      [](sph_particle_data *p,FileFloatOrDouble *b) { p->Chemistry.gamma=*b; });
  }
  void onCalculateNonStandardPhysicsEndOfTimestep(simparticles *particles) override {
    std::vector<sph_particle_data *> owners;
    std::vector<chemistry::ArepoNgWork> work;
    for(gas_iterator_up_to_timebin p(particles,All.HighestSynchronizedTimeBin);p();p++) {
      if(p->getTimeBin()==0 || p->getMass()<=0) continue;
      chemistry::ArepoNgWork w{};
      w.environment.density_cgs=p->Density*All.UnitDensity_in_cgs;
      w.environment.metallicity=p->Chemistry.metallicity;
      std::copy(std::begin(p->Chemistry.abundance),std::end(p->Chemistry.abundance),w.abundance);
      w.internal_energy_cgs=p->Utherm*unit_energy(); w.gamma=p->Chemistry.gamma;
      // NG advances cells on the actual particle bin (min of gravity/hydro),
      // also rescaled by restart::readjust_timebase; TimeBinHydro is a limit.
      w.duration_s=std::ldexp(All.Timebase_interval,p->getTimeBin())*All.UnitTime_in_s;
      w.id=p->ID.get();
      owners.push_back(p()); work.push_back(w);
    }
    ChemistryGlobalEnvironment env{};
    env.scale_factor=1; env.cmb_temperature=cmb_temperature; env.simulation_time=All.Time;
    const auto status=adapter.advance(env,work,distribute!=0);
    if(status!=CHEMISTRY_STATUS_OK) Terminate("CHEMISTRY: advance failed (%d)",int(status));
    for(size_t i=0;i<work.size();++i) {
      auto *p=owners[i]; const auto &w=work[i];
      std::copy(std::begin(w.abundance),std::end(w.abundance),p->Chemistry.abundance);
      const double floor=std::max(All.MinEgySpec,p->Chemistry.energy_at_temperature(All.MinGasTemp)/unit_energy());
      const double u=std::max(w.internal_energy_cgs/unit_energy(),floor);
      p->Energy+=p->getMass()*(u-p->Utherm);
      p->Utherm=u;
      p->Chemistry.store_conserved(p->getMass());
      if(!p->Chemistry.update_thermodynamics(p->Utherm*unit_energy())) Terminate("CHEMISTRY: invalid output EOS");
      p->Pressure=p->get_pressure(); p->Csnd=p->get_sound_speed();
    }
  }
};
} // namespace

BaseModule *create_arepo_ng_chemistry_module() { return new ChemistryModule; }
void chemistry_ng_initialize_particle(sph_particle_data *p, bool initialize_energy) {
  double hydrogen=p->Chemistry.abundance[sHI]+p->Chemistry.abundance[sHII]+2*p->Chemistry.abundance[sH2I];
  if(hydrogen==0 && p->Chemistry.network_signature==0) {
    p->Chemistry.initialize(p->getMass(),initial_metallicity);
    const double h=CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON;
    p->Chemistry.abundance[sHII]=h*initial_ionized_fraction;
    p->Chemistry.abundance[sELEC]=p->Chemistry.abundance[sHII];
    p->Chemistry.abundance[sH2I]=h*initial_molecular_fraction;
    p->Chemistry.abundance[sHI]=h*(1-initial_ionized_fraction-2*initial_molecular_fraction);
  }
  if(p->Chemistry.network_signature!=chemistry::arepo_ng_network_signature())
    Terminate("CHEMISTRY: IC chemistry network signature is absent or incompatible");
  if(initialize_energy) p->Utherm=p->Chemistry.energy_at_temperature(All.InitGasTemp)/unit_energy();
  p->Utherm=std::max<double>(p->Utherm,std::max(All.MinEgySpec,p->Chemistry.energy_at_temperature(All.MinGasTemp)/unit_energy()));
  if(!p->Chemistry.update_thermodynamics(p->Utherm*unit_energy())) Terminate("CHEMISTRY: invalid IC composition or thermal energy");
  p->Chemistry.store_conserved(p->getMass());
}
