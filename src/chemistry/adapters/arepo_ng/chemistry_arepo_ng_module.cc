// NG lifecycle/units only. Numerical chemistry and MPI scheduling live in the adapter.
#include "arepoconfig.h"
#include "modules/basemodule.h"
#include "chemistry_arepo_ng.h"
#include "thermodynamics/thermodynamics.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

#if !defined(RIEMANN_HLLD) || defined(ENABLE_ENTROPY_TREATMENT) || defined(COLLAPSE_EOS) || defined(ISOTHERM_EQS) || defined(COOLING) || defined(STARFORMATION)
#error "CHEMISTRY requires HLLD and the thermal-energy equation, without another cooling/EOS module"
#endif

namespace {
double initial_metallicity=0;
// Bootstrap seeds for a fresh IC/no-chemistry-data particle only (see
// hydrogen==0 below), not a physical initial condition most runs need to tune.
// Optional parameters (see onRegisterParameters): these defaults apply unless
// overridden in the parameter file. INITUNIFORM, when enabled, bypasses this
// path entirely with its own per-species parameters.
double initial_ionized_fraction=1.e-4;
double initial_molecular_fraction=1.e-6;
double unit_energy() { return All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s; }
#ifdef STIFF_EOS
// Cells at or above this hydrogen number density (cm^-3) keep their composition and
// thermal energy frozen for the step: chemistry heating/cooling is skipped so the gas
// responds adiabatically, which stiffens the effective EOS above the threshold. Mirrors
// the legacy "arepo" adapter's STIFF_EOS branch (adapters/arepo/chem_balance.cc) and the
// original AREPO/chemistry_ONECHEM chem_balance.c; same parameter name ("nH_adiabatic")
// for parameter-file compatibility. INCREASE_STIFF_DENSITY / REFINEMENT_AROUND_BH_SC
// (variable threshold near sink-forming regions) are not ported.
double nH_stiff_eos=0;
#endif

#ifdef INITUNIFORM
// Per-species uniform initial abundance, applied unconditionally on a fresh IC start
// (All.RestartFlag==RST_BEGIN), overriding whatever the IC file supplied -- mirrors
// AREPO's INITUNIFORM (init.c: SphP[i].Abund[j] = XH * All.InitAbund[j]; read_ic.c
// there additionally skips reading the ChemicalAbundances dataset in this case, which
// here is unnecessary since it's READ_IF_PRESENT and gets overwritten below anyway).
// Indexed by the species enum from species.h, so it automatically tracks whichever
// species set is actually compiled in (CHEMISTRY_SPECIES_HELIUM/H2_PLUS/DEUTERIUM/
// CARBON_OXYGEN and their sub-options) rather than assuming AREPO's CHEMISTRY_SC list.
double init_abund[TGCHEM_NUM_ABUNDANCES]={};
#endif

class ChemistryModule final: public BaseModule {
  chemistry::ArepoNgAdapter adapter;
  sim *simulation_=nullptr;
  ArepoNgEnvironmentProvider environment_provider_=nullptr;
  char table_root[200]{};
  double cmb_temperature=2.728;
  int distribute=0;
public:
  explicit ChemistryModule(ArepoNgEnvironmentProvider provider) : environment_provider_(provider) { moduleName="Chemistry"; }
  void onRegisterParameters(tParameters *params) override {
    params->add_param("ChemistryTableRoot",table_root,PARAM_STRING,PARAM_FIXED);
    params->add_param("ChemistryInitialMetallicity",&initial_metallicity,PARAM_DOUBLE,PARAM_FIXED);
    // Optional: keep the compiled-in default (see the comment above) unless overridden.
    params->add_param("ChemistryInitialIonizedFraction",&initial_ionized_fraction,PARAM_DOUBLE,PARAM_FIXED,false);
    params->add_param("ChemistryInitialMolecularFraction",&initial_molecular_fraction,PARAM_DOUBLE,PARAM_FIXED,false);
    params->add_param("ChemistryCMBTemperature",&cmb_temperature,PARAM_DOUBLE,PARAM_FIXED);
    params->add_param("ChemistryBalanceWork",&distribute,PARAM_INT,PARAM_FIXED);
#ifdef STIFF_EOS
    params->add_param("nH_adiabatic",&nH_stiff_eos,PARAM_DOUBLE,PARAM_FIXED);
#endif
#ifdef INITUNIFORM
    // One "Init<SpeciesName>" parameter per currently-compiled stored species
    // (same tag names as AREPO's Init* parameters, e.g. InitElec/InitHI/InitHII).
    // species.h's species names already match AREPO's Init<Name> capitalization for
    // every species except the electron, which AREPO spells "Elec" (species.h: "ELEC").
#define CHEM_INITUNIFORM_REGISTER(symbol,name,charge,h,he,d,c,o) \
    { char tag[40]; std::snprintf(tag,sizeof(tag),"Init%s",std::strcmp(name,"ELEC")==0?"Elec":name); params->add_param(tag,&init_abund[s##symbol],PARAM_DOUBLE,PARAM_FIXED); }
    CHEM_SPECIES_STORED_LIST(CHEM_INITUNIFORM_REGISTER)
#undef CHEM_INITUNIFORM_REGISTER
#endif
  }
  void onBegrunBeforeInit(sim *simulation) override {
    simulation_=simulation;
    // Module parameter storage is outside All's automatically broadcast block.
    MPI_Bcast(table_root,sizeof(table_root),MPI_CHAR,0,simulation->Communicator);
    MPI_Bcast(&initial_metallicity,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&initial_ionized_fraction,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&initial_molecular_fraction,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&cmb_temperature,1,MPI_DOUBLE,0,simulation->Communicator);
    MPI_Bcast(&distribute,1,MPI_INT,0,simulation->Communicator);
#ifdef STIFF_EOS
    MPI_Bcast(&nH_stiff_eos,1,MPI_DOUBLE,0,simulation->Communicator);
    if(!std::isfinite(nH_stiff_eos) || nH_stiff_eos<=0)
      Terminate("CHEMISTRY: nH_adiabatic (STIFF_EOS threshold) must be finite and positive");
#endif
#ifdef INITUNIFORM
    MPI_Bcast(init_abund,TGCHEM_NUM_ABUNDANCES,MPI_DOUBLE,0,simulation->Communicator);
    for(double a: init_abund)
      if(!std::isfinite(a) || a<0) Terminate("CHEMISTRY: Init<species> (INITUNIFORM) values must be finite and non-negative");
#endif
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
  void onStartOfTimestep(simparticles *, mymesh *mesh) override {
    simulation_->NgbTree->update_data_of_ghost_cells<ngbtree::primitives>();
    mesh->calculate_gradients();
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
#ifdef STIFF_EOS
      // Comoving integration is rejected above, so density needs no a^-3/h^2 factors here.
      if(chemistry_hydrogen_number_density(p->Density*All.UnitDensity_in_cgs) > nH_stiff_eos) continue;
#endif
      chemistry::ArepoNgWork w{};
      w.environment.density_cgs=p->Density*All.UnitDensity_in_cgs;
      w.environment.metallicity=p->Chemistry.metallicity;
      if(environment_provider_) environment_provider_(p(),w.environment);
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
      // The solver already converges gamma self-consistently with the returned
      // abundances/energy (chemistry_solver.cc's energy residual uses the same O(1)
      // chemistry_temperature()/chemistry_gamma_from_abundances() pair to do so).
      // AREPO's own adapter (adapters/arepo/chem_balance.cc: compute_temp(),
      // particle_chemistry_state_commit()) likewise just stores/reuses the solver's
      // gamma directly, never re-solving it. Re-deriving it here via
      // update_thermodynamics()'s 80-iteration bracketing was redundant; only the MinGasTemp
      // floor clamp above can make the stored u differ from the solver's own u, in
      // which case gamma is very slightly stale until the next chemistry step --
      // exactly the same staleness AREPO already accepted.
      if(!(w.gamma>1) || !(w.gamma<2) || !std::isfinite(w.gamma)) Terminate("CHEMISTRY: invalid output gamma");
      p->Chemistry.gamma=w.gamma;
      p->Pressure=p->get_pressure(); p->Csnd=p->get_sound_speed();
    }
  }
};
} // namespace

BaseModule *create_arepo_ng_chemistry_module(ArepoNgEnvironmentProvider provider) { return new ChemistryModule(provider); }
void chemistry_ng_initialize_particle(sph_particle_data *p, bool initialize_energy) {
#ifdef INITUNIFORM
  if(All.RestartFlag==RST_BEGIN) {
    // Unconditional override from Init<Species> parameters, matching AREPO's
    // INITUNIFORM: SphP[i].Abund[j] = XH * All.InitAbund[j] for every stored species,
    // regardless of whatever the IC file supplied.
    p->Chemistry={};
    p->Chemistry.metallicity=initial_metallicity;
    p->Chemistry.network_signature=chemistry::arepo_ng_network_signature();
    const double h=CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON;
#define CHEM_INITUNIFORM_APPLY(symbol,name,charge,h_,he,d,c,o) p->Chemistry.abundance[s##symbol]=h*init_abund[s##symbol];
    CHEM_SPECIES_STORED_LIST(CHEM_INITUNIFORM_APPLY)
#undef CHEM_INITUNIFORM_APPLY
  } else
#endif
  {
    double hydrogen=p->Chemistry.abundance[sHI]+p->Chemistry.abundance[sHII]+2*p->Chemistry.abundance[sH2I];
    if(hydrogen==0 && p->Chemistry.network_signature==0) {
      p->Chemistry.initialize(p->getMass(),initial_metallicity);
      const double h=CHEMISTRY_HYDROGEN_NUCLEI_PER_BARYON;
      p->Chemistry.abundance[sHII]=h*initial_ionized_fraction;
      p->Chemistry.abundance[sELEC]=p->Chemistry.abundance[sHII];
      p->Chemistry.abundance[sH2I]=h*initial_molecular_fraction;
      p->Chemistry.abundance[sHI]=h*(1-initial_ionized_fraction-2*initial_molecular_fraction);
    }
  }
  if(p->Chemistry.network_signature!=chemistry::arepo_ng_network_signature())
    Terminate("CHEMISTRY: IC chemistry network signature is absent or incompatible");
  // InitGasTemp is a known target temperature: set u and gamma together (O(1)) rather
  // than setting u from it and immediately re-solving that same u back to a
  // temperature via update_thermodynamics()'s bisection.
  if(initialize_energy) p->Utherm=p->Chemistry.initialize_at_temperature(All.InitGasTemp)/unit_energy();
  const double min_u=std::max(All.MinEgySpec,p->Chemistry.energy_at_temperature(All.MinGasTemp)/unit_energy());
  if(p->Utherm<min_u) {
    // Floor engaged: MinGasTemp is also a known temperature, so set directly.
    p->Utherm=p->Chemistry.initialize_at_temperature(All.MinGasTemp)/unit_energy();
  } else if(!initialize_energy) {
    // Utherm came from the IC/restart file at an a priori unknown temperature.
    if(!p->Chemistry.update_thermodynamics(p->Utherm*unit_energy()))
      Terminate("CHEMISTRY: invalid IC composition or thermal energy");
  }
  // else: initialize_energy path, not floored -- gamma was already set above.
  if(!(p->Chemistry.gamma>1) || !(p->Chemistry.gamma<2) || !std::isfinite(p->Chemistry.gamma))
    Terminate("CHEMISTRY: invalid IC composition or thermal energy");
  p->Chemistry.store_conserved(p->getMass());
}
