#include "adapters/arepo_ng/chemistry_arepo_ng.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  int rank,size; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  chemistry::ArepoNgAdapter adapter;
  ChemistryConfig config{}; chemistry_config_defaults(&config); config.table_root=argv[1]; config.print_configuration=0;
  assert(adapter.initialize(config,MPI_COMM_WORLD)==CHEMISTRY_STATUS_OK);
  chemistry::ArepoNgCell cell{}; cell.initialize(2,0);
  cell.abundance[sHI]=0.56-1e-4; cell.abundance[sH2I]=0.1;
  cell.abundance[sHII]=1e-4; cell.abundance[sELEC]=1e-4; cell.abundance[sHM]=1e-20;
  assert(cell.update_thermodynamics(1e11));
  assert(cell.gamma>1.3 && cell.gamma<5.0/3.0);
  cell.store_conserved(2);
  auto recovered=cell;
  assert(recovered.recover_conserved(2,1e11));
  assert(std::abs(recovered.gamma-cell.gamma)<1e-14);
  ChemistryGlobalEnvironment global{}; global.scale_factor=1; global.cmb_temperature=2.728;
  std::vector<chemistry::ArepoNgWork> local(rank==0 ? 5 : (rank==1 ? 0 : 2));
  for(size_t i=0;i<local.size();++i) {
    auto &w=local[i]; w={};
    w.environment.density_cgs=1e-20; w.environment.metallicity=0;
    std::memcpy(w.abundance,cell.abundance,sizeof(w.abundance));
    w.internal_energy_cgs=1e11; w.gamma=cell.gamma; w.duration_s=1e3;
    w.id=(1ULL<<40)+rank*100+i;
  }
  auto balanced=local;
  assert(adapter.advance(global,local,false)==CHEMISTRY_STATUS_OK);
  assert(adapter.advance(global,balanced,true)==CHEMISTRY_STATUS_OK);
  assert(local.size()==balanced.size());
  for(size_t i=0;i<local.size();++i) {
    assert(local[i].id==balanced[i].id);
    assert(local[i].internal_energy_cgs==balanced[i].internal_energy_cgs);
    assert(local[i].gamma==balanced[i].gamma);
    for(int j=0;j<TGCHEM_NUM_ABUNDANCES;++j) assert(local[i].abundance[j]==balanced[i].abundance[j]);
  }
  // Cross the bounded communication batch boundary with empty remote owners.
  std::vector<chemistry::ArepoNgWork> many;
  if(rank==0) {
    many.assign(4101,local.front());
    for(size_t i=0;i<many.size();++i) { many[i].id=(1ULL<<42)+i; many[i].duration_s=0; }
  }
  assert(adapter.advance(global,many,true)==CHEMISTRY_STATUS_OK);
  for(size_t i=0;i<many.size();++i) assert(many[i].id==(1ULL<<42)+i);
  auto failed=balanced;
  if(rank==0) failed[0].environment.density_cgs=-1;
  const auto before=failed;
  assert(adapter.advance(global,failed,true)!=CHEMISTRY_STATUS_OK);
  assert(failed.size()==before.size());
  if(!failed.empty()) assert(std::memcmp(failed.data(),before.data(),failed.size()*sizeof(failed[0]))==0);
  std::vector<chemistry::ArepoNgWork> empty;
  assert(adapter.advance(global,empty,true)==CHEMISTRY_STATUS_OK);
  if(rank==0) std::printf("NG chemistry MPI: local/balanced, empty owner, 64-bit ID and atomic failure passed (%d ranks)\n",size);
  MPI_Finalize();
}
