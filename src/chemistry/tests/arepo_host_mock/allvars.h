#ifndef CHEMISTRY_AREPO_TEST_ALLVARS_H
#define CHEMISTRY_AREPO_TEST_ALLVARS_H

#define TGCHEM_NUM_ABUNDANCES 8

struct global_data_all_processes {
  char ChemistryTableDir[1024];
  double Epsilon;
  double UnitEnergy_in_cgs;
  double UnitMass_in_g;
  double T_CMB;
  double cf_atime;
  double Time;
  int ComovingIntegrationOn;
  double J21;
  double rate_Hminus_pd;
  double rate_H2_diss;
  double G0;
};

extern struct global_data_all_processes All;

#endif
