#ifndef CHEMISTRY_AREPO_TEST_ALLVARS_H
#define CHEMISTRY_AREPO_TEST_ALLVARS_H

#include "species.h"
#include <stddef.h>

struct global_data_all_processes {
  char ChemistryTableDir[1024];
  double Epsilon;
  double UnitDensity_in_cgs, HubbleParam, cf_a3inv, LimitUBelowThisDensity;
  double Timebase_interval, cf_time_hubble_a, UnitTime_in_s, Metallicity, nH_StiffEos;
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

typedef long long integertime;
typedef long long MyIDType;
#define TIMER_STORE ((void)0)
#define N_METAL_SPECIES 3
struct sph_particle_data {
  double Metallicity_spec[N_METAL_SPECIES];
  double Abund[TGCHEM_NUM_ABUNDANCES], AbundMass[TGCHEM_NUM_ABUNDANCES];
  double Gamma, Utherm, Energy, Density, t_chem, t_cool, Tdust, Trad;
  double rate_Hminus_pd, rate_H2_diss, rate_HD_diss, rate_CO_diss, J21, G0;
  double k22, piHI, photoionization_rate_HeI, photoheating_rate_HeI;
  double photoionization_rate_HeII, photoheating_rate_HeII;
  double photoionization_rate_C[3], photoheating_rate_C[3];
  double photoionization_rate_O[3], photoheating_rate_O[3];
};
struct particle_data { double Mass; long long ID; int TimeBinHydro; };
extern struct sph_particle_data *SphP;
extern struct particle_data *P;
struct mock_timebins { int NActiveParticles; int *ActiveParticleList; };
extern struct mock_timebins TimeBinsHydro;
enum {CPU_MISC, CPU_CHEM, CPU_CHEMCOMP, CPU_CHEMWAIT, CPU_CHEMCOMM};
extern double CPU_Step[5];
extern int ThisTask, NTask, PTask;
extern int *Send_count, *Recv_count, *Send_offset, *Recv_offset;
#define TIMER_START(x) ((void)(x))
#define TIMER_STOP(x) ((void)(x))
#define TIMER_STOPSTART(x,y) ((void)(x), (void)(y))
#define TAG_CHEM_A 123
#define TAG_CHEM_B 124
double measure_time(void);
double second(void);
double timediff(double, double);
void mpi_printf(const char *, ...);
size_t FreeBytes(void);
void *mymalloc(const char *, size_t);
void myfree(void *);
int data_index_compare(const void *, const void *);
#endif
