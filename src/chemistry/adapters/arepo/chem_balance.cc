#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

extern "C" {
#include "allvars.h"
}
#include "chemistry_constants.h"
#include "chemistry_arepo.h"
#include "chemistry_model.hpp"
#include "dust_cooling.h"
#include "thermodynamics.h"

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT

#ifdef CHEMISTRY_HYBRID_IMPLICIT
int Nimplicit  = 0;
int Nimplicit0 = 0;
#endif
#ifdef CHEMISTRY_ENERGY_IMPLICIT
int Ndiv  = 0;
int Ndiv0 = 0;
#endif
#ifdef CHEM_BALANCE
struct chemcomm_in
{
  MyIDType ID;

  double Utherm_cgs;                /*!< current value of entropy (actually entropic function) of particle */
  double Density;		/*!< current baryonic mass density of particle */

  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;

#ifdef THERMAL_FEEDBACK
  double Trad;
#endif

#ifdef CHEMISTRY_LW_RADIATION
#ifdef LW_RADIATION_LOCAL
  double rate_Hminus_pd;
  double rate_H2_diss;
#ifdef RYTR_NCO
  double rate_CO_diss;
#endif
  double J21;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  double rate_HD_diss;
#endif
#if CHEMISTRY_DUST_COOLING
  double G0;
#endif
#endif
#endif

#ifdef IONIZING_PHOTON
  double k22, piHI;
#ifdef RYTR_HELIUM_IONIZATION
  double photoionization_rate_HeI, photoheating_rate_HeI;
  double photoionization_rate_HeII, photoheating_rate_HeII;
#endif
#ifdef RYTR_CARBON_OXYGEN_IONIZATION
  double photoionization_rate_C[3], photoheating_rate_C[3];
  double photoionization_rate_O[3], photoheating_rate_O[3];
#endif
#endif

#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  double Metallicity;
#endif
  double dtime;

}
  *ChemDataIn,
  *ChemDataGet;

struct chemcomm_out
{
  double Utherm_cgs;                /*!< current value of entropy (actually entropic function) of particle */

  double Abund[TGCHEM_NUM_ABUNDANCES];

  double Gamma;
  double t_chem, t_cool;

#ifdef METALLICITY
  double Tdust;
#endif

#ifdef CHEMISTRY_ENERGY_IMPLICIT
  char flag_imp;
  char flag_conv;
#endif
}
  *ChemDataResult,
  *ChemDataOut;

void chem_balance(const chemistry::Model *context,
                  const ChemistryStepContext *global_step);
void chem_evaluate(int target,
                   const chemistry::Model *context,
                   const ChemistryStepContext *global_step);
static inline void particle2in_chem(struct chemcomm_in *in, int i);
static inline void out2particle_chem(struct chemcomm_out *out, int i);
#endif

// Host storage/transport policy ends here. Core always receives solar units.
static double particle_metallicity_solar(int i)
{
#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  double metallicity = 0.;
  for(int k = 0; k < N_METAL_SPECIES; ++k)
    metallicity += SphP[i].Metallicity_spec[k];
  return metallicity / Z_SUN;
#elif defined(METALLICITY)
  (void)i;
  return All.Metallicity;
#else
  (void)i;
  return 0.;
#endif
}

static inline void particle2in_arepo(int i, ChemistryCellEnvironment *cell,
                                    ChemistryState *state);
static inline void out2particle_arepo(int i, const ChemistryState *state,
                                     const ChemistryResult *result);

static inline void particle_chemistry_state_read(int i, double *abund, double *gamma)
{
  /* The host stores chemistry's public abundance representation unchanged. */
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    abund[j] = SphP[i].Abund[j];

  *gamma = SphP[i].Gamma;
}

static inline void particle_chemistry_state_commit(int i, double utherm_cgs,
                                                   const double *abund, double gamma)
{
  double du = utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) - SphP[i].Utherm;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
  SphP[i].Utherm = utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);

  /* No normalization conversion occurs at the host boundary. */
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    {
      SphP[i].Abund[j] = abund[j];
      SphP[i].AbundMass[j] = abund[j] * P[i].Mass;
    }

  SphP[i].Gamma = gamma;
}

static inline void particle2in_arepo(int i, ChemistryCellEnvironment *cell,
                                    ChemistryState *state)
{
  memset(cell, 0, sizeof(*cell));
  cell->diagnostic_id = (int) P[i].ID;
  cell->density_cgs = All.UnitDensity_in_cgs * SphP[i].Density *
                  All.HubbleParam * All.HubbleParam * All.cf_a3inv;

  cell->metallicity = particle_metallicity_solar(i);

#if defined(CHEMISTRY_LW_RADIATION) && defined(LW_RADIATION_LOCAL)
  cell->local_rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  cell->local_rate_H2_diss = SphP[i].rate_H2_diss;
#ifdef RYTR_NCO
  cell->local_rate_CO_diss = SphP[i].rate_CO_diss;
#endif
  cell->local_J21 = SphP[i].J21;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  cell->local_rate_HD_diss = SphP[i].rate_HD_diss;
#endif
#if CHEMISTRY_DUST_COOLING
  cell->local_G0 = SphP[i].G0;
#endif
#endif

#ifdef THERMAL_FEEDBACK
  cell->stellar_radiation_flux = 4. * STEFAN * pow(SphP[i].Trad, 4.);
#endif

#ifdef IONIZING_PHOTON
  cell->photoionization_rate_HI = SphP[i].k22;
  cell->photoheating_rate_HI = SphP[i].piHI;
#ifdef RYTR_HELIUM_IONIZATION
  cell->photoionization_rate_HeI = SphP[i].photoionization_rate_HeI;
  cell->photoheating_rate_HeI = SphP[i].photoheating_rate_HeI;
  cell->photoionization_rate_HeII = SphP[i].photoionization_rate_HeII;
  cell->photoheating_rate_HeII = SphP[i].photoheating_rate_HeII;
#endif
#ifdef RYTR_CARBON_OXYGEN_IONIZATION
  for(int stage = 0; stage < 3; stage++) {
    cell->photoionization_rate_C[stage] =
        SphP[i].photoionization_rate_C[stage];
    cell->photoheating_rate_C[stage] = SphP[i].photoheating_rate_C[stage];
    cell->photoionization_rate_O[stage] =
        SphP[i].photoionization_rate_O[stage];
    cell->photoheating_rate_O[stage] = SphP[i].photoheating_rate_O[stage];
  }
#endif
#endif

  state->specific_thermal_energy_cgs = SphP[i].Utherm * All.UnitEnergy_in_cgs/All.UnitMass_in_g;
  particle_chemistry_state_read(i, state->abundances, &state->gamma);
}

static inline void out2particle_arepo(int i, const ChemistryState *state,
                                     const ChemistryResult *result)
{
  particle_chemistry_state_commit(i, state->specific_thermal_energy_cgs, state->abundances,
                                  state->gamma);

#ifndef CHEMISTRY_IMPLICIT
  SphP[i].t_chem = result->chemistry_time_s;
#endif
#ifndef CHEMISTRY_ENERGY_IMPLICIT
  SphP[i].t_cool = result->cooling_time_s;
#endif

#ifdef METALLICITY
  SphP[i].Tdust = result->dust_temperature;
#endif
}

extern "C" double compute_temp(int ithis)
{
  /* Preserve the operation order of the historical AREPO implementation. */
  const double utherm_cgs =
      (All.UnitEnergy_in_cgs / All.UnitMass_in_g) * SphP[ithis].Utherm;
  double abund[TGCHEM_NUM_ABUNDANCES];
  double gamma;
  particle_chemistry_state_read(ithis, abund, &gamma);
  return chemistry_temperature(abund, gamma, utherm_cgs);
}

extern "C" double compute_xmu(struct sph_particle_data sph_p)
{
  double abund[TGCHEM_NUM_ABUNDANCES];
  for(int species = 0; species < TGCHEM_NUM_ABUNDANCES; species++)
    abund[species] = sph_p.Abund[species];
  return chemistry_molecular_weight(abund);
}

extern "C" void cooling_only(void){
  CPU_Step[CPU_MISC] += measure_time();

#ifdef CHEMISTRY_HYBRID_IMPLICIT
  Nimplicit  = 0;
#endif
#ifdef CHEMISTRY_ENERGY_IMPLICIT
  Ndiv  = 0;
#endif
  int diagnostic_retry_used = 0;

  const chemistry::Model *chemistry_context = chemistry_arepo_host_context();
  if(chemistry_context == NULL)
    {
      mpi_printf("CHEMISTRY: AREPO chemistry context is not initialized.\n");
      return;
    }

  ChemistryStepContext global_step;
  if(!chemistry_arepo_host_prepare_step(&global_step))
    {
      mpi_printf("CHEMISTRY: failed to prepare the shared chemistry step.\n");
      return;
    }

#ifdef AVOID_ADIABATIC_COOLING
  if(All.Time < All.Time_avoid_adiabatic_cooling)
    return;
#endif

  mpi_printf("CHEM: init...\n");
  double t0 = second();

#ifndef CHEM_BALANCE
  /* Now assign new timesteps  */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;  // skip cells that have been swallowed or eliminated
      if(SphP[i].Density < All.LimitUBelowThisDensity)
        continue;

#ifdef STOP_REFINMENT_DISTANCE
      if(SphP[i].AllowRefinement == 0)
	continue;
#endif

#if defined(CHEMISTRY_IGNORE_LOWRES) && defined(REFINEMENT_HIGH_RES_GAS)
      if(SphP[i].AllowRefinement == 0)
	continue;

      if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
        continue;
#endif

#ifdef STIFF_EOS
      const double rho_cgs = All.UnitDensity_in_cgs * SphP[i].Density *
          All.HubbleParam * All.HubbleParam * All.cf_a3inv;
      double nH = chemistry_hydrogen_number_density(rho_cgs);

      double nH_StiffEos = All.nH_StiffEos;
#ifdef REFINEMENT_AROUND_BH_SC
      nH_StiffEos = All.nH_StiffEosInit * pow(SphP[i].Dist_BH / All.DistReference, -2.);
#ifdef REDUCE_RESOLUTION_ABOVE_RADIUS
      if(nH_StiffEos < All.nH_StiffEos)
#else
	if(nH_StiffEos > All.nH_StiffEos)
#endif
	  nH_StiffEos = All.nH_StiffEos;
#endif

      if(nH > nH_StiffEos)
	continue;
#endif

      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      if(dt == 0)
        continue;
      double dtime = All.cf_atime * dt / All.cf_time_hubble_a;
      const double duration_s = dtime * All.UnitTime_in_s / All.HubbleParam;

      ChemistryCellEnvironment cell;
      double abundance_storage[TGCHEM_NUM_ABUNDANCES];
      ChemistryState state = {abundance_storage, TGCHEM_NUM_ABUNDANCES, 0., 0.};
      ChemistryResult result = {};
      const ChemistryAdvanceOptions options = {CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};
      ChemWorkspace workspace{};
      const auto workspace_status = CHEMISTRY_STATUS_OK;
      particle2in_arepo(i, &cell, &state);

      ChemistryStatus ifunc = workspace_status;
      if(ifunc == CHEMISTRY_STATUS_OK)
        ifunc = chemistry::Driver{}.advance(*chemistry_context, global_step, workspace,
                                  cell, duration_s, options, state, result);
      if(ifunc != CHEMISTRY_STATUS_OK && !diagnostic_retry_used){
        cell.diagnostics_enabled = 1;
        ifunc = chemistry::Driver{}.advance(*chemistry_context, global_step, workspace,
                                  cell, duration_s, options, state, result);
        diagnostic_retry_used = 1;
      }

      out2particle_arepo(i, &state, &result);

#ifdef CHEMISTRY_HYBRID_IMPLICIT
      if(result.implicit_flag)
        Nimplicit++;
#endif

#ifdef CHEMISTRY_ENERGY_IMPLICIT
      if(result.implicit_flag && ifunc){
        Ndiv++;
      }
#endif
    }
#else
  chem_balance(chemistry_context, &global_step);
#endif

#ifdef CHEMISTRY_HYBRID_IMPLICIT
  int ierr;
  ierr = MPI_Reduce(&Nimplicit, &Nimplicit0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef CHEMISTRY_ENERGY_IMPLICIT
  ierr = MPI_Reduce(&Ndiv, &Ndiv0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  if(ThisTask == 0){
    printf("chemistries of %d particles are calculated by implicit method.\n",Nimplicit0);
#ifdef CHEMISTRY_ENERGY_IMPLICIT
    if(Ndiv0)
      printf("energy-loop of %d particles are not converged\n",Ndiv0);
#endif
  }
#endif

  double t1 = second();
  mpi_printf("CHEMISTRY: took %g sec.\n", timediff(t0, t1));

  CPU_Step[CPU_CHEM] += measure_time();
}



#ifdef CHEM_BALANCE
struct data_index
{
  int Task;
  int Index;
  int IndexGet;
} *DataIndexTable;

size_t sizemax(size_t a, size_t b)
{
  if(a < b)   return b;
  else    return a;
}

void chem_balance(const chemistry::Model *context,
                  const ChemistryStepContext *global_step)
{
  int j, ndone;
  MPI_Status status;

  mpi_printf("\n  [chem_balance.c] start --> \n");

  int BunchSize = (int) ((0.3 * FreeBytes()) / (sizeof(struct data_index) +
						sizeof(struct chemcomm_in) +
						sizeof(struct chemcomm_out) +
						sizemax(sizeof(struct chemcomm_in),
							sizeof(struct chemcomm_out))));
  DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", BunchSize * sizeof(struct data_index));

  TIMER_STORE;

  int idx = 0;
  int n_exported = 0;
  do{
    TIMER_START(CPU_CHEMCOMP);

    int nexport = 0;
    for(int j = 0; j < NTask; j++)
      Send_count[j] = 0;

    /* do local particles and prepare export list */
    for(;idx < TimeBinsHydro.NActiveParticles; idx++){
      int i = TimeBinsHydro.ActiveParticleList[idx];

      if(nexport >= BunchSize)
	break;

      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;  // skip cells that have been swallowed or eliminated
#if defined(CHEMISTRY_IGNORE_LOWRES) && defined(REFINEMENT_HIGH_RES_GAS)
      if(SphP[i].AllowRefinement == 0)
	continue;

      if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
        continue;
#endif

#ifdef STIFF_EOS
      const double rho_cgs = All.UnitDensity_in_cgs * SphP[i].Density *
          All.HubbleParam * All.HubbleParam * All.cf_a3inv;
      double nH = chemistry_hydrogen_number_density(rho_cgs);
      if(nH > All.nH_StiffEos)
	continue;
#endif

      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      if(dt == 0)
	continue;

      int stask = nexport % (NTask - 1);
      if(stask >= ThisTask) stask++;

      DataIndexTable[nexport].Task  = stask;
      DataIndexTable[nexport].Index = i;
      DataIndexTable[nexport].IndexGet = nexport;

      nexport++;
    }

    n_exported += nexport;
    for(int j = 0; j < NTask; j++)
      Send_count[j] = 0;
    for(int j = 0; j < nexport; j++)
      Send_count[DataIndexTable[j].Task]++;

    qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

    TIMER_STOPSTART(CPU_CHEMCOMP,CPU_CHEMWAIT);

    MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

    TIMER_STOPSTART(CPU_CHEMWAIT,CPU_CHEMCOMM);

    int nimport = 0;
    for(j = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
      {
	nimport += Recv_count[j];
	if(j > 0)
	  {
	    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	    Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	  }
      }


    ChemDataGet = (struct chemcomm_in *) mymalloc("ChemDataGet", nimport * sizeof(struct chemcomm_in));
    ChemDataIn = (struct chemcomm_in *) mymalloc("ChemDataIn", nexport * sizeof(struct chemcomm_in));

    /* prepare particle data for export */
    for(int j = 0; j < nexport; j++)
      {
	int place = DataIndexTable[j].Index;
	particle2in_chem(&ChemDataIn[j], place);		// MADE D_IND CHANGE IN HERE
      }

    for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
      {
	int recvTask = ThisTask ^ ngrp;

	if(recvTask < NTask)
	  {
	    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	      {
		/* get the particles */
		MPI_Sendrecv(&ChemDataIn[Send_offset[recvTask]],
			     Send_count[recvTask] * sizeof(struct chemcomm_in), MPI_BYTE,
			     recvTask, TAG_CHEM_A,
			     &ChemDataGet[Recv_offset[recvTask]],
			     Recv_count[recvTask] * sizeof(struct chemcomm_in), MPI_BYTE,
			     recvTask, TAG_CHEM_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	  }
      }

    TIMER_STOPSTART(CPU_CHEMCOMM,CPU_CHEMCOMP);

    myfree(ChemDataIn);
    ChemDataResult = (struct chemcomm_out *) mymalloc("ChemDataResult", nimport * sizeof(struct chemcomm_out));
    ChemDataOut = (struct chemcomm_out *) mymalloc("ChemDataOut", nexport * sizeof(struct chemcomm_out));
    //report_memory_usage(&HighMark_sphhydro, "SPH_CHEM");

    /* Each evaluation owns its diagnostics flag and workspace. */
    for(int j=0;j<nimport;j++)
      chem_evaluate(j, context, global_step);

    int ndone_flag = 0;
    if(idx == TimeBinsHydro.NActiveParticles)
      ndone_flag = 1;

    TIMER_STOPSTART(CPU_CHEMCOMP, CPU_CHEMWAIT);

    MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    TIMER_STOPSTART(CPU_CHEMWAIT, CPU_CHEMCOMM);

    /* get the result */
    for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
      {
	int recvTask = ThisTask ^ ngrp;
	if(recvTask < NTask)
	  {
	    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	      {
		/* send the results */
		MPI_Sendrecv(&ChemDataResult[Recv_offset[recvTask]],
			     Recv_count[recvTask] * sizeof(struct chemcomm_out),
			     MPI_BYTE, recvTask, TAG_CHEM_B,
			     &ChemDataOut[Send_offset[recvTask]],
			     Send_count[recvTask] * sizeof(struct chemcomm_out),
			     MPI_BYTE, recvTask, TAG_CHEM_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	  }
      }
    TIMER_STOP(CPU_CHEMCOMM);

    /* add the result to the local particles */
    for(int j = 0; j < nexport; j++)
      {
	int place = DataIndexTable[j].Index;
	out2particle_chem(&ChemDataOut[j], place);
      }

    myfree(ChemDataOut);
    myfree(ChemDataResult);
    myfree(ChemDataGet);
  }while(ndone < NTask);

  myfree(DataIndexTable);
}




void chem_evaluate(int target,
                   const chemistry::Model *context,
                   const ChemistryStepContext *global_step)
{
  ChemistryCellEnvironment cell;
  memset(&cell, 0, sizeof(cell));
  cell.diagnostic_id = (int) ChemDataGet[target].ID;
  cell.density_cgs = All.UnitDensity_in_cgs * ChemDataGet[target].Density *
                 All.HubbleParam * All.HubbleParam * All.cf_a3inv;

#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  cell.metallicity = ChemDataGet[target].Metallicity;
#elif defined(METALLICITY)
  // Uniform metallicity is reconstructed locally, not sent per particle.
  cell.metallicity = All.Metallicity;
#else
  cell.metallicity = 0.;
#endif

#ifdef CHEMISTRY_LW_RADIATION
#ifdef LW_RADIATION_LOCAL
  cell.local_rate_Hminus_pd = ChemDataGet[target].rate_Hminus_pd;
  cell.local_rate_H2_diss = ChemDataGet[target].rate_H2_diss;
#ifdef RYTR_NCO
  cell.local_rate_CO_diss = ChemDataGet[target].rate_CO_diss;
#endif
  cell.local_J21 = ChemDataGet[target].J21;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  cell.local_rate_HD_diss = ChemDataGet[target].rate_HD_diss;
#endif
#if CHEMISTRY_DUST_COOLING
  cell.local_G0 = ChemDataGet[target].G0;
#endif
#endif
#endif

#ifdef THERMAL_FEEDBACK
  cell.stellar_radiation_flux = 4. * STEFAN * pow(ChemDataGet[target].Trad, 4.);
#endif

#ifdef IONIZING_PHOTON
  cell.photoionization_rate_HI = ChemDataGet[target].k22;
  cell.photoheating_rate_HI = ChemDataGet[target].piHI;
#ifdef RYTR_HELIUM_IONIZATION
  cell.photoionization_rate_HeI = ChemDataGet[target].photoionization_rate_HeI;
  cell.photoheating_rate_HeI = ChemDataGet[target].photoheating_rate_HeI;
  cell.photoionization_rate_HeII = ChemDataGet[target].photoionization_rate_HeII;
  cell.photoheating_rate_HeII = ChemDataGet[target].photoheating_rate_HeII;
#endif
#ifdef RYTR_CARBON_OXYGEN_IONIZATION
  for(int stage = 0; stage < 3; stage++) {
    cell.photoionization_rate_C[stage] =
        ChemDataGet[target].photoionization_rate_C[stage];
    cell.photoheating_rate_C[stage] =
        ChemDataGet[target].photoheating_rate_C[stage];
    cell.photoionization_rate_O[stage] =
        ChemDataGet[target].photoionization_rate_O[stage];
    cell.photoheating_rate_O[stage] =
        ChemDataGet[target].photoheating_rate_O[stage];
  }
#endif
#endif

  double abundance_storage[TGCHEM_NUM_ABUNDANCES];
  ChemistryState state = {abundance_storage, TGCHEM_NUM_ABUNDANCES, 0., 0.};
  state.specific_thermal_energy_cgs = ChemDataGet[target].Utherm_cgs;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    state.abundances[j] = ChemDataGet[target].Abund[j];
  state.gamma = ChemDataGet[target].Gamma;

  ChemistryResult result = {};
  const ChemistryAdvanceOptions options = {CHEMISTRY_MODE_DEFAULT, CHEMISTRY_SOLVER_AUTO};
  ChemWorkspace workspace{};
  const auto workspace_status = CHEMISTRY_STATUS_OK;
  const double dtime = ChemDataGet[target].dtime;
  ChemistryStatus ifunc = workspace_status;
  if(ifunc == CHEMISTRY_STATUS_OK)
    ifunc = chemistry::Driver{}.advance(*context, *global_step, workspace,
                              cell, dtime, options, state, result);
  if(ifunc != CHEMISTRY_STATUS_OK){
    cell.diagnostics_enabled = 1;
    ifunc = chemistry::Driver{}.advance(*context, *global_step, workspace, cell, dtime,
                              options, state, result);
  }

#ifdef CHEMISTRY_ENERGY_IMPLICIT
#ifdef CHEMISTRY_HYBRID_IMPLICIT
  if(result.implicit_flag){
    Nimplicit++;
  }
#endif
  if(result.implicit_flag && ifunc){
    Ndiv++;
  }
#endif

  ChemDataResult[target].Utherm_cgs = state.specific_thermal_energy_cgs;

  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    ChemDataResult[target].Abund[j] = state.abundances[j];

  ChemDataResult[target].Gamma  = state.gamma;
  ChemDataResult[target].t_chem = result.chemistry_time_s;
  ChemDataResult[target].t_cool = result.cooling_time_s;

#ifdef METALLICITY
  ChemDataResult[target].Tdust = result.dust_temperature;
#endif

#ifdef CHEMISTRY_ENERGY_IMPLICIT
  ChemDataResult[target].flag_imp = result.implicit_flag;
  ChemDataResult[target].flag_conv = result.convergence_flag;
#endif
}


static inline void particle2in_chem(struct chemcomm_in *in, int i)
{
  in->ID    = P[i].ID;

  in->Utherm_cgs = All.UnitEnergy_in_cgs/All.UnitMass_in_g * SphP[i].Utherm;
  in->Density = SphP[i].Density;

  particle_chemistry_state_read(i, in->Abund, &in->Gamma);

#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  in->Metallicity = particle_metallicity_solar(i);
#endif

#ifdef CHEMISTRY_LW_RADIATION
#ifdef LW_RADIATION_LOCAL
  in->rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  in->rate_H2_diss = SphP[i].rate_H2_diss;
#ifdef RYTR_NCO
  in->rate_CO_diss = SphP[i].rate_CO_diss;
#endif
  in->J21 = SphP[i].J21;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  in->rate_HD_diss = SphP[i].rate_HD_diss;
#endif
#if CHEMISTRY_DUST_COOLING
  in->G0 = SphP[i].G0;
#endif
#endif
#endif

#ifdef THERMAL_FEEDBACK
  in->Trad = SphP[i].Trad;
#endif

#ifdef IONIZING_PHOTON
  in->k22  = SphP[i].k22;
  in->piHI = SphP[i].piHI;
#ifdef RYTR_HELIUM_IONIZATION
  in->photoionization_rate_HeI = SphP[i].photoionization_rate_HeI;
  in->photoheating_rate_HeI = SphP[i].photoheating_rate_HeI;
  in->photoionization_rate_HeII = SphP[i].photoionization_rate_HeII;
  in->photoheating_rate_HeII = SphP[i].photoheating_rate_HeII;
#endif
#ifdef RYTR_CARBON_OXYGEN_IONIZATION
  for(int stage = 0; stage < 3; stage++) {
    in->photoionization_rate_C[stage] =
        SphP[i].photoionization_rate_C[stage];
    in->photoheating_rate_C[stage] = SphP[i].photoheating_rate_C[stage];
    in->photoionization_rate_O[stage] =
        SphP[i].photoionization_rate_O[stage];
    in->photoheating_rate_O[stage] = SphP[i].photoheating_rate_O[stage];
  }
#endif
#endif

  double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
  double dtime = All.cf_atime * dt / All.cf_time_hubble_a;
  const double duration_s = dtime * All.UnitTime_in_s / All.HubbleParam;
  in->dtime = duration_s;

}


/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_chem(struct chemcomm_out *out, int i)
{
  particle_chemistry_state_commit(i, out->Utherm_cgs, out->Abund, out->Gamma);
#ifndef CHEMISTRY_IMPLICIT
  SphP[i].t_chem = out->t_chem;
#endif
#ifndef CHEMISTRY_ENERGY_IMPLICIT
  SphP[i].t_cool = out->t_cool;
#endif

#ifdef METALLICITY
  SphP[i].Tdust  = out->Tdust;
#endif
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU. */
extern "C" int data_index_compare(const void *a, const void *b)
{
  if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
    return -1;

  if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
    return +1;

  if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
    return -1;

  if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
    return +1;

  if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
    return -1;

  if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
    return +1;

  return 0;
}
#endif

/* AREPO run lifecycle and shared radiation input. */
#include <string.h>

#include "allvars.h"

#include "chemistry_arepo.h"
#include "radiation/radiation_rates.h"

/*
 * Existing AREPO calls cooling_only() without an argument, so the thin host
 * layer must retain one run-lifetime handle. Tables and numerical state remain
 * inside the explicitly passed context; the chemistry core has no accessor for
 * this object.
 */
static chemistry::Model::Handle arepo_chemistry_context;

#ifdef CHEMISTRY_LW_RADIATION
static RadiationRateContext arepo_radiation_context;
#ifdef CHEMISTRY_BACKGROUND_LW_RADIATION
static RadiationBackgroundState arepo_background_state;

#endif
#endif

int chemistry_arepo_host_init(void)
{
  ChemistryConfig config;
  chemistry_config_defaults(&config);
  config.table_root = All.ChemistryTableDir;
  config.timescale_safety_factor = All.Epsilon;
  config.print_configuration = (ThisTask == 0);
  auto created = chemistry::Model::create(config);
  if(created.status != CHEMISTRY_STATUS_OK) return 0;
  arepo_chemistry_context = std::move(created.model);

#ifdef CHEMISTRY_LW_RADIATION
  if(!radiation_rate_context_init(&arepo_radiation_context)) {
    arepo_chemistry_context.reset();

    return 0;
  }
#ifdef CHEMISTRY_BACKGROUND_LW_RADIATION
  RadiationBlackbodyInput background = {
      All.J21,
      All.Teff};
  if(!radiation_background_state_update(
         &arepo_radiation_context, &background,
         &arepo_background_state)) {
    radiation_rate_context_destroy(&arepo_radiation_context);
    radiation_background_state_reset(&arepo_background_state);
    arepo_chemistry_context.reset();

    return 0;
  }
#endif
#endif
  return 1;
}

void chemistry_arepo_host_destroy(void)
{
#ifdef CHEMISTRY_LW_RADIATION
#ifdef CHEMISTRY_BACKGROUND_LW_RADIATION
  radiation_background_state_reset(&arepo_background_state);
#endif
  radiation_rate_context_destroy(&arepo_radiation_context);
#endif
  arepo_chemistry_context.reset();

}

const chemistry::Model *chemistry_arepo_host_context(void)
{
  return arepo_chemistry_context.get();
}

int chemistry_arepo_host_prepare_step(ChemistryStepContext *step)
{
  if(step == NULL)
    return 0;

  ChemistryGlobalEnvironment global = {};
  global.cmb_temperature = All.T_CMB;
  global.scale_factor = All.cf_atime;
  global.simulation_time = All.Time;
  global.comoving_integration = All.ComovingIntegrationOn;

#if defined(CHEMISTRY_LW_RADIATION) && defined(CHEMISTRY_BACKGROUND_LW_RADIATION)
  const RadiationBlackbodyInput background = {
      All.J21,
      All.Teff};
  if(!radiation_background_state_update(
         &arepo_radiation_context, &background,
         &arepo_background_state))
    return 0;

  global.background_J21 = arepo_background_state.input.J21;
  global.background_rate_Hminus_pd =
      arepo_background_state.rates.hminus_photodetachment_rate;
  global.background_rate_H2_diss =
      arepo_background_state.rates.h2_photodissociation_rate;
#ifdef CHEMISTRY_SPECIES_DEUTERIUM
  global.background_rate_HD_diss =
      arepo_background_state.rates.hd_photodissociation_rate;
#endif
#ifdef METALLICITY
  global.background_G0 = arepo_background_state.rates.g0;
#endif
#endif
  return arepo_chemistry_context &&
      chemistry::Driver::prepare(global, *step) == CHEMISTRY_STATUS_OK;
}

int chemistry_arepo_blackbody_rates(double effective_temperature,
                                    double *photodetachment_rate, double *g0)
{
  RadiationBlackbodySpectrumRates rates = {};
  if(photodetachment_rate == NULL || g0 == NULL)
    return 0;

#ifdef CHEMISTRY_LW_RADIATION
  if(effective_temperature < arepo_radiation_context.effective_temperature_min)
    effective_temperature = arepo_radiation_context.effective_temperature_min;
  if(effective_temperature > arepo_radiation_context.effective_temperature_max)
    effective_temperature = arepo_radiation_context.effective_temperature_max;
  if(!radiation_rate_blackbody_spectrum(&arepo_radiation_context,
                                        effective_temperature, &rates))
    return 0;
#else
  return 0;
#endif
  *photodetachment_rate = rates.hminus_photodetachment_rate;
  *g0 = rates.g0;
  return 1;
}

int chemistry_arepo_dust_opacities(double temperature,
                                   double *kappa_rosseland,
                                   double *kappa_planck)
{
  if(kappa_rosseland == NULL || kappa_planck == NULL)
    return 0;

  *kappa_rosseland = 0.0;
  *kappa_planck = 0.0;
#if CHEMISTRY_DUST_COOLING
  if(!arepo_chemistry_context || !isfinite(temperature) || temperature < 0.) return 0;
  chemistry_dust_opacities_from_tables(arepo_chemistry_context->dust_opacity(),
                                      temperature, kappa_rosseland, kappa_planck);
  return 1;
#else
  (void)temperature;
  return 0;
#endif
}
