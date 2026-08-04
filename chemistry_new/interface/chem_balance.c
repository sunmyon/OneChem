#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"
#include "chemistry_constants.h"

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT

#ifndef ID_DEBUG_PRINT
#define ID_DEBUG_PRINT -1 //1005794452 //1006065322
#endif

//#define DEBUG_CHEM 1

#ifdef HYBRID_IMPLICIT
int Nimplicit  = 0;
int Nimplicit0 = 0;
#endif 
#ifdef ENERGY_IMPLICIT
int Ndiv  = 0;
int Ndiv0 = 0;
#endif
extern int count_print;

#ifdef CHEM_BALANCE
struct chemcomm_in
{
  double Mass;
  MyIDType ID;

  double Utherm_cgs;                /*!< current value of entropy (actually entropic function) of particle */
  double Density;		/*!< current baryonic mass density of particle */

  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;

#ifdef THERMAL_FEEDBACK
  double stellar_rad_flux;
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  double rate_Hminus_pd;
  double rate_H2_diss;
  double J21;
#ifdef HD
  double rate_HD_diss;
#endif
#ifdef DUST_COOLING
  double G0;
#endif
#endif
#endif

#ifdef IONIZATION_PHOTON
  double k22, piHI;
#endif

#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  double Metallicity;
#endif
  double dtime;

#ifdef DEBUG_CHEM
  int Task, Index;
#endif
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

#ifdef RYTR_TAU_DUST
  double kappa_dust;
#endif

#ifdef ENERGY_IMPLICIT
  char flag_imp;
  char flag_conv;
#endif
}
  *ChemDataResult,
  *ChemDataOut;

void chem_balance(void);
void chem_evaluate(int target, int thread_id);
static inline void particle2in_chem(struct chemcomm_in *in, int i);
static inline void out2particle_chem(struct chemcomm_out *out, int i);
void *chem_evaluate_secondary(void *p);
#endif

int compute_abundances_d(int mode, struct chemdata_in *sph_in, struct chemdata_out *sph_out, double dtime);

void cooling_only(void){
  CPU_Step[CPU_MISC] += measure_time();

#ifdef HYBRID_IMPLICIT
  Nimplicit  = 0;
#endif 
#ifdef ENERGY_IMPLICIT
  Ndiv  = 0;
#endif
  count_print = 0;

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
      double nH = All.UnitDensity_in_cgs * SphP[i].Density * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS * XH;

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
      
      struct chemdata_in sph_in;
      struct chemdata_out sph_out;
      particle2in_chemdata(i, &sph_in, &sph_out);
      
      int ifunc = compute_abundances_d(1, &sph_in, &sph_out, duration_s);
      if(ifunc == 1 && count_print == 0){
        ifunc = compute_abundances_d(1, &sph_in, &sph_out, duration_s);
        count_print++;
      }

      out2particle_chemdata(i, &sph_out);
      
#ifdef HYBRID_IMPLICIT
      if(sph_out.flag_implicit)
        Nimplicit++;
#endif

#ifdef ENERGY_IMPLICIT
      if(sph_out.flag_implicit && ifunc){
        Ndiv++;
      }
#endif
    }
#else
  chem_balance();
#endif

#ifdef HYBRID_IMPLICIT
  int ierr;
  ierr = MPI_Reduce(&Nimplicit, &Nimplicit0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef ENERGY_IMPLICIT
  ierr = MPI_Reduce(&Ndiv, &Ndiv0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  if(ThisTask == 0){
    printf("chemistries of %d particles are calculated by implicit method.\n",Nimplicit0);
#ifdef ENERGY_IMPLICIT
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

void chem_balance(void)
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
      double nH = All.UnitDensity_in_cgs * SphP[i].Density * All.HubbleParam * All.HubbleParam * All.cf_a3inv / PROTONMASS * XH;
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
      int mainthreadid = omp_get_thread_num();
#else
      int mainthreadid = 0;
#endif
      for(int j=0;j<nimport;j++)
	chem_evaluate(j, mainthreadid);
    }

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




void chem_evaluate(int target, int thread_id)
{
  struct chemdata_in sph_in;
  sph_in.Mass   = ChemDataGet[target].Mass;  
  sph_in.ID     = ChemDataGet[target].ID;
  sph_in.rho_cgs = All.UnitDensity_in_cgs * ChemDataGet[target].Density * All.HubbleParam * All.HubbleParam * All.cf_a3inv; /* physical gram / cm3 */

#ifdef METALLICITY
#ifdef CONSTANT_METALLICITY
  sph_in.Metallicity = All.Metallicity;
#elif defined(METAL_ADVECTION_SC)
  sph_in.Metallicity = ChemDataGet[target].Metallicity;
#endif
#endif
  
#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  sph_in.rate_Hminus_pd = ChemDataGet[target].rate_Hminus_pd;
  sph_in.rate_H2_diss = ChemDataGet[target].rate_H2_diss;  
#ifdef HD
  sph_in.rate_HD_diss = ChemDataGet[target].rate_HD_diss;
#endif
#ifdef DUST_COOLING
  sph_in.G0 = ChemDataGet[target].G0;
#endif
#endif
#endif

#ifdef THERMAL_FEEDBACK
  sph_in.stellar_rad_flux = ChemDataGet[target].stellar_rad_flux;
#endif

#ifdef IONIZATION_PHOTON
  sph_in.k22   = ChemDataGet[target].k22;
  sph_in.piHI  = ChemDataGet[target].piHI;
#endif

  sph_in.flag_print = 0;
#ifdef ID_DEBUG_PRINT
  if(sph_in.ID == ID_DEBUG_PRINT) sph_in.flag_print=1;
#endif
  
  struct chemdata_out sph_out;
  sph_out.Utherm_cgs = ChemDataGet[target].Utherm_cgs;

  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    sph_out.Abund[j] = ChemDataGet[target].Abund[j];
  
  sph_out.Gamma = ChemDataGet[target].Gamma;
  
  sph_out.t_chem = MAX_REAL_NUMBER;
  sph_out.t_cool = MAX_REAL_NUMBER;
#ifdef METALLICITY
  sph_out.Tdust = 0.;
#endif
#ifdef RYTR_TAU_DUST
  sph_out.kappa_dust = 0.;
#endif
  sph_out.flag_implicit = 0;
  sph_out.conv_flag = 0;
  
  double dtime= ChemDataGet[target].dtime;
  
  int ifunc = compute_abundances_d(1, &sph_in, &sph_out, dtime);  
  if(ifunc == 1 && count_print == 0){
    sph_in.flag_print = 1;
    ifunc = compute_abundances_d(1, &sph_in, &sph_out, dtime);    
    count_print++;
  }
  
#ifdef ENERGY_IMPLICIT
#ifdef HYBRID_IMPLICIT
  if(sph_out.flag_implicit){
    Nimplicit++;
  }
#endif
  if(sph_out.flag_implicit && ifunc){
    Ndiv++;
  }
#endif
  
  ChemDataResult[target].Utherm_cgs = sph_out.Utherm_cgs;

  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    ChemDataResult[target].Abund[j] = sph_out.Abund[j];

  ChemDataResult[target].Gamma  = sph_out.Gamma;
  ChemDataResult[target].t_chem = sph_out.t_chem;
  ChemDataResult[target].t_cool = sph_out.t_cool;

#ifdef METALLICITY
  ChemDataResult[target].Tdust = sph_out.Tdust;
#endif

#ifdef RYTR_TAU_DUST
  ChemDataResult[target].kappa_dust = sph_out.kappa_dust;
#endif

#ifdef ENERGY_IMPLICIT
  ChemDataResult[target].flag_imp = sph_out.flag_implicit;
  ChemDataResult[target].flag_conv = sph_out.conv_flag;
#endif
}


static inline void particle2in_chemdata(int i, struct chemdata_in *in, struct chemdata_out *out)
{
  in->flag_print = 0;
#ifdef ID_DEBUG_PRINT
  if(P[i].ID == ID_DEBUG_PRINT)
    in->flag_print=1;
#endif   

  in->Mass  = P[i].Mass;  
  in->ID    = P[i].ID;
  in->Density = SphP[i].Density;
  in->rho_cgs = All.UnitDensity_in_cgs * SphP[i].Density * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  
#ifdef METALLICITY
#ifdef CONSTANT_METALLICITY
  in->Metallicity = All.Metallicity;
#elif defined(METAL_ADVECTION_SC)
  double metallicity = 0.;
  for(int k=0;k<N_METAL_SPECIES;k++)
    metallicity += SphP[i].Metallicity_spec[k];
  
  in->Metallicity = metallicity / Z_SUN;    
#endif
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  in->rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  in->rate_H2_diss = SphP[i].rate_H2_diss;
  in->J21 = SphP[i].J21;
#ifdef HD
  in->rate_HD_diss = SphP[i].rate_HD_diss;
#endif
#ifdef DUST_COOLING
  in->G0 = SphP[i].G0;
#endif
#endif
#endif

#ifdef THERMAL_FEEDBACK
  in->stellar_rad_flux = SphP[i].stellar_rad_flux;
#endif

#ifdef IONIZATION_PHOTON
  in->k22  = SphP[i].k22;
  in->piHI = SphP[i].piHI;
#endif

  out->Utherm_cgs = SphP[i].Utherm * All.UnitEnergy_in_cgs/All.UnitMass_in_g;
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    out-> Abund[j] = SphP[i].Abund[j];

  out->Gamma = SphP[i].Gamma;  
  out->t_chem = MAX_REAL_NUMBER;
  out->t_cool = MAX_REAL_NUMBER;
  out->Tdust = 0.;
  out->kappa_dust = 0.;
  out->flag_implicit = 0;
  out->conv_flag = 0;
}


/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_chemdata(int i, struct chemdata_out *out)
{
  double du = out->Utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) - SphP[i].Utherm;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
  SphP[i].Utherm = out->Utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);

  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++){
    SphP[i].Abund[j] = out->Abund[j];
    SphP[i].AbundMass[j] = out->Abund[j] * P[i].Mass;
  }

  SphP[i].Gamma  = out->Gamma;

#ifndef CHEMISTRY_IMPLICIT
  SphP[i].t_chem = out->t_chem;
#endif
#ifndef ENERGY_IMPLICIT
  SphP[i].t_cool = out->t_cool;
#endif
  
#ifdef METALLICITY
  SphP[i].Tdust  = out->Tdust;
#endif
#ifdef RYTR_TAU_DUST
  SphP[i].kappa_dust = out->kappa_dust;
#endif
}


static inline void particle2in_chem(struct chemcomm_in *in, int i)
{
  in->Mass  = P[i].Mass;  
  in->ID    = P[i].ID;

  in->Utherm_cgs = All.UnitEnergy_in_cgs/All.UnitMass_in_g * SphP[i].Utherm;
  in->Density = SphP[i].Density;
	    
  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    in -> Abund[j] = SphP[i].Abund[j];

  in->Gamma = SphP[i].Gamma;

#if defined(METALLICITY) && defined(METAL_ADVECTION_SC) && !defined(CONSTANT_METALLICITY)
  double metallicity = 0.;
  for(int k=0;k<N_METAL_SPECIES;k++)
    metallicity += SphP[i].Metallicity_spec[k];
  
  in->Metallicity = metallicity / Z_SUN;    
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  in->rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  in->rate_H2_diss = SphP[i].rate_H2_diss;
  in->J21 = SphP[i].J21;
#ifdef HD
  in->rate_HD_diss = SphP[i].rate_HD_diss;
#endif
#ifdef DUST_COOLING
  in->G0 = SphP[i].G0;
#endif
#endif
#endif

#ifdef THERMAL_FEEDBACK
  in->stellar_rad_flux = SphP[i].stellar_rad_flux;
#endif

#ifdef IONIZATION_PHOTON
  in->k22  = SphP[i].k22;
  in->piHI = SphP[i].piHI;
#endif

  double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
  double dtime = All.cf_atime * dt / All.cf_time_hubble_a;
  const double duration_s = dtime * All.UnitTime_in_s / All.HubbleParam;
  in->dtime = duration_s;

#ifdef DEBUG_CHEM
  in->Task = ThisTask;
  in->Index = i;
#endif
}


/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_chem(struct chemcomm_out *out, int i)
{
  double du = out->Utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g) - SphP[i].Utherm;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;    
  SphP[i].Utherm = out->Utherm_cgs / (All.UnitEnergy_in_cgs/All.UnitMass_in_g);

  for(int j = 0; j < TGCHEM_NUM_ABUNDANCES; j++){
    SphP[i].Abund[j] = out->Abund[j];
    SphP[i].AbundMass[j] = out->Abund[j] * P[i].Mass;
  }

  SphP[i].Gamma  = out->Gamma;
#ifndef CHEMISTRY_IMPLICIT
  SphP[i].t_chem = out->t_chem;
#endif
#ifndef ENERGY_IMPLICIT  
  SphP[i].t_cool = out->t_cool;
#endif
  
#ifdef METALLICITY
  SphP[i].Tdust  = out->Tdust;
#endif
#ifdef RYTR_TAU_DUST
  SphP[i].kappa_dust = out->kappa_dust;
#endif
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU. */
int data_index_compare(const void *a, const void *b)
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
