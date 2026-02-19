#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "../allvars.h"
#include "../proto.h"
#include "../tags.h"

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT

#define ID_CHEM_DEBUG 56918776

#ifdef HYBRID_IMPLICIT
int Nimplicit  = 0;
int Nimplicit0 = 0;
#endif 
#ifdef ENERGY_IMPLICIT
int Ndiv  = 0;
int Ndiv0 = 0;
#endif


#ifdef CHEM_BALANCE
struct chemdata_in
{
  FLOAT Pos[3];
  int   ID;

  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */
  FLOAT Density;		/*!< current baryonic mass density of particle */

  FLOAT elec;
  FLOAT HI;
  FLOAT HII;
#ifdef HELIUM
  FLOAT HeI;
  FLOAT HeII;
  FLOAT HeIII;
#endif
  FLOAT H2I;
#ifdef H2_PLUS
  FLOAT H2II;
#endif
  FLOAT HM;

#ifdef HD
  FLOAT DI;
  FLOAT DII;
  FLOAT HDI;
#ifdef HD_PLUS
  FLOAT HDII;
#endif
#ifdef DMINUS
  FLOAT DM;
#endif
#endif

#ifdef METALLICITY
  FLOAT Metallicity;
#endif

  FLOAT Gamma;
  FLOAT t_chem, t_cool;

#ifdef SOBOLEV_APP
  FLOAT DelV[3];
#endif

#ifdef SHIELD_6RAY
  FLOAT ColumnDens_H2[6];  /* [cm^-2] */
  FLOAT ColumnDens_HI[6];
#ifdef HD
  FLOAT ColumnDens_HD[6];
#endif
#endif

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  FLOAT rate_Hminus_pd;
  FLOAT rate_H2_diss;
#endif
#ifdef DUST_COOLING
  FLOAT G0;
#endif
#endif

#ifdef IONIZATION_PHOTON
  FLOAT k22, piHI;
  FLOAT tau_ion;
#endif

#ifdef THERMAL_FEEDBACK
  FLOAT stellar_rad_flux;
#endif

#ifdef DENSITY_GRAD
  FLOAT l_core;
#endif
  
#ifdef FLAG_DEBUG_CHEMISTRY
  //FLOAT Entropy_old;
  //FLOAT Gamma_old;
  //FLOAT Density_old;
#endif

  integertime   Ti_begstep;
  integertime   Ti_endstep;
  integertime   ti_step;
}
  *ChemDataIn,
  *ChemDataGet;

struct chemdata_out
{
  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */

  FLOAT elec;
  FLOAT HI;
  FLOAT HII;
#ifdef HELIUM
  FLOAT HeI;
  FLOAT HeII;
  FLOAT HeIII;
#endif
  FLOAT H2I;
#ifdef H2_PLUS
  FLOAT H2II;
#endif
  FLOAT HM;

#ifdef HD
  FLOAT DI;
  FLOAT DII;
  FLOAT HDI;
#ifdef HD_PLUS
  FLOAT HDII;
#endif
#ifdef DMINUS
  FLOAT DM;
#endif
#endif

  FLOAT Gamma;
  FLOAT t_chem, t_cool;

#ifdef ENERGY_IMPLICIT
  bool flag_imp;
  bool flag_conv;
#endif

#ifdef METALLICITY
  FLOAT Tdust;
#endif
  
#ifdef FLAG_DEBUG_CHEMISTRY
  FLOAT edot_rad;
  FLOAT edot_cont;
  FLOAT edot_H2;
  FLOAT edot_ciHI;
  FLOAT edot_ceHI;
#ifdef HD
  FLOAT edot_HD;
#endif
#ifdef METALLICITY
  FLOAT edot_CI;
  FLOAT edot_CII;
  FLOAT edot_OI;
  FLOAT edot_dust;
#endif
  FLOAT edot_chem;
  FLOAT edot_H2form;
  FLOAT edot_HMform;

  FLOAT edot_H2form_3body;
  FLOAT edot_H2form_HM_loss;
  FLOAT edot_H2form_dust_loss;

  //FLOAT t_edot_chem;
  //FLOAT t_edot_rad;

  //FLOAT ThermalEnergy_new; 
  //FLOAT ThermalEnergy_old; 
  //FLOAT drhoU_comp; 
  //FLOAT drhoU_shock;
  //FLOAT drhoU_chem; 
  //FLOAT drhoU_rad;

  FLOAT f_reduce;
  FLOAT tau_cont;
#endif
}
  *ChemDataResult,
  *ChemDataOut;


void chem_balance(void);
void chem_evaluate(int target, int thread_id);
static inline void particle2in_chem(struct chemdata_in *in, int i);
static inline void out2particle_chem(struct chemdata_out *out, int i);
void *chem_evaluate_secondary(void *p);
#endif


void do_chemistry(void){
#ifdef HYBRID_IMPLICIT
  Nimplicit  = 0;
#endif 
#ifdef ENERGY_IMPLICIT
  Ndiv  = 0;
#endif

  MPI_Barrier(MPI_COMM_WORLD);

#ifndef CHEM_BALANCE
  int i;
  double t0, t1, timeall;

  t0 = my_second();
  /* Now assign new timesteps  */    
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type != 0)
	continue;

      if(SphP[i].flag_chem == false)
	continue;

      int bin;
      integertime ti_step;
  
      bin = P[i].TimeBin;              
      ti_step = bin ? (((integertime) 1) << bin) : 0;

      integertime tstart, tend;

      tstart = 0.5 * (P[i].Ti_begstep + P[i].Ti_endstep);
      tend   = P[i].Ti_endstep + 0.5 * ti_step;

      double a_start, a_end;
      if(All.ComovingIntegrationOn){  /* time in cosmic expansion parameter */ 
	a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
	a_end   = All.TimeBegin * exp(tend   * All.Timebase_interval);
      }else{
	a_start = tstart * All.Timebase_interval;
	a_end   = tend   * All.Timebase_interval;
      }



      int flag_implicit, flag_print, ifunc, mode;        
      ifunc = compute_abundances_d(mode = 1, i, a_start, a_end, &flag_implicit, flag_print=0);  
      if(ifunc == 1){
	ifunc = compute_abundances_d(mode = 1, i, a_start, a_end, &flag_implicit, flag_print=1);  
      }
      
#ifdef HYBRID_IMPLICIT
      if(flag_implicit) 
	Nimplicit++;    
#endif

#ifdef ENERGY_IMPLICIT
      if(flag_implicit)
	SphP[i].flag_imp = true;
      else
	SphP[i].flag_imp = false;

      if(flag_implicit && ifunc){
	SphP[i].flag_conv = true;
	Ndiv++;
      }else{
	SphP[i].flag_conv = false;
      }
#endif    
    }

  t1 = WallclockTime = my_second();
  timeall  = timediff(t0, t1);

  CPU_Step[CPU_CHEM]    += timeall;
#else
  chem_balance();
#endif


#ifdef HYBRID_IMPLICIT
  int ierr;

  ierr = MPI_Reduce(&Nimplicit, &Nimplicit0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  ierr = MPI_Reduce(&Ndiv, &Ndiv0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0){
    printf("chemistries of %d particles are calculated by implicit method.\n",Nimplicit0);
    if(Ndiv0)
      printf("energy-loop of %d particles are not converged\n",Ndiv0);
  }
#endif
}



#ifdef CHEM_BALANCE
void chem_balance(void)
{
  int i, j, k, m, n, ngrp, ndone, ndone_flag;
  int n_exported;
  int level, sendTask, recvTask, nexport, place;
  double t0, t1, timeall;
  double tstart, tend;
  double timewait1=0., timewait2=0., timewait;
  double timecomp2=0., timecomp;
  double timecommsumm1=0., timecommsumm2=0., timecomm;
  MPI_Status status;

  struct particle_data p_temp;

  
  if(ThisTask == 0) printf("\n  [timestep.c/chem_balance] start --> \n");  
 
  All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + 
							   sizeof(struct chemdata_in) +
							   sizeof(struct chemdata_out) +
							   sizemax(sizeof(struct chemdata_in),
								   sizeof(struct chemdata_out))));
  DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  NextParticle = FirstActiveParticle;	/* beginn with this index */

  t0 = my_second();
  do{ 

    Nexport = 0;
    for(j = 0; j < NTask; j++)      
      Send_count[j] = 0;
          
    /* do local particles and prepare export list */

    int exitFlag;

    while(1){
      exitFlag = 0;
      if(NextParticle < 0)
	{
	  exitFlag = 1;
	}
      else
	{
	  i = NextParticle;
	  NextParticle = NextActiveParticle[NextParticle];
	}

      if(Nexport >= All.BunchSize){
	BufferFullFlag = 1;
	exitFlag = 1;
      }

      if(exitFlag)
	break;

      if(P[i].Type != 0)
	continue;

      if(SphP[i].flag_chem == false)
	continue;

      int stask;
      stask = Nexport % (NTask - 1);
      if(stask >= ThisTask) stask++;

      DataIndexTable[Nexport].Task  = stask;
      DataIndexTable[Nexport].Index = i;
      DataIndexTable[Nexport].IndexGet = Nexport;
      
      Nexport++;
    }

    n_exported += Nexport;
    for(j = 0; j < NTask; j++)
      Send_count[j] = 0;
    for(j = 0; j < Nexport; j++)
      Send_count[DataIndexTable[j].Task]++;
        
    MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
    tstart = my_second();
    MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
    tend = my_second();
    timewait1 += timediff(tstart, tend);


    for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
      {
	Nimport += Recv_count[j];
	if(j > 0)
	  {
	    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	    Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	  }
      }
    ChemDataGet = (struct chemdata_in *) mymalloc("ChemDataGet", Nimport * sizeof(struct chemdata_in));
    ChemDataIn = (struct chemdata_in *) mymalloc("ChemDataIn", Nexport * sizeof(struct chemdata_in));

    /* prepare particle data for export */
    for(j = 0; j < Nexport; j++)
      {
	place = DataIndexTable[j].Index;
	particle2in_chem(&ChemDataIn[j], place);		// MADE D_IND CHANGE IN HERE
      }


    tstart = my_second();
    for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
      {
	recvTask = ThisTask ^ ngrp;
            
	if(recvTask < NTask)
	  {
	    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	      {
		/* get the particles */
		MPI_Sendrecv(&ChemDataIn[Send_offset[recvTask]],
			     Send_count[recvTask] * sizeof(struct chemdata_in), MPI_BYTE,
			     recvTask, TAG_CHEM_A,
			     &ChemDataGet[Recv_offset[recvTask]],
			     Recv_count[recvTask] * sizeof(struct chemdata_in), MPI_BYTE,
			     recvTask, TAG_CHEM_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	  }
      }
    tend = my_second();
    timecommsumm1 += timediff(tstart, tend);


    myfree(ChemDataIn);
    ChemDataResult = (struct chemdata_out *) mymalloc("ChemDataResult", Nimport * sizeof(struct chemdata_out));
    ChemDataOut = (struct chemdata_out *) mymalloc("ChemDataOut", Nexport * sizeof(struct chemdata_out));
    //report_memory_usage(&HighMark_sphhydro, "SPH_CHEM");

    tstart = my_second();
    NextJ = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
      int mainthreadid = omp_get_thread_num();
#else
      int mainthreadid = 0;
#endif
      chem_evaluate_secondary(&mainthreadid);
    }

    tend = my_second();
    timecomp2 += timediff(tstart, tend);
        
    if(NextParticle < 0)
      ndone_flag = 1;
    else
      ndone_flag = 0;
        

    tstart = my_second();
    MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    tend = my_second();
    timewait2 += timediff(tstart, tend);

    /* get the result */
    tstart = my_second();
    for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
      {
	recvTask = ThisTask ^ ngrp;
	if(recvTask < NTask)
	  {
	    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	      {
		/* send the results */
		MPI_Sendrecv(&ChemDataResult[Recv_offset[recvTask]],
			     Recv_count[recvTask] * sizeof(struct chemdata_out),
			     MPI_BYTE, recvTask, TAG_CHEM_B,
			     &ChemDataOut[Send_offset[recvTask]],
			     Send_count[recvTask] * sizeof(struct chemdata_out),
			     MPI_BYTE, recvTask, TAG_CHEM_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	  }
      }
    tend = my_second();
    timecommsumm2 += timediff(tstart, tend);
        
    /* add the result to the local particles */
    for(j = 0; j < Nexport; j++)
      {
	place = DataIndexTable[j].Index;
	out2particle_chem(&ChemDataOut[j], place);
      }

    myfree(ChemDataOut);
    myfree(ChemDataResult);
    myfree(ChemDataGet);
  }while(ndone < NTask);
  
  myfree(DataIndexTable);
  
  t1 = WallclockTime = my_second();
  
  timeall  = timediff(t0, t1);
  timecomp = timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;
  
  CPU_Step[CPU_CHEMCOMPUTE] += timecomp;
  CPU_Step[CPU_CHEMWAIT]    += timewait;
  CPU_Step[CPU_CHEMCOMM]    += timecomm;
  CPU_Step[CPU_CHEM]        += timeall;
}




void chem_evaluate(int target, int thread_id)
{
  int k, index;
  struct particle_data p_save;

  index = N_gas + thread_id;
  p_save = P[index];

  for(k=0;k<3;k++){
    P[index].Pos[k]           = ChemDataGet[target].Pos[k];
#ifdef SOBOLEV_APP
    SphP[index].DelV[k]       = ChemDataGet[target].DelV[k];
#endif
  }

#ifdef SHIELD_6RAY
  for(k=0; k<6; k++){
    SphP[index].ColumnDens_HI  =  ChemDataGet[target].ColumnDens_HI ;
    SphP[index].ColumnDens_H2I =  ChemDataGet[target].ColumnDens_H2I ;
#ifdef HD
    SphP[index].ColumnDens_HD  =  ChemDataGet[target].ColumnDens_HD ;
#endif
 }
#endif

  P[index].ID     = ChemDataGet[target].ID;

  SphP[index].Entropy = ChemDataGet[target].Entropy;
  SphP[index].Density = ChemDataGet[target].Density;

  SphP[index].elec  = ChemDataGet[target].elec;
  SphP[index].HI    = ChemDataGet[target].HI;
  SphP[index].HII   = ChemDataGet[target].HII;
#ifdef HELIUM
  SphP[index].HeI   = ChemDataGet[target].HeI;
  SphP[index].HeII  = ChemDataGet[target].HeII;
  SphP[index].HeIII = ChemDataGet[target].HeIII;
#endif
  SphP[index].H2I   = ChemDataGet[target].H2I;
#ifdef H2_PLUS
  SphP[index].H2II  = ChemDataGet[target].H2II;
#endif
  SphP[index].HM    = ChemDataGet[target].HM;
#ifdef HD
  SphP[index].DI   = ChemDataGet[target].DI;
  SphP[index].DII  = ChemDataGet[target].DII;
  SphP[index].HDI  = ChemDataGet[target].HDI;
#ifdef HD_PLUS
  SphP[index].HDII = ChemDataGet[target].HDII;
#endif
#ifdef DMINUS
  SphP[index].DM   = ChemDataGet[target].DM;
#endif
#endif
#ifdef METALLICITY
  SphP[index].Metallicity = ChemDataGet[target].Metallicity;
#endif

  SphP[index].Gamma = ChemDataGet[target].Gamma;
  SphP[index].t_chem = ChemDataGet[target].t_chem;
  SphP[index].t_cool = ChemDataGet[target].t_cool;

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  SphP[index].rate_Hminus_pd = ChemDataGet[target].rate_Hminus_pd;
  SphP[index].rate_H2_diss = ChemDataGet[target].rate_H2_diss;  
#endif

#ifdef DUST_COOLING
  SphP[index].G0 = ChemDataGet[target].G0;
#endif
#endif

#ifdef IONIZATION_PHOTON
  SphP[index].k22   = ChemDataGet[target].k22;
  SphP[index].piHI  = ChemDataGet[target].piHI;
  SphP[index].tau_ion  = ChemDataGet[target].tau_ion;
#endif

#ifdef THERMAL_FEEDBACK
  SphP[index].stellar_rad_flux = ChemDataGet[target].stellar_rad_flux;
#endif

#ifdef DENSITY_GRAD
  SphP[index].l_core  = ChemDataGet[target].l_core;  
#endif
  
#ifdef FLAG_DEBUG_CHEMISTRY
  //SphP[index].Entropy_old = ChemDataGet[target].Entropy_old;
  //SphP[index].Density_old = ChemDataGet[target].Density_old;
  //SphP[index].Gamma_old = ChemDataGet[target].Gamma_old;
#endif

  int ifunc, mode;
  integertime ti_step;
  integertime tstart, tend;
  double a_start, a_end;

  ti_step = ChemDataGet[target].ti_step;

  tstart = (ChemDataGet[target].Ti_begstep + ChemDataGet[target].Ti_endstep) / 2.;  /* midpoint of old step */
  tend   =  ChemDataGet[target].Ti_endstep + ti_step / 2.;            /* midpoint of new step */
  
  if(All.ComovingIntegrationOn){  /* time in cosmic expansion parameter */ 
    a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
    a_end   = All.TimeBegin * exp(tend   * All.Timebase_interval);
  }else{
    a_start = tstart * All.Timebase_interval;
    a_end   = tend   * All.Timebase_interval;
  }

  int flag_implicit, flag_print=0;

  ifunc = compute_abundances_d(mode = 1, index, a_start, a_end, &flag_implicit, flag_print=0);  
  if(ifunc == 1){
    //ifunc = compute_abundances_d(mode = 1, index, a_start, a_end, &flag_implicit, flag_print=1);  
  }
  

#ifdef ENERGY_IMPLICIT
#ifdef HYBRID_IMPLICIT
  if(flag_implicit){
    SphP[index].flag_imp = true;
    Nimplicit++;
  }else
    SphP[index].flag_imp = false;
#endif

  if(flag_implicit && ifunc){
    SphP[index].flag_conv = true;
    Ndiv++;
  }else{
    SphP[index].flag_conv = false;
  }
#endif

  P[index] = p_save;

  ChemDataResult[target].Entropy = SphP[index].Entropy;

  ChemDataResult[target].elec  = SphP[index].elec;
  ChemDataResult[target].HI    = SphP[index].HI;
  ChemDataResult[target].HII   = SphP[index].HII;
#ifdef HELIUM
  ChemDataResult[target].HeI   = SphP[index].HeI;
  ChemDataResult[target].HeII  = SphP[index].HeII;
  ChemDataResult[target].HeIII = SphP[index].HeIII;
#endif
  ChemDataResult[target].H2I   = SphP[index].H2I;
#ifdef H2_PLUS
  ChemDataResult[target].H2II  = SphP[index].H2II;
#endif
  ChemDataResult[target].HM    = SphP[index].HM;
#ifdef HD
  ChemDataResult[target].DI   = SphP[index].DI;
  ChemDataResult[target].DII  = SphP[index].DII;
  ChemDataResult[target].HDI  = SphP[index].HDI;
#ifdef HD_PLUS
  ChemDataResult[target].HDII = SphP[index].HDII;
#endif
#ifdef DMINUS
  ChemDataResult[target].DM   = SphP[index].DM;
#endif
#endif
  ChemDataResult[target].Gamma  = SphP[index].Gamma;
  ChemDataResult[target].t_chem = SphP[index].t_chem;
  ChemDataResult[target].t_cool = SphP[index].t_cool;

#ifdef METALLICITY
  ChemDataResult[target].Tdust = SphP[index].Tdust;
#endif
  
#ifdef ENERGY_IMPLICIT
  ChemDataResult[target].flag_imp  = SphP[index].flag_imp;
  ChemDataResult[target].flag_conv = SphP[index].flag_conv;
#endif

#ifdef FLAG_DEBUG_CHEMISTRY
  ChemDataResult[target].edot_rad  = SphP[index].edot_rad;
  ChemDataResult[target].edot_chem = SphP[index].edot_chem;
  ChemDataResult[target].edot_cont = SphP[index].edot_cont;
  ChemDataResult[target].edot_H2   = SphP[index].edot_H2;
  ChemDataResult[target].edot_ceHI = SphP[index].edot_ceHI;
  ChemDataResult[target].edot_ciHI = SphP[index].edot_ciHI;
#ifdef HD
  ChemDataResult[target].edot_HD   = SphP[index].edot_HD;
#endif
#ifdef METALLICITY
  ChemDataResult[target].edot_CI  = SphP[index].edot_CI;
  ChemDataResult[target].edot_CII = SphP[index].edot_CII;
  ChemDataResult[target].edot_OI  = SphP[index].edot_OI;
  ChemDataResult[target].edot_dust= SphP[index].edot_dust;
#endif
  ChemDataResult[target].edot_H2form = SphP[index].edot_H2form;
  ChemDataResult[target].edot_HMform = SphP[index].edot_HMform;
  
  ChemDataResult[target].edot_H2form_3body = SphP[index].edot_H2form_3body;
  ChemDataResult[target].edot_H2form_HM_loss = SphP[index].edot_H2form_HM_loss;
  ChemDataResult[target].edot_H2form_dust_loss = SphP[index].edot_H2form_dust_loss;

  //ChemDataResult[target].ThermalEnergy_new  = SphP[index].ThermalEnergy_new;
  //ChemDataResult[target].ThermalEnergy_old  = SphP[index].ThermalEnergy_old;
  //ChemDataResult[target].drhoU_rad  = SphP[index].drhoU_rad;
  //ChemDataResult[target].drhoU_chem = SphP[index].drhoU_chem;
  //ChemDataResult[target].drhoU_comp = SphP[index].drhoU_comp;
  //ChemDataResult[target].drhoU_shock = SphP[index].drhoU_shock;

  ChemDataResult[target].f_reduce = SphP[index].f_reduce;
  ChemDataResult[target].tau_cont = SphP[index].tau_cont;

  //ChemDataResult[target].t_edot_rad  = SphP[index].t_edot_rad;
  //ChemDataResult[target].t_edot_chem = SphP[index].t_edot_chem;
#endif
}




static inline void particle2in_chem(struct chemdata_in *in, int i)
{
  int k;
  for(k = 0; k < 3; k++){
    in->Pos[k] = P[i].Pos[k];
#ifdef SOBOLEV_APP
    in->DelV[k] = SphP[i].DelV[k];
#endif
  }
  
#ifdef SHIELD_6RAY
  for(k=0; k<6; k++){
    in->ColumnDens_HI  = SphP[i].ColumnDens_HI;
    in->ColumnDens_H2I = SphP[i].ColumnDens_H2I;
#ifdef HD
    in->ColumnDens_HD  = SphP[i].ColumnDens_HD;
#endif
  }
#endif
  in->ID    = P[i].ID;

  in->Entropy = SphP[i].Entropy;
  in->Density = SphP[i].Density;
	    
  in->elec  = SphP[i].elec;
  in->HI    = SphP[i].HI;
  in->HII   = SphP[i].HII;
#ifdef HELIUM
  in->HeI   = SphP[i].HeI;
  in->HeII  = SphP[i].HeII;
  in->HeIII = SphP[i].HeIII;
#endif
  in->H2I   = SphP[i].H2I;
#ifdef H2_PLUS
  in->H2II  = SphP[i].H2II;
#endif
  in->HM    = SphP[i].HM;
#ifdef HD
  in->DI   = SphP[i].DI;
  in->DII  = SphP[i].DII;
  in->HDI  = SphP[i].HDI;
#ifdef HD_PLUS
  in->HDII = SphP[i].HDII;
#endif
#ifdef DMINUS
  in->DM   = SphP[i].DM;
#endif
#endif
#ifdef METALLICITY
  in->Metallicity = SphP[i].Metallicity;
#endif

  in->Gamma = SphP[i].Gamma;
  in->t_chem = SphP[i].t_chem;
  in->t_cool = SphP[i].t_cool;

#ifdef LYMANWERNER
#ifdef LW_RADIATION_LOCAL
  in->rate_Hminus_pd = SphP[i].rate_Hminus_pd;
  in->rate_H2_diss = SphP[i].rate_H2_diss;
#endif

#ifdef DUST_COOLING
  in->G0 = SphP[i].G0;
#endif
#endif

#ifdef IONIZATION_PHOTON
  in->k22  = SphP[i].k22;
  in->piHI = SphP[i].piHI;
  in->tau_ion = SphP[i].tau_ion;
#endif

#ifdef THERMAL_FEEDBACK
  in->stellar_rad_flux = SphP[i].stellar_rad_flux;
#endif

#ifdef DENSITY_GRAD
  in->l_core = SphP[i].l_core;
#endif

#if defined(FLAG_FORBID_COLLAPSE) &&  defined(FLAG_FORBID_RAD)
  in->flag_chem = SphP[i].flag_chem;
#endif

#ifdef FLAG_DEBUG_CHEMISTRY
  //in->Entropy_old = SphP[i].Entropy_old;
  //in->Density_old = SphP[i].Density_old;
  //in->Gamma_old   = SphP[i].Gamma_old;
#endif

  in->Ti_begstep = P[i].Ti_begstep;
  in->Ti_endstep = P[i].Ti_endstep;

  int bin;
  integertime ti_step;
  
  bin = P[i].TimeBin;              
  ti_step = bin ? (((integertime) 1) << bin) : 0;

  in->ti_step = ti_step;
}



/* --------------------------------------------------------------------------------- */
/* this subroutine adds the output variables back to the particle values */
/* --------------------------------------------------------------------------------- */
static inline void out2particle_chem(struct chemdata_out *out, int i)
{
  SphP[i].Entropy = out->Entropy;
	      
  SphP[i].elec  = out->elec;
  SphP[i].HI    = out->HI;
  SphP[i].HII   = out->HII;
#ifdef HELIUM
  SphP[i].HeI   = out->HeI;
  SphP[i].HeII  = out->HeII;
  SphP[i].HeIII = out->HeIII;
#endif
  SphP[i].H2I   = out->H2I;
#ifdef H2_PLUS
  SphP[i].H2II  = out->H2II;
#endif
  SphP[i].HM    = out->HM;
#ifdef HD
  SphP[i].DI    = out->DI;
  SphP[i].DII   = out->DII;
  SphP[i].HDI   = out->HDI;
#ifdef HD_PLUS
  SphP[i].HDII  = out->HDII;
#endif
#ifdef DMINUS
  SphP[i].DM    = out->DM;     
#endif
#endif
  SphP[i].Gamma  = out->Gamma;

  SphP[i].t_chem = out->t_chem;
  SphP[i].t_cool = out->t_cool;
#ifdef METALLICITY
  SphP[i].Tdust  = out->Tdust;
#endif
#ifdef ENERGY_IMPLICIT
  SphP[i].flag_imp = out->flag_imp;
  SphP[i].flag_conv = out->flag_conv;
#endif

#ifdef FLAG_DEBUG_CHEMISTRY
  SphP[i].edot_rad  = out->edot_rad;
  SphP[i].edot_chem = out->edot_chem;
  SphP[i].edot_cont = out->edot_cont;
  SphP[i].edot_H2   = out->edot_H2;
  SphP[i].edot_ceHI = out->edot_ceHI;
  SphP[i].edot_ciHI = out->edot_ciHI;
#ifdef HD
  SphP[i].edot_HD   = out->edot_HD;  
#endif
#ifdef METALLICITY
  SphP[i].edot_CI   = out->edot_CI;
  SphP[i].edot_CII  = out->edot_CII;
  SphP[i].edot_OI   = out->edot_OI;
  SphP[i].edot_dust = out->edot_dust;
#endif
  SphP[i].edot_H2form = out->edot_H2form;
  SphP[i].edot_HMform = out->edot_HMform;

  SphP[i].edot_H2form_3body = out->edot_H2form_3body;
  SphP[i].edot_H2form_HM_loss = out->edot_H2form_HM_loss;
  SphP[i].edot_H2form_dust_loss = out->edot_H2form_dust_loss;
#endif

}

/* --------------------------------------------------------------------------------- */
/* one of the core sub-routines used to do the MPI version of the hydro evaluation
 (don't put actual operations here!!!) */
/* --------------------------------------------------------------------------------- */
void *chem_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;
  int j;
  while(1)
    {
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }
      UNLOCK_NEXPORT;
        
      if(j >= Nimport)
	break;
        
      chem_evaluate(j, thread_id);
    }
  return NULL;
}
#endif
