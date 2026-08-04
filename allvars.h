#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef CHEMISTRY_OLD
#ifdef HD
#define N_SPECIES 14
#else
#define N_SPECIES 9
#endif
#define TGCHEM_NUM_ABUNDANCES N_SPECIES
#endif

#ifdef CHEMISTRY_OLD
#ifndef METAL
#define k_MAX 60
#else
#define k_MAX 200
#endif
#else
#include "chemistry_new/chemistry_config.h"
#define k_MAX CHEM_RATE_COUNT
#endif

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif

#ifndef T_CMB0
#define  T_CMB0 2.73
#endif

#ifdef CHEMISTRY_OLD
#include "chemistry_old/chemistry.h"
#else
#include "chemistry_new/chemistry.h"
#endif

#define  mymalloc(x, y)            malloc(y)
#define  myfree(x)                 free(x)
#define  imax(x, y) ((x)>(y)?(x):(y))
#define  imin(x, y) ((x)<(y)?(x):(y))

#define  CMB_BG 0
#define  Ngrid 1

#define  NREACTION  200
/* ----- Declare some physical constants (cgs). */
#define  yr_to_sec    3.155815e7
#define  kms_to_cms   1.e5

#define  eV_to_K      11606.0
#define  eV_to_erg    1.60184e-12
#define  eV_to_Hz     2.41838e14

#define  GRAVITY     6.672e-8 //unit in cgs
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  BOLTZMANN   1.3806e-16

#define  STEFAN      5.6692e-5  /* erg cm-2 K-4 s-1 */ 
#define  GAMMA_MINUS1  0.6666666
#define  C_LIGHT      2.9979e10
#define  PLANCK      6.6262e-27
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28

#define  HYDROGEN_MASSFRAC 0.76  /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  MAX_REAL_NUMBER   1e30
#define  MAX_DBL_NUMBER   1e80

#define  HUBBLE      3.2407789e-18   /* in h/sec */

/* cosmology (WMAP 7 LCDM model) */
/*
#define OMEGA0      0.266
#define OMEGALAMBDA 0.734
#define OMEGABARYON 0.0449
#define H_PARAM 0.710 */
/* cosmology (PLANCK13 LCDM model) */
#define RHOCRITICAL 1.8783e-29
#define OMEGA0      0.3086
#define OMEGALAMBDA 0.6914
#define OMEGABARYON 0.04825
//#define H_PARAM     0.6777
#define H_PARAM      1.


struct io_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};


extern struct sph_particle_data
{
  double DMdensity;
  double Density;
  double DtDensity;
  double Entropy;
  double DtEntropy;
  double n_nuclei;
  double temperature;
  double temp_hold;
  double Gamma, molecular_weight;
  double pressure,ThermalEnergy;

  double elec;
  double HI,HII,HM;
#ifdef HELIUM
  double HeI,HeII,HeIII;
#endif
  double H2I;
#ifdef H2_PLUS
  double H2II;
#endif
#ifdef HD
  double DI,DII,HDI;
#ifdef HD_PLUS
  double HDII;
#endif
#ifdef DMINUS
  double DM;
#endif
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
  double CI, CII, CO, CO2;
  double OI, OH, H2O, O2;
#ifndef CHEM_INTERMEDIATE_EQUILIBRIUM
  double COII, CH, CH2;
  double OII, OHII, H2OII, H3OII, O2II;
#endif
#endif

  double t_col, t_cool, t_chem;   
  double tcool_rad, tcool_chem;
  
  double edot_ceHI,edot_ceHeI,edot_ceHeII;
  double edot_ciHI,edot_ciHeI,edot_ciHeIS,edot_ciHeII;
  double edot_ciHI0;
  double edot_reHII,edot_reHeII1,edot_reHeII2,edot_reHeIII;
  double edot_comp, edot_brem, edot_CIE, edot_H2, edot_HD;
  double edot_H2form_3body, edot_H2form_gas, edot_chem;
  double edot_H2form_HM, edot_H2form_H2II, edot_H2form_dust;
  double edot_adia, edot_adm, edot_grav;
  double edot_radHI, edot_radHeI, edot_radHeII;
  double edot_cont, edot_2ph;
  double edot_CI, edot_CII, edot_OI;
#ifdef CARBON_OXYGEN_CHEMISTRY
  double edot_CO, edot_OH, edot_H2O;
#endif
  double edot_dust, edot_PE;
  double edot_rad;

  double edot_H2form, edot_HMform;
  double edot_H2form_HM_loss, edot_H2form_dust_loss;
  
  double f_reduce;

  double k[k_MAX], k_inv[k_MAX];

  double n_crit;
  double Jeans_Length;

#ifdef SHIELD_6RAY
  double ColumnDens_H2[6], ColumnDens_HI[6], ColumnDens_HD[6];
#endif

#ifdef H2_OPACITY
  double DelV[3];
#endif
  double tau_cont;
  double fesc_H2;
#ifdef HD
  double fesc_HD;
#endif

#ifdef LYMANWERNER
  FLOAT rate_Hminus_pd;
  FLOAT rate_H2_diss;
#ifdef DUST_COOLING
  FLOAT G0;
#endif
#endif

#ifdef DENSITY_GRAD
  FLOAT GradRho[3];
  FLOAT l_core;
#endif

#if defined(IONIZATION_PHOTON) || defined(X_RAY)
  FLOAT  d_tau_0, d_r, dist_source;
  FLOAT  tau_ion, d_tau;
  double k22, piHI;
#ifdef X_RAY
  double k24, piHeI;
#ifdef X_RAY_2nd
  double k22_2nd,k24_2nd,piHI_2nd;
#endif
#endif

  FLOAT n_nuclei_ch;
  int   flag_ch;
#endif

  double Hsml;
  double SHR;
  double Entropy_old;
  
  int flag_div;

  double Tdust;
  double Metallicity;
} SphP[Ngrid];


extern double z_max, DMdensity_virial;

extern double g_weight[20];
extern double n_pop[3][20];

extern double A[3][20][3][20], B[3][20][3][20];
extern double E_vJ[3][20], E_vJ_erg[3][20];
extern double E_vJ_eV[3][20], E_vJ_K[3][20];
extern double del_E[3][20][3][20],del_E_K[3][20][3][20],del_nu[3][20][3][20];

extern double Lambda_frac[3][20][3][20], Lambda_frac_thick[3][20][3][20];
extern double alpha[3][20][3][20], n_pop_dummy[3][20][3][20], alpha_h2[3][20][3][20];


#ifdef CHEMISTRY_OLD
extern int N_SP, N_RE;
#else
extern int N_RE;
#endif

/* ----- Tables for the old chemistry module ------- */
#ifdef CHEMISTRY_OLD
extern double  T[N_T],J0_nu[N_nu],J_nu[N_nu],nu[N_nu];
extern double  k_table[k_MAX][N_T];
extern double  k_table_inverse[k_MAX][N_T];

extern double  ciHIa[N_T],ciHeIa[N_T],ciHeIIa[N_T],ciHeISa[N_T],reHIIa[N_T],brema[N_T],H2_CIE[N_T];
extern double  ceHIa[N_T],ceHeIa[N_T],ceHeIIa[N_T],reHeII1a[N_T],reHeII2a[N_T],reHeIIIa[N_T];

#if defined(BG_RADIATION) || defined(IONIZATION_PHOTON) || defined(X_RAY) || defined(CALC_CONT_OPACITY)
#define s_MAX 30

extern double sigma[s_MAX][N_nu];
extern double SIGMA_e22, SIGMA_e23, SIGMA_e24;
#endif
#endif


/* ----- Declare some control parameters. */
#define Y_C_sun  0.927e-4
#define Y_O_sun  3.568e-4

#define A_C_sun  3.16e-4
#define A_O_sun  5.75e-4

extern const int ThisTask;
extern const int NTask;

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))

#ifndef dmax
#define dmax(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef dmin
#define dmin(a,b) ((a)>(b)?(b):(a))
#endif

/*! This structure contains data which is the SAME for all tasks 
 *    (mostly code parameters read from the parameter file).  
 *  Holding this data in a structure is convenient for writing/reading the restart file, 
 *    and it allows the introduction of new global variables in a simple way. 
 *  The only thing to do is to introduce them into this structure.
 */
extern struct global_data_all_processes
{
  int TotNumPart;		/*!< total     particle numbers (global value) */
  int TotN_gas;	                /*!< total gas particle numbers (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of     particles that can be stored on one processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */

  double BoxSize;               /*!< Boxsize in case periodic boundary conditions are used */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;/*!< maximum number of files that may be written simultaneously when
				  writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				  NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				  the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				  the maximum(!) number of particles.  Note: A typical local tree for N
				  particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  double DesNumNgb;             /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


  /* some force counters  */

  long long TotNumOfForces;	             /*!< counts total number of force computations  */
  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */


  /* system of units  */

  double G;                        /*!< Gravity-constant in internal units */
  double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
  double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
  double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
  double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
  double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
  double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
  double UnitTime_in_Megayears;    /*!< factor to convert internal time to megayears/h */
  double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
				     gravitational constant is set to the Newtonian value based on the system of
				     units specified. Otherwise the value provided is taken as internal gravity constant G. */


  /* Cosmological parameters */

  double Hubble;       /*!< Hubble-constant in internal units */
  double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
  double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
  double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
  double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */
  

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */

  char ChemistryTableDir[1024];

  /* Parameters determining output frequency */
  int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */

#ifdef VARIABLEGAMMA
  double GM1;
#endif 


  double LenColumn;
  double ScaleFactor;



  /* Current time of the simulation, global step, and end of simulation */

  double Time;                  /*!< current time of the simulation */
  double TimeBegin;             /*!< time of initial conditions of the simulation */
  double TimeStep;              /*!< difference between current times of previous and current timestep */
  double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */


  /* variables for organizing discrete timeline */

  int Ti_Current;               /*!< current time on integer timeline */ 

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				  timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
				  timestep determined by the timestep criteria falls below this limit. */ 
  double MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
				  in comoving coordinates.  To this end, the code computes the rms velocity
				  of all particles, and limits the timestep such that the rms displacement
				  is a fraction of the mean particle separation (determined from the
				  particle mass and the cosmological parameters). This parameter specifies
				  this fraction. */

  double CourantFac;		/*!< SPH-Courant factor */

  /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
   * Five groups of particles are supported 0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"*/
  double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of SPH gravitational softening length */
  double MinGasHsml;            /*!< minimum allowed SPH smoothing length */

  double SofteningGas;          /*!< comoving gravitational softening lengths for type 0 */ 
  double SofteningHalo;         /*!< comoving gravitational softening lengths for type 1 */ 
  double SofteningDisk;         /*!< comoving gravitational softening lengths for type 2 */ 
  double SofteningBulge;        /*!< comoving gravitational softening lengths for type 3 */ 
  double SofteningStars;        /*!< comoving gravitational softening lengths for type 4 */ 
  double SofteningBndry;        /*!< comoving gravitational softening lengths for type 5 */ 

  double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */ 
  double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */ 
  double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */ 
  double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */ 
  double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */ 
  double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */ 

  double SofteningTable[6];     /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];     /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */

  double MassTable[6];          /*!< Table with particle masses for particle types with equal mass.
				  If particle masses are all equal for one type, the corresponding entry in MassTable 
				  is set to this value, allowing the size of the snapshot files to be reduced. */


#if defined(STATICNFW) || defined(PHYS_COSMO)
  double ExpansionParam;  /* only needed for runs with physical coordinate */
#endif

  double Epsilon;


#ifdef INITUNIFORM
  double InitElec, InitHeI, InitHeII, InitHeIII;
  double InitHI, InitHII, InitHM;
  double InitDI, InitDII, InitDM;
  double InitH2I, InitH2II, InitHDI, InitHDII;
  double InitLevel;
#endif


#ifdef BG_RADIATION
  double J21;
#else
#ifdef LYMANWERNER
  double J21, LymanWerner;

  FLOAT rate_Hminus_pd;
  FLOAT rate_H2_diss;
  //double J_LW;
#ifdef DUST_COOLING
  double G0;
#endif
#endif
#endif

#ifdef X_RAY
  double X_J21;
#endif

#ifdef SHIELD_LW
  double LenColumnDens;
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double SoftFraction;
#endif


  /* for particle splitting */
#ifdef FLAG_SPLIT
  double array[13][3];
#ifdef FLAG_SPLIT_DM
  FLOAT  pos_cen_split[VARI_MAX_SPLV][3];
  int    n_splv_sph[VARI_MAX_SPLV];
  int    n_splv_nbody[VARI_MAX_SPLV];
  int    n_split_sph[VARI_MAX_SPLV];
  int    n_split_nbody[VARI_MAX_SPLV];
#endif
#endif

  double T_CMB;
  FILE *fp_debug;

  FILE *fp_REAC_rate;
  char fname_REAC_rate[255];

  double ReacRate[60];
  double cf_atime, cf_hubble_a, cf_a3inv;
  double Metallicity;

} All;  /*!< a container variable for global variables that are equal on all processors */

extern struct particle_data
{
  FLOAT Pos[3];                 /*!< particle position at its current time */
  FLOAT Vel[3];                 /*!< particle velocity at its current time */
  FLOAT Mass;                   /*!< particle mass */
  FLOAT GravAccel[3];   /*!< particle acceleration due to gravity */

#ifndef LONGIDS
  unsigned int ID;              /*!< particle identifier */
#else
  unsigned long long ID;        /*!< particle identifier */
#endif

  int Type;                     /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */

  int Ti_endstep;               /*!< marks start of current timestep of particle on integer timeline */
  int Ti_begstep;               /*!< marks end of current timestep of particle on integer timeline */
}P[Ngrid];


extern double eta_a[8], kappa_a[8];
extern int FlagPrint;
