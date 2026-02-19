EXEC   = 1ZONE

#.KEEP_STATE:

OPT1  += -DOMUKAI05
#OPT1  += -DSHR_ITER
#OPT1  += -DHIGH_TEMP
#OPT1  += -DFULLY_MOLECULAR
#OPT1  += -DHIGH_Z

#OPT1  +=  -DDEBUG_PRINT_FLAG
#OPT1  +=  -DID_CHEM_DEBUG=0

OPT1  += -DSINK
OPT1  += -DSPLIT
OPT1  += -DGAMMA_NEW

#OPT1  +=  -DFLOAT=float
OPT1  +=  -DFLOAT=double
#OPT1  +=  -DDEBUG_1ZONE
#OPT1  +=  -DMODE_TEST
#OPT1  +=  -DMODE_REPEAT
#OPT1  +=  -DSLOW_COLLAPSE=0.5
#OPT1  +=  -DRANDOM_EXPANSION=0.01
#OPT1  +=  -DADIABATIC_EXPANSION
#OPT1  +=  -DSTELLAR_FUV
#OPT1  +=  -DSURVEY_NT_REGION
#OPT1  += -DSURVEY_METALLICITY
#OPT1  += -DSURVEY_CMB_REDSHIFT
OPT1  += -DSURVEY_LYMANWERNER

#OPT1  += -DDENSITY_GRAD
#OPT1  +=  -DJEANS_LENGTH_TMAX=1.e4

OPT1  +=  -DVARIABLE_GAMMA
OPT1  +=  -DPHYS_COSMO
OPT1  +=  -DCMB_HEAT_PROPER
OPT1  +=  -DCMB_HEATING

OPT1  +=  -DCHEMISTRY

OPT1  +=  -DENERGY_IMPLICIT
OPT1  +=  -DCHEMISTRY_IMPLICIT
#OPT1  +=  -DHYBRID_IMPLICIT
OPT1  +=  -DFLAG_BRENT
OPT1  +=  -DCHEM_DIFF_ANALYTIC

#OPT1  +=  -DLIMIT_DT_IMPLICIT=1.e11
#OPT1  +=  -DHELIUM
#OPT1  +=  -DH2_PLUS

OPT1  +=  -DHD
#OPT1  +=  -DHD_PLUS
#OPT1  +=  -DDMINUS

OPT1  +=  -DINVERSE_REACTIONS
OPT1  +=  -DSAHA

#OPT1  +=  -DCALC_CONTINUUM_OPACITY
OPT1  +=  -DTHREEBODY 
OPT1  +=  -DMARTIN97
OPT1  +=  -DGA08
OPT1  +=  -DH2_COOLING_G15
#OPT1  +=  -DH2_COOLING_O00
#OPT1  +=  -DTBR_FLOW07
OPT1  +=  -DTBR_FORR13 
#OPT1  +=  -DENZO_CHEM
#OPT1  +=  -DO05_CHEM
OPT1  +=  -DCONTINUUM_COOLING_MATSUKOBA

OPT1  +=  -DCASE_B

#OPT1  +=  -DMETALLICITY=-4.
#OPT1  +=  -DMETAL_COOLING_LW
#OPT1  +=  -DDUST_COOLING
#OPT1  +=  -DH2_FORMATION_ON_DUST
#OPT1  +=  -DNO_H2_DESTRUCTION

OPT1  +=  -DH2_OPACITY 
OPT1  +=  -DH2_OPACITY_FIT
OPT1  +=  -DCONT_OPACITY

OPT1  +=  -DLYMANWERNER
OPT1  +=  -DBACKGROUND_LW_RADIATION
#OPT1  +=  -DTYPE_A
OPT1  +=  -DTYPE_B   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_3E4   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_POW   #type of the spectrum of LW intensity
OPT1  +=  -DSHIELD_LW 
#OPT1  +=  -DSHIELD_DENS=1.e9

#OPT1  +=  -DIONIZATION_PHOTON
#OPT1  +=  -DSUBGRID_CHAMPAGNE_FLOW
#OPT1  +=  -DX_RAY=1.e-3
#OPT1  +=  -DX_RAY_2nd

GSL_INCL =  -I${GSL_DIR}/include -I/home/chon/usr/local/include
GSL_LIBS =  -L${GSL_DIR}/lib -L/home/chon/usr/local/lib

CFLAGS =  $(OPT1) -g -O0 -I/home/chon/usr/local/include -L/home/chon/usr/local/lib -lblas -llapack -lgslcblas -Wformat

LIBS   =  $(GSL_INCL) $(GSL_LIBS) -lm -llapack -lblas -lgsl -lgslcblas -lgfortran

CC     =  cc

OBJS   = main.o allvars.o timescales.o \
	 chemistry/H2cooling.o chemistry/dust_cooling.o \
	 chemistry/chemistry_200601.o chemistry/cooling_func_200601.o \
         chemistry/saha.o chemistry/partition.o chemistry/stellar_spectrum.o

INCL   = allvars.h proto.h tags.h chemistry/chemistry.h

BUILDDIR := build
OBJDIR_OBJS := $(addprefix $(BUILDDIR)/,$(OBJS))

$(EXEC): $(OBJDIR_OBJS)
	$(CC) $(CFLAGS) $(OBJDIR_OBJS) $(LIBS) -o $(EXEC)

$(BUILDDIR)/%.o: %.c $(INCL)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(EXEC)
