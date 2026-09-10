EXEC   = 1ZONE

#.KEEP_STATE:

OPT1  += -DOMUKAI05
OPT1  += -DCHEMISTRY_STANDALONE_HOST
#OPT1  += -DSHR_ITER
#OPT1  += -DHIGH_TEMP
#OPT1  += -DFULLY_MOLECULAR
#OPT1  += -DHIGH_Z

#OPT1  +=  -DDEBUG_PRINT_FLAG
#OPT1  +=  -DID_CHEM_DEBUG=0
USE_CHEMISTRY_OLD ?= 0
ifeq ($(USE_CHEMISTRY_OLD),1)
USE_CHEMISTRY_PORTABLE ?= 0
else
USE_CHEMISTRY_PORTABLE ?= 1
endif
ifeq ($(USE_CHEMISTRY_OLD),1)
OPT1  +=  -DCHEMISTRY_OLD
endif
ifeq ($(USE_CHEMISTRY_OLD)$(USE_CHEMISTRY_PORTABLE),11)
$(error USE_CHEMISTRY_OLD and USE_CHEMISTRY_PORTABLE cannot both be 1)
endif
ifeq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 += -DCHEMISTRY_PORTABLE
endif
ENABLE_LW ?= 0
WRITE_CHEMISTRY_DIAGNOSTIC_TABLES ?= 0
UPDATE_L_CORE ?= 0
EXPLICIT_CHEMISTRY_UPDATE ?= SEQUENTIAL
# Intermediate treatment: equilibrium guess on/off, then FULL or REDUCED Newton.
INTERMEDIATE_EQUILIBRIUM ?= 1
EQUILIBRIUM_NEWTON_SYSTEM ?= FULL
NEWTON_LU_REUSE ?= 1
EQUILIBRIUM_SPECIES_CONVERGENCE ?= 0
CARBON_OXYGEN_CHEMISTRY ?= 1
CARBON_OXYGEN_CONSERVATION ?= 1
MOLECULAR_LINE_COOLING ?= 1
EXTRA_OPT ?=
OPTIMIZATION ?= -O2
OPT1 += $(EXTRA_OPT)
ifeq ($(EXPLICIT_CHEMISTRY_UPDATE),SIMULTANEOUS)
OPT1 += -DCHEM_EXPLICIT_SIMULTANEOUS
else ifneq ($(EXPLICIT_CHEMISTRY_UPDATE),SEQUENTIAL)
$(error EXPLICIT_CHEMISTRY_UPDATE must be SEQUENTIAL or SIMULTANEOUS)
endif
ifeq ($(CARBON_OXYGEN_CHEMISTRY),1)
OPT1 += -DCARBON_OXYGEN_CHEMISTRY
ifeq ($(CARBON_OXYGEN_CONSERVATION),0)
OPT1 += -DCHEM_DISABLE_CARBON_OXYGEN_CONSERVATION
else ifneq ($(CARBON_OXYGEN_CONSERVATION),1)
$(error CARBON_OXYGEN_CONSERVATION must be 0 or 1)
endif
ifeq ($(MOLECULAR_LINE_COOLING),1)
OPT1 += -DMOLECULAR_LINE_COOLING
else ifneq ($(MOLECULAR_LINE_COOLING),0)
$(error MOLECULAR_LINE_COOLING must be 0 or 1)
endif
else ifneq ($(CARBON_OXYGEN_CHEMISTRY),0)
$(error CARBON_OXYGEN_CHEMISTRY must be 0 or 1)
endif
ifeq ($(INTERMEDIATE_EQUILIBRIUM),1)
ifeq ($(EQUILIBRIUM_NEWTON_SYSTEM),FULL)
OPT1 += -DCHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
else ifeq ($(EQUILIBRIUM_NEWTON_SYSTEM),REDUCED)
OPT1 += -DCHEM_INTERMEDIATE_EQUILIBRIUM
else
$(error EQUILIBRIUM_NEWTON_SYSTEM must be FULL or REDUCED)
endif
else ifeq ($(INTERMEDIATE_EQUILIBRIUM),0)
ifneq ($(EQUILIBRIUM_NEWTON_SYSTEM),FULL)
$(error EQUILIBRIUM_NEWTON_SYSTEM=REDUCED requires INTERMEDIATE_EQUILIBRIUM=1)
endif
else
$(error INTERMEDIATE_EQUILIBRIUM must be 0 or 1)
endif
ifeq ($(NEWTON_LU_REUSE),1)
ifeq ($(EQUILIBRIUM_NEWTON_SYSTEM),FULL)
OPT1 += -DCHEM_NEWTON_LU_REUSE
endif
else ifneq ($(NEWTON_LU_REUSE),0)
$(error NEWTON_LU_REUSE must be 0 or 1)
endif
ifeq ($(EQUILIBRIUM_SPECIES_CONVERGENCE),1)
ifeq ($(INTERMEDIATE_EQUILIBRIUM)$(EQUILIBRIUM_NEWTON_SYSTEM),1FULL)
OPT1 += -DCHEM_EQUILIBRIUM_SPECIES_CONVERGENCE
endif
else ifneq ($(EQUILIBRIUM_SPECIES_CONVERGENCE),0)
$(error EQUILIBRIUM_SPECIES_CONVERGENCE must be 0 or 1)
endif
ifeq ($(UPDATE_L_CORE),1)
OPT1 += -DUPDATE_L_CORE_IN_SOLVER
endif

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
OPT1  += -DSURVEY_METALLICITY
OPT1  += -DSURVEY_CMB_REDSHIFT
ifeq ($(ENABLE_LW),1)
OPT1  += -DSURVEY_LYMANWERNER
endif

ifeq ($(WRITE_CHEMISTRY_DIAGNOSTIC_TABLES),1)
OPT1 += -DCHEMISTRY_WRITE_DIAGNOSTIC_TABLES
else ifneq ($(WRITE_CHEMISTRY_DIAGNOSTIC_TABLES),0)
$(error WRITE_CHEMISTRY_DIAGNOSTIC_TABLES must be 0 or 1)
endif

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

OPT1  +=  -DMETALLICITY=-4.
OPT1  +=  -DMETAL_COOLING_LW
OPT1  +=  -DDUST_COOLING
OPT1  +=  -DH2_FORMATION_ON_DUST
#OPT1  +=  -DNO_H2_DESTRUCTION

OPT1  +=  -DH2_OPACITY 
OPT1  +=  -DH2_OPACITY_FIT
OPT1  +=  -DCONT_OPACITY

ifeq ($(ENABLE_LW),1)
OPT1  +=  -DLYMANWERNER
OPT1  +=  -DBACKGROUND_LW_RADIATION
#OPT1  +=  -DTYPE_A
OPT1  +=  -DTYPE_B   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_3E4   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_POW   #type of the spectrum of LW intensity
OPT1  +=  -DSHIELD_LW 
#OPT1  +=  -DSHIELD_DENS=1.e9
endif

#OPT1  +=  -DIONIZATION_PHOTON
#OPT1  +=  -DSUBGRID_CHAMPAGNE_FLOW
#OPT1  +=  -DX_RAY=1.e-3
#OPT1  +=  -DX_RAY_2nd

GSL_DIR ?= $(shell gsl-config --prefix 2>/dev/null || echo /opt/homebrew)
LOCAL_PREFIX ?= /home/chon/usr/local
LOCAL_INCL = $(if $(wildcard $(LOCAL_PREFIX)/include),-I$(LOCAL_PREFIX)/include)
LOCAL_LIBS = $(if $(wildcard $(LOCAL_PREFIX)/lib),-L$(LOCAL_PREFIX)/lib)
GSL_INCL =  -I$(GSL_DIR)/include $(LOCAL_INCL)
GSL_LIBS =  -L$(GSL_DIR)/lib $(LOCAL_LIBS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LAPACK_LIBS = -framework Accelerate
else
LAPACK_LIBS = -llapack -lblas -lgfortran
endif

CFLAGS =  $(OPT1) -g $(OPTIMIZATION) $(GSL_INCL) -I. -Isrc -Isrc/standalone $(CHEMISTRY_INCLUDE_FLAGS) -Wformat

LIBS   =  $(GSL_INCL) $(GSL_LIBS) -lm $(LAPACK_LIBS) -lgsl -lgslcblas

CC     =  mpicc

ifneq (,$(findstring -DCHEMISTRY_OLD,$(OPT1)))
include legacy/chemistry_old/standalone.mk
else ifeq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 += -DCHEMISTRY_SC
CHEMISTRY_DIR = src/chemistry
include src/standalone/standalone.mk
CHEMISTRY_OBJS = $(CHEMISTRY_PORTABLE_CORE_SOURCES:.c=.o) \
                 $(CHEMISTRY_STANDALONE_ADAPTER_SOURCES:.c=.o)
CHEMISTRY_SUPPORT_OBJS =
CHEMISTRY_INCL = $(CHEMISTRY_DIR)/chemistry_core.h \
                 $(CHEMISTRY_DIR)/chemistry_build_config.h \
                 $(CHEMISTRY_DIR)/chemistry_types.h \
                 $(CHEMISTRY_DIR)/initialization/chemistry_init.h \
                 $(CHEMISTRY_DIR)/solver/chemistry_newton_diagnostics.h \
                 src/standalone/chemistry_standalone_adapter.h \
                 src/standalone/chemistry_standalone_types.h
CHEMISTRY_INCLUDE_FLAGS = -I$(CHEMISTRY_DIR) \
                          -I$(CHEMISTRY_DIR)/include \
                          -I$(CHEMISTRY_DIR)/common \
                          -I$(CHEMISTRY_DIR)/cooling \
                          -I$(CHEMISTRY_DIR)/driver \
                          -I$(CHEMISTRY_DIR)/initialization \
                          -I$(CHEMISTRY_DIR)/reaction_network \
                          -I$(CHEMISTRY_DIR)/solver \
                          -I$(CHEMISTRY_DIR)/thermodynamics \
                          -I$(CHEMISTRY_DIR)/tables \
                          -Isrc/standalone
else
OPT1 += -DCHEMISTRY_SC
CHEMISTRY_DIR  = legacy/chemistry_new
CHEMISTRY_OBJS = $(CHEMISTRY_DIR)/driver/chemistry_driver.o \
                 $(CHEMISTRY_DIR)/species.o \
                 $(CHEMISTRY_DIR)/initialization/chemistry_init.o \
                 $(CHEMISTRY_DIR)/solver/chemistry_solver.o \
                 $(CHEMISTRY_DIR)/solver/saha.o \
                 $(CHEMISTRY_DIR)/common/chemistry_models.o \
                 $(CHEMISTRY_DIR)/common/chemistry_profile.o \
                 $(CHEMISTRY_DIR)/reaction_network/forward_rates.o \
                 $(CHEMISTRY_DIR)/reaction_network/reaction_rates.o \
                 $(CHEMISTRY_DIR)/reaction_network/reaction_network.o \
                 $(CHEMISTRY_DIR)/reaction_network/shielding.o \
                 $(CHEMISTRY_DIR)/cooling/H2cooling.o \
                 $(CHEMISTRY_DIR)/cooling/chemistry_cooling.o \
                 $(CHEMISTRY_DIR)/cooling/cooling_func.o \
                 $(CHEMISTRY_DIR)/cooling/molecular_cooling.o \
                 $(CHEMISTRY_DIR)/cooling/dust_cooling.o \
                 $(CHEMISTRY_DIR)/thermodynamics/partition.o \
                 $(CHEMISTRY_DIR)/thermodynamics/thermodynamics.o
CHEMISTRY_SUPPORT_OBJS =
CHEMISTRY_INCL = $(CHEMISTRY_DIR)/chemistry.h \
                 $(CHEMISTRY_DIR)/chemistry_api.h \
                 $(CHEMISTRY_DIR)/species.h \
                 $(CHEMISTRY_DIR)/species_definitions.h \
                 $(CHEMISTRY_DIR)/chemistry_build_config.h \
                 $(CHEMISTRY_DIR)/chemistry_types.h \
                 $(CHEMISTRY_DIR)/common/chemistry_constants.h \
                 $(CHEMISTRY_DIR)/common/chemistry_diagnostics.h \
                 $(CHEMISTRY_DIR)/common/chemistry_models.h \
                 $(CHEMISTRY_DIR)/common/chemistry_profile.h \
                 $(CHEMISTRY_DIR)/cooling/chemistry_cooling.h \
                 $(CHEMISTRY_DIR)/cooling/cooling_func.h \
                 $(CHEMISTRY_DIR)/cooling/molecular_cooling.h \
                 $(CHEMISTRY_DIR)/cooling/dust_cooling.h \
                 $(CHEMISTRY_DIR)/cooling/H2cooling.h \
                 $(CHEMISTRY_DIR)/thermodynamics/partition.h \
                 $(CHEMISTRY_DIR)/thermodynamics/thermodynamics.h \
                 $(CHEMISTRY_DIR)/reaction_network/forward_rates.h \
                 $(CHEMISTRY_DIR)/reaction_network/reaction_network.h \
                 $(CHEMISTRY_DIR)/reaction_network/reaction_rates.h \
                 $(CHEMISTRY_DIR)/reaction_network/shielding.h \
                 $(CHEMISTRY_DIR)/solver/chemistry_solver.h \
                 $(CHEMISTRY_DIR)/solver/saha.h \
                 $(CHEMISTRY_DIR)/initialization/chemistry_init.h \
                 $(CHEMISTRY_DIR)/tables/chemistry_tables.h
CHEMISTRY_INCLUDE_FLAGS = -I$(CHEMISTRY_DIR) \
                          -I$(CHEMISTRY_DIR)/common \
                          -I$(CHEMISTRY_DIR)/cooling \
                          -I$(CHEMISTRY_DIR)/driver \
                          -I$(CHEMISTRY_DIR)/initialization \
                          -I$(CHEMISTRY_DIR)/reaction_network \
                          -I$(CHEMISTRY_DIR)/solver \
                          -I$(CHEMISTRY_DIR)/thermodynamics \
                          -I$(CHEMISTRY_DIR)/tables
endif

OBJS   = src/standalone/main.o src/standalone/standalone_run.o \
         src/standalone/standalone_collapse.o src/standalone/standalone_host.o \
         src/standalone/timescales.o src/radiation/radiation_rates.o \
	 $(CHEMISTRY_SUPPORT_OBJS) \
	 $(CHEMISTRY_OBJS)

INCL   = src/standalone/allvars.h src/standalone/standalone_chemistry_compat.h \
         src/standalone/standalone_run.h src/standalone/standalone_host.h \
         $(CHEMISTRY_INCL)

BUILDDIR ?= build
OBJDIR_OBJS := $(addprefix $(BUILDDIR)/,$(OBJS))

$(EXEC): $(OBJDIR_OBJS)
	$(CC) $(CFLAGS) $(OBJDIR_OBJS) $(LIBS) -o $(EXEC)

$(BUILDDIR)/%.o: %.c $(INCL)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(EXEC)
