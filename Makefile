EXEC   = 1ZONE

#.KEEP_STATE:

OPT1  += -DOMUKAI05
OPT1  += -DCHEMISTRY_STANDALONE_HOST
# Every standalone entry point (--chem-step, --chem-replay-log, the SphP[]
# particle loop, ...) supplies photoionization/photoheating rates as
# external input (CLI record fields or loaded particle data), never by
# computing them internally.  The legacy internal ionization-table
# initializers main.c calls when IONIZING_PHOTON is set and this is NOT
# defined (init_ionization_rad(), make_ionization_rate_table()) only exist
# in legacy/chemistry_old/, which is not linked into this (CHEMISTRY_PORTABLE)
# build, so without this the link fails whenever IONIZING_PHOTON is on.
OPT1  += -DSTANDALONE_EXTERNAL_RADIATION
#OPT1  += -DSHR_ITER
#OPT1  += -DHIGH_TEMP
#OPT1  += -DFULLY_MOLECULAR
#OPT1  += -DHIGH_Z

#OPT1  +=  -DDEBUG_PRINT_FLAG
#OPT1  +=  -DCHEMISTRY_DEBUG=1
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
ifneq ($(origin EXPLICIT_CHEMISTRY_UPDATE),undefined)
CHEMISTRY_EXPLICIT_UPDATE ?= $(EXPLICIT_CHEMISTRY_UPDATE)
endif
ifneq ($(origin INTERMEDIATE_EQUILIBRIUM),undefined)
CHEMISTRY_EQUILIBRIUM_INITIALIZATION ?= $(INTERMEDIATE_EQUILIBRIUM)
endif
ifneq ($(origin EQUILIBRIUM_NEWTON_SYSTEM),undefined)
CHEMISTRY_NEWTON_SYSTEM ?= $(EQUILIBRIUM_NEWTON_SYSTEM)
endif
ifneq ($(origin NEWTON_LU_REUSE),undefined)
CHEMISTRY_NEWTON_LU_REUSE ?= $(NEWTON_LU_REUSE)
endif
ifneq ($(origin EQUILIBRIUM_SPECIES_CONVERGENCE),undefined)
CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE ?= $(EQUILIBRIUM_SPECIES_CONVERGENCE)
endif
ifneq ($(origin CARBON_OXYGEN_CONSERVATION),undefined)
CHEMISTRY_CARBON_OXYGEN_CONSERVATION ?= $(CARBON_OXYGEN_CONSERVATION)
endif
ifneq ($(origin UPDATE_L_CORE),undefined)
CHEMISTRY_UPDATE_CORE_LENGTH ?= $(UPDATE_L_CORE)
endif
ifneq ($(origin WRITE_CHEMISTRY_DIAGNOSTIC_TABLES),undefined)
CHEMISTRY_WRITE_DIAGNOSTIC_TABLES ?= $(WRITE_CHEMISTRY_DIAGNOSTIC_TABLES)
endif
ENABLE_LW ?= 0
CHEMISTRY_WRITE_DIAGNOSTIC_TABLES ?= 0
CHEMISTRY_UPDATE_CORE_LENGTH ?=
# Preserve the current portable default (simultaneous); legacy used sequential.
ifeq ($(USE_CHEMISTRY_PORTABLE),1)
CHEMISTRY_EXPLICIT_UPDATE ?= SIMULTANEOUS
else
CHEMISTRY_EXPLICIT_UPDATE ?= SEQUENTIAL
endif
# Intermediate treatment: equilibrium guess on/off, then FULL or REDUCED Newton.
CHEMISTRY_EQUILIBRIUM_INITIALIZATION ?= 1
CHEMISTRY_NEWTON_SYSTEM ?= FULL
CHEMISTRY_NEWTON_LU_REUSE ?= 1
CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE ?= 0
# Accept the former Make variable at the host configuration boundary.
CHEMISTRY_SPECIES_CARBON_OXYGEN ?= $(if $(filter undefined,$(origin CARBON_OXYGEN_CHEMISTRY)),1,$(CARBON_OXYGEN_CHEMISTRY))
CHEMISTRY_CARBON_OXYGEN_CONSERVATION ?= 1
CHEMISTRY_MOLECULAR_LINE_COOLING ?= $(if $(filter undefined,$(origin MOLECULAR_LINE_COOLING)),,$(MOLECULAR_LINE_COOLING))
EXTRA_OPT ?=
# Derive chemistry-relevant compile options directly from a real AREPO
# Config.sh instead of hand-copying flags into EXTRA_OPT (see
# tools/arepo_config_to_extra_opt.py for the exact rule: only option names
# that ONEZONE's own source actually #ifdefs on are passed through, so
# AREPO-only physics in Config.sh -- gravity, mesh, sinks, MPI, ... -- is
# never blindly forwarded).
#   Usage: make AREPO_CONFIG=/path/to/Config.sh
AREPO_CONFIG ?=
ifneq ($(AREPO_CONFIG),)
AREPO_DERIVED_OPT := $(shell python3 $(CURDIR)/tools/arepo_config_to_extra_opt.py $(AREPO_CONFIG))
$(info Derived from $(AREPO_CONFIG): $(AREPO_DERIVED_OPT))
OPT1 += $(AREPO_DERIVED_OPT)
endif
OPTIMIZATION ?= -O2
OPT1 += $(EXTRA_OPT)
ifeq ($(CHEMISTRY_EXPLICIT_UPDATE),SIMULTANEOUS)
OPT1 += -DCHEMISTRY_EXPLICIT_SIMULTANEOUS
else ifeq ($(CHEMISTRY_EXPLICIT_UPDATE),SEQUENTIAL)
ifeq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 += -DCHEMISTRY_EXPLICIT_SEQUENTIAL
endif
else
$(error CHEMISTRY_EXPLICIT_UPDATE must be SEQUENTIAL or SIMULTANEOUS)
endif
ifeq ($(CHEMISTRY_SPECIES_CARBON_OXYGEN),1)
OPT1 += -DCHEMISTRY_SPECIES_CARBON_OXYGEN
ifeq ($(CHEMISTRY_CARBON_OXYGEN_CONSERVATION),0)
OPT1 += -DCHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION
else ifneq ($(CHEMISTRY_CARBON_OXYGEN_CONSERVATION),1)
$(error CHEMISTRY_CARBON_OXYGEN_CONSERVATION must be 0 or 1)
endif
else ifneq ($(CHEMISTRY_SPECIES_CARBON_OXYGEN),0)
$(error CHEMISTRY_SPECIES_CARBON_OXYGEN must be 0 or 1)
endif
ifneq ($(CHEMISTRY_MOLECULAR_LINE_COOLING),)
ifneq ($(filter $(CHEMISTRY_MOLECULAR_LINE_COOLING),0 1),$(CHEMISTRY_MOLECULAR_LINE_COOLING))
$(error CHEMISTRY_MOLECULAR_LINE_COOLING must be 0 or 1)
endif
OPT1 += -DCHEMISTRY_MOLECULAR_LINE_COOLING=$(CHEMISTRY_MOLECULAR_LINE_COOLING)
endif
ifneq ($(USE_CHEMISTRY_PORTABLE),1)
ifeq ($(CHEMISTRY_SPECIES_CARBON_OXYGEN),1)
ifneq ($(CHEMISTRY_MOLECULAR_LINE_COOLING),0)
OPT1 += -DMOLECULAR_LINE_COOLING
endif
endif
endif
ifeq ($(CHEMISTRY_EQUILIBRIUM_INITIALIZATION),1)
ifeq ($(CHEMISTRY_NEWTON_SYSTEM),FULL)
OPT1 += -DCHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION
else ifeq ($(CHEMISTRY_NEWTON_SYSTEM),REDUCED)
OPT1 += -DCHEMISTRY_INTERMEDIATE_EQUILIBRIUM
else
$(error CHEMISTRY_NEWTON_SYSTEM must be FULL or REDUCED)
endif
else ifeq ($(CHEMISTRY_EQUILIBRIUM_INITIALIZATION),0)
ifneq ($(CHEMISTRY_NEWTON_SYSTEM),FULL)
$(error CHEMISTRY_NEWTON_SYSTEM=REDUCED requires CHEMISTRY_EQUILIBRIUM_INITIALIZATION=1)
endif
else
$(error CHEMISTRY_EQUILIBRIUM_INITIALIZATION must be 0 or 1)
endif
ifeq ($(CHEMISTRY_NEWTON_LU_REUSE),1)
ifeq ($(CHEMISTRY_NEWTON_SYSTEM),FULL)
OPT1 += -DCHEMISTRY_NEWTON_LU_REUSE
endif
else ifneq ($(CHEMISTRY_NEWTON_LU_REUSE),0)
$(error CHEMISTRY_NEWTON_LU_REUSE must be 0 or 1)
endif
ifeq ($(CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE),1)
ifeq ($(CHEMISTRY_EQUILIBRIUM_INITIALIZATION)$(CHEMISTRY_NEWTON_SYSTEM),1FULL)
OPT1 += -DCHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE
endif
else ifneq ($(CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE),0)
$(error CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE must be 0 or 1)
endif
ifeq ($(USE_CHEMISTRY_PORTABLE),1)
ifneq ($(CHEMISTRY_UPDATE_CORE_LENGTH),)
ifneq ($(filter $(CHEMISTRY_UPDATE_CORE_LENGTH),0 1),$(CHEMISTRY_UPDATE_CORE_LENGTH))
$(error CHEMISTRY_UPDATE_CORE_LENGTH must be 0 or 1)
endif
OPT1 += -DCHEMISTRY_UPDATE_CORE_LENGTH=$(CHEMISTRY_UPDATE_CORE_LENGTH)
endif
else
ifeq ($(CHEMISTRY_UPDATE_CORE_LENGTH),1)
OPT1 += -DCHEMISTRY_UPDATE_CORE_LENGTH
endif
OPT1 += -DGAMMA_NEW
endif

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

ifeq ($(CHEMISTRY_WRITE_DIAGNOSTIC_TABLES),1)
OPT1 += -DCHEMISTRY_WRITE_DIAGNOSTIC_TABLES
else ifneq ($(CHEMISTRY_WRITE_DIAGNOSTIC_TABLES),0)
$(error CHEMISTRY_WRITE_DIAGNOSTIC_TABLES must be 0 or 1)
endif

#OPT1  += -DDENSITY_GRAD
#OPT1  +=  -DCHEMISTRY_JEANS_TEMPERATURE_MAX=1.e4

OPT1  +=  -DPHYS_COSMO
OPT1  +=  -DCHEMISTRY_CMB_HEATING

OPT1  +=  -DCHEMISTRY

OPT1  +=  -DCHEMISTRY_ENERGY_IMPLICIT
# Define this directly only when chemistry is implicit but energy is explicit.
#OPT1  +=  -DCHEMISTRY_IMPLICIT
#OPT1  +=  -DCHEMISTRY_HYBRID_IMPLICIT
ifeq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 += -DCHEMISTRY_ENERGY_BRACKET_FALLBACK
else
OPT1 += -DFLAG_BRENT
endif
OPT1  +=  -DCHEMISTRY_DIFF_ANALYTIC

#OPT1  +=  -DLIMIT_DT_IMPLICIT=1.e11
#OPT1  +=  -DCHEMISTRY_SPECIES_HELIUM
#OPT1  +=  -DCHEMISTRY_SPECIES_H2_PLUS

OPT1  +=  -DCHEMISTRY_SPECIES_DEUTERIUM
#OPT1  +=  -DCHEMISTRY_SPECIES_HD_PLUS
#OPT1  +=  -DCHEMISTRY_SPECIES_D_MINUS

OPT1  +=  -DCHEMISTRY_INVERSE_REACTIONS
OPT1  +=  -DCHEMISTRY_SAHA

# Physics defaults live in chemistry_build_config.h for the portable core.
ifneq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1  +=  -DTHREEBODY
OPT1  +=  -DMARTIN97
OPT1  +=  -DH2_COOLING_G15
#OPT1  +=  -DTBR_FLOW07
OPT1  +=  -DTBR_FORR13
#OPT1  +=  -DENZO_CHEM
#OPT1  +=  -DO05_CHEM
OPT1  +=  -DCONTINUUM_COOLING_MATSUKOBA
OPT1  +=  -DCASE_B
OPT1  +=  -DDUST_COOLING
OPT1  +=  -DH2_FORMATION_ON_DUST
OPT1  +=  -DH2_OPACITY
OPT1  +=  -DH2_OPACITY_FIT
OPT1  +=  -DCONT_OPACITY
endif
OPT1 += -DGA08
OPT1 += -DMETALLICITY=-4.
OPT1 += -DMETAL_COOLING_LW

ifeq ($(ENABLE_LW),1)
OPT1  +=  -DCHEMISTRY_LW_RADIATION
OPT1  +=  -DCHEMISTRY_BACKGROUND_LW_RADIATION
#OPT1  +=  -DTYPE_A
OPT1  +=  -DTYPE_B   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_3E4   #type of the spectrum of LW intensity
#OPT1  +=  -DTYPE_POW   #type of the spectrum of LW intensity
OPT1  +=  -DCHEMISTRY_LW_SELF_SHIELDING
#OPT1  +=  -DSHIELD_DENS=1.e9
endif

#OPT1  +=  -DIONIZING_PHOTON
#OPT1  +=  -DSUBGRID_CHAMPAGNE_FLOW

# The archived implementations retain their historical option names.
ifneq ($(USE_CHEMISTRY_PORTABLE),1)
ifneq (,$(filter -DCHEMISTRY_SPECIES_HELIUM,$(OPT1)))
OPT1 += -DHELIUM
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_DEUTERIUM,$(OPT1)))
OPT1 += -DHD
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_H2_PLUS,$(OPT1)))
OPT1 += -DH2_PLUS
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_HD_PLUS,$(OPT1)))
OPT1 += -DHD_PLUS
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_D_MINUS,$(OPT1)))
OPT1 += -DDMINUS
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_CARBON_OXYGEN,$(OPT1)))
OPT1 += -DCARBON_OXYGEN_CHEMISTRY
endif
ifneq (,$(filter -DCHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION,$(OPT1)))
OPT1 += -DCARBON_OXYGEN_HIGH_IONIZATION
endif
endif

ifneq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_DIFF_ANALYTIC%,-DCHEM_DIFF_ANALYTIC%,$(filter -DCHEMISTRY_DIFF_ANALYTIC -DCHEMISTRY_DIFF_ANALYTIC=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_DIFF_NUMERICAL%,-DCHEM_DIFF_NUMERICAL%,$(filter -DCHEMISTRY_DIFF_NUMERICAL -DCHEMISTRY_DIFF_NUMERICAL=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION%,-DCHEM_DISABLE_CARBON_OXYGEN_CONSERVATION%,$(filter -DCHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION -DCHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_ENERGY_BRACKET_FALLBACK%,-DCHEM_ENERGY_BRACKET_FALLBACK%,$(filter -DCHEMISTRY_ENERGY_BRACKET_FALLBACK -DCHEMISTRY_ENERGY_BRACKET_FALLBACK=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE%,-DCHEM_EQUILIBRIUM_SPECIES_CONVERGENCE%,$(filter -DCHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE -DCHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_EXPLICIT_SEQUENTIAL%,-DCHEM_EXPLICIT_SEQUENTIAL%,$(filter -DCHEMISTRY_EXPLICIT_SEQUENTIAL -DCHEMISTRY_EXPLICIT_SEQUENTIAL=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_EXPLICIT_SIMULTANEOUS%,-DCHEM_EXPLICIT_SIMULTANEOUS%,$(filter -DCHEMISTRY_EXPLICIT_SIMULTANEOUS -DCHEMISTRY_EXPLICIT_SIMULTANEOUS=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION%,-DCHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION%,$(filter -DCHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION -DCHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_INTERMEDIATE_EQUILIBRIUM%,-DCHEM_INTERMEDIATE_EQUILIBRIUM%,$(filter -DCHEMISTRY_INTERMEDIATE_EQUILIBRIUM -DCHEMISTRY_INTERMEDIATE_EQUILIBRIUM=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_NEWTON_LU_REUSE%,-DCHEM_NEWTON_LU_REUSE%,$(filter -DCHEMISTRY_NEWTON_LU_REUSE -DCHEMISTRY_NEWTON_LU_REUSE=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT%,-DCHEM_NEWTON_RESIDUAL_GATE_LIMIT%,$(filter -DCHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT -DCHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_SAHA_DIAGNOSTICS%,-DCHEM_SAHA_DIAGNOSTICS%,$(filter -DCHEMISTRY_SAHA_DIAGNOSTICS -DCHEMISTRY_SAHA_DIAGNOSTICS=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_CO_O2_DEPENDENT_DENSITY%,-DCHEM_CO_O2_DEPENDENT_DENSITY%,$(filter -DCHEMISTRY_CO_O2_DEPENDENT_DENSITY -DCHEMISTRY_CO_O2_DEPENDENT_DENSITY=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_ENERGY_IMPLICIT%,-DENERGY_IMPLICIT%,$(filter -DCHEMISTRY_ENERGY_IMPLICIT -DCHEMISTRY_ENERGY_IMPLICIT=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_HYBRID_IMPLICIT%,-DHYBRID_IMPLICIT%,$(filter -DCHEMISTRY_HYBRID_IMPLICIT -DCHEMISTRY_HYBRID_IMPLICIT=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_SAHA%,-DSAHA%,$(filter -DCHEMISTRY_SAHA -DCHEMISTRY_SAHA=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_INVERSE_REACTIONS%,-DINVERSE_REACTIONS%,$(filter -DCHEMISTRY_INVERSE_REACTIONS -DCHEMISTRY_INVERSE_REACTIONS=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_CMB_HEATING%,-DCMB_HEATING%,$(filter -DCHEMISTRY_CMB_HEATING -DCHEMISTRY_CMB_HEATING=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_UPDATE_CORE_LENGTH%,-DUPDATE_L_CORE_IN_SOLVER%,$(filter -DCHEMISTRY_UPDATE_CORE_LENGTH -DCHEMISTRY_UPDATE_CORE_LENGTH=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_JEANS_TEMPERATURE_MAX%,-DJEANS_LENGTH_TMAX%,$(filter -DCHEMISTRY_JEANS_TEMPERATURE_MAX -DCHEMISTRY_JEANS_TEMPERATURE_MAX=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_MAX_ITERATIONS%,-DMAX_ITER_CHEM%,$(filter -DCHEMISTRY_MAX_ITERATIONS -DCHEMISTRY_MAX_ITERATIONS=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_FIXED_TEMPERATURE%,-DFIX_TEMPERATURE%,$(filter -DCHEMISTRY_FIXED_TEMPERATURE -DCHEMISTRY_FIXED_TEMPERATURE=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_TEMPERATURE_MAX%,-DSET_TMAX%,$(filter -DCHEMISTRY_TEMPERATURE_MAX -DCHEMISTRY_TEMPERATURE_MAX=%,$(OPT1)))
endif

ifneq ($(USE_CHEMISTRY_PORTABLE),1)
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_LW_RADIATION%,-DLYMANWERNER%,$(filter -DCHEMISTRY_LW_RADIATION -DCHEMISTRY_LW_RADIATION=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_BACKGROUND_LW_RADIATION%,-DBACKGROUND_LW_RADIATION%,$(filter -DCHEMISTRY_BACKGROUND_LW_RADIATION -DCHEMISTRY_BACKGROUND_LW_RADIATION=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DCHEMISTRY_LW_SELF_SHIELDING%,-DSHIELD_LW%,$(filter -DCHEMISTRY_LW_SELF_SHIELDING -DCHEMISTRY_LW_SELF_SHIELDING=%,$(OPT1)))
OPT1 := $(OPT1) $(patsubst -DIONIZING_PHOTON%,-DIONIZATION_PHOTON%,$(filter -DIONIZING_PHOTON -DIONIZING_PHOTON=%,$(OPT1)))
endif

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
CXX    =  mpicxx
CXXFLAGS = $(filter-out -std=%,$(CFLAGS)) -std=c++17 -pthread
LINKER = $(CC)
LINKFLAGS = $(CFLAGS)

ifneq (,$(findstring -DCHEMISTRY_OLD,$(OPT1)))
include legacy/chemistry_old/standalone.mk
else ifeq ($(USE_CHEMISTRY_PORTABLE),1)
CHEMISTRY_DIR = src/chemistry
include src/standalone/standalone.mk
LINKER = $(CXX)
LINKFLAGS = $(CXXFLAGS)
CHEMISTRY_OBJS = $(patsubst %.cc,%.o,$(CHEMISTRY_PORTABLE_CORE_SOURCES:.c=.o)) \
                 $(patsubst %.cc,%.o,$(CHEMISTRY_STANDALONE_ADAPTER_SOURCES:.c=.o))
CHEMISTRY_SUPPORT_OBJS =
CHEMISTRY_INCL = $(wildcard $(CHEMISTRY_DIR)/*.h $(CHEMISTRY_DIR)/*/*.h) \
                 $(wildcard $(CHEMISTRY_DIR)/include/chemistry/*.h) \
                 $(CHEMISTRY_DIR)/tables/chemistry_tables.h \
                 $(CHEMISTRY_DIR)/chemistry_core.h \
                 $(CHEMISTRY_DIR)/chemistry_build_config.h \
                 $(CHEMISTRY_DIR)/chemistry_types.h \
                 $(wildcard $(CHEMISTRY_DIR)/*.hpp) \
                 $(wildcard $(CHEMISTRY_DIR)/solver/*.hpp) \
                 $(CHEMISTRY_DIR)/solver/chemistry_newton_internal.h \
                 $(CHEMISTRY_DIR)/chemistry_core.h \
                 $(CHEMISTRY_DIR)/solver/chemistry_newton_diagnostics.h \
                 src/standalone/chemistry_standalone_adapter.h \
                 src/standalone/chemistry_standalone_types.h
CHEMISTRY_INCLUDE_FLAGS = -I$(CHEMISTRY_DIR) \
                          -I$(CHEMISTRY_DIR)/include \
                          -I$(CHEMISTRY_DIR)/cooling \
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
	$(LINKER) $(LINKFLAGS) $(LDFLAGS) $(OBJDIR_OBJS) $(LIBS) -o $(EXEC)

$(BUILDDIR)/%.o: %.c $(INCL)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILDDIR)/%.o: %.cc $(INCL)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: clean check-chemistry-initialization check-energy-replay check-chemistry-cpp

# Replays captured AREPO energy-step failures (tests/replay) through the
# standalone build; builds its own binaries in build-replay*.
check-energy-replay:
	sh tests/replay/run_energy_step_replay.sh

$(BUILDDIR)/initialization_regression: $(BUILDDIR)/src/chemistry/tests/initialization_regression.o $(filter-out $(BUILDDIR)/src/standalone/main.o,$(OBJDIR_OBJS))
	$(LINKER) $(LINKFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

check-chemistry-initialization: $(BUILDDIR)/initialization_regression
	$(abspath $(BUILDDIR))/initialization_regression src/chemistry

$(BUILDDIR)/cpp_api_link: $(BUILDDIR)/src/chemistry/tests/cpp_api_link.o $(BUILDDIR)/src/chemistry/tests/c_status_link.o $(filter-out $(BUILDDIR)/src/standalone/main.o,$(OBJDIR_OBJS))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

check-chemistry-cpp: $(BUILDDIR)/cpp_api_link
	$(abspath $(BUILDDIR))/cpp_api_link src/chemistry

.PHONY: check-chemistry-model
$(BUILDDIR)/model_lifecycle: $(BUILDDIR)/src/chemistry/tests/model_lifecycle.o $(filter-out $(BUILDDIR)/src/standalone/main.o,$(OBJDIR_OBJS))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

check-chemistry-model: $(BUILDDIR)/model_lifecycle
	sh src/chemistry/tests/run_model_lifecycle.sh $(abspath $(BUILDDIR))/model_lifecycle src/chemistry

.PHONY: check-chemistry-fixed-rate
$(BUILDDIR)/fixed_rate_problem: $(BUILDDIR)/src/chemistry/tests/fixed_rate_problem.o $(filter-out $(BUILDDIR)/src/standalone/main.o,$(OBJDIR_OBJS))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

check-chemistry-fixed-rate: $(BUILDDIR)/fixed_rate_problem
	$(abspath $(BUILDDIR))/fixed_rate_problem src/chemistry

clean:
	rm -rf $(BUILDDIR) $(EXEC)

.PHONY: check-chemistry-numerics
$(BUILDDIR)/numerics_dense: $(BUILDDIR)/src/chemistry/tests/numerics_dense.o $(BUILDDIR)/src/chemistry/solver/linear_algebra.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

check-chemistry-numerics: $(BUILDDIR)/numerics_dense
	$(abspath $(BUILDDIR))/numerics_dense

-include $(OBJDIR_OBJS:.o=.d) $(wildcard $(BUILDDIR)/src/chemistry/tests/*.d)

.PHONY: check-build-options
check-build-options:
	python3 tests/test_build_options.py
