# Include this fragment from AREPO's makefiles/modules.make after chemistry/
# and radiation/ have been copied below AREPO's src/ directory.
CHEMISTRY_ROOT = chemistry
include $(SRC_DIR)/chemistry/mk/core.mk

CHEMISTRY_AREPO_ADAPTER_SOURCES += \
  radiation/radiation_rates.c \
  chemistry/adapters/arepo/chem_balance.cc

CHEMISTRY_AREPO_INCLUDE_FLAGS += \
  -Ichemistry/include \
  -Ichemistry \
  -Ichemistry/cooling \
  -Ichemistry/reaction_network \
  -Ichemistry/solver \
  -Ichemistry/thermodynamics \
  -Ichemistry/tables \
  -Ichemistry/adapters/arepo

CHEMISTRY_PORTABLE_OBJECTS := \
  $(patsubst %.cc,%.o,$(CHEMISTRY_PORTABLE_CORE_SOURCES:.c=.o)) \
  $(patsubst %.cc,%.o,$(CHEMISTRY_AREPO_ADAPTER_SOURCES:.c=.o))

CHEMISTRY_PORTABLE_BUILD_OBJECTS := \
  $(addprefix $(BUILD_DIR)/,$(CHEMISTRY_PORTABLE_OBJECTS))

OBJS += $(CHEMISTRY_PORTABLE_OBJECTS)
SUBDIRS += $(sort $(dir $(CHEMISTRY_PORTABLE_OBJECTS)))
INCL += \
  chemistry/tables/chemistry_tables.h \
  chemistry/chemistry_core.h \
  chemistry/chemistry_types.h \
  chemistry/chemistry_model.hpp \
  chemistry/solver/linear_algebra.hpp \
  chemistry/solver/fixed_rate_solver.hpp \
  chemistry/solver/chemistry_newton_internal.h \
  chemistry/chemistry_core.h \
  radiation/radiation_rates.h \
  chemistry/adapters/arepo/chemistry_arepo.h

# Keep the portable core's configuration and private include paths target-local.
# In particular, do not force AREPO's generated configuration header into every
# unrelated translation unit merely because portable chemistry is enabled.
$(CHEMISTRY_PORTABLE_BUILD_OBJECTS): CFLAGS += \
  -include $(BUILD_DIR)/arepoconfig.h \
  -I$(SRC_DIR) \
  $(addprefix -I$(SRC_DIR)/,$(patsubst -I%,%,$(CHEMISTRY_AREPO_INCLUDE_FLAGS)))

# Keep the host's C compiler and flags intact. Some AREPO systypes do not
# define CPPC; the portable module can select its matching MPI C++ wrapper.
CHEMISTRY_CXX ?= $(if $(strip $(CPPC)),$(CPPC),mpicxx)
CHEMISTRY_PORTABLE_CPP_OBJECTS := $(addprefix $(BUILD_DIR)/,$(patsubst %.cc,%.o,$(filter %.cc,$(CHEMISTRY_PORTABLE_CORE_SOURCES))))
# A pattern rule must not replace the host Makefile's default goal.
$(BUILD_DIR)/$(CHEMISTRY_ROOT)/%.o: $(SRC_DIR)/$(CHEMISTRY_ROOT)/%.cc
	$(CHEMISTRY_CXX) $(filter-out -std=%,$(CFLAGS)) -std=c++17 -c $< -o $@
