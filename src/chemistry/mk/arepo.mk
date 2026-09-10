# Include this fragment from AREPO's makefiles/modules.make after chemistry/
# and radiation/ have been copied below AREPO's src/ directory.
CHEMISTRY_ROOT = chemistry
include $(SRC_DIR)/chemistry/mk/core.mk

CHEMISTRY_AREPO_ADAPTER_SOURCES += \
  radiation/radiation_rates.c \
  chemistry/adapters/arepo/chemistry_arepo_adapter.c \
  chemistry/adapters/arepo/chemistry_arepo_host.c \
  chemistry/interface/chem_balance.c

CHEMISTRY_AREPO_INCLUDE_FLAGS += \
  -Ichemistry/include \
  -Ichemistry \
  -Ichemistry/common \
  -Ichemistry/cooling \
  -Ichemistry/driver \
  -Ichemistry/initialization \
  -Ichemistry/reaction_network \
  -Ichemistry/solver \
  -Ichemistry/thermodynamics \
  -Ichemistry/tables \
  -Ichemistry/adapters/arepo

CHEMISTRY_PORTABLE_OBJECTS := \
  $(CHEMISTRY_PORTABLE_CORE_SOURCES:.c=.o) \
  $(CHEMISTRY_AREPO_ADAPTER_SOURCES:.c=.o)

CHEMISTRY_PORTABLE_BUILD_OBJECTS := \
  $(addprefix $(BUILD_DIR)/,$(CHEMISTRY_PORTABLE_OBJECTS))

OBJS += $(CHEMISTRY_PORTABLE_OBJECTS)
SUBDIRS += $(sort $(dir $(CHEMISTRY_PORTABLE_OBJECTS)))
INCL += \
  chemistry/chemistry_core.h \
  chemistry/chemistry_types.h \
  radiation/radiation_rates.h \
  chemistry/adapters/arepo/chemistry_arepo_adapter.h \
  chemistry/adapters/arepo/chemistry_arepo_host.h

# Keep the portable core's configuration and private include paths target-local.
# In particular, do not force AREPO's generated configuration header into every
# unrelated translation unit merely because portable chemistry is enabled.
$(CHEMISTRY_PORTABLE_BUILD_OBJECTS): CFLAGS += \
  -include $(BUILD_DIR)/arepoconfig.h \
  -I$(SRC_DIR) \
  $(addprefix -I$(SRC_DIR)/,$(patsubst -I%,%,$(CHEMISTRY_AREPO_INCLUDE_FLAGS)))
