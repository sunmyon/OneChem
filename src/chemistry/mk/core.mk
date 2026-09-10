# Host-independent chemistry implementation. Host build fragments include this
# file and append exactly one adapter source list.
CHEMISTRY_ROOT ?= src/chemistry
CHEMISTRY_PORTABLE_CORE_SOURCES += \
  $(CHEMISTRY_ROOT)/interface/chemistry_api.c \
  $(CHEMISTRY_ROOT)/driver/chemistry_driver.c \
  $(CHEMISTRY_ROOT)/species.c \
  $(CHEMISTRY_ROOT)/initialization/chemistry_init.c \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_control.c \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_diagnostics.c \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_reset.c \
  $(CHEMISTRY_ROOT)/solver/chemistry_solver.c \
  $(CHEMISTRY_ROOT)/solver/saha.c \
  $(CHEMISTRY_ROOT)/common/chemistry_models.c \
  $(CHEMISTRY_ROOT)/common/chemistry_profile.c \
  $(CHEMISTRY_ROOT)/common/chemistry_status.c \
  $(CHEMISTRY_ROOT)/reaction_network/forward_rates.c \
  $(CHEMISTRY_ROOT)/reaction_network/reaction_rates.c \
  $(CHEMISTRY_ROOT)/reaction_network/reaction_network.c \
  $(CHEMISTRY_ROOT)/reaction_network/shielding.c \
  $(CHEMISTRY_ROOT)/cooling/H2cooling.c \
  $(CHEMISTRY_ROOT)/cooling/chemistry_cooling.c \
  $(CHEMISTRY_ROOT)/cooling/cooling_func.c \
  $(CHEMISTRY_ROOT)/cooling/molecular_cooling.c \
  $(CHEMISTRY_ROOT)/cooling/dust_cooling.c \
  $(CHEMISTRY_ROOT)/thermodynamics/partition.c \
  $(CHEMISTRY_ROOT)/thermodynamics/thermodynamics.c
