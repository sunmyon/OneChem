# Host-independent chemistry implementation. Host build fragments include this
# file and append exactly one adapter source list.
CHEMISTRY_ROOT ?= src/chemistry
CHEMISTRY_PORTABLE_CORE_SOURCES += \
  $(CHEMISTRY_ROOT)/chemistry_driver.cc \
  $(CHEMISTRY_ROOT)/chemistry_model.cc \
  $(CHEMISTRY_ROOT)/solver/fixed_rate_solver.cc \
  $(CHEMISTRY_ROOT)/solver/linear_algebra.cc \
  $(CHEMISTRY_ROOT)/species.cc \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_control.cc \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_diagnostics.cc \
  $(CHEMISTRY_ROOT)/solver/chemistry_newton_reset.cc \
  $(CHEMISTRY_ROOT)/solver/chemistry_solver.cc \
  $(CHEMISTRY_ROOT)/solver/saha.cc \
  $(CHEMISTRY_ROOT)/chemistry_profile.cc \
  $(CHEMISTRY_ROOT)/reaction_network/forward_rates.cc \
  $(CHEMISTRY_ROOT)/reaction_network/reaction_rates.cc \
  $(CHEMISTRY_ROOT)/reaction_network/reaction_network.cc \
  $(CHEMISTRY_ROOT)/reaction_network/shielding.cc \
  $(CHEMISTRY_ROOT)/cooling/H2cooling.cc \
  $(CHEMISTRY_ROOT)/cooling/chemistry_cooling.cc \
  $(CHEMISTRY_ROOT)/cooling/cooling_func.cc \
  $(CHEMISTRY_ROOT)/cooling/molecular_cooling.cc \
  $(CHEMISTRY_ROOT)/cooling/dust_cooling.cc \
  $(CHEMISTRY_ROOT)/thermodynamics/partition.cc \
  $(CHEMISTRY_ROOT)/thermodynamics/thermodynamics.cc
