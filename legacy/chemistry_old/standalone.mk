# Standalone build integration for the detachable historical implementation.
CHEMISTRY_DIR = legacy/chemistry_old

CHEMISTRY_OBJS = $(CHEMISTRY_DIR)/chemistry_200601.o \
                 $(CHEMISTRY_DIR)/cooling_func_200601.o \
                 $(CHEMISTRY_DIR)/stellar_spectrum.o

CHEMISTRY_SUPPORT_OBJS = $(CHEMISTRY_DIR)/H2cooling.o \
                         $(CHEMISTRY_DIR)/dust_cooling.o \
                         $(CHEMISTRY_DIR)/saha.o \
                         $(CHEMISTRY_DIR)/partition.o \
                         $(CHEMISTRY_DIR)/standalone_compat.o

CHEMISTRY_INCL = $(CHEMISTRY_DIR)/chemistry.h \
                 $(CHEMISTRY_DIR)/chemistry_old_api.h \
                 $(CHEMISTRY_DIR)/standalone_compat.h

CHEMISTRY_INCLUDE_FLAGS = -I$(CHEMISTRY_DIR)
