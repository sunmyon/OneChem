# Include before NG prefixes its native object paths with BUILD_DIR.
include $(CHEMISTRY_ROOT)/mk/core.mk
CHEMISTRY_SOURCES = $(CHEMISTRY_PORTABLE_CORE_SOURCES) $(wildcard $(CHEMISTRY_ROOT)/adapters/arepo_ng/*.cc)
CHEMISTRY_OBJECTS = $(patsubst $(CHEMISTRY_ROOT)/%.cc,$(BUILD_DIR)/chemistry/%.o,$(CHEMISTRY_SOURCES))
OPT += -I$(CHEMISTRY_ROOT)/.. $(addprefix -I$(CHEMISTRY_ROOT)/, . include cooling reaction_network solver thermodynamics tables adapters/arepo_ng)
# No default: plain reference LAPACK/BLAS is often not what's wanted (MKL, qMKL,
# etc. are common instead), and every other third-party lib here (GSL_LIBS,
# HDF5_LIBS, FFTW_LIBS) is likewise supplied by the caller, not defaulted in this
# Makefile. Pass CHEMISTRY_LAPACK_LIBS explicitly (command line or a systype file).
CHEMISTRY_LAPACK_LIBS ?=
LIBS += $(CHEMISTRY_LAPACK_LIBS)
CHEMISTRY_HEADER_DIRS = $(addprefix $(CHEMISTRY_ROOT)/, . include/chemistry cooling reaction_network solver thermodynamics tables adapters/arepo_ng)
CHEMISTRY_HEADERS = $(wildcard $(addsuffix /*.h,$(CHEMISTRY_HEADER_DIRS)) $(addsuffix /*.hpp,$(CHEMISTRY_HEADER_DIRS)))
CHEMISTRY_CONFIG_CHECK = $(BUILD_DIR)/chemistry.options.check
# NG-side headers the adapter includes (e.g. logs/timer.h for TIMER_START/
# TIMER_CREATE) live outside CHEMISTRY_ROOT and so aren't covered by
# $(CHEMISTRY_HEADERS); without this, an incremental build can relink stale
# chemistry .o files against a changed logs::timers enum layout. This is included
# before the Makefile's own $(INCL) := $(addprefix $(SRC_DIR)/,$(INCL)) ... rewrite,
# so apply the same prefix here rather than relying on the (still-unprefixed,
# nonexistent-as-named) bare $(INCL).
$(BUILD_DIR)/chemistry/%.o: $(CHEMISTRY_ROOT)/%.cc $(BUILD_DIR)/arepoconfig.h $(CHEMISTRY_HEADERS) $(addprefix $(SRC_DIR)/,$(INCL)) $(AREPO_MAKEFILES)
	@mkdir -p $(dir $@)
	$(CPP) $(CFLAGS) -include $(BUILD_DIR)/arepoconfig.h -c $< -o $@
# Core owns its internal option checks; expose used options to NG's Config check.
$(CHEMISTRY_CONFIG_CHECK): $(CHEMISTRY_SOURCES) $(CHEMISTRY_HEADERS) buildsystem/check.py
	@$(PYTHON) -c 'import sys,io; sys.path.insert(0,"buildsystem"); from check import filter_code,write; used=set(); [used.update(filter_code(io.StringIO("#include <arepoconfig.h>\n"+open(p).read()))) for p in sys.argv[2:]]; write(used,sys.argv[1])' $@ $(CHEMISTRY_SOURCES) $(CHEMISTRY_HEADERS)
