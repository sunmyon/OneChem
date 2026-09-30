# Include before NG prefixes its native object paths with BUILD_DIR.
include $(CHEMISTRY_ROOT)/mk/core.mk
CHEMISTRY_SOURCES = $(CHEMISTRY_PORTABLE_CORE_SOURCES) $(wildcard $(CHEMISTRY_ROOT)/adapters/arepo_ng/*.cc)
CHEMISTRY_OBJECTS = $(patsubst $(CHEMISTRY_ROOT)/%.cc,$(BUILD_DIR)/chemistry/%.o,$(CHEMISTRY_SOURCES))
OPT += -I$(CHEMISTRY_ROOT)/.. $(addprefix -I$(CHEMISTRY_ROOT)/, . include cooling reaction_network solver thermodynamics tables adapters/arepo_ng)
CHEMISTRY_LAPACK_LIBS ?= -llapack -lblas
LIBS += $(CHEMISTRY_LAPACK_LIBS)
CHEMISTRY_HEADER_DIRS = $(addprefix $(CHEMISTRY_ROOT)/, . include/chemistry cooling reaction_network solver thermodynamics tables adapters/arepo_ng)
CHEMISTRY_HEADERS = $(wildcard $(addsuffix /*.h,$(CHEMISTRY_HEADER_DIRS)) $(addsuffix /*.hpp,$(CHEMISTRY_HEADER_DIRS)))
CHEMISTRY_CONFIG_CHECK = $(BUILD_DIR)/chemistry.options.check
$(BUILD_DIR)/chemistry/%.o: $(CHEMISTRY_ROOT)/%.cc $(BUILD_DIR)/arepoconfig.h $(CHEMISTRY_HEADERS)
	@mkdir -p $(dir $@)
	$(CPP) $(CFLAGS) -include $(BUILD_DIR)/arepoconfig.h -c $< -o $@
# Core owns its internal option checks; expose used options to NG's Config check.
$(CHEMISTRY_CONFIG_CHECK): $(CHEMISTRY_SOURCES) $(CHEMISTRY_HEADERS) buildsystem/check.py
	@$(PYTHON) -c 'import sys,io; sys.path.insert(0,"buildsystem"); from check import filter_code,write; used=set(); [used.update(filter_code(io.StringIO("#include <arepoconfig.h>\n"+open(p).read()))) for p in sys.argv[2:]]; write(used,sys.argv[1])' $@ $(CHEMISTRY_SOURCES) $(CHEMISTRY_HEADERS)
