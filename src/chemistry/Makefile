CC ?= cc
CXX ?= c++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -Werror -pedantic
CPPFLAGS ?=
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic

PUBLIC_INCLUDE := -Iinclude


.PHONY: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check check
api-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PUBLIC_INCLUDE) -fsyntax-only tests/public_api_compile.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PUBLIC_INCLUDE) -I. -fsyntax-only tests/cpp_api_link.cc

arepo-adapter-check:
	sh tests/check_arepo_adapter_contract.sh

arepo-adapter-compile-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iinclude -Iadapters/arepo -fsyntax-only tests/arepo_public_api_compile.c
	mpicxx $(CPPFLAGS) -std=c++17 -DCHEMISTRY_SPECIES_DEUTERIUM -DCHEMISTRY_LW_RADIATION \
	  -DCHEMISTRY_BACKGROUND_LW_RADIATION -DMETALLICITY=-4. \
	  -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.. \
	  -fsyntax-only adapters/arepo/chem_balance.cc
	mpicxx $(CPPFLAGS) -std=c++17 -DCHEMISTRY_SPECIES_DEUTERIUM -DCHEM_BALANCE \
	  -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.. \
	  -fsyntax-only adapters/arepo/chem_balance.cc

arepo-background-check:
	@set -e; for metal in "" "-DCONSTANT_METALLICITY" "-DMETAL_ADVECTION_SC" "-DMETAL_ADVECTION_SC -DCONSTANT_METALLICITY"; do for mode in local imported; do \
	  sources="adapters/arepo/chem_balance.cc"; options=""; \
	  if [ "$$mode" = imported ]; then sources=""; options="-DCHEM_BALANCE"; fi; \
	  CC=mpicc CXX=mpicxx CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-DCHEMISTRY_SPECIES_DEUTERIUM -DCHEMISTRY_LW_RADIATION -DCHEMISTRY_BACKGROUND_LW_RADIATION -DMETALLICITY=-4. $$metal $$options -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh \
	      ../radiation/radiation_rates.cc $$sources tests/arepo_background_bridge_test.cc; \
	done; done

core-dependency-check:
	CC="$(CC)" CXX="$(CXX)" sh tests/check_core_dependencies.sh

no-lazy-init-check:
	sh tests/check_no_lazy_initialization.sh

no-process-termination-check:
	sh tests/check_no_process_termination.sh

global-state-check:
	sh tests/check_global_state_migration.sh

host-boundary-check:
	sh tests/check_host_boundary.sh

radiation-check:
	CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" CHEM_TEST_FLAGS="-I.." CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh ../radiation/radiation_rates.cc tests/radiation_rates_test.c

.PHONY: high-ionization-option-check metal-ionization-check photoionization-check co-photodissociation-check carbon-oxygen-normalization-check
high-ionization-option-check:
	@set -e; for options in "" "-DCHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION -DEXPECT_HIGH_IONIZATION"; do \
	  CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-Wall -Wextra -Werror -pedantic -DCHEMISTRY_SPECIES_CARBON_OXYGEN -DMETALLICITY=-4. $$options -I." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh species.cc tests/high_ionization_option_test.cc; \
	done

metal-ionization-check:
	CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	  CHEM_TEST_FLAGS="-Wall -Wextra -Wno-unused-function -DCHEMISTRY_SPECIES_CARBON_OXYGEN -DCHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION -DRYTR_CARBON_OXYGEN_IONIZATION -DRYTR_NCO -DMETALLICITY=-4. -I. -Ireaction_network -Itables" \
	  CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh \
	  species.cc reaction_network/forward_rates.cc reaction_network/reaction_network.cc tests/metal_ionization_network_test.cc

photoionization-check:
	CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" CHEM_TEST_FLAGS="-I../radiation" CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh ../radiation/photoionization_bins.cc tests/photoionization_bins_test.c

co-photodissociation-check:
	CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" CHEM_TEST_FLAGS="-I../radiation" CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh ../radiation/co_photodissociation.cc tests/co_photodissociation_test.c

carbon-oxygen-normalization-check:
	@set -e; for options in "" "-DCHEMISTRY_SPECIES_CARBON_OXYGEN" \
	  "-DCHEMISTRY_SPECIES_CARBON_OXYGEN -DCHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION"; do \
	  CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-Wall -Wextra -Werror -pedantic -DMETALLICITY -DCHEMISTRY_SPECIES_HELIUM -DCHEMISTRY_SPECIES_DEUTERIUM -DCHEMISTRY_SPECIES_HD_PLUS -DCHEMISTRY_SPECIES_D_MINUS -DCHEMISTRY_SPECIES_H2_PLUS $$options -Iinclude -I." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh species.cc tests/carbon_oxygen_normalization_test.cc; \
	done

check: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check

# Reuse precisely the core objects/configuration linked into the NG executable.
# Build NG first; this target is intentionally separate from portable checks.
NG_BUILD_DIR ?= /tmp/arepo_ng_chemistry_build
NG_LINALG_LIBS ?= -llapack -lblas
NG_GSL_LIBS ?= -lgsl -lgslcblas
NG_CHEMISTRY_OBJECTS = $(filter-out %/chemistry_arepo_ng_module.o,$(wildcard $(NG_BUILD_DIR)/chemistry/*.o $(NG_BUILD_DIR)/chemistry/*/*.o $(NG_BUILD_DIR)/chemistry/adapters/arepo_ng/*.o))
.PHONY: arepo-ng-mpi-check
arepo-ng-mpi-check:
	@test -n "$(NG_CHEMISTRY_OBJECTS)" || { echo "Build NG with CHEMISTRY first and set NG_BUILD_DIR"; exit 1; }
	mpicxx -std=c++20 -O2 -I. -include $(NG_BUILD_DIR)/arepoconfig.h tests/arepo_ng_mpi_test.cc $(NG_CHEMISTRY_OBJECTS) $(NG_GSL_LIBS) $(NG_LINALG_LIBS) -o $(NG_BUILD_DIR)/chemistry_ng_mpi_test
	@set -e; for ranks in 1 2 3; do mpiexec -n $$ranks $(NG_BUILD_DIR)/chemistry_ng_mpi_test "$(CURDIR)"; done
