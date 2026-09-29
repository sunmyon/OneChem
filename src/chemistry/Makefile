CC ?= cc
CXX ?= c++
CXXFLAGS ?= -std=c++17 -Wall -Wextra -Werror -pedantic
CPPFLAGS ?=
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic

PUBLIC_INCLUDE := -Iinclude


.PHONY: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check check
api-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PUBLIC_INCLUDE) -fsyntax-only tests/public_api_compile.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PUBLIC_INCLUDE) -I. -DCHEMISTRY_SC -fsyntax-only tests/cpp_api_link.cc

arepo-adapter-check:
	sh tests/check_arepo_adapter_contract.sh

arepo-adapter-compile-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -Iinclude -Iadapters/arepo -fsyntax-only tests/arepo_public_api_compile.c
	mpicxx $(CPPFLAGS) -std=c++17 -DCHEMISTRY_SC -DHD -DLYMANWERNER \
	  -DBACKGROUND_LW_RADIATION -DMETALLICITY=-4. \
	  -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.. \
	  -fsyntax-only adapters/arepo/chem_balance.cc
	mpicxx $(CPPFLAGS) -std=c++17 -DCHEMISTRY_SC -DHD -DCHEM_BALANCE \
	  -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.. \
	  -fsyntax-only adapters/arepo/chem_balance.cc

arepo-background-check:
	@set -e; for mode in local imported; do \
	  sources="adapters/arepo/chem_balance.cc"; options=""; \
	  if [ "$$mode" = imported ]; then sources=""; options="-DCHEM_BALANCE"; fi; \
	  CC=mpicc CXX=mpicxx CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-DCHEMISTRY_SC -DHD -DLYMANWERNER -DBACKGROUND_LW_RADIATION -DTYPE_B -DMETALLICITY=-4. $$options -Itests/arepo_host_mock -Iinclude -I. -Ithermodynamics -Icooling -Iadapters/arepo -I.." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh \
	      ../radiation/radiation_rates.c $$sources tests/arepo_background_bridge_test.cc; \
	done

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
	$(CC) $(CPPFLAGS) $(CFLAGS) -I.. ../radiation/radiation_rates.c tests/radiation_rates_test.c -lm -o tests/radiation_rates_test
	./tests/radiation_rates_test
	$(RM) tests/radiation_rates_test

.PHONY: high-ionization-option-check metal-ionization-check photoionization-check co-photodissociation-check carbon-oxygen-normalization-check
high-ionization-option-check:
	@set -e; for options in "" "-DCARBON_OXYGEN_HIGH_IONIZATION -DEXPECT_HIGH_IONIZATION"; do \
	  CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-Wall -Wextra -Werror -pedantic -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DMETALLICITY=-4. $$options -I." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh species.cc tests/high_ionization_option_test.cc; \
	done

metal-ionization-check:
	CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	  CHEM_TEST_FLAGS="-Wall -Wextra -Wno-unused-function -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DCARBON_OXYGEN_HIGH_IONIZATION -DRYTR_CARBON_OXYGEN_IONIZATION -DRYTR_NCO -DMETALLICITY=-4. -I. -Ireaction_network -Itables" \
	  CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh \
	  species.cc reaction_network/forward_rates.cc reaction_network/reaction_network.cc tests/metal_ionization_network_test.cc

photoionization-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -I../radiation \
	  ../radiation/photoionization_bins.c tests/photoionization_bins_test.c -lm \
	  -o tests/photoionization_bins_test
	./tests/photoionization_bins_test
	$(RM) tests/photoionization_bins_test

co-photodissociation-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -I../radiation \
	  ../radiation/co_photodissociation.c tests/co_photodissociation_test.c -lm \
	  -o tests/co_photodissociation_test
	./tests/co_photodissociation_test
	$(RM) tests/co_photodissociation_test

carbon-oxygen-normalization-check:
	@set -e; for options in "" "-DCARBON_OXYGEN_CHEMISTRY" \
	  "-DCARBON_OXYGEN_CHEMISTRY -DCARBON_OXYGEN_HIGH_IONIZATION"; do \
	  CC="$(CC)" CXX="$(CXX)" CPPFLAGS="$(CPPFLAGS)" \
	    CHEM_TEST_FLAGS="-Wall -Wextra -Werror -pedantic -DCHEMISTRY_SC -DMETALLICITY -DHELIUM -DHD -DHD_PLUS -DDMINUS -DH2_PLUS $$options -Iinclude -I." \
	    CHEM_TEST_LIBS="-lm" sh tests/run_mixed_core_test.sh species.cc tests/carbon_oxygen_normalization_test.cc; \
	done

check: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check
