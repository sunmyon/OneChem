CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic

PUBLIC_INCLUDE := -Iinclude

.PHONY: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check check
api-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(PUBLIC_INCLUDE) -fsyntax-only tests/public_api_compile.c

arepo-adapter-check:
	sh tests/check_arepo_adapter_contract.sh

arepo-adapter-compile-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DHD -DMETALLICITY=-4. \
	  -Iinclude -I. -Icommon -Iadapters/arepo \
	  -fsyntax-only tests/arepo_public_api_compile.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DHD -DMETALLICITY=-4. \
	  -Iinclude -I. -Icommon -Icooling -Idriver -Iinitialization -Ireaction_network \
	  -Isolver -Ithermodynamics -Itables -Iadapters/arepo \
	  -fsyntax-only adapters/arepo/chemistry_arepo_adapter.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DHD -DLYMANWERNER \
	  -DBACKGROUND_LW_RADIATION -DMETALLICITY=-4. \
	  -Itests/arepo_host_mock -Iinclude -I. -Icommon -Iadapters/arepo -I.. \
	  -fsyntax-only adapters/arepo/chemistry_arepo_host.c

arepo-background-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DHD -DLYMANWERNER \
	  -DBACKGROUND_LW_RADIATION -DTYPE_B -DMETALLICITY=-4. \
	  -Itests/arepo_host_mock -Iinclude -I. -Icommon -Iadapters/arepo -I.. \
	  ../radiation/radiation_rates.c \
	  adapters/arepo/chemistry_arepo_host.c \
	  tests/arepo_background_bridge_test.c -lm \
	  -o tests/arepo_background_bridge_test
	./tests/arepo_background_bridge_test
	$(RM) tests/arepo_background_bridge_test

core-dependency-check:
	CC="$(CC)" sh tests/check_core_dependencies.sh

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
	  $(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DMETALLICITY=-4. $$options \
	    -I. -Icommon species.c tests/high_ionization_option_test.c -lm -o tests/high_ionization_option_test; \
	  ./tests/high_ionization_option_test; \
	done
	$(RM) tests/high_ionization_option_test

metal-ionization-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -Wno-unused-function -DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY \
	  -DCARBON_OXYGEN_HIGH_IONIZATION -DRYTR_CARBON_OXYGEN_IONIZATION -DRYTR_NCO -DMETALLICITY=-4. \
	  -I. -Icommon -Ireaction_network -Itables \
	  species.c reaction_network/forward_rates.c reaction_network/reaction_network.c \
	  tests/metal_ionization_network_test.c -lm -o tests/metal_ionization_network_test
	./tests/metal_ionization_network_test
	$(RM) tests/metal_ionization_network_test

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
	  $(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DMETALLICITY -DHELIUM -DHD -DHD_PLUS -DDMINUS -DH2_PLUS $$options \
	    -Iinclude -I. -Icommon species.c interface/chemistry_normalize.c tests/carbon_oxygen_normalization_test.c \
	    -lm -o tests/carbon_oxygen_normalization_test; \
	  ./tests/carbon_oxygen_normalization_test; \
	done
	$(RM) tests/carbon_oxygen_normalization_test

check: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check
