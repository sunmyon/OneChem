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

check: api-check arepo-adapter-check arepo-adapter-compile-check arepo-background-check core-dependency-check global-state-check host-boundary-check no-lazy-init-check no-process-termination-check radiation-check
