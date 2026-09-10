CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic

.PHONY: check background-bridge-check

background-bridge-check:
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEMISTRY_SC -DHD \
	  -I../chemistry -I../chemistry/common \
	  -I../chemistry/include -I../.. -I.. -I. \
	  ../radiation/radiation_rates.c \
	  chemistry_standalone_background.c \
	  tests/standalone_background_bridge_test.c -lm \
	  -o tests/standalone_background_bridge_test
	./tests/standalone_background_bridge_test
	$(RM) tests/standalone_background_bridge_test

check: background-bridge-check
