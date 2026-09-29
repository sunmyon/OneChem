CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -std=c++17 -Wall -Wextra -Werror -pedantic

.PHONY: check background-bridge-check

background-bridge-check:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCHEMISTRY_SPECIES_DEUTERIUM \
	  -I../chemistry -I../chemistry/common \
	  -I../chemistry/include -I../.. -I.. -I. \
	  ../radiation/radiation_rates.cc \
	  chemistry_standalone_background.cc \
	  tests/standalone_background_bridge_test.cc -lm \
	  -o tests/standalone_background_bridge_test
	./tests/standalone_background_bridge_test
	$(RM) tests/standalone_background_bridge_test

check: background-bridge-check
