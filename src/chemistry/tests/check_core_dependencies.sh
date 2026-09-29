#!/bin/sh
set -eu

cc_command=${CC:-cc}
include_flags="-Iinclude -I.  -Icooling   -Ireaction_network -Isolver -Ithermodynamics -Itables"
feature_flags="-DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DHD -DMETALLICITY=-4."

modules="
chemistry_driver.cc
species.cc
chemistry_profile.cc
chemistry_model.cc
solver/fixed_rate_solver.cc
solver/linear_algebra.cc
cooling/H2cooling.cc
cooling/chemistry_cooling.cc
cooling/cooling_func.cc
cooling/molecular_cooling.cc
cooling/dust_cooling.cc
reaction_network/forward_rates.cc
reaction_network/reaction_network.cc
reaction_network/reaction_rates.cc
reaction_network/shielding.cc
solver/chemistry_solver.cc
solver/saha.cc
thermodynamics/partition.cc
thermodynamics/thermodynamics.cc
"

for module in $modules; do
  compiler=$cc_command
  case "$module" in *.cc) compiler="${CXX:-c++} -std=c++17" ;; esac
  dependencies=$($compiler -MM $feature_flags $include_flags "$module")
  if printf '%s\n' "$dependencies" | grep -Eq '(^|[ /])(allvars|proto|tags)\.h'; then
    printf 'host dependency reached from %s:\n%s\n' "$module" "$dependencies" >&2
    exit 1
  fi
done

printf 'core dependency check passed for migrated modules\n'
