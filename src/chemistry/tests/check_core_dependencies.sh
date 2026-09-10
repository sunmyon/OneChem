#!/bin/sh
set -eu

cc_command=${CC:-cc}
include_flags="-Iinclude -I. -Icommon -Icooling -Idriver -Iinitialization -Ireaction_network -Isolver -Ithermodynamics -Itables"
feature_flags="-DCHEMISTRY_SC -DCARBON_OXYGEN_CHEMISTRY -DHD -DMETALLICITY=-4."

modules="
interface/chemistry_api.c
driver/chemistry_driver.c
species.c
common/chemistry_profile.c
cooling/H2cooling.c
cooling/chemistry_cooling.c
cooling/cooling_func.c
cooling/molecular_cooling.c
cooling/dust_cooling.c
initialization/chemistry_init.c
reaction_network/forward_rates.c
reaction_network/reaction_network.c
reaction_network/reaction_rates.c
reaction_network/shielding.c
solver/chemistry_solver.c
solver/saha.c
thermodynamics/partition.c
thermodynamics/thermodynamics.c
"

for module in $modules; do
  dependencies=$($cc_command -MM $feature_flags $include_flags "$module")
  if printf '%s\n' "$dependencies" | grep -Eq '(^|[ /])(allvars|proto|tags)\.h'; then
    printf 'host dependency reached from %s:\n%s\n' "$module" "$dependencies" >&2
    exit 1
  fi
done

printf 'core dependency check passed for migrated modules\n'
