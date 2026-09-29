#!/bin/sh
set -eu
binary=$1
table_root=$(cd "$2" && pwd)
fixture=$(mktemp -d "${TMPDIR:-/tmp}/chemistry-model.XXXXXX")
trap 'rm -rf "$fixture"' EXIT HUP INT TERM
mkdir "$fixture/tables"
for table in "$table_root"/tables/*; do
  case "$table" in */kP.out) continue ;; esac
  ln -s "$table" "$fixture/tables/$(basename "$table")"
done
mkdir "$fixture/alternate" "$fixture/alternate/tables"
for table in "$table_root"/tables/*; do
  case "$table" in
    */reaction_network.tsv)
      sed -e 's/CHEMISTRY_LW_RADIATION/LYMANWERNER/g' -e 's/IONIZING_PHOTON/IONIZATION_PHOTON/g' -e 's/CHEMISTRY_INVERSE_REACTIONS/INVERSE_REACTIONS/g' -e 's/CHEMISTRY_THREE_BODY/THREEBODY/g' -e 's/CHEMISTRY_H2_FORMATION_ON_DUST/H2_FORMATION_ON_DUST/g' -e 's/CHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION/CARBON_OXYGEN_HIGH_IONIZATION/g' -e 's/CHEMISTRY_SPECIES_CARBON_OXYGEN/CARBON_OXYGEN_CHEMISTRY/g' -e 's/CHEMISTRY_SPECIES_DEUTERIUM/HD/g' -e 's/CHEMISTRY_SPECIES_H2_PLUS/H2_PLUS/g' -e 's/CHEMISTRY_SPECIES_HD_PLUS/HD_PLUS/g' -e 's/CHEMISTRY_SPECIES_D_MINUS/DMINUS/g' -e 's/CHEMISTRY_SPECIES_HELIUM/HELIUM/g' "$table" > "$fixture/alternate/tables/reaction_network.tsv" ;;
    */continuum_metal.dat)
      awk '{ printf "%.17g %.17g\n", $1, 2 * $2 }' "$table" > "$fixture/alternate/tables/continuum_metal.dat" ;;
    *) ln -s "$table" "$fixture/alternate/tables/$(basename "$table")" ;;
  esac
done
if ! "$binary" "$table_root" "$fixture" "$fixture/alternate" > "$fixture/report"; then
  cat "$fixture/report"
  exit 1
fi
cat "$fixture/report"
# Only the first successful creation requests a report; later ones are silent.
test "$(grep -c '^Chemistry configuration$' "$fixture/report")" -eq 1
grep -Eq '^  recombination +: CASE [AB]$' "$fixture/report"
grep -Eq '^  chemistry solver +: (hybrid implicit|implicit|explicit)$' "$fixture/report"
if grep -Eq '^  .*: .*CHEMISTRY_' "$fixture/report"; then
  echo 'configuration report exposed a macro name' >&2
  exit 1
fi
