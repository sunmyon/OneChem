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
    */continuum_metal.dat)
      awk '{ printf "%.17g %.17g\n", $1, 2 * $2 }' "$table" > "$fixture/alternate/tables/continuum_metal.dat" ;;
    *) ln -s "$table" "$fixture/alternate/tables/$(basename "$table")" ;;
  esac
done
"$binary" "$table_root" "$fixture" "$fixture/alternate"
