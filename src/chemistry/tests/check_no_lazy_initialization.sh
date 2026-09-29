#!/bin/sh
set -eu

# These process-wide flags previously caused table loading from cooling
# evaluation paths.  Initialization must remain an explicit lifecycle phase.
forbidden='flag_dust_init|flag_read_SD93|flag_read_opacities|flag_tabulate'

if grep -R -n -E "$forbidden" cooling --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: a cooling-time lazy-initialization flag was reintroduced" >&2
  exit 1
fi

echo "cooling evaluation contains no legacy lazy-initialization flags"
