#!/bin/sh
set -eu

adapter=interface/chem_balance.c
api_adapter=adapters/arepo/chemistry_arepo_adapter.c

if ! grep -q 'return chemistry_temperature(abund, gamma, utherm_cgs);' "$adapter"; then
  echo "error: AREPO compute_temp no longer delegates to chemistry_temperature" >&2
  exit 1
fi

if grep -n 'compute_abundances_d' "$adapter"; then
  echo "error: standalone compatibility entry point leaked into chem_balance" >&2
  exit 1
fi

if ! grep -q 'chemistry_arepo_advance' "$adapter"; then
  echo "error: chem_balance does not call the explicit AREPO adapter" >&2
  exit 1
fi

if grep -n 'chemistry_core_advance' "$api_adapter"; then
  echo "error: AREPO adapter bypasses the public chemistry API" >&2
  exit 1
fi

if ! grep -q 'status = chemistry_advance(' "$api_adapter"; then
  echo "error: AREPO adapter does not call chemistry_advance" >&2
  exit 1
fi

if grep -n -E 'Abund[^;]*(/|\*)[[:space:]]*(XH|XHe|H_number_fraction|He_number_fraction)' "$adapter"; then
  echo "error: abundance normalization leaked into the AREPO adapter" >&2
  exit 1
fi

if grep -n -E 'sph_in\.Mass|in->Mass|double[[:space:]]+Mass;' "$adapter"; then
  echo "error: unused particle mass was added to the chemistry input/communication contract" >&2
  exit 1
fi

commit_count=$(grep -c 'SphP\[i\]\.Abund\[j\] = abund\[j\];' "$adapter")
if [ "$commit_count" -ne 1 ]; then
  echo "error: AREPO abundance commit must have one source of truth" >&2
  exit 1
fi

echo "AREPO storage and thermodynamic-boundary contract check passed"
