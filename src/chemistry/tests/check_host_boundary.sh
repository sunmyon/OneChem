#!/bin/sh
set -eu

if find adapters -type d -name standalone -print | grep -q .; then
  echo "error: standalone host adapter directory is inside the chemistry core" >&2
  exit 1
fi

if rg -n 'CHEMISTRY_STANDALONE_HOST|chemdata_(in|out)|chemistry_standalone' \
    --glob '*.c' --glob '*.h' --glob '*.mk' \
    --glob '!tests/**' .; then
  echo "error: standalone host contract leaked into portable production code" >&2
  exit 1
fi

printf 'portable host-boundary check passed\n'
