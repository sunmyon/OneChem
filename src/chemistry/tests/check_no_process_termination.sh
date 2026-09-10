#!/bin/sh
set -eu

if rg -n '\b(endrun|exit|_Exit|quick_exit|abort|raise|MPI_Abort)[[:space:]]*\(' \
    --glob '*.c' --glob '*.h' --glob '!tests/**' .; then
  echo "error: chemistry core must return status instead of terminating the host process" >&2
  exit 1
fi

printf 'portable chemistry contains no process-termination calls\n'
