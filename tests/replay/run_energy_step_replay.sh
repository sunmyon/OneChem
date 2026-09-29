#!/bin/sh
# Energy-step replay regressions (run from anywhere; `make check-energy-replay`).
#
# Builds the standalone 1ZONE twice with the chemistry options of the AREPO run
# that produced the fixtures (tests/replay/arepo_chemistry_flags.txt):
#
#   default  - the shipped configuration (CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT = 1e3)
#   gate1    - CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT = 1.0, the pre-relaxation
#              certification.  Chemistry trials of ID 1333015846 call 2 fail
#              again, which exercises the energy solver's handling of failed
#              trial temperatures (solve_bisection).
#
# and replays `--energy-step` on each fixture.  See tests/replay/README.md.

root=$(cd "$(dirname "$0")/../.." && pwd) || exit 2
cd "$root" || exit 2

flags=$(cat tests/replay/arepo_chemistry_flags.txt)
failures=0

fail() { echo "FAIL: $*"; failures=$((failures + 1)); }

# The standalone locates its tables relative to its own executable, so the
# binaries live in the repository root (both names are git-ignored).
build() {
  dir=$1
  exe=$2
  extra=$3
  echo "== building $exe ($extra)"
  if ! make -s BUILDDIR="$dir" EXEC="$exe" \
       EXTRA_OPT="$flags $extra" > "$dir.log" 2>&1; then
    echo "build of $exe failed; see $root/$dir.log"
    exit 2
  fi
  mkdir -p "$dir"
}

# run_case <binary> <label> <fixture> <expected T_final [K]> <tolerance [K]> <mode>
#   mode = clean   : no failed chemistry trials may be reported
#   mode = failing : failed trials MUST occur (the case must stay a real test
#                    of the failure path) yet the answer must still be physical
run_case() {
  bin=$1; label=$2; fixture=$3; expect=$4; tol=$5; mode=$6
  out="build-replay-$label.log"
  "$bin" --energy-step "tests/replay/$fixture" > "$out" 2>&1
  status=$?

  final=$(sed -n 's/^CHEM_STEP converged=\([01]\).*T_final_K=\([^ ]*\).*/\1 \2/p' "$out" | tail -1)
  converged=${final%% *}
  temperature=${final##* }
  nan_msgs=$(grep -c 'dT=nan' "$out")
  guard_msgs=$(grep -c 'solve_bisection: non-finite' "$out")
  not_conv=$(grep -c 'CHEM_NOT_CONVERGED' "$out")
  # normalized energy residual of the last committed evaluation:
  #   |e_new - e_old - de| / max(|e_new|, |e_old|)
  residual=$(grep '^resid(e_new-e_old-de)=' "$out" | tail -1 | awk '{
      split($1, a, "="); g = a[2] < 0 ? -a[2] : a[2];
      split($2, b, "="); e = b[2] < 0 ? -b[2] : b[2];
      split($3, c, "="); o = c[2] < 0 ? -c[2] : c[2];
      m = e > o ? e : o; if (m < 1) m = 1; printf "%.3e", g / m }')

  printf "   %-34s converged=%s T_final=%s (expect %s +/- %s) residual=%s dT=nan:%s bisect-guard:%s\n" \
         "$label" "${converged:-?}" "${temperature:-?}" "$expect" "$tol" "${residual:-?}" "$nan_msgs" "$guard_msgs"

  [ "$status" -eq 0 ] && [ "${converged:-0}" = "1" ] || { fail "$label: solve did not converge (exit $status)"; return; }
  awk -v t="$temperature" -v e="$expect" -v d="$tol" 'BEGIN { x = t - e; if (x < 0) x = -x; exit !(x <= d) }' \
    || fail "$label: T_final=$temperature is not within $tol K of $expect K"
  awk -v r="${residual:-1}" 'BEGIN { exit !(r < 1.0e-3) }' \
    || fail "$label: energy residual $residual is not small; the returned T is not a root of the energy equation"
  [ "$not_conv" -eq 0 ] || fail "$label: CHEM_NOT_CONVERGED reported"

  case $mode in
    clean)
      [ "$nan_msgs" -eq 0 ] || fail "$label: $nan_msgs failed-trial (dT=nan) message(s); expected none"
      [ "$guard_msgs" -eq 0 ] || fail "$label: bisection guard fired; expected none"
      ;;
    failing)
      [ "$nan_msgs" -ge 1 ] || fail "$label: no failed chemistry trials occurred, so this case no longer exercises the failure path; regenerate the fixture"
      [ "$guard_msgs" -ge 1 ] || fail "$label: solve_bisection was never given a non-finite residual, so the guard is not exercised"
      ;;
  esac
}

build build-replay 1ZONE_replay ""
build build-replay-gate1 1ZONE_replay_gate1 "-DCHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT=1.0"

echo "== replaying"
d=./1ZONE_replay
g=./1ZONE_replay_gate1
# ID 1333015846 call 2: 23 of 38 chemistry trials failed with the old gate and
# the solver returned 138196 K (residual ~3.9); the root is ~19630.6 K.
run_case $d call2_default energy_step_1333015846_call2.txt 19630.6 1.0 clean
run_case $g call2_gate1   energy_step_1333015846_call2.txt 19630.6 1.0 failing
# other cases that were problematic earlier
run_case $d call1_default energy_step_1333015846_call1.txt 18663.8 0.5 clean
run_case $d id1339855177  energy_step_1339855177.txt          7995.06 0.2 clean

if [ "$failures" -ne 0 ]; then
  echo "energy-step replay regressions: $failures failure(s)"
  exit 1
fi
echo "energy-step replay regressions: all passed"
