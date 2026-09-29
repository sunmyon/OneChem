# Energy-step replay regressions

`make check-energy-replay` replays captured AREPO chemistry/energy steps through
the standalone `1ZONE --energy-step` and checks that the energy solver returns a
physical answer.  It builds two binaries (`build-replay`, `build-replay-gate1`)
with the chemistry options in `arepo_chemistry_flags.txt`.

## Fixtures

Each `energy_step_*.txt` is one `/*CHEMISTRY_RERUN_V1*/` line printed by the AREPO
run (state entering chemistry, timestep, density, radiation rates) for the
particle ID in its name.

| fixture | what it is | expected |
|---|---|---|
| `energy_step_1333015846_call2.txt` | T_old = 18874 K, strong local radiation (n ~ 1.9e3 cm^-3). With the old certification 23 of 38 trial chemistry solves failed and the solver returned 138196 K. | T_final = 19630.6 K |
| `energy_step_1333015846_call1.txt` | earlier call of the same particle (T_old = 19475 K) | T_final = 18663.8 K |
| `energy_step_1339855177.txt` | earlier problem case; chemistry trials at 2.7e5-4.4e5 K stalled | T_final = 7995.06 K |

`arepo_chemistry_flags.txt` is the output of
`tools/arepo_config_to_extra_opt.py <Config.sh of the AREPO run>` (only options
that 1ZONE's own source tests with `#ifdef`).  Without them the fixtures do not
parse (unknown radiation keys).

## Assertions

For every case: `converged=1`, `T_final` within the stated tolerance, no
`CHEM_NOT_CONVERGED`, and the normalized energy residual
`|e_new - e_old - de| / max(e_new, e_old)` below 1e-3 (the 138196 K answer had
about 3.9, so the returned temperature must be a root of the energy equation).

* `call2_default` - default build (`CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT = 1e3`): no
  failed trial (`dT=nan`) messages and the bisection guard silent.
* `call2_gate1` - built with `-DCHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT=1.0` so chemistry
  trials fail as before.  Failed trials **must** occur and `solve_bisection` **must**
  report a non-finite residual, and the answer must still be 19630.6 K.  This
  guards the bisection fix independently of the residual-gate relaxation; if a
  future change stops the trials failing, the test says so instead of passing
  vacuously.
* `call1_default`, `id1339855177` - default build, clean.

## Tolerances

Different builds/compilers land on slightly different roots of the energy
equation (the energy solver stops at a relative step of ~1e-5), so temperatures
are compared to a few 1e-5 relative, not bit-exactly.

The captured cases predate the species-dependent helium continuum cooling
default. Their flags explicitly set `CHEMISTRY_HELIUM_COOLING=0` to reproduce
the original physics; new helium-enabled builds default to `1`.

The captured cases also retain `CHEMISTRY_UPDATE_CORE_LENGTH=0`: their original
solver used the host-supplied length rather than recomputing it at each trial.
