# ONEZONE / AREPO chemistry merge (2026-09-18)

## Shared-source integration, 2026-09-28

Compared ONEZONE `b67a6f5` with AREPO `feature/1zone-chemistry`
at `ce5680c273` (the branch worktree, not the sibling `arepo` checkout,
which currently runs `sink-wind-feedback`). The AREPO worktree is:

```text
/Users/j/Documents/Codex/2026-09-11/referenced-chatgpt-conversation-this-is-an/work/arepo-1zone
```

The Newton and energy-solver numerical changes were already present in both
trees. Shared-core differences were diagnostic output and API formatting.
Imported AREPO's diagnostics, cached photoionization coefficients and their
reference regression data, radiation compatibility function, and host adapter
blackbody/dust-opacity entry points into ONEZONE. Fixed chemical-heating
diagnostics in both trees to print the newly computed `rhoU_dot_chem`, rather
than the previous workspace value; added the missing shielding-log newline.

All 65 shared runtime C/header files now match byte for byte. Check again with:

```sh
python3 tools/check_chemistry_parity.py /path/to/arepo-1zone
```

The command is read-only and exits nonzero for differences or missing files.
It compares the core directories, public API/normalization implementations,
and radiation module. It deliberately excludes host adapters, `chem_balance.c`,
snapshot integration, build files, tests, docs, and editor backup files.
Matching source does not establish matching host inputs or compiler options.

Validation in ONEZONE: `make -C src/chemistry check`, high-ionization option,
metal-ionization, photoionization (including cached-coefficient references),
CO-photodissociation, and elemental normalization checks all passed.
`make check-energy-replay` passed all four cases, including forced chemistry
trial failures exercising the non-finite bisection guard.

The additional initialization regression reports `short-step solve failed`
three times, both in the integrated tree and in a clean archive of pre-merge
`b67a6f5`; this is an existing failure, not a passing validation. The root
Makefile's initialization target also prefixes its executable with `./`, so
an absolute `BUILDDIR` cannot be used to run that target directly; the
integrated test binary was invoked by its absolute path after building.

This is the first integration slice. AREPO-specific snapshot handling,
`chem_balance.c` (including a hard-coded particle diagnostic trigger), build
fragments, and host integration tests still differ. AREPO's full executable
was not rebuilt here. The two repositories still carry separate source
copies; the parity check detects drift but does not automatically synchronize
them. Changes remain uncommitted in both working trees.

## Previous integration (2026-09-18)

Source of truth for imported changes: AREPO `feature/1zone-chemistry`, commit
`121b27cb2e`, under `src/chemistry_ONECHEM`. The comparison base was the older
local AREPO `src/chemistry_ONECHEM`; the destination was the existing, dirty
ONEZONE `src/chemistry`. A three-way file merge was followed by a full shared-
core comparison. The final runtime chemistry source files, including the
solver, species metadata, reactions, cooling, thermodynamics, and public API,
are byte-identical between the two trees. The `mk/core.mk` default path and
AREPO/ONEZONE host adapters remain host-specific. No changes were committed.
In particular, the AREPO solver version was selected for the shared solver;
ONEZONE's former uncommitted solver variant is not retained as a separate
runtime implementation. Its previous content can be recovered from the
pre-merge local AREPO copy used as the merge base.

Imported shared features include optional `CARBON_OXYGEN_HIGH_IONIZATION`
species and reactions, He/C/O photoionization and photoheating inputs, CO
photodissociation, and elemental normalization. The ONEZONE standalone
`--chem-step`/`--energy-step` parser and adapter now pass these rates through;
the AREPO energy diagnostic prints them using the same keys. The replay line
also records `timescale_safety_factor`; ONEZONE applies that value to its
standalone chemistry run context rather than silently using its default.

Example high-ionization replay build from ONEZONE root:

```sh
make -B -j4 ENABLE_LW=1 EXTRA_OPT='-DHELIUM -DCARBON_OXYGEN_HIGH_IONIZATION -DIONIZATION_PHOTON -DRYTR_CARBON_OXYGEN_IONIZATION -DRYTR_NCO -DLW_RADIATION_LOCAL -DTHERMAL_FEEDBACK -DCHEMISTRY_NORMALIZE_ABUNDANCES -DCHEM_REPLAY_EXTERNAL_RADIATION'
```

`CHEM_REPLAY_EXTERNAL_RADIATION` skips the standalone's unavailable local
ionizing-source table initialization. It is for supplied-rate replay, not for
normal radiative transfer. Absent new rate keys default to zero.

Verified in the staging copy: standard and high-ionization builds. Both
repositories' chemistry regression targets pass, including public API/host
checks, high-ionization option, metal-ionization network,
photoionization-bin, CO-photodissociation, and normalization tests. A fixed-T
chemistry smoke test changed He/C/O/CO abundances when nonzero rates were
supplied, and an energy-step diagnostic echoed the exact supplied rates.

The historical particle 1336657807 energy replay is **not an equal-input
comparison**: the old AREPO log does not contain the new rates or a verified
matching build configuration. The merged high-ionization replay with the
missing rates set to zero returned a NaN final temperature at the original
long timestep; this cannot establish a numerical mismatch against AREPO.
Its short-timestep smoke test completed. Rerun AREPO with the expanded
diagnostic before comparing the failed chemistry-loop trial.
