# `5d30197` full-survey baseline

This baseline was generated from a clean archive of commit `5d30197` on
2026-08-06 with the default Makefile configuration:

```sh
make -j2
./1ZONE
```

The run completed successfully, produced 36 data files (six metallicities by
six CMB redshifts), and occupied approximately 35 MB.  The SHA-256 manifest is
stored in `full_survey.sha256`; paths are relative to the standalone run root.

An initial explicit-value propagation prototype was built from a separate clean
tree and all 36 files compared byte-for-byte equal with `diff -qr`.  That
prototype was then superseded by the lifetime-correct run/step context design.
The standalone adapter was configured once after unit initialization and its step
context refreshed at each host-wide chemistry step.  On 2026-08-06 that version
also completed the full survey and all 36 files compared byte-for-byte equal to
this baseline with `diff -qr`.

The same full comparison was repeated after changing the transitional eager
initializer to `InitChem(const ChemistryInitConfig *, ChemRunContext *)` and passing `table_root`
explicitly.  All 36 outputs again matched byte-for-byte.

The comparison was repeated after removing the dust, continuum-opacity,
continuum-cooling, and SD93 runtime lazy-initialization fallbacks.  Explicit
`InitChem` initialization produced all 36 outputs byte-for-byte equal to this
baseline.

It was repeated once more after removing host headers from `InitChem` and the
migrated aggregate/H2/molecular/LW-opacity cooling kernels.  Background-rate
and abundance-floor mutation was also removed from `InitChem`; all 36 outputs
remained byte-for-byte equal.

The full comparison was repeated after moving the parsed reaction array,
reaction count, and species reaction indices from mutable file-scope globals
into an explicitly owned `ChemReactionNetwork`.  All network evaluation paths
received that object through the run context, and all 36 outputs remained
byte-for-byte equal.

The comparison was then repeated after moving the common temperature grid,
forward/inverse reaction rates, atomic cooling coefficients, and Saha formation
energies into a heap-owned `ChemistryRateTables` object.  All 36 outputs again
matched byte-for-byte.

The CO, OH, and H2O molecular cooling tables were subsequently moved from
mutable module-static objects into the run-owned `ChemistryCoolingTables`.
Eager loading, interpolation, and destruction completed successfully, and all
36 outputs remained byte-for-byte equal.

The dust temperature grid, Planck/Rosseland opacities, and interpolation
metadata were then moved from mutable globals into a rank-local
`ChemistryDustTables`.  The existing rank-0 read and MPI broadcast ordering was
retained.  All 36 outputs again matched byte-for-byte.

The gas temperature/density grids and Planck/Rosseland continuum-opacity
matrices were then moved into a rank-local
`ChemistryContinuumOpacityTables`.  Evaluation now follows an explicit
run-context pointer, while the existing rank-0 read and MPI broadcast ordering
was retained.  The full 36-case survey again matched byte-for-byte.

The blackbody effective-temperature, H-minus photodetachment, and G0 tables
were subsequently moved through an intermediate `ChemistryRadiationTables`
state.  The ordinary 36-case configuration remained byte-for-byte equal.

That context-owned table was an intermediate migration state. It now lives
outside chemistry in the sibling module's `RadiationRateContext`; the fixed
J21=1, Teff=10^4 K output and J21 scaling are covered by
`make radiation-check`.  A full-density selective LW survey covering all five
J21 values produced 45 files; the active and portable paths both matched the
clean `5d30197` output byte-for-byte after H2, HD, H-minus, and G0 were all
passed explicitly through the host-owned shared-step boundary.

The generic radiation module was subsequently separated from the standalone
adapter and given explicit initialization/destruction plus independent-context
tests. The ordinary portable 36-case survey completed in about 5.7 seconds and
all files again matched this manifest byte-for-byte.

The standalone adapter, build fragment, symbols, and documentation were then
renamed to remove host-specific naming from chemistry. The duplicated
standalone molecular-weight calculation was routed through
`chemistry_molecular_weight()`, and its Jeans-radius wrapper through
`chemistry_jeans_radius()` with the historical `1.0` thermal-mass factor kept
explicit. The 36-case survey again matched every manifest entry byte-for-byte.

The data files themselves are intentionally not committed here.  A later test
tool will compare selected columns/fixtures for portability across compilers
while retaining the full bytewise comparison as the same-toolchain gate.

The 36-case comparison was repeated after defining the AREPO abundance-storage
contract, routing the host `compute_temp(i)` closure through
`chemistry_temperature()`, and consolidating particle-state commit logic in the
adapter.  All manifest checks remained byte-for-byte equal.

The comparison was repeated after replacing both `chem_balance` evaluation
paths with `chemistry_arepo_advance()` and removing MPI table broadcasts from
the core in favor of rank-local context initialization.  The standalone
standalone path again matched all 36 manifest entries byte-for-byte.

The table-root singleton and process-wide `FlagPrint/count_print` diagnostics
were then removed. Every loader receives the initialization root explicitly,
and diagnostic enablement is per `ChemWorkspace`. The standalone portable
survey completed in about 6.0 seconds and all 36 outputs again matched the
manifest byte-for-byte.

The remaining chemistry-owned process termination paths were subsequently
replaced by `ChemistryStatus` propagation. Initialization failures now unwind
owned tables, cooling failures survive implicit-energy residual probes through
the per-call workspace, and host adapters decide how to handle the returned
status. The portable survey completed in about 5.9 seconds and all 36 outputs
again matched this manifest byte-for-byte.
