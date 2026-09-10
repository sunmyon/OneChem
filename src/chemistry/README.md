# Chemistry core

This directory is the actively maintained, host-independent chemistry and
cooling implementation. The 1-zone compatibility layer lives in the sibling
`../standalone/` directory; historical implementations live under
`../../legacy/` and are not selected by the default build.

The host-independent public API is implemented in `interface/chemistry_api.c`.
Public headers live under `include/chemistry/` and compile without `allvars.h`,
`proto.h`, `tags.h`, MPI, or GSL headers. AREPO now reaches the numerical core
through this API; its storage, MPI, time-bin, and unit conversion remain in the
AREPO adapter.

Run the public-header check with:

```sh
make -C src/chemistry api-check
```

Diagnostic table dumps are disabled by default. The root build enables
`partition_functions_H2.dat`, and the legacy LW implementation's H-minus
lookup dump when that implementation is selected, only with
`WRITE_CHEMISTRY_DIAGNOSTIC_TABLES=1`.

Build the standalone integration path without sharing object files with the
default build:

```sh
make -j2 USE_CHEMISTRY_PORTABLE=1 \
  BUILDDIR=/tmp/standalone-portable-build EXEC=1ZONE_portable
```

The physical sources started as an exact snapshot of `legacy/chemistry_new/` from
`5d30197`. Migration is now occurring only inside this directory. The first
completed slice introduces run/step/cell lifetimes; legacy state/unit conversion
is outside this tree. The full default survey remains byte-for-byte equal to
`5d30197`.

Uniform-background spectrum conversion is deliberately outside this tree.
The sibling `../radiation/` module owns an explicit `RadiationRateContext`
and converts J21 and effective temperature into H-minus photodetachment, H2
and HD photodissociation, and G0 values. This directory consumes those final
values but contains no background spectrum model or radiation lookup table.

See:

- `docs/ARCHITECTURE.md` for object lifetimes and dependency direction;
- `docs/PUBLIC_API.md` for the external call and data contract;
- `docs/RESPONSIBILITIES.md` for conservation, charge, and energy ownership;
- `../standalone/LEGACY_API.md` for the legacy
  `compute_abundances_d()` compatibility boundary.
- `docs/CONTEXT_LIFETIMES.md` for the exact update frequency of every context.
- `docs/INITIALIZATION.md` for eager initialization and table ownership.
- `docs/AREPO_DATA_CONTRACT.md` for `SphP.Abund`, units, and temperature ownership.
- `docs/AREPO_ADAPTER.md` for the division between `chem_balance`, the AREPO
  adapter, and the chemistry core.
- `docs/AREPO_INTEGRATION_AUDIT.md` for the tested host-side change surface and
  the separate RayTrace/radiation migration.
