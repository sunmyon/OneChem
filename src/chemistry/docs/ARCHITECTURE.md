# Architecture contract

## Dependency direction

The chemistry implementation depends only on its public/internal headers and
the C runtime.  Host adapters depend on chemistry and their host; chemistry
never depends on an adapter.

```text
new host ---------------------------> public API --> chemistry core
AREPO host ------> AREPO adapter -----^

standalone host --> legacy compatibility adapter --> chemistry core
```

The public entry point converts the public CGS structures into the private
solver structures and calls the unchanged numerical core. Host headers are
confined to adapter translation units.

## Lifetime model

`All` is not replaced by another undifferentiated global structure.  Values are
classified by lifetime:

| Lifetime | Owner | Examples |
|---|---|---|
| run | `ChemistryContext` | tables, network, subcycle factor |
| global chemistry sweep | `ChemistryStepContext` | CMB temperature, scale factor, uniform radiation |
| cell call | `ChemistryCellEnvironment` | density, metallicity, local radiation rates |
| evolving cell | `ChemistryState` | abundances, thermal energy, gamma |
| temporary call | `ChemistryWorkspace` | rates, Jacobian, trial states |

The context becomes read-only after creation.  A step context is prepared once
and passed by pointer to all cells.  A workspace belongs to one concurrent
call.  This removes per-cell copies of global-time values without introducing
hidden mutable process state.

## Units

The core accepts physical CGS quantities.  Unit conversion remains entirely in
each host adapter.  Consequently the public chemistry headers do not expose a
unit structure. The standalone adapter owns `ChemStandaloneUnits`; an AREPO adapter
will own the factors required by AREPO's code-unit convention.

The standalone `compute_abundances_d()` entry remains a compatibility boundary
and currently calls the internal core entry. New hosts call the public API
directly; a host-specific adapter is needed only for storage, unit conversion,
scheduling, or communication.
