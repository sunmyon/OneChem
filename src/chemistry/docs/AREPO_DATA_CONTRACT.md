# AREPO data contract

`SphP[i].Abund` is host storage for a chemistry-owned state.  AREPO does not
define its normalization or species ordering.  Chemistry defines both through
its species configuration and `chemistry_composition.h`.

## Stored state

- `SphP[i].Abund[j]` stores `n_j / (rho / proton_mass)` in chemistry order.
- `SphP[i].AbundMass[j]` is the corresponding AREPO advected quantity,
  `SphP[i].Abund[j] * P[i].Mass`.
- `SphP[i].Gamma` is the gamma accepted by the chemistry thermodynamic closure.
- `SphP[i].Utherm` remains in AREPO units; the adapter converts it to and from
  CGS specific energy.

The adapter copies abundances without multiplying or dividing by `XH`.
Hydrogen-relative values used by reaction kernels are a chemistry-internal
representation and never cross the host boundary.

## Temperature

The chemistry source of truth is:

```c
const double utherm_cgs =
    (All.UnitEnergy_in_cgs / All.UnitMass_in_g) * SphP[i].Utherm;
const double temperature = chemistry_temperature(
    SphP[i].Abund, SphP[i].Gamma, utherm_cgs);
```

An AREPO `compute_temp(i)` function is therefore a thin wrapper that reads the
particle and performs the unit conversion.  It must not maintain a second
molecular-weight or abundance-normalization formula.

The same rule applies to standalone host helpers. A host may import its stored
species into a temporary chemistry abundance array, but molecular weight is
always evaluated by `chemistry_molecular_weight()`. Jeans-radius evaluation is
provided by `chemistry_jeans_radius(rho_cgs, temperature,
thermal_mass_in_proton_masses)`. The last argument is explicit because the
historical chemistry shielding closure uses `1.2`, while the standalone
collapse helper uses `1.0`; hiding that difference in duplicate formulas would
make later unification unsafe.

The multiplication order above intentionally matches the historical AREPO
`compute_temp(i)` implementation.  Other adapter paths retain their existing
unit-conversion order until an AREPO-level bitwise comparison proves that it is
safe to unify them.

## Adapter-owned operations

`/Users/j/Desktop/codex/1ZONE/src/chemistry/adapters/arepo/chem_balance.cc` owns AREPO communication buffers, time-bin conversion,
physical-density conversion, abundance storage/advection, and committing the
accepted thermal-energy difference to `SphP[i].Energy`.  None of those types or
operations belongs in the chemistry core.
