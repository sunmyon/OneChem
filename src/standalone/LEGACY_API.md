# Legacy `compute_abundances_d()` boundary

The existing signature remains available only through a compatibility adapter:

```c
int compute_abundances_d(int mode,
                         const struct chemdata_in *input,
                         struct chemdata_out *output,
                         double duration_s);
```

The adapter preserves the current observable contract:

- the integer `mode` meaning and return-code mapping;
- all input radiation, density, metallicity, and diagnostic fields;
- the feature-dependent stored abundance order and count;
- host-unit versus CGS conversions exactly as in `5d30197`;
- output gamma, timescales, dust state, rates, diagnostics, and solver flags;
- state commit/rollback behavior on convergence failure.

Internally the adapter performs only these operations:

1. validate and convert legacy input/output into global-step, cell, and state
   objects;
2. invoke `chemistry_advance()`;
3. translate status and accepted outputs back to the legacy structures.

The standalone host configures the adapter once with a borrowed `ChemRunView` from `chemistry::Model` and
`ChemStandaloneUnits`, then prepares `ChemGlobalStep` only when shared conditions
change.  CMB temperature, scale factor, and uniform backgrounds are therefore
not rebuilt or copied from `All` for every cell call.

It does not evaluate rates, enforce charge or element constraints, select an
ODE method, or modify thermal energy beyond unit conversion.  Those operations
belong to the core and therefore remain identical for standalone and AREPO.

Before replacing an existing call site, a dual-call regression must feed both
paths identical CGS state/environment values and compare every abundance,
thermal energy, gamma, timescale, convergence flag, and enabled diagnostic.
