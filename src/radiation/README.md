# Radiation-to-chemistry rate provider

This sibling module owns spectrum-dependent radiation data that does not
belong to the chemistry reaction/cooling core.  It has no dependency on
the sibling chemistry core, AREPO, standalone globals, MPI, or table I/O.

## Lifetime

Each process or MPI rank owns one `RadiationRateContext`.  The host calls
`radiation_rate_context_init()` once, shares the initialized context as a
read-only object, and calls `radiation_rate_context_destroy()` at shutdown.
No context pointer or lookup array is stored in file-scope mutable state.

## Outputs

- `radiation_rate_uniform_background()` converts J21 plus blackbody effective
  temperature into final H-minus, H2, HD, and G0 values.  A host copies those
  values into `ChemGlobalStep` only when its uniform background changes.
- `RadiationBackgroundState` caches those final values. Calling
  `radiation_background_state_update()` with unchanged J21 and effective
  temperature performs no spectral conversion and leaves its revision fixed.
- `radiation_rate_blackbody_spectrum()` returns the undiluted H-minus and G0
  quantities needed by point-source transport.  RayTrace remains responsible
  for geometry, attenuation, and writing cell-local final rates.
- `radiation_hminus_photodetachment_cross_section()` provides the same cross
  section to LW opacity/transport without making radiation include a chemistry
  implementation header.

The numerical constants, 200-bin temperature lookup, 500-bin spectral
integrals, interpolation, and operation order preserve the former standalone and
`chemistry_SC/react_rate_FUV.c` implementation.

For the uniform legacy spectrum, radiation owns the explicit unshielded HD
rate `1.55e-12 J21 s^-1`; chemistry never reconstructs it from the H2 rate.
Spatially varying transport similarly supplies its final, independently
shielded HD rate per cell.
