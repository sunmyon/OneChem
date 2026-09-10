# Standalone chemistry host layer

This directory owns the compatibility boundary between the 1-zone host and
`../chemistry/`. It contains the legacy `chemdata_in` / `chemdata_out`
transfer structures, host-unit conversion, compatibility entry points, and the
optional uniform-background bridge.

None of these files belongs to the portable chemistry core or its AREPO source
list. The repository root includes `standalone.mk` only when
`USE_CHEMISTRY_PORTABLE=1`. Isolation is therefore structural and does not rely
on `CHEMISTRY_STANDALONE_HOST` conditionals inside the chemistry core.

Run the boundary test with:

```sh
make -C src/standalone check
```

When `SURVEY_LYMANWERNER` is enabled, the standalone executable evaluates
`J21 = 0, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000`. Set
`ONEZONE_J21_INDEX` to an index from 0 through 8 to run only one member of this
sequence. Leaving the variable unset runs the complete survey.
