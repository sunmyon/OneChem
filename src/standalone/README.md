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
