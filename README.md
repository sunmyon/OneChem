# 1ZONE chemistry layout

The actively maintained implementation lives under `src/`:

- `src/chemistry/`: host-independent chemistry and cooling core;
- `src/radiation/`: spectrum-to-rate calculations outside the chemistry core;
- `src/standalone/`: the 1-zone executable, host state, and compatibility adapter.

Historical implementations are retained under `legacy/` for regression and
comparison only:

- `legacy/chemistry_new/`;
- `legacy/chemistry_old/`;
- `legacy/chemistry_original/`.

The default build uses `src/chemistry/`:

```sh
make -j2
```

The historical implementations remain explicitly selectable:

```sh
make -j2 USE_CHEMISTRY_PORTABLE=0
make -j2 USE_CHEMISTRY_OLD=1
```

Core and standalone boundary checks are available with:

```sh
make -C src/chemistry check
make -C src/standalone check
```
