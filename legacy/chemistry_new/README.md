# Portable chemistry module

`chemistry_new/` is being reorganized into a chemistry and cooling library that
can be used by both the standalone driver and AREPO. The current
source tree is an intermediate state: the physical implementation is already
grouped by responsibility, but a number of files still include the host
program's `allvars.h`, `proto.h`, and `tags.h` directly.

The target architecture and the staged migration plan are specified in
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).  In particular, the document
defines:

- the boundary between chemistry physics and host-code policy;
- the public API and ownership of configuration, units, state, and workspaces;
- the ODE solver interface;
- the roles of the AREPO and standalone adapters;
- the numerical-regression gates required before changing an integration path.

The architecture document is a design contract, not a claim that all of the
described directories and interfaces have already been implemented.
