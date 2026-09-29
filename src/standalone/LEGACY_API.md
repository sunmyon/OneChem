# Standalone C++ adapter

ChemistryStandaloneAdapter replaces the old standalone compute_abundances_d
function and its hidden singleton. The host creates an adapter, calls configure
with a borrowed Model run view and unit factors, prepares its step, then calls
advance or solve_fixed. Configuration invalidates any previously prepared step.
Each instance is independent; its Model must remain alive until it is detached
or destroyed. Callers synchronize configuration/step changes against evaluation.

The host input/output value structures, mode semantics, status mapping, CGS unit
conversion and commit/rollback behavior are preserved. The adapter calls native
Driver and solver methods. Rate evaluation and physical constraints remain in
the shared core. Radiation background conversion is a member operation defined
in the existing background implementation file.
