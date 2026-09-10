# Physics responsibility boundaries

## Abundance representation and composition

Chemistry defines the species ordering, abundance count, and normalization.
For every stored species `i`, `abundance[i] = n_i / (rho / proton_mass)`.
The host provides storage and advection for this array but does not reinterpret
or renormalize it.  Hydrogen, helium, and deuterium reference abundances are
immutable chemistry ABI metadata declared in
`include/chemistry/chemistry_composition.h`; they are neither `All` fields nor
AREPO hydrodynamic inputs.

The adapter copies `SphP[i].Abund` to and from `ChemistryState` unchanged.  It
alone converts host internal energy units to CGS and commits the accepted
energy change to the host conserved variables.

## Charge neutrality

Charge neutrality is a chemistry-core constraint.  Species charge comes from
the network metadata.  The solver/projection layer reconstructs the dependent
electron abundance according to the selected compatibility profile and reports
whether projection changed the input state.  An adapter must not independently
repair electrons; doing so would make standalone and AREPO results diverge.

The initial compatibility profile preserves the exact `5d30197` set of terms,
including its deliberate omission of trace charged intermediates.  A fully
general charge sum is a later physics/numerics change with separate regression.

## Element conservation

Reaction-table validation checks exact stoichiometric element and charge
conservation at context creation.  During integration, the core owns the
existing reduced conservation system, dependent-species reconstruction,
positivity limiting, and rollback.

Adapters validate array size and species layout but do not impose H/He/D/C/O
conservation.  Changing which trace species enter a reduced conservation
constraint is explicitly a solver-policy change, not an integration change.

## Energy coupling

The host supplies physical specific thermal energy at the beginning of a local
chemistry interval.  The core advances thermal energy using the selected
chemical and radiative terms and returns the accepted energy plus its change.

The core owns:

- chemical formation/destruction heating and cooling;
- radiative and continuum cooling enabled by its profile;
- coupled energy/chemistry iteration;
- gamma/temperature closure used inside that iteration;
- convergence, accepted-step rollback, and chemistry subcycling.

The host adapter owns:

- conversion from entropy/conserved variables to CGS specific energy;
- hydrodynamic work and all non-chemistry source terms;
- conversion of the accepted energy change back to host variables;
- scheduling the cell timestep.

On a hard failure the public call leaves `ChemistryState` unchanged.  A
non-converged result must never be silently committed by an adapter.
