# LW integration regression

The default full-survey gate does not compile `LYMANWERNER`.  After migrating
the FUV/LW blackbody lookup table, both the `5d30197` implementation and the
portable implementation were therefore built with `ENABLE_LW=1` and run at
one metallicity and one CMB redshift across all five survey LW strengths.

The table has moved completely into the sibling radiation-rate module.
`radiation-check` fixes the legacy J21=1, Teff=10^4 K H-minus, H2, HD, and G0
rates as a direct regression gate.  A full-density selective survey was then
run for all five LW strengths with one metallicity and one CMB redshift.  It
produced 45 files. Both the historical `legacy/chemistry_new` path and the standalone
`USE_CHEMISTRY_PORTABLE=1` path matched the clean `5d30197` results
byte-for-byte for all 45 files.

This comparison caught two missing shared-step inputs before the gate passed:
the HD photodissociation rate and the G0 value consumed by photoelectric
heating.  Both are now produced outside chemistry and passed through
`ChemGlobalStep`; neither is reconstructed per cell.

The comparison also exposed a pre-existing missing declaration for
`dust_PEheating` in LW builds.  Its existing implementation is now declared in
`dust_cooling.h`; no physical expression was changed.
