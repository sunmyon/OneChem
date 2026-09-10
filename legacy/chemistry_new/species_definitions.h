#pragma once

#define CHEM_SPECIES_BASE(X) \
  X(ELEC, "ELEC", -1, 0, 0, 0, 0, 0) \
  X(HI,   "HI",    0, 1, 0, 0, 0, 0) \
  X(HII,  "HII",   1, 1, 0, 0, 0, 0)

#ifdef HELIUM
#define CHEM_SPECIES_HELIUM(X) \
  X(HeI,   "HeI",   0, 0, 1, 0, 0, 0) \
  X(HeII,  "HeII",  1, 0, 1, 0, 0, 0) \
  X(HeIII, "HeIII", 2, 0, 1, 0, 0, 0)
#else
#define CHEM_SPECIES_HELIUM(X)
#endif

#ifdef H2_PLUS
#define CHEM_SPECIES_H2_PLUS(X) X(H2II, "H2II", 1, 2, 0, 0, 0, 0)
#else
#define CHEM_SPECIES_H2_PLUS(X)
#endif

#ifdef HD
#define CHEM_SPECIES_DEUTERIUM_BASE(X) \
  X(DI,  "DI",  0, 0, 0, 1, 0, 0) \
  X(DII, "DII", 1, 0, 0, 1, 0, 0) \
  X(HDI, "HDI", 0, 1, 0, 1, 0, 0)
#ifdef HD_PLUS
#define CHEM_SPECIES_HD_PLUS(X) X(HDII, "HDII", 1, 1, 0, 1, 0, 0)
#else
#define CHEM_SPECIES_HD_PLUS(X)
#endif
#ifdef DMINUS
#define CHEM_SPECIES_D_MINUS(X) X(DM, "DM", -1, 0, 0, 1, 0, 0)
#else
#define CHEM_SPECIES_D_MINUS(X)
#endif
#define CHEM_SPECIES_DEUTERIUM(X) \
  CHEM_SPECIES_DEUTERIUM_BASE(X) \
  CHEM_SPECIES_HD_PLUS(X) \
  CHEM_SPECIES_D_MINUS(X)
#else
#define CHEM_SPECIES_DEUTERIUM(X)
#endif

#ifdef CARBON_OXYGEN_CHEMISTRY
#define CHEM_SPECIES_CARBON_OXYGEN_STORED(X) \
  X(CI,   "CI",    0, 0, 0, 0, 1, 0) \
  X(CII,  "CII",   1, 0, 0, 0, 1, 0) \
  X(CO,   "CO",    0, 0, 0, 0, 1, 1) \
  X(CO2,  "CO2",   0, 0, 0, 0, 1, 2) \
  X(OI,   "OI",    0, 0, 0, 0, 0, 1) \
  X(OH,   "OH",    0, 1, 0, 0, 0, 1) \
  X(H2O,  "H2O",   0, 2, 0, 0, 0, 1) \
  X(O2,   "O2",    0, 0, 0, 0, 0, 2)

#define CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X) \
  X(COII,  "COII",  1, 0, 0, 0, 1, 1) \
  X(CH,    "CH",    0, 1, 0, 0, 1, 0) \
  X(CH2,   "CH2",   0, 2, 0, 0, 1, 0) \
  X(OII,   "OII",   1, 0, 0, 0, 0, 1) \
  X(OHII,  "OHII",  1, 1, 0, 0, 0, 1) \
  X(H2OII, "H2OII", 1, 2, 0, 0, 0, 1) \
  X(H3OII, "H3OII", 1, 3, 0, 0, 0, 1) \
  X(O2II,  "O2II",  1, 0, 0, 0, 0, 2)
#else
#define CHEM_SPECIES_CARBON_OXYGEN_STORED(X)
#define CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X)
#endif

#define CHEM_SPECIES_STORED_LIST(X) \
  CHEM_SPECIES_BASE(X) \
  CHEM_SPECIES_HELIUM(X) \
  X(H2I, "H2I", 0, 2, 0, 0, 0, 0) \
  CHEM_SPECIES_H2_PLUS(X) \
  X(HM, "HM", -1, 1, 0, 0, 0, 0) \
  CHEM_SPECIES_DEUTERIUM(X) \
  CHEM_SPECIES_CARBON_OXYGEN_STORED(X)

#define CHEM_SPECIES_INTERNAL_LIST(X) \
  CHEM_SPECIES_CARBON_OXYGEN_INTERMEDIATES(X)

#define CHEM_SPECIES_LIST(X) \
  CHEM_SPECIES_STORED_LIST(X) \
  CHEM_SPECIES_INTERNAL_LIST(X)
