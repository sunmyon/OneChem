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

#define CHEM_SPECIES_LIST(X) \
  CHEM_SPECIES_BASE(X) \
  CHEM_SPECIES_HELIUM(X) \
  X(H2I, "H2I", 0, 2, 0, 0, 0, 0) \
  CHEM_SPECIES_H2_PLUS(X) \
  X(HM, "HM", -1, 1, 0, 0, 0, 0) \
  CHEM_SPECIES_DEUTERIUM(X)
