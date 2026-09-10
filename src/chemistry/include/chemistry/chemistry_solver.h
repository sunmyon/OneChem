#ifndef CHEMISTRY_SOLVER_H
#define CHEMISTRY_SOLVER_H

typedef enum ChemistryMode {
  /* Values intentionally preserve the legacy integer mode contract. */
  CHEMISTRY_MODE_DEFAULT = 1
} ChemistryMode;

typedef enum ChemistrySolverKind {
  /* The compiled chemistry configuration selects explicit/implicit/Saha. */
  CHEMISTRY_SOLVER_AUTO = 0
} ChemistrySolverKind;

typedef struct ChemistryAdvanceOptions {
  ChemistryMode mode;
  ChemistrySolverKind solver;
} ChemistryAdvanceOptions;

#endif
