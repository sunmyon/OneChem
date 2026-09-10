#ifndef CHEMISTRY_AREPO_ADAPTER_H
#define CHEMISTRY_AREPO_ADAPTER_H

#include "../../chemistry_core.h"
#include "../../include/chemistry/chemistry.h"

typedef struct {
  const char *table_root;
  double subcycle_safety_factor;
} ChemArepoConfig;

/* One explicitly owned instance per MPI rank; immutable after initialization. */
typedef struct {
  ChemistryContext *chemistry;
} ChemArepoContext;

/* Physical per-cell inputs. Shared sweep values remain in the public step. */
typedef struct {
  double rho_cgs;
  double metallicity;

  double J21;
  double rate_Hminus_pd;
  double rate_H2_diss;
#ifdef HD
  double rate_HD_diss;
#endif
  double G0;
  double stellar_radiation_flux;
  double photoionization_rate_HI;
  double photoheating_rate_HI;

  /* Reserved until AREPO supplies these radiation rates. */
  double photoionization_rate_HeI;
  double photoheating_rate_HeI;
  double photoionization_rate_HeII;
  double photoheating_rate_HeII;
  double rate_CO_diss;

  int diagnostics_enabled;
  int diagnostic_id;
} ChemArepoCell;

/* Physical evaluation state. chem_balance converts AREPO energy exactly once. */
typedef struct {
  double abund[TGCHEM_NUM_ABUNDANCES];
  double utherm_cgs;
  double gamma;
} ChemArepoState;

typedef struct {
  double chemistry_time_s;
  double cooling_time_s;
  double dust_temperature;
  double dust_kappa_planck;
  int implicit_flag;
  int convergence_flag;
} ChemArepoResult;

ChemistryStatus chemistry_arepo_context_init(const ChemArepoConfig *config,
                                             ChemArepoContext *context);
void chemistry_arepo_context_destroy(ChemArepoContext *context);

ChemistryStatus chemistry_arepo_advance(const ChemArepoContext *context,
                                        const ChemistryStepContext *step,
                                        const ChemArepoCell *cell,
                                        double duration_s,
                                        ChemArepoState *state,
                                        ChemArepoResult *result);

#endif
