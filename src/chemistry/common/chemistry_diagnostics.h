#pragma once

typedef struct {
  double H2_escape_fraction;
#ifdef HD
  double HD_escape_fraction;
#endif
  double edot_rad, edot_chem;
  double edot_H2, edot_comp, edot_ceHI, edot_radHI;
  double edot_ciHI, edot_H2form_gas, edot_HMform;
  double edot_H2form_HM, edot_H2form_dust, edot_H2form_3body;
  double n_crit;
#ifdef HELIUM
  double edot_ceHeI, edot_ceHeII;
  double edot_ciHeI, edot_ciHeII, edot_ciHeIS;
  double edot_reHeII1, edot_reHeII2, edot_reHeIII;
#endif
#ifdef HD
  double edot_HD;
#endif
#ifdef METALLICITY
  double edot_CI, edot_CII, edot_OI;
#endif
#ifdef CARBON_OXYGEN_CHEMISTRY
  double edot_CO, edot_OH, edot_H2O;
#endif
#ifdef DUST_COOLING
  double edot_dust;
#ifdef LYMANWERNER
  double edot_PE;
#endif
#endif
#ifdef CONTINUUM_COOLING_MATSUKOBA
  double edot_cont;
#else
  double edot_reHI, edot_brem, edot_CIE, edot_cont;
#endif
} ChemDiagnostics;
