#ifndef CHEMISTRY_CONFIG_H
#define CHEMISTRY_CONFIG_H

#include <stdint.h>

#define CHEMISTRY_ABI_VERSION 1u

typedef struct ChemistryConfig {
  uint32_t abi_version;
  /* Used only during chemistry_context_create; caller ownership may then end. */
  const char *table_root;

  double subcycle_safety_factor;
} ChemistryConfig;

void chemistry_config_defaults(ChemistryConfig *config);

#endif
