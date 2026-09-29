#ifndef CHEMISTRY_CONFIG_H
#define CHEMISTRY_CONFIG_H


#include <stdint.h>

#define CHEMISTRY_ABI_VERSION 2u

typedef struct ChemistryConfig {
  uint32_t abi_version;
  /* Used only during Model::create (or legacy context creation). */
  const char *table_root;

  double timescale_safety_factor;
  /* Print resolved build settings after successful creation; default: enabled. */
  int print_configuration;
} ChemistryConfig;


#ifdef __cplusplus
void chemistry_config_defaults(ChemistryConfig *config);
#endif


#endif
