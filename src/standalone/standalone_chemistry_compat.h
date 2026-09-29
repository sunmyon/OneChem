#ifndef STANDALONE_CHEMISTRY_COMPAT_H
#define STANDALONE_CHEMISTRY_COMPAT_H

#if defined(CHEM_REPLAY_EXTERNAL_RADIATION) && !defined(STANDALONE_EXTERNAL_RADIATION)
#define STANDALONE_EXTERNAL_RADIATION
#endif

/* Build-time compatibility for the three chemistry implementations. */
#ifdef CHEMISTRY_OLD
#include "legacy/chemistry_old/standalone_compat.h"
#elif defined(CHEMISTRY_PORTABLE)
#include "chemistry/chemistry_build_config.h"
#define k_MAX CHEM_RATE_COUNT
#else
#include "legacy/chemistry_new/chemistry_config.h"
#define k_MAX CHEM_RATE_COUNT
#endif

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif

#ifdef CHEMISTRY_PORTABLE
#include "chemistry/chemistry_core.h"
#elif !defined(CHEMISTRY_OLD)
#include "legacy/chemistry_new/chemistry.h"
#endif

#endif
