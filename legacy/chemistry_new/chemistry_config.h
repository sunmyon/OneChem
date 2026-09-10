#pragma once

#ifndef CHEM_CO_O2_DEPENDENT_DENSITY
#define CHEM_CO_O2_DEPENDENT_DENSITY 1.e9
#endif

#if defined(METAL) || defined(CARBON_OXYGEN_CHEMISTRY)
#define CHEM_RATE_COUNT 200
#else
#define CHEM_RATE_COUNT 60
#endif
