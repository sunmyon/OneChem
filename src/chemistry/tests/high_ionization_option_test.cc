#include <assert.h>
#include <stdio.h>

#include "species.h"

int main(void)
{
  enum SPECIES species;
  assert(chemistry_species_from_name("CI", &species));
  assert(chemistry_species_from_name("CII", &species));
  assert(chemistry_species_from_name("OI", &species));

#ifdef EXPECT_HIGH_IONIZATION
  assert(chemistry_species_from_name("CIII", &species));
  assert(chemistry_species_from_name("CIV", &species));
  assert(chemistry_species_from_name("OIII", &species));
  assert(chemistry_species_from_name("OIV", &species));
#else
  assert(!chemistry_species_from_name("CIII", &species));
  assert(!chemistry_species_from_name("CIV", &species));
  assert(!chemistry_species_from_name("OIII", &species));
  assert(!chemistry_species_from_name("OIV", &species));
#endif

  puts("high-ionization option boundary check passed");
  return 0;
}
