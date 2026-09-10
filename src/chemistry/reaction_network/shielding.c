/* Chemistry implementation split by responsibility. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "chemistry_core.h"
#include "chemistry_models.h"
#include "chemistry_types.h"
#include "cooling_func.h"
#include "forward_rates.h"
#include "partition.h"
#include "reaction_network.h"
#include "shielding.h"
#include "thermodynamics.h"

#ifdef CHEMISTRY_SC

#ifdef SHIELD_LW
double SelfShield(double ColumnDens, double T_gas, double b)
{
  double x = ColumnDens/5.e14;
  double b5 = b/1.;  /* Since "b[km/sec]"*/

  return 0.965/pow(1.+x/b5, 1.1) + 0.035/pow(1.+x, 0.5)*exp(-8.5e-4*pow(1.+x, 0.5));
}

double EachShield(double ColumnDens, double alpha, double delta, double eta)
{
  double x = ColumnDens/eta;

  return exp(-alpha*x)/pow(1.+x, delta);
}
#endif

#endif /* CHEMISTRY_SC */
