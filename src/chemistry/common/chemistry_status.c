#include "include/chemistry/chemistry_status.h"

const char *chemistry_status_string(ChemistryStatus status)
{
  switch(status) {
    case CHEMISTRY_STATUS_OK:
      return "ok";
    case CHEMISTRY_STATUS_INVALID_ARGUMENT:
      return "invalid argument";
    case CHEMISTRY_STATUS_INCOMPATIBLE_STATE:
      return "incompatible state";
    case CHEMISTRY_STATUS_IO_ERROR:
      return "I/O error";
    case CHEMISTRY_STATUS_OUT_OF_MEMORY:
      return "out of memory";
    case CHEMISTRY_STATUS_NOT_CONVERGED:
      return "not converged";
    case CHEMISTRY_STATUS_INTERNAL_ERROR:
      return "internal error";
  }
  return "unknown chemistry status";
}
