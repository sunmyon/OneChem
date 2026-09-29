#include "linear_algebra.hpp"

extern "C" {
void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
void dgetrf_(int *, int *, double *, int *, int *, int *);
void dgetrs_(char *, int *, int *, const double *, int *, const int *,
             double *, int *, int *);
}

namespace numerics {

int LapackDense::solve(int size, double *matrix, int *pivots, double *rhs) {
  int count = 1, info = 1;
  dgesv_(&size, &count, matrix, &size, pivots, rhs, &size, &info);
  return info;
}

int LapackDense::factorize(int size, double *matrix, int *pivots) {
  int info = 1;
  dgetrf_(&size, &size, matrix, &size, pivots, &info);
  return info;
}

int LapackDense::solve_factored(int size, const double *factors,
                                const int *pivots, double *rhs) {
  char transpose = 'N';
  int count = 1, info = 1;
  dgetrs_(&transpose, &size, &count, factors, &size, pivots, rhs, &size, &info);
  return info;
}

} // namespace numerics
