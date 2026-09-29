#include "solver/linear_algebra.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

static void near(double actual, double expected) {
  assert(std::fabs(actual - expected) < 1.e-13);
}

int main() {
  // Non-symmetric Jacobian detects transposition; BE gives [[2,3],[4,7]].
  const double jacobian[2][2] = {{-2., -6.}, {-8., -12.}};
  double matrix[4];
  numerics::assemble_backward_euler_matrix(2, 0.5,
      [&](int row, int col) { return jacobian[row][col]; }, matrix);
  const double expected[4] = {2., 4., 3., 7.};
  for(int i = 0; i < 4; ++i) near(matrix[i], expected[i]);
  double rhs[2];
  numerics::assemble_backward_euler_rhs(2, 0.5,
      [](int row) { return row == 0 ? 18. : 40.; },
      [](int row) { return row == 0 ? 1. : 2.; }, rhs);
  near(rhs[0], 8.); near(rhs[1], 18.);
  int pivots[2];
  assert(numerics::LapackDense::solve(2, matrix, pivots, rhs) == 0);
  near(rhs[0], 1.); near(rhs[1], 2.);
  for(int i = 0; i < 4; ++i) matrix[i] = expected[i];
  assert(numerics::LapackDense::factorize(2, matrix, pivots) == 0);
  for(int repeat = 1; repeat <= 2; ++repeat) {
    rhs[0] = 8. * repeat; rhs[1] = 18. * repeat;
    assert(numerics::LapackDense::solve_factored(2, matrix, pivots, rhs) == 0);
    near(rhs[0], repeat); near(rhs[1], 2. * repeat);
  }
  double singular[4] = {1., 2., 2., 4.};
  assert(numerics::LapackDense::factorize(2, singular, pivots) > 0);
  double singular_again[4] = {1., 2., 2., 4.};
  assert(numerics::LapackDense::solve(2, singular_again, pivots, rhs) > 0);
  // Independent dimension, no retained factor cache or chemistry headers/config.
  double scalar = 4., scalar_rhs = 12.;
  assert(numerics::LapackDense::solve(1, &scalar, pivots, &scalar_rhs) == 0);
  near(scalar_rhs, 3.);
  std::puts("numerics dense checks passed");
}
