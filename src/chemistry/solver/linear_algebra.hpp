#pragma once

namespace numerics {

// Stateless LP64 double-precision LAPACK backend. size must be positive.
// Caller owns contiguous column-major size*size matrix, size pivots and size RHS.
// All methods return LAPACK INFO (0 success); fresh solve and factorize overwrite
// the matrix. solve_factored requires successful factors of the SAME dimension.
// Factor lifetime, invalidation and reuse policy belong to the calling session.
class LapackDense final {
public:
  static int solve(int size, double *matrix, int *pivots, double *rhs);
  static int factorize(int size, double *matrix, int *pivots);
  static int solve_factored(int size, const double *factors,
                            const int *pivots, double *rhs);
};

} // namespace numerics


#include <cmath>

namespace numerics {

// Caller supplies reduced-coordinate accessors. Output is contiguous column-major.
// derivative(row, column) returns dR_row/dy_column; buffers contain size^2/size doubles.
template<class Derivative>
void assemble_backward_euler_matrix(int size, double dt, Derivative derivative,
                                   double *matrix) {
  for(int column = 0; column < size; ++column)
    for(int row = 0; row < size; ++row)
      matrix[column * size + row] = column == row
          ? 1.0 - dt * derivative(row, column)
          : -dt * derivative(row, column);
}

template<class Rate, class Change>
void assemble_backward_euler_rhs(int size, double dt, Rate rate, Change change,
                                double *rhs) {
  for(int row = 0; row < size; ++row)
    rhs[row] = std::fma(rate(row), dt, -change(row));
}

} // namespace numerics
