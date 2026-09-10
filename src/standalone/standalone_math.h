#ifndef STANDALONE_MATH_H
#define STANDALONE_MATH_H

static inline double standalone_min(double left, double right)
{
  return left < right ? left : right;
}

static inline double standalone_max(double left, double right)
{
  return left > right ? left : right;
}

#endif
