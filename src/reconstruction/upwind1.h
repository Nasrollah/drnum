#ifndef UPWIND1_H
#define UPWIND1_H

#include "cartesianpatch.h"

struct Upwind1
{
  void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
               size_t i1, size_t j1, size_t k1,
               size_t,    size_t,    size_t,
               real = 0, real = 0, real = 0,
               real = 0, real = 0, real = 0);
};

inline void Upwind1::project(CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
                             size_t i1, size_t j1, size_t k1,
                             size_t   , size_t   , size_t,
                             real, real, real,
                             real, real, real)
{
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = patch->f(i_field, i_var, i1, j1, k1);
  }
}

#endif // UPWIND1_H
