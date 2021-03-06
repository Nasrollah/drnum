// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

#include "drnum.h"

namespace GlobalDebug
{

using namespace std;

extern size_t i;
extern size_t j;
extern size_t k;
extern real   x;
extern real   y;
extern real   z;

#ifndef __CUDACC__
  inline void CUDA_DH print()
  {
    printf("last variables:\n");
    printf("---------------\n");
    printf("i = %d\n", int(i));
    printf("j = %d\n", int(j));
    printf("k = %d\n", int(k));
    printf("x = %f\n", x);
    printf("y = %f\n", y);
    printf("z = %f\n", z);
    printf("\n");
  }
#else
  inline void CUDA_DH print() {}
#endif


#ifdef DEBUG
  #ifndef __CUDACC__
    inline CUDA_DH void xyz(real X, real Y, real Z)
    {
      x = X;
      y = Y;
      z = Z;
    }
  #else
    inline CUDA_DH void xyz(real, real, real) {}
  #endif
#else
  inline CUDA_DH void xyz(real, real, real) {}
#endif

#ifdef DEBUG
  #ifndef __CUDACC__
    inline CUDA_DH void ijk(size_t I, size_t J, size_t K)
    {
      i = I;
      j = J;
      k = K;
    }
  #else
    inline CUDA_DH void ijk(size_t, size_t, size_t) {}
  #endif
#else
  inline CUDA_DH void ijk(size_t, size_t, size_t) {}
#endif

}

#define DBG_PRT_INT(X) printf("variable %s = %d\n", #X, X);
#define DBG_PRT_REAL(X) printf("variable %s = %f\n", #X, X);

#endif // DEBUG_H
