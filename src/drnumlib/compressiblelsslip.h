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

#ifndef COMPRESSIBLELSSLIP_H
#define COMPRESSIBLELSSLIP_H

#include "drnum.h"

template <typename TPatch, typename TGas>
class CompressibleLsSlip
{

public:

  CUDA_DH real extrapolate(real h0, real h1, real h2, real f1, real f2)
  {
    real dh = min(h0-h1, h1-h2);
    return f1 + (h0-h1)*(f1-f2)/dh;

    /*
    real H = h1*h1 - h2*h2;
    real a = (f2*h1*h1 - f1*h2*h2)/H;
    real b = (f1 - f2)/H;
    //printf("h0=%f, h1=%f, h2=%f, f0=%f, f1=%f, f2=%f\n", h0, h1, h2, a + b*h0, f1, f2);
    return a + b*h0*h0;
    */
  }

  CUDA_DH void operate(real* var1, real* var2, real* var, real h0, real h1, real h2, real wgt, real gx, real gy, real gz)
  {
    real p1, T1, u1, v1, w1;
    real p2, T2, u2, v2, w2;
    TGas::conservativeToPrimitive(var1, p1, T1, u1, v1, w1);
    TGas::conservativeToPrimitive(var2, p2, T2, u2, v2, w2);
    real u = wgt*u1 + (1-wgt)*u2;
    real v = wgt*v1 + (1-wgt)*v2;
    real w = wgt*w1 + (1-wgt)*w2;
    real h = wgt*h1 + (1-wgt)*h2;
    real s = 1.0;
    if (h2 >= 0) {
      s = (h-h0)/h*(u*gx + v*gy + w*gz);
    }
    u -= s*gx;
    v -= s*gy;
    w -= s*gz;
    real p = extrapolate(h0, h1, h2, p1, p2);
    real T = extrapolate(h0, h1, h2, T1, T2);
    TGas::primitiveToConservative(p, T, u, v, w, var);

    bool bad = false;
    for (size_t i_var = 0; i_var < 5; ++i_var) {
      if (isnan(var[i_var]) || isinf(var[i_var])) {
        bad = true;
      }
    }
    if (bad) {
      printf("p=%f, T=%f, u=%f, v=%f, w=%f, rho=%f, h0=%f, h1=%f, h2=%f\n", p, T, u, v, w, var[0], h0, h1, h2);
      asm("trap;");
    }

  }

};

#endif // COMPRESSIBLELSSLIP_H
