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
#ifndef AUSMDV_H
#define AUSMDV_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <unsigned int DIM, typename TReconstruction, typename TGas>
class AusmDV : public CompressibleFlux<DIM, TGas>
{

  TReconstruction m_Reconstruction;

public: // methods

  using CompressibleFlux<DIM, TGas>::M1;
  using CompressibleFlux<DIM, TGas>::M2;
  using CompressibleFlux<DIM, TGas>::P5;

  template <typename PATCH> CUDA_DH void xField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJX;
    COMPRESSIBLE_RIGHT_PROJX;

    real a  = 0.5*(a_l + a_r);
    real fl = p_l/r_l;
    real fr = p_r/r_r;
    real wp = 2*fl/(fl+fr);
    real wm = 2*fr/(fl+fr);
    real Mp = wp*M2(u_l/a, 1) + (1-wp)*M1(u_l/a, 1);
    real Mm = wm*M2(u_r/a,-1) + (1-wm)*M1(u_r/a,-1);
    real p  = P5(u_l/a,1)*p_l + P5(u_r/a,-1)*p_r;
    countFlops(25);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += 0.5*flux[0]*(u_l + u_r) + A*p - 0.5*fabs(flux[0])*(u_r - u_l);
    flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
    flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
    flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }

  template <typename PATCH> CUDA_DH void yField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJY;
    COMPRESSIBLE_RIGHT_PROJY;

    real a  = 0.5*(a_l + a_r);
    real fl = p_l/r_l;
    real fr = p_r/r_r;
    real wp = 2*fl/(fl+fr);
    real wm = 2*fr/(fl+fr);
    real Mp = wp*M2(v_l/a, 1) + (1-wp)*M1(v_l/a, 1);
    real Mm = wm*M2(v_r/a,-1) + (1-wm)*M1(v_r/a,-1);
    real p  = P5(v_l/a,1)*p_l + P5(v_r/a,-1)*p_r;
    countFlops(25);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
    flux[2] += 0.5*flux[0]*(v_l + v_r) + A*p - 0.5*fabs(flux[0])*(v_r - v_l);
    flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
    flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }

  template <typename PATCH> CUDA_DH void zField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJZ;
    COMPRESSIBLE_RIGHT_PROJZ;

    real a  = 0.5*(a_l + a_r);
    real fl = p_l/r_l;
    real fr = p_r/r_r;
    real wp = 2*fl/(fl+fr);
    real wm = 2*fr/(fl+fr);
    real Mp = wp*M2(w_l/a, 1) + (1-wp)*M1(w_l/a, 1);
    real Mm = wm*M2(w_r/a,-1) + (1-wm)*M1(w_r/a,-1);
    real p  = P5(w_l/a,1)*p_l + P5(w_r/a,-1)*p_r;
    countFlops(25);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
    flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
    flux[3] += 0.5*flux[0]*(w_l + w_r) + A*p - 0.5*fabs(flux[0])*(w_r - w_l);
    flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }

};

#endif // AUSMDV_H
