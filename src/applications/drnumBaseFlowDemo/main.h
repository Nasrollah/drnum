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

#ifndef BASEFLOWDEMO_MAIN_H
#define BASEFLOWDEMO_MAIN_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
//#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
//#include "fluxes/ausmdv.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "fluxes/compressibleslipflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"
#include "externalexchangelist.h"

#ifdef GPU
#include "iterators/gpu_cartesianiterator.h"
#else
#include "iterators/cartesianiterator.h"
#endif

#include "rungekutta.h"
#include "iteratorfeeder.h"
#include "cubeincartisianpatch.h"

#include "configmap.h"
#include "timeaverage.h"

typedef Upwind2<5, VanAlbada> reconstruction_t;
//typedef Upwind1<5> reconstruction_t;

#include "baseflowflux.h"

void run()
{
  dim_t<5> dim;

  // control files
  ConfigMap config;
  config.addDirectory("control");

  real Ma             = config.getValue<real>("Mach-number");
  real p              = config.getValue<real>("pressure");
  real T              = config.getValue<real>("temperature");
  real u              = Ma*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
  real L              = 2*config.getValue<real>("radius");
  real time           = L/u;
  real cfl_target     = config.getValue<real>("CFL-number");
  real write_interval = config.getValue<real>("write-interval")*time;
  real total_time     = config.getValue<real>("total-time")*time;
  bool mesh_preview   = config.getValue<bool>("mesh-preview");
  real scale          = config.getValue<real>("scale");
  real sample_rate    = config.getValue<real>("sample-rate");
  int  thread_limit   = 0;
  int  base_patch_id  = config.getValue<int>("base-patch");


  if (config.exists("thread-limit")) {
    thread_limit = config.getValue<int>("thread-limit");
  }

#ifdef GPU
  int  cuda_device    = config.getValue<int>("cuda-device");
#endif

  bool start_from_zero = false;
  if (config.exists("start-from-zero")) {
    start_from_zero = config.getValue<bool>("start-from-zero");
  }

  // Patch grid
  PatchGrid patch_grid;
  //.. general settings (apply to all subsequent patches)
  patch_grid.setNumberOfFields(3);
  patch_grid.setNumberOfVariables(5);
  patch_grid.defineVectorVar(1);
  patch_grid.setInterpolateData();
  patch_grid.setNumSeekLayers(2);  /// @todo check default = 2
  patch_grid.setTransferType("padded_direct");
  patch_grid.readGrid("patches/standard.grid", scale);
  //patch_grid.readGrid("patches/V1");
  patch_grid.computeDependencies(true);

  // Time step
  real ch_speed = u + sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
  real dt       = cfl_target*patch_grid.computeMinChLength()/ch_speed;

  cout << " patch_grid.computeMinChLength() = " << patch_grid.computeMinChLength() << endl;
  cout << " dt  =  " << dt << endl;

  // Initialize
  real init_var[5];
  PerfectGas::primitiveToConservative(0.1*p, T, 0, 1e-3*u, 1e-3*u, init_var);
  patch_grid.setFieldToConst(0, init_var);

  patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), 0);

  if (mesh_preview) {
    exit(EXIT_SUCCESS);
  }

  BaseFlowFlux flux;
  flux.setup(u, p, T, config.getValue<real>("radius"), base_patch_id);

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

#ifdef GPU
  GPU_CartesianIterator<5, BaseFlowFlux>  iterator_std(flux, cuda_device, thread_limit);
#else
  BUG;
  CartesianIterator<5, BaseFlowFlux>  iterator_std(flux);
#endif
  iterator_std.setCodeString(CodeString("fx fy fz far far far far far far 0"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator_std);
  iterator_feeder.feed(patch_grid);

  runge_kutta.addIterator(&iterator_std);

  int write_counter = 0;
  int iter = 0;
  real t = 0;

  QString restart_file = config.getValue<QString>("restart-file");
  if (restart_file.toLower() != "none") {
    write_counter = restart_file.right(6).toInt();
    t = patch_grid.readData(0, "data/" + restart_file);
    if (start_from_zero) {
      t = 0;
    }
  }

#ifdef GPU
  iterator_std.updateDevice();
#endif

  QVector<TimeAverage> centre_averages, vertical_averages;

  startTiming();
  real next_write_time = t + write_interval;

  while (t < total_time) {

#ifdef GPU

    {
      CartesianPatch *patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(base_patch_id));
      if (!patch) {
        BUG;
      }
      GPU_CartesianPatch gpu_patch(patch);

      size_t j_centre = patch->sizeJ()/2;
      size_t k_centre = patch->sizeK()/2;
      //gpu_patch.partialCopyFromDevice<5>(patch, 0, 0, patch->sizeI(), j_centre, j_centre + 1, k_centre, k_centre + 1); // centre line for velocity
      //gpu_patch.partialCopyFromDevice<5>(patch, 0, 0, 1, j_centre, j_centre + 1, 0, patch->sizeK());                   // vertical line for pressure
      gpu_patch.copyFromDevice(patch);

      if (vertical_averages.size() == 0) {
        vertical_averages.fill(TimeAverage(time/sample_rate), patch->sizeK());
      }
      if (centre_averages.size() == 0) {
        centre_averages.fill(TimeAverage(time/sample_rate), patch->sizeI());
      }

      real t_old = t;

      omp_set_num_threads(2);
      #pragma omp parallel
      {
        if (omp_get_num_threads() != 2) {
          BUG;
        }
        if (omp_get_thread_num() == 0) {
          runge_kutta(dt);
          t += dt;
        } else {
          for (size_t i = 0; i < patch->sizeI(); ++i) {
            real var[5], p, T, u, v, w;
            patch->getVarDim(0, i, j_centre, k_centre, var);
            PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
            centre_averages[i].add(t_old, u);
          }
          for (size_t k = 0; k < patch->sizeK(); ++k) {
            real var[5], p_base, T;
            patch->getVarDim(0, 0, j_centre, k, var);
            PerfectGas::conservativeToPrimitive(var, p_base, T);
            real cp = 2*(p_base/p - 1)/(PerfectGas::gamma()*Ma*Ma);
            vertical_averages[k].add(t_old, cp);
          }
          if (vertical_averages[0].valid()) {
            ofstream f("cp.csv", std::ios_base::out);
            f << "r/R, cp\n";
            for (size_t k = 0; k < patch->sizeK(); ++k) {
              real kr = real(k);
              real krc = real(k_centre);
              real r = (kr - krc)*patch->dz();
              f << 2*r/L << ", " << vertical_averages[k].value() << "\n";
            }
            f.flush();
          }
          if (centre_averages[0].valid()) {
            ofstream f("u.csv", std::ios_base::out);
            f << "x/R, u\n";
            for (size_t i = 0; i < patch->sizeI(); ++i) {
              real x = i*patch->dx();
              f << 2*x/L << ", " << centre_averages[i].value()/u << "\n";
            }
            f.flush();
          }
        }
      }
    }

#else

    runge_kutta(dt);
    t += dt;

#endif

    if (t >= next_write_time) {

      // Do some diagnose on patches
      real CFL_max = 0;
      real rho_min = 1000;
      real rho_max = 0;
      real max_norm_allpatches = 0.;
      real l2_norm_allpatches;
      real ql2_norm_allpatches = 0.;

#ifdef GPU
      runge_kutta.copyDonorData(0);
      iterator_std.updateHost();
#endif

      for (size_t i_p = 0; i_p < patch_grid.getNumPatches(); i_p++) {
        CartesianPatch& patch = *(dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_p)));
        size_t NI = patch.sizeI();
        size_t NJ = patch.sizeJ();
        size_t NK = patch.sizeK();

        for (size_t i = 0; i < NI; ++i) {
          for (size_t j = 0; j < NJ; ++j) {
            for (size_t k = 0; k < NK; ++k) {
              real p, u, v, w, T, var[5];
              patch.getVar(dim, 0, i, j, k, var);
              rho_min = min(var[0], rho_min);
              rho_max = max(var[0], rho_max);
              PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
              real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
              CFL_max = max(CFL_max, fabs(u)*dt/patch.dx());
              CFL_max = max(CFL_max, fabs(u+a)*dt/patch.dx());
              CFL_max = max(CFL_max, fabs(u-a)*dt/patch.dx());
              countSqrts(1);
              countFlops(10);
            }
          }
        }
        real max_norm, l2_norm;
        patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
        if (max_norm > max_norm_allpatches) {
          max_norm_allpatches = max_norm;
        }
        ql2_norm_allpatches += l2_norm * l2_norm;
      }
      ql2_norm_allpatches /= patch_grid.getNumPatches();
      l2_norm_allpatches = sqrt(ql2_norm_allpatches);
      ++iter;
      cout << iter << " iterations,  t=" << t/time << "*L/u_oo,  dt: " << dt;
      cout << "  CFL: " << CFL_max;
      cout << "  max: " << max_norm_allpatches << "  L2: " << l2_norm_allpatches;
      cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
      printTiming();

      ConfigMap config;
      config.addDirectory("control");
      cfl_target     = config.getValue<real>("CFL-number");
      write_interval = config.getValue<real>("write-interval")*time;
      total_time     = config.getValue<real>("total-time")*time;

      next_write_time = (int(t/write_interval) + 1)*write_interval;
      cout << "next output is at t=" << next_write_time/time << "*L/u_oo" << endl;

      dt *= cfl_target/CFL_max;

      ++write_counter;
      patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), write_counter);
      patch_grid.writeData(0, "data/step", t, write_counter);
    } else {
      ++iter;
      cout << iter << " iterations,  t=" << t << " = " << t/time << "*L/u_oo,  dt: " << dt << endl;
    }
  }

  stopTiming();
  cout << iter << " iterations" << endl;

#ifdef GPU
  runge_kutta.copyDonorData(0);
  iterator_std.updateHost();
#endif

  patch_grid.writeToVtk(0, "VTK/final", CompressibleVariables<PerfectGas>(), -1);
}

#endif // EXTERNAL_AERO_H
