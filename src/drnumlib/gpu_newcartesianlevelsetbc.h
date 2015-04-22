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

#ifndef GPU_CARTESIANLEVELSETBC_H
#define GPU_CARTESIANLEVELSETBC_H

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
class GPU_CartesianLevelSetBC;

#include "gpu_levelsetbc.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

#include "perfectgas.h"

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
class GPU_CartesianLevelSetBC : public GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>
{

public: // data types

  struct cell_t
  {
    size_t index;
    size_t src_n_weight[8];
    size_t src_index[8];
  };

protected: // attributes

  QVector<QList<cell_t> >   m_InsideCells;
  QVector<cell_t*>          m_GpuInsideCells;
  bool                      m_UpdateRequired;
  LS                        m_Ls;
  BC                        m_Bc;


protected: // methods

  void update();


public: // methods

  GPU_CartesianLevelSetBC(PatchGrid* patch_grid, LS ls, BC bc, int cuda_device = 0, size_t thread_limit = 0);

  template <typename T_Patch>
  CUDA_DH static void grad(T_Patch& patch, LS& ls,
                           size_t i, size_t j, size_t k,
                           real& gx, real& gy, real& gz)
  {
    if      (i == 0)                 gx = patch.idx()*(ls.G(patch, i+1, j, k) - ls.G(patch, i,   j, k));
    else if (i == patch.sizeI() - 1) gx = patch.idx()*(ls.G(patch, i  , j, k) - ls.G(patch, i-1, j, k));
    else                             gx = 0.5*patch.idx()*(ls.G(patch, i+1, j, k) - ls.G(patch, i-1, j, k));

    if      (j == 0)                 gy = patch.idy()*(ls.G(patch, i, j+1, k) - ls.G(patch, i, j  , k));
    else if (j == patch.sizeJ() - 1) gy = patch.idy()*(ls.G(patch, i, j  , k) - ls.G(patch, i, j-1, k));
    else                             gy = 0.5*patch.idy()*(ls.G(patch, i, j+1, k) - ls.G(patch, i, j-1, k));

    if      (k == 0)                 gz = patch.idz()*(ls.G(patch, i, j, k+1) - ls.G(patch, i, j, k  ));
    else if (k == patch.sizeK() - 1) gz = patch.idz()*(ls.G(patch, i, j, k  ) - ls.G(patch, i, j, k-1));
    else                             gz = 0.5*patch.idz()*(ls.G(patch, i, j, k+1) - ls.G(patch, i, j, k-1));
  }

  CUDA_HO virtual void operator()();

};


template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::GPU_CartesianLevelSetBC(PatchGrid* patch_grid, LS ls, BC bc, int cuda_device, size_t thread_limit)
  : GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>(patch_grid, cuda_device, thread_limit)
{
  m_Bc = bc;
  m_Ls = ls;
  m_UpdateRequired = true;
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    CartesianPatch& patch = *(this->m_Patches[i_patch]);
    for (size_t i = 0; i < patch.sizeI(); ++i) {
      for (size_t j = 0; j < patch.sizeJ(); ++j) {
        for (size_t k = 0; k < patch.sizeK(); ++k) {
          if (ls.G(patch, i, j, k) < 0) {
            patch.deactivate(i,j,k);
          }
        }
      }
    }
  }
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::update()
{
  if (!m_UpdateRequired) {
    return;
  }
  m_InsideCells.resize(this->m_Patches.size());
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    m_InsideCells[i_patch].clear();
    CartesianPatch* patch = this->m_Patches[i_patch];
    int imax = patch->sizeI();
    int jmax = patch->sizeJ();
    int kmax = patch->sizeK();

    for (size_t i = 0; i < imax; ++i) {
      for (size_t j = 0; j < jmax; ++j) {
        for (size_t k = 0; k < kmax; ++k) {
          if (m_Ls.G(*patch, i, j, k) < 0) {
            size_t index = patch->index(i, j, k);
            if (notInBcList(index)) {
              cell_t cell_i;
              cell_i.index = index;
              real gx, gy, gz;
              int count = 0;
              real total_weight = 0;
              this->grad(patch, ls, i, j, k, gx, gy, gz);
              int i_stt = -1;
              int i_end =  1;
              int j_stt = -1;
              int j_end =  1;
              int k_stt = -1;
              int k_end =  1;
              if (i == 0)         i_stt = 0;
              if (i == imax - 1)  i_end = 0;
              if (j == 0)         j_stt = 0;
              if (j == jmax - 1)  j_end = 0;
              if (k == 0)         k_stt = 0;
              if (k == kmax - 1)  k_end = 0;
              for (int di = i_stt; di <= i_end; ++di) {
                for (int dj = j_stt; dj <= j_end; ++dj) {
                  for (int dk = k_stt; dk <= k_end; ++dk) {
                    if (di != 0 || dj != 0 || dk != 0) {
                      if (ls.G(patch, i + di, j + dj, k + dk) > ls.G(patch, i, j, k)) {
                        int cell_n = patch->index(i + di, j + dj, k +dk);
                        real dx = di*patch->dx();
                        real dy = dj*patch->dy();
                        real dz = dk*patch->dz();
                        real weight = fabs(dx*gx + dy*gy + dz*gz)/sqrt(dx*dx + dy*dy + dz*dz);
                        cell_i.src_index[count]    = cell_n;
                        cell_i.src_n_weight[count] = weight;
                        total_weight += weight;
                        count++;
                      }
                    }
                  }
                }
              }
              if (count > 0) {
                for (size_t m = 0; m != count; ++m ) {
                  cell_i.src_n_weight[count] /= total_weight;
                }
              }
              m_InsideCells << cell_i;
            }
          }
        }
      }
    }
  }

  // delete old GPU arrays
  foreach (cell_t* cells, m_GpuInsideCells) {
    cudaFree(cells);
    CUDA_CHECK_ERROR;
  }
  m_GpuInsideCells.clear();

  // allocate new arrays
  m_GpuInsideCells.resize(this->m_Patches.size());
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    if (m_InsideCells[i_patch].size() > 0) {
      cudaMalloc(&m_GpuInsideCells[i_patch], m_InsideCells[i_patch].size()*sizeof(cell_t));
      CUDA_CHECK_ERROR;
      cell_t* cells = new cell_t[m_InsideCells[i_patch].size()];
      for (int i = 0; i < m_InsideCells[i_patch].size(); ++i) {
        cells[i] = m_InsideCells[i_patch][i];
      }
      cudaMemcpy(m_GpuInsideCells[i_patch], cells, m_InsideCells[i_patch].size()*sizeof(cell_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;
      delete [] cells;
    } else {
      m_GpuInsideCells[i_patch] = NULL;
    }
  }

  m_UpdateRequired = false;
}


template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernel(GPU_CartesianPatch patch, typename GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::bccell_t* cells, size_t num_cells, LS ls, BC bc)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= num_cells) {
    return;
  }

  size_t  cell_i   = cells[idx].index;
  real*   weight_d = cells[idx].octet_wt;
  int*    cell_d   = cells[idx].octet_id;


  dim_t<DIM> dim;
  real var_1[DIM];
  real p = 0;
  real T = 0;
  real u = 0;
  real v = 0;
  real w = 0;
  real p_d, T_d, u_d, v_d, w_d;
  for(size_t i = 0; i != 8; ++i) {
    patch.getVar(dim, 0, cell_d[i], var_1);
    PerfectGas::conservativeToPrimitive(var_1, p_d, T_d, u_d, v_d, w_d);

    u += u_d * weight_d[i];
    v += v_d * weight_d[i];
    w += w_d * weight_d[i];
    p += p_d * weight_d[i];
    T += T_d * weight_d[i];
  }

  // Needs mirroring for u,v,w.  What to do with p and T?
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetInsideCellsBC_kernel(GPU_CartesianPatch patch, typename GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::cell_t* cells, size_t num_cells)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= num_cells) {
    return;
  }

  size_t index   = cells[idx].index;
  int* src_index = cells[idx].src_index;
  int* n_weight  = cells[idx].src_n_weight;

  dim_t<DIM> dim;
  real var[DIM], var_bc[DIM];
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    var[i_var] = 0.0;
  }
  for (size_t n_bc = 0; n_bc != 8; ++n_bc) {
    patch.getVar(dim, 0, src_index[n_bc], var_bc);
    for (size_t i_var = 0; i_var < DIM - NUM_LS; ++i_var) {
      var[i_var] += n_weight[n_bc]*var_bc[i_var];
    }
  }
  // Euh.. Check this..
  patch.setVar(dim, 2, i, j, k, var);

}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::operator()()
{
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CUDA_CHECK_ERROR;

  update();
  size_t max_num_threads = this->m_MaxNumThreads;
  cudaDeviceSynchronize();

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    CUDA_CHECK_ERROR;


    if (m_InsideCells[i_patch].size() > 0) {
      int num_cells   = m_InsideCells[i_patch].size();
      int num_blocks  = max(int(16), int(num_cells/max_num_threads) + 1);
      int num_threads = num_cells/num_blocks + 1;
      if (num_cells > num_blocks*num_threads) BUG;
      cudaMemcpy(this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
      GPU_CartesianLevelSetInsideCellsBC_kernel<DIM,NUM_LS,LS,BC> <<<num_blocks, num_threads>>>(this->m_GpuPatches[i_patch], m_GpuInsideCells[i_patch], num_cells, m_Ls, m_Bc);
      CUDA_CHECK_ERROR;
      cudaDeviceSynchronize();
      cudaMemcpy(this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
    }

  }
}

#endif // GPU_CARTESIANLEVELSETBC_H
