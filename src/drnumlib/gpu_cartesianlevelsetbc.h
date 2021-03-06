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

template <unsigned int SIZE, unsigned int DIM>
class GPU_CartesianLevelSetBC;

#include "gpu_levelsetbc.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"
#include "lsobject.h"
#include "lsbc_cell_t.h"
#include "discretelevelset.h"

#include "perfectgas.h"


template <unsigned int SIZE, unsigned int DIM>
class GPU_CartesianLevelSetBC : public GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>
{

public: // data types

  struct lsextp_cell_t
  {
    int  index;
    //real src_weight[SIZE];
    //int  src_index[SIZE];
    real src_weight[26];
    int  src_index[26];
    real tmp_data[DIM];

  };


protected: // attributes

  QVector<QList<lsextp_cell_t> > m_InsideCells;
  QVector<lsextp_cell_t*>        m_GpuInsideCells;
  QSet<size_t>                   m_BcIndex;

  lsbc_list_t<SIZE,DIM>  m_BcCells;
  lsbc_cell_t<SIZE,DIM>* m_GpuBcCells;

  bool m_UpdateRequired;
  DiscreteLevelSet& m_Dls;


protected: // methods

  void update();
  void grad(CartesianPatch* patch, DiscreteLevelSet& dls, size_t i, size_t j, size_t k, real& gx, real& gy, real& gz);


public: // methods

  GPU_CartesianLevelSetBC(PatchGrid* patch_grid, DiscreteLevelSet& dls, lsbc_list_t<SIZE,DIM>* lsbc_list, int cuda_device = 0, size_t thread_limit = 0);


  CUDA_HO virtual void operator()();

};

template <unsigned int SIZE, unsigned int DIM>
void GPU_CartesianLevelSetBC<SIZE,DIM>::grad(CartesianPatch* patch, DiscreteLevelSet& dls, size_t i, size_t j, size_t k, real& gx, real& gy, real& gz)
{
  if      (i == 0)                  gx =     patch->idx()*(dls.G(patch, i+1, j, k) - dls.G(patch, i,   j, k));
  else if (i == patch->sizeI() - 1) gx =     patch->idx()*(dls.G(patch, i  , j, k) - dls.G(patch, i-1, j, k));
  else                              gx = 0.5*patch->idx()*(dls.G(patch, i+1, j, k) - dls.G(patch, i-1, j, k));

  if      (j == 0)                  gy =     patch->idy()*(dls.G(patch, i, j+1, k) - dls.G(patch, i, j  , k));
  else if (j == patch->sizeJ() - 1) gy =     patch->idy()*(dls.G(patch, i, j  , k) - dls.G(patch, i, j-1, k));
  else                              gy = 0.5*patch->idy()*(dls.G(patch, i, j+1, k) - dls.G(patch, i, j-1, k));

  if      (k == 0)                  gz =     patch->idz()*(dls.G(patch, i, j, k+1) - dls.G(patch, i, j, k  ));
  else if (k == patch->sizeK() - 1) gz =     patch->idz()*(dls.G(patch, i, j, k  ) - dls.G(patch, i, j, k-1));
  else                              gz = 0.5*patch->idz()*(dls.G(patch, i, j, k+1) - dls.G(patch, i, j, k-1));
}

template <unsigned int SIZE, unsigned int DIM>
GPU_CartesianLevelSetBC<SIZE,DIM>::GPU_CartesianLevelSetBC(PatchGrid* patch_grid, DiscreteLevelSet& dls, lsbc_list_t<SIZE,DIM>* lsbc_list, int cuda_device, size_t thread_limit)
  : GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>(patch_grid, cuda_device, thread_limit), m_Dls(dls)
{
  m_UpdateRequired = true;
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    CartesianPatch& patch = *(this->m_Patches[i_patch]);
    for (size_t i = 0; i < patch.sizeI(); ++i) {
      for (size_t j = 0; j < patch.sizeJ(); ++j) {
        for (size_t k = 0; k < patch.sizeK(); ++k) {
          if (dls.G(&patch, i, j, k) < 0) {
            patch.deactivate(i,j,k);
          }
        }
      }
    }
  }
  cout << "Number of BC cells->" << lsbc_list->dst_cells_size << endl;

  // Create map of cpu pointers to gpu pointers
  QMap<real*,real*> cpu2gpu;
  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    if (!this->m_Patches[i_patch]->gpuDataSet()) BUG;
    cpu2gpu.insert( this->m_Patches[i_patch]->getData(), this->m_Patches[i_patch]->getGpuData() );
  }

  // Deep copy elements of lsbc_list
  m_BcCells.dst_cells_size = lsbc_list->dst_cells_size;
  m_BcCells.num_groups     = lsbc_list->num_groups;
  m_BcCells.dst_cells      = new lsbc_cell_t<SIZE,DIM>[lsbc_list->dst_cells_size];
  m_BcCells.group_start    = new int[lsbc_list->num_groups];
  m_BcCells.group_size     = new int[lsbc_list->num_groups];
  for(int i_grp = 0; i_grp < m_BcCells.num_groups; ++i_grp) {
    m_BcCells.group_start[i_grp] = lsbc_list->group_start[i_grp];
    m_BcCells.group_size[i_grp]  = lsbc_list->group_size[i_grp];
  }


  // 1. Map cpu pointers to gpu pointers
  // 2. Deep copy dst_cells data
  // 3. Append index of cell to m_BcIndex list
  for (int c_i = 0; c_i < m_BcCells.dst_cells_size; ++c_i) {
    lsbc_cell_t<SIZE,DIM>& cell_i = m_BcCells.dst_cells[c_i];

    cell_i.src_data = cpu2gpu.value( lsbc_list->dst_cells[c_i].src_data );
    cell_i.dst_data = cpu2gpu.value( lsbc_list->dst_cells[c_i].dst_data );
    for (int i = 0; i < SIZE; ++i) {
      cell_i.src_index[i]  = lsbc_list->dst_cells[c_i].src_index[i];
      cell_i.src_weight[i] = lsbc_list->dst_cells[c_i].src_weight[i];
    }
    cell_i.src_step  = lsbc_list->dst_cells[c_i].src_step;
    cell_i.dst_index = lsbc_list->dst_cells[c_i].dst_index;
    cell_i.dst_step  = lsbc_list->dst_cells[c_i].dst_step;
    cell_i.gx        = lsbc_list->dst_cells[c_i].gx;
    cell_i.gy        = lsbc_list->dst_cells[c_i].gy;
    cell_i.gz        = lsbc_list->dst_cells[c_i].gz;
    cell_i.h         = lsbc_list->dst_cells[c_i].h;
    cell_i.layer     = lsbc_list->dst_cells[c_i].layer;
    cell_i.transform_index = lsbc_list->dst_cells[c_i].transform_index;

    m_BcIndex.insert( cell_i.dst_index );
  }

  // Create new gpu pointer storage variable
  m_GpuBcCells = new lsbc_cell_t<SIZE,DIM>[m_BcCells.dst_cells_size];

  // Allocate new arrays
  cudaMalloc(&m_GpuBcCells, m_BcCells.dst_cells_size*sizeof(lsbc_cell_t<SIZE,DIM>));
  CUDA_CHECK_ERROR;

  // Copy to GPU
  cudaMemcpy(m_GpuBcCells,  m_BcCells.dst_cells, m_BcCells.dst_cells_size*sizeof(lsbc_cell_t<SIZE,DIM>), cudaMemcpyHostToDevice);
  CUDA_CHECK_ERROR;
}

template <unsigned int SIZE, unsigned int DIM>
void GPU_CartesianLevelSetBC<SIZE,DIM>::update()
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
          if (m_Dls.G(patch, i, j, k) < 0) {
            size_t index = patch->index(i, j, k);
            if (!m_BcIndex.contains(index)) {
              lsextp_cell_t cell_i;
              cell_i.index = index;
              //To Improve!!
              for(int c_n = 0; c_n < 26; ++c_n) {
                cell_i.src_weight[c_n] = 0;
              }
              real gx, gy, gz;
              int count = 0;
              real total_weight = 0;
              this->grad(patch, m_Dls, i, j, k, gx, gy, gz);
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
                      if (m_Dls.G(patch, i + di, j + dj, k + dk) > m_Dls.G(patch, i, j, k)) {
                        if (count >= 26) {
                          cout << " cout->" << count << endl;
                          cout << " SIZE->" << SIZE  << endl;
                          BUG;
                        }
                        int cell_n = patch->index(i + di, j + dj, k +dk);
                        real dx = di*patch->dx();
                        real dy = dj*patch->dy();
                        real dz = dk*patch->dz();
                        real weight = fabs(dx*gx + dy*gy + dz*gz)/sqrt(dx*dx + dy*dy + dz*dz);
                        cell_i.src_index[count]  = cell_n;
                        cell_i.src_weight[count] = weight;
                        total_weight += weight;
                        count++;
                      }
                    }
                  }
                }
              }
              if (count > 0) {
                for (size_t m = 0; m != count; ++m ) {
                  cell_i.src_weight[count] /= total_weight;
                }
                m_InsideCells[i_patch] << cell_i;
              }
            }
          }
        }
      }
    }
  }

  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    cout << "Patch->" << i_patch << ", number of Inside cells->" << m_InsideCells[i_patch].size() << endl;
  }

  // delete old GPU arrays
  foreach (lsextp_cell_t* cells, m_GpuInsideCells) {
    cudaFree(cells);
    CUDA_CHECK_ERROR;
  }
  m_GpuInsideCells.clear();

  // allocate new arrays
  m_GpuInsideCells.resize(this->m_Patches.size());
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    if (m_InsideCells[i_patch].size() > 0) {
      cudaMalloc(&m_GpuInsideCells[i_patch], m_InsideCells[i_patch].size()*sizeof(lsextp_cell_t));
      CUDA_CHECK_ERROR;
      lsextp_cell_t* cells = new lsextp_cell_t[m_InsideCells[i_patch].size()];
      for (int i = 0; i < m_InsideCells[i_patch].size(); ++i) {
        cells[i] = m_InsideCells[i_patch][i];
      }
      cudaMemcpy(m_GpuInsideCells[i_patch], cells, m_InsideCells[i_patch].size()*sizeof(lsextp_cell_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;
      delete [] cells;
    } else {
      m_GpuInsideCells[i_patch] = NULL;
    }
  }

  m_UpdateRequired = false;
}

CUDA_DH inline real norm(real x, real y, real z)
{
  return sqrt(x*x + y*y + z*z);
}

CUDA_DH inline real dot(real x1, real y1, real z1, real x2, real y2, real z2)
{
  return (x1*x2 + y1*y2 + z1*z2);
}

template <unsigned int SIZE, unsigned int DIM>
__global__ void GPU_CartesianLevelSetComputeInsideCells_kernel(GPU_CartesianPatch patch, typename GPU_CartesianLevelSetBC<SIZE,DIM>::lsextp_cell_t* cells, size_t num_cells)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= num_cells) {
    return;
  }

  typename GPU_CartesianLevelSetBC<SIZE,DIM>::lsextp_cell_t& cell_i = cells[idx];
  dim_t<DIM> dim;
  real var[DIM], src_var[DIM];
  // Allocate
  for (size_t i_var = 0; i_var != DIM; ++i_var) {
    var[i_var] = 0;
  }

  // Compute
  // Warning!  The weight should have been normalize by the total weight beforehand
  // for (size_t n_src = 0; n_src < SIZE; ++n_src) {
  // Hardcoded to 26, the number of total useable points of the rubik's cube
  for (size_t n_src = 0; n_src < 26; ++n_src) {
    if (cell_i.src_weight[n_src] > 0) {
      patch.getVar(dim, 0, cell_i.src_index[n_src], src_var);
      for (size_t i_var = 0; i_var < DIM; ++i_var) {
        var[i_var] += cell_i.src_weight[n_src]*src_var[i_var];
      }
    }
  }

  // 1. Store
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    cell_i.tmp_data[i_var] = var[i_var];
  }
}

template <unsigned int SIZE, unsigned int DIM>
__global__ void GPU_CartesianLevelSetWriteInsideCells_kernel(GPU_CartesianPatch patch, typename GPU_CartesianLevelSetBC<SIZE,DIM>::lsextp_cell_t* cells, size_t num_cells)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= num_cells) {
    return;
  }

  int index = cells[idx].index;
  real var[DIM];
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    var[i_var] = cells[idx].tmp_data[i_var];
  }
  dim_t<DIM> dim;
  patch.GPU_Patch::setVar(dim, 0, index, var);
}

template <unsigned int SIZE, unsigned int DIM>
__global__ void GPU_CartesianLevelSetComputeBC_kernel(lsbc_cell_t<SIZE,DIM>* cells, int dst_cells_size)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= dst_cells_size) {
    return;
  }
  lsbc_cell_t<SIZE,DIM>& cell_i = cells[idx];

  for(size_t i_var = 0; i_var < DIM; ++i_var) {
    cell_i.tmp_data[i_var] = cell_i.src_weight[0] * cell_i.src_data[i_var*cell_i.src_step + cell_i.src_index[0]];
    for(size_t cell_n = 1; cell_n < SIZE; ++cell_n) {
      cell_i.tmp_data[i_var] += cell_i.src_weight[cell_n] * cell_i.src_data[i_var*cell_i.src_step + cell_i.src_index[cell_n]];
    }
  }
}

template <unsigned int SIZE, unsigned int DIM>
__global__ void GPU_CartesianLevelSetWriteBCCells_kernel(lsbc_cell_t<SIZE,DIM>* cells, int grp_nb, int grp_start, int grp_size)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= grp_size) {
    return;
  }

  lsbc_cell_t<SIZE,DIM>& cell_i = cells[idx + grp_start];

  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    if (grp_nb == 0) {
      cell_i.dst_data[i_var*cell_i.dst_step + cell_i.dst_index] = cell_i.tmp_data[i_var];
    }
    else {
      cell_i.dst_data[i_var*cell_i.dst_step + cell_i.dst_index] += cell_i.tmp_data[i_var];
    }
  }
}

template <unsigned int SIZE, unsigned int DIM>
__global__ void GPU_CartesianLevelSetMirrorBCCells_kernel(lsbc_cell_t<SIZE,DIM>* cells, int grp_size)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= grp_size) {
    return;
  }
  lsbc_cell_t<SIZE,DIM>& cell_i = cells[idx];

  real var[DIM];
  for(size_t i_var = 0; i_var < DIM; ++i_var) {
    var[i_var] = cell_i.dst_data[i_var*cell_i.dst_step + cell_i.dst_index];
  }

  real p, T, u, v, w;
  PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
  real n[3];
  n[0] = cell_i.gx;
  n[1] = cell_i.gy;
  n[2] = cell_i.gz;
  real norm_g = norm(n[0], n[1], n[2]);
  n[0] /= norm_g;
  n[1] /= norm_g;
  n[2] /= norm_g;

  real dot_n = dot(u, v, w, n[0], n[1], n[2]);
  //u -= 2*n[0]*dot_n;
  //v -= 2*n[1]*dot_n;
  //w -= 2*n[2]*dot_n;
  u -= 1.8*n[0]*dot_n;
  v -= 1.8*n[1]*dot_n;
  w -= 1.8*n[2]*dot_n;
  PerfectGas::primitiveToConservative(p, T, u, v, w, var);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    cell_i.dst_data[i_var*cell_i.dst_step + cell_i.dst_index] = var[i_var];
  }
}


template <unsigned int SIZE, unsigned int DIM>
void GPU_CartesianLevelSetBC<SIZE,DIM>::operator()()
{
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CUDA_CHECK_ERROR;

  update();
  size_t max_num_threads = this->m_MaxNumThreads;
  cudaDeviceSynchronize();

  CUDA_CHECK_ERROR;

  if (m_BcCells.dst_cells_size > 0) {
    int num_cells   = m_BcCells.dst_cells_size;
    int num_blocks  = max(int(16), int(num_cells/max_num_threads) + 1);
    int num_threads = num_cells/num_blocks + 1;
    if (num_cells > num_blocks*num_threads) BUG;

    GPU_CartesianLevelSetComputeBC_kernel<SIZE,DIM> <<<num_blocks, num_threads>>>(m_GpuBcCells, m_BcCells.dst_cells_size);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR;

    for(int i_grp = 0; i_grp < m_BcCells.num_groups; ++i_grp) {
      GPU_CartesianLevelSetWriteBCCells_kernel<SIZE,DIM> <<<num_blocks, num_threads>>>(m_GpuBcCells, i_grp, m_BcCells.group_start[i_grp], m_BcCells.group_size[i_grp]);
      cudaDeviceSynchronize();
      CUDA_CHECK_ERROR;
    }

    GPU_CartesianLevelSetMirrorBCCells_kernel<SIZE,DIM> <<<num_blocks, num_threads>>>(m_GpuBcCells, m_BcCells.group_size[0]);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR;

  }

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    CUDA_CHECK_ERROR;

    if (m_InsideCells[i_patch].size() > 0) {
      int num_cells   = m_InsideCells[i_patch].size();
      int num_blocks  = max(int(16), int(num_cells/max_num_threads) + 1);
      int num_threads = num_cells/num_blocks + 1;
      if (num_cells > num_blocks*num_threads) BUG;

      //GPU_CartesianLevelSetComputeInsideCells_kernel<SIZE,DIM> <<<num_blocks, num_threads>>>(this->m_GpuPatches[i_patch], m_GpuInsideCells[i_patch], num_cells);
      //cudaDeviceSynchronize();
      //CUDA_CHECK_ERROR;

      //GPU_CartesianLevelSetWriteInsideCells_kernel<SIZE,DIM> <<<num_blocks, num_threads>>>(this->m_GpuPatches[i_patch], m_GpuInsideCells[i_patch], num_cells);
      //cudaDeviceSynchronize();
      //CUDA_CHECK_ERROR;
    }

  }
}


#endif // GPU_CARTESIANLEVELSETBC_H
