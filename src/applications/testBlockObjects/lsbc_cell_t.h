#ifndef LSBC_CELL_T_H
#define LSBC_CELL_T_H

template <unsigned int SIZE, unsigned int DIM>
struct lsbc_cell_t;

#include "drnum.h"

/** Struct to define access points per receiving cell of an LSObject boundary simulation.
  * One or more entries of type "lsbc_cell_t" will be stored per receiving cell.
  *  - One entry:
  *     Single unique donor patch (no donor overlaps) and sufficiently large SIZE for all
  *     donor cell contributions.
  *  - More than one entry:
  *     a) access point in overlap of donor zones of multiple patches and/or
  *     b) insufficient SIZE to cope with all contributing donor cells
  * SIZE: Number of contributing interpolation donors.
  *       If SIZE is not sufficient, annother
  * DIM:  Actual number of variables minus the number of variables which are used for
  *       level sets themselves.
  * For SIZE=8 and DIM=5 it consumes 176 bytes for real=float and 244 bytes for
  * real=double.
*/
template <unsigned int SIZE, unsigned int DIM>
struct lsbc_cell_t
{
  real*  src_data;          ///< pointer to the source data field (m_Data + i_field * fieldsize)
  int    src_index[SIZE];   ///< contributing cell indices
  // int    src_step[SIZE];    ///< step size of the source data (var[i_var] = src_data[...][src_step[...]*i_var])
  int    src_step;          ///< step size of the source data (var_ptr[i_var] = src_data + src_step*i_var)
  real   src_weight[SIZE];  ///< weights of the contributing cells
  real*  dst_data;          ///< pointer to the destination data origin
  int    dst_index;
  int    dst_step;          ///< step size of the destination data
  real   tmp_data[DIM];     ///< field to intermediately store the interpolation result
  real   gx, gy, gz;        ///< gradient of G
  real   h;                 ///< wall distance
  int    transform_index;   ///< index in the field containing the transformations (-1 -> no transform)
  int    layer;             ///< cell layer relative to boundary. Inside:-1,-2,... , outside:1,2,... , false marker: 0
};

#endif // LSBC_CELL_T_H
