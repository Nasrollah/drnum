#ifndef LSBC_LIST_T_H
#define LSBC_LIST_T_H

template <unsigned int SIZE, unsigned int DIM>
struct lsbc_list_t;

#include "drnum.h"
#include "lsbc_cell_t.h"

/** @todo Check, if it wouldnt be smarter to have it layer by layer. */

/** The following data structure represents a list which contains all required information for
  * LSObject based boundary simulation and exists once per computation. */
template <unsigned int SIZE, unsigned int DIM>
struct lsbc_list_t
{
 int                     dst_cells_size;  ///< total number of destination cells
 int                     num_groups;
 lsbc_cell_t<SIZE,DIM>*  dst_cells;       ///< array with the interpolation details (length of dst_cells_size)
 int*                    group_start;
 int*                    group_size;
};

#endif // LSBC_LIST_T_H
