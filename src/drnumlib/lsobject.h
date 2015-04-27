#ifndef LSOBJECT_H
#define LSOBJECT_H

template <unsigned int SIZE, unsigned int DIM>
class LSObject;

#include "patchgrid.h"
#include "lsbc_list_t.h"
#include "weightedset.h"

template <unsigned int SIZE, unsigned int DIM>
class LSObject
{
  PatchGrid* m_PatchGrid;
  size_t m_GExtraIndex;     ///< index of extra data, starting ptr m_ExtraCPUData + variableSize() * m_GExtraIndex
  size_t m_NumLayers;       ///< total number of layers to create: prefer inner layers, but may deviate to outwards
  real m_MinInnerGDist;     ///< minimum distance of center of an inner cell to surface allowed
  real m_MinOuterGDist;     ///< minimum distance of center of an outer cell to surface allowed (needed??)
  size_t m_TransFieldIndex; ///< index of variable field for which to access/extrapolte data for BC construction

  /** Interpolation list */
  lsbc_list_t<SIZE, DIM> m_LSBCList;   /// @todo make "this" a template class too to pass over number of vars

  bool m_LSBCListOK;        ///< bool indicating, m_LSBCList is available.

  /// @todo m_MinInnerGDist may better rely on relative G-values: thresholds definedrelative to cell sizes.

  /**
   * Intermediate struct for per-patch layer data */
  struct layer_cell_entry_pp
  {
    size_t i_patch;  // patch index in m_PatchGrid
    size_t l_cell;   // 1D cell index in patch
    real   g;        // g-value (equiv. to wall distance, inside negative)
    int    layer;    // layer index
  };

public:
  /**
   * @brief LSObject
   * @param patch_grid PatchGrid to work on
   * @param g_extra_index Index of extra field, on which levelset info is stored
   * @param num_layers numer of layers to create
   * @param min_inner_g_dist minimum g-distance of inner layer cells from surface (try 0.)
   * @param trans_field_index field index, on which later access/mirror ops are performed
   */
  LSObject(PatchGrid *patch_grid,
           size_t g_extra_index,
           size_t num_layers,
           real min_inner_g_dist,
           size_t trans_field_index);

  /**
   * @brief Build up layer cell access/mirror info for the present lavel set.
   */
  void extractBCellLayers();

  /**
   * @brief getLSBCList Access lsbc_list_t produced by this::extractBCellLayers() for external storage
   * @return lsbc_list_t ...
   */
  lsbc_list_t<SIZE, DIM> getLSBCList() ;  /// @todo check: will the std-copy constructoro do the newing job??


  /**
   * @brief getLSBCListPtr Access pointer to lsbc_list_t produced by this::extractBCellLayers()
   * and stored permanently.
   * @return pointer to lsbc_list_t ...
   */
  lsbc_list_t<SIZE, DIM>* getLSBCListPtr();
};

#include "lsobject.cpp"

#endif // LSOBJECT_H
