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


// ============================

template <unsigned int SIZE, unsigned int DIM>
LSObject<SIZE, DIM>::LSObject(PatchGrid *patch_grid,
                              size_t g_extra_index,
                              size_t num_layers,
                              real min_inner_g_dist,
                              size_t trans_field_index)
{
  m_PatchGrid       = patch_grid;
  m_GExtraIndex     = g_extra_index;
  m_NumLayers       = num_layers;
  m_MinInnerGDist   = min_inner_g_dist;
  m_TransFieldIndex = trans_field_index;

  m_LSBCListOK      = false;
}


template <unsigned int SIZE, unsigned int DIM>
void LSObject<SIZE, DIM>::extractBCellLayers()
{
  // Find boundary cells affected and store in m_LSBCList

  // Do the following:
  // 0) Set up
  // 1) Loop over patches and cells
  //   1a) Build 1st inner layer (layer marker: -1)
  //   1b) Build further layers
  //   1c) Compute mirror points, interpolate and transfer cells found in this patch to layer_cells
  // 2) Transfer to array type data structure m_LSBCList, containing all access/miorror infos for this
  //    level set. Scatter into non-recursive vector groups, while ensuring, that the first group 0
  //    contains ALL dst-cells (all receiving cells) once. Reason: This avoids need to set dst-data
  //    to zero before gathering contributions (in hot zone!). Vectorization: No group is allowed to
  //    contain any dst-address more than once.
  // 3) Perform some simple checks on integrity and vectorization.

  // 0) Set up
  // Intermediate data storage as vector. Reason: Size unknown.
  vector<lsbc_cell_t<SIZE,DIM> > layer_cells;
  // Count cells with 1,2,3,... contributions
  //   example: count_contribs[1] will later be the number of dst-cells with 2 contributors
  //            count_contribs[i_cc] will later be the number of dst-cells with i_cc+1 contributors
  vector<size_t> count_contribs;

  // 1) Loop over patches and cells
  // Loop for all patches. Note 10: direct indexing on patches!
  /// @todo midterm: apply earlier search reduction to exclude patches w/o boundaries in early stage
  for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++) {
    Patch* patch = m_PatchGrid->getPatch(i_p);
    real* g_array = patch->getExtraCPUVarset(m_GExtraIndex);
    vector<size_t> ind_cell_neighbours;
    vector<int> cell_marker;  /// @todo there was a nice class scratch list in MOUSE: rewrite in stl ??
    cell_marker.resize(patch->variableSize(), 0);

    // Need an Eps for g to prevent precision hazards
    /// @todo better Eps handling needed (efficiency)
    real g_max = g_array[0];
    real g_min = g_array[0];
    for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
      if (g_min > g_array[l_c]) g_min = g_array[l_c];
      if (g_max < g_array[l_c]) g_max = g_array[l_c];
    }
    real g_eps = (g_max - g_min) * 1.e-7;

    // Exclude patches with G-values that are or only greater or only lower than threshold -m_MinInnerGDist
    /// @todo need a detection system for baffles later: no negative G-values on these, but local minima in direction of nabla(G)
    real test = (g_max + m_MinInnerGDist + g_eps) * (g_min + m_MinInnerGDist - g_eps);
    if(test < 0.) {  // patch containes a boundary

      // Find cells, that have at least one face neighbour with different sign of
      // the levelset function. These are 0th layer cells.
      /// @todo need a detection system for baffles ...
      //.. Loop for cells of patch
      //Postponed: building outer layers. Only inside at present
      //      for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
      //        real g = var[l_c];
      //        //.... check face neighbours
      //        patch->cellOverFaceNeighbours(l_c,
      //                                      ind_cell_neighbours);
      //        bool any_other_sign = false;
      //        for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
      //          size_t l_cn = ind_cell_neighbours[ll_cn];
      //          real test_1 = (var[l_cn] + g_eps) * (g - g_eps);
      //          real test_2 = (var[l_cn] - g_eps) * (g + g_eps);
      //          if (test_1 < 0. || test_2 < 0.) {  // boundary between cell neighbourship with
      //            any_other_sign = true;
      //            cell_marker[l_c] = true;
      //            break;
      //          }
      //        }
      //        if (any_other_sign) {
      //          LSLayerDataExtrapol lslde_h;
      //          lslde_h.m_Data = patch->getData();
      //          lslde_h.m_FieldSize = patch->fieldSize();
      //          lslde_h.m_VariableSize = patch->variableSize();
      //          lslde_h.m_Cell = l_c;
      //          lslde_h.m_G = g;
      //          if (g < 0.) {
      //            if(m_NumInnerLayers > 0) { // only if al least a 0th layer is requested
      //              // m_InnerCellsLayers[i_p][0].push_back(LSLayerDataExtrapol(l_c, g));
      //              m_InnerCellsLayers[i_p][0].push_back(lslde_h);
      //            }
      //          } else {
      //            if(m_NumOuterLayers > 0) { // only if al least a 0th layer is requested
      //              m_OuterCellsLayers[i_p][0].push_back(lslde_h);
      //            }
      //          }
      //        }
      //      }

      // layer_cells_pp[layer][cell]: layer cells per patch
      vector<vector<layer_cell_entry_pp> > layer_cells_pp;
      layer_cells_pp.resize(m_NumLayers);

      // 1a) Build 1st inner layer (layer marker: -1)
      if (m_NumLayers > 0) {
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          real g = g_array[l_c];
          if (g < -m_MinInnerGDist) {  // consider inside. Note 20: make it more prone to get found.
            //.... check face neighbours
            patch->cellOverFaceNeighbours(l_c,
                                          ind_cell_neighbours);
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cn = ind_cell_neighbours[ll_cn];
              real gn = g_array[l_cn];
              real test = (gn +m_MinInnerGDist + g_eps) * (g +m_MinInnerGDist - g_eps);
              if (test < 0.) {  // boundary between cell l_c and l_cn
                cell_marker[l_c] = -1;
                //                pair<size_t, int> layer_cell;
                //                layer_cell.first  = l_c; // cell index in patch
                //                layer_cell.second = -1;  // first inner cell layer
                //                layer_cells.push_back(layer_cell);
                layer_cell_entry_pp entry;
                entry.i_patch = i_p;
                entry.l_cell  = l_c;
                entry.g       = g;
                entry.layer   = -1;
                layer_cells_pp[0].push_back(entry);
                break;
              }
            }
          }
        }
      }

      // 1b) Build further layers
      /// @todo at present only towards inside region of objects
      for (size_t i_layer = 1; i_layer < m_NumLayers; i_layer++) {
        size_t i_layer_last = i_layer - 1;
        for (size_t ll_c = 0; ll_c < layer_cells_pp[i_layer_last].size(); ll_c++) {
          size_t l_c = layer_cells_pp[i_layer_last][ll_c].l_cell;
          real g     = layer_cells_pp[i_layer_last][ll_c].g;
          // check neighbours of these cells to find non marked ones, that have a lower g-value
          patch->cellOverFaceNeighbours(l_c,
                                        ind_cell_neighbours);
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cn = ind_cell_neighbours[ll_cn];
            if (cell_marker[l_cn] == 0) {  // cell has not been considered yet
              real gn = g_array[l_cn];
              if (gn < g) {  // cell l_cn is on the next inner cell layer
                cell_marker[l_cn] = -(i_layer + 1);
                layer_cell_entry_pp entry;
                entry.i_patch = i_p;
                entry.l_cell  = l_cn;
                entry.g       = gn;
                entry.layer   = -(i_layer + 1); // first inner layer: -1, second inner layer: -2 , ...
                layer_cells_pp[i_layer].push_back(entry);
              }
            }
          }
        }
      }

      // 1c) Compute mirror points, interpolate and transfer cells found in this patch to layer_cells
      // NOTE 25: Search restricted to own patch and its donor patches.
      // NOTE 30: If multiple patches contribute to a single dst-cell, this single dst-cell
      //          will receive a corresponding number of entries in layer_cells. Destination related
      //          data will be the same for these multiple entries (dst_data, dst_index, dst_step,
      //          layer, h, gx, gy, gz).
      // NOTE 40: Multiple contributing patches will receive consecutive positions in layer_cells.
      // NOTE 50: This intermediate step is required, as m_LSBCList.dst_cells (see below) must
      //          be allocated at correct size and the size is only known after processing all patches.

      // Max number of contribs for any dst-cell: number of vector groups required later.
      //size_t max_contribs = 0;


      for (size_t i_layer = 0; i_layer < m_NumLayers; i_layer++) {
        for (size_t ll_c = 0; ll_c < layer_cells_pp[i_layer].size(); ll_c++) {
          layer_cell_entry_pp lce_pp = layer_cells_pp[i_layer][ll_c];
          //.. Create a first new entry for layer_cells and fill in yet known qunatities
          lsbc_cell_t<SIZE,DIM> lc;
          lc.dst_data = patch->getField(m_TransFieldIndex);
          lc.dst_index = lce_pp.l_cell;
          lc.dst_step  = patch->variableSize();
          lc.layer     = lce_pp.layer;
          lc.h         = lce_pp.g;
          lc.transform_index = -1;  /// @todo Implement/debug transformations
          //.... Nabla G. Note 60: Coordinate system is own system
          patch->computeNablaArray(g_array, lc.dst_index,
                                   lc.gx, lc.gy, lc.gz);
          //.. Construct mirror point (on oposite side)
          vec3_t nablaxyz_g(lc.gx, lc.gy, lc.gz);
          real shift_len = -2. * lce_pp.g;
          vec3_t mirror_shift = shift_len * nablaxyz_g;
          real xc, yc, zc;
          patch->xyzCell(lce_pp.l_cell,
                         xc, yc, zc);
          vec3_t xyz_c(xc, yc, zc);
          vec3_t xyz_mirror = xyz_c + mirror_shift;
          //.. Interpolate
          vector<pair<Patch*, WeightedSet<real> > > all_ws;
          if (patch->interpolCoeffs4XYZRestrict(xyz_mirror[0], xyz_mirror[1], xyz_mirror[2],
                                                all_ws)) {
            for (size_t i_c = 0; i_c < all_ws.size(); i_c++) {
              lc.src_data = all_ws[i_c].first->getField(m_TransFieldIndex);
              all_ws[i_c].second.tranferToFixedArrays(SIZE,
                                                      &(lc.src_index[0]), &(lc.src_weight[0]));
              lc.src_step = all_ws[i_c].first->variableSize();
              layer_cells.push_back(lc);
            }
            if (count_contribs.size() < all_ws.size()) {
              count_contribs.resize(all_ws.size(),0);
            }
            count_contribs[all_ws.size()-1]++;  // it is sure that all_ws.size()>=1
          }
        }
      }
    } // if(test < 0.)
  }   // for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++)

  // 2) Transfer to array type data structure m_LSBCList, containing all access/miorror infos for this
  //    level set. Scatter into non-recursive vector groups, while ensuring, that the first group 0
  //    contains ALL dst-cells (all receiving cells) once. Reason: This avoids need to set dst-data
  //    to zero before gathering contributions (in hot zone!). Vectorization: No group is allowed to
  //    contain any dst-address twice.
  //    Consider NOTE 40: multiple contribs in direct consecution.
  //.. Size and internal list dst_cells
  m_LSBCList.dst_cells_size = layer_cells.size();
  m_LSBCList.dst_cells = new lsbc_cell_t<SIZE,DIM>[layer_cells.size()]; // make same size as layer_cells
  //.. Vector groups and insection indices
  size_t num_groups = count_contribs.size();
  m_LSBCList.num_groups = num_groups;
  m_LSBCList.group_start = new int[num_groups];
  m_LSBCList.group_size  = new int[num_groups];
  //.... Compute vector group sizes:
  for (size_t i_g = 0; i_g < num_groups; i_g++) {
    size_t group_size = 0;
    for (size_t i_g2 = i_g; i_g2 < num_groups; i_g2++) {
      group_size += count_contribs[i_g2];
    }
    m_LSBCList.group_size[i_g] = group_size;
  }
  //.... Compute start indices of vector groups (insection indices for m_LSBCList.dst_cells)
  size_t running_index = 0;
  for (size_t i_g = 0; i_g < num_groups; i_g++) {
    m_LSBCList.group_start[i_g] = running_index;
    running_index += m_LSBCList.group_size[i_g];
  }
  vector<size_t> group_counter;  // count entries in the vector groups to know insertion positions
  group_counter.resize(num_groups, 0);
  //.. Fill in data
  //   prime with very first entry. No recurrence check needed.
  if (layer_cells.size() > 0) {
    m_LSBCList.dst_cells[0] = layer_cells[0];  // put in first group 0
    group_counter[0]++;
  }
  //   all others
  size_t i_group = 0;
  for (size_t i_lc = 1; i_lc < layer_cells.size(); i_lc++) {
    //.... Check for recursion. On recursion: shift up to next vector group
    // if (layer_cells[i_lc].dst_data + layer_cells[i_lc].dst_index ==
    //     layer_cells[i_lc-1].dst_data + layer_cells[i_lc-1].dst_index) {
    if (layer_cells[i_lc].dst_data == layer_cells[i_lc-1].dst_data &&
        layer_cells[i_lc].dst_index == layer_cells[i_lc-1].dst_index) {  // same destination
      i_group++; // put this in the next vector group
    }
    else {  // other destination than previous
      i_group = 0;  // reset i_group to again start at first group 0
    }
    size_t entry_index = m_LSBCList.group_start[i_group] + group_counter[i_group];
    m_LSBCList.dst_cells[entry_index] = layer_cells[i_lc];
    group_counter[i_group]++;
  }

  // 3) Perform some simple checks on integrity and vectorization.
  //#ifdef DEBUG
  // 3a) Check all groups are filled:
  bool error_1 = false;
  for (size_t i_g = 0; i_g < num_groups; i_g++) {
    if (group_counter[i_g] != m_LSBCList.group_size[i_g]) {
      error_1 = true;
      //BUG;
    }
  }

  // 3b) Check vectorization (uniqueness) and if group 0 contains all dst-cells
  vector<real*> addresses;
  bool no_error = true;
  for (size_t i_g = 0; i_g < m_LSBCList.num_groups; i_g++) {  // loop groups
    addresses.clear();
    // insert all left side mem addresses og group i_g (maps 1:1 to device)
    for (size_t i_c = m_LSBCList.group_start[i_g]; i_c < m_LSBCList.group_size[i_g]; i_c++) {
      lsbc_cell_t<SIZE,DIM> contrib = m_LSBCList.dst_cells[i_c];
      real* address = contrib.dst_data + size_t(contrib.dst_index);  // ptr-shift
      addresses.push_back(address);
    }
    // Compare raw size and unique size: if unequal, multiple same left side addresses were present
    // (equivalent to a vectoriztation error)
    sort(addresses.begin(), addresses.end());
    size_t raw_size = addresses.size();
    vector<real*>::iterator it = unique(addresses.begin(), addresses.end());
    addresses.resize(it - addresses.begin());
    size_t unique_size = addresses.size();
    if(raw_size != unique_size) {
      cout << "Error, vector group " << i_g << " with multiple left side addresses." << endl;
      no_error = false;
    }
    // Continue checking, if first group 0 contains all entries
    if (i_g == 0) {  // only for first group
      if (m_LSBCList.num_groups > 1) {  // only meaningful if more than a single group exists
        size_t unique_size_0 = unique_size;
        // append all entries of all subsequent groups
        for (size_t i_c = m_LSBCList.group_start[0]+m_LSBCList.group_size[0]; i_c < m_LSBCList.dst_cells_size; i_c++) {
        //for (size_t i_c = m_LSBCList.group_start[1]; i_c < m_LSBCList.dst_cells_size; i_c++) {
          lsbc_cell_t<SIZE,DIM> contrib = m_LSBCList.dst_cells[i_c];
          real* address = contrib.dst_data + size_t(contrib.dst_index);  // ptr-shift
          addresses.push_back(address);
        }
        // sort and make unique, then check size again: If size is bigger than unique_size_0, there
        // were left side addresses in groups of higher index, that did not exist in the first group 0
        // This is a concat initializer error
        sort(addresses.begin(), addresses.end());
        it = unique(addresses.begin(), addresses.end());
        addresses.resize(it - addresses.begin());
        unique_size = addresses.size();
        if (unique_size > unique_size_0) {
          cout << "Concat-init-error: vector group 0 did not contain all left side addresses" << endl;
          no_error = false;
        }
      }
    }
  }
  if(!no_error) {
    BUG;
  }
  //#endif

  m_LSBCListOK = true;  // assume that getting here means its fine
}


template <unsigned int SIZE, unsigned int DIM>
lsbc_list_t<SIZE, DIM> LSObject<SIZE, DIM>::getLSBCList()
{
  if (!m_LSBCListOK) {
    extractBCellLayers();
  }
  return m_LSBCList;
}


template <unsigned int SIZE, unsigned int DIM>
lsbc_list_t<SIZE, DIM>* LSObject<SIZE, DIM>::getLSBCListPtr()
{
  if (!m_LSBCListOK) {
    extractBCellLayers();
  }
  return &m_LSBCList;
}

#endif // LSOBJECT_H
