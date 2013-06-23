#include "patch.h"

Patch::Patch(size_t num_seeklayers, size_t num_addprotectlayers)
{
  m_Data = NULL;
  m_NumFields = 0;
  m_NumVariables = 0;
  m_VariableSize = 0;
  m_FieldSize = 0;
  m_SeekExceptions = false;
  m_NumAddProtectLayers = num_addprotectlayers;
  m_NumSeekLayers = num_seeklayers;
  m_InterpolateData = false;
  //  m_InterpolateGrad1N = false;
  m_receiveCells_OK = false;
  m_bbox_OK = false;

  m_NumDonorPatches = 0;
  m_NumReceivingCellsConcat = 0;
  m_NumReceivingCellsUnique = 0;
  m_NumDonorWIConcat = 0;
}

Patch::~Patch()
{
  deleteData();
}

bool Patch::readFromFile(istringstream& iss_input)
{
  vec3_t xyzoref;        // reference point in parental coords
  vec3_t base_i, base_j; // base vectors of bloc orientation in parental coords.
  iss_input >> xyzoref[0];
  iss_input >> xyzoref[1];
  iss_input >> xyzoref[2];
  iss_input >> base_i[0];
  iss_input >> base_i[1];
  iss_input >> base_i[2];
  iss_input >> base_j[0];
  iss_input >> base_j[1];
  iss_input >> base_j[2];
  iss_input >> m_ioscale;
  setupTransformation(xyzoref,
                      base_i, base_j);
  /// @todo check before returning "true"
  return true;
}


void Patch::buildDonorTransferData()
{
  // Number of donor patches influencing this (receiver) patch
  m_NumDonorPatches = m_InterCoeffData.size();

  // Number of concatenated cells receiving data from any of all donors (multiple indexing)
  m_NumReceivingCellsConcat = 0;
  for (size_t i_donor = 0; i_donor < m_NumDonorPatches; i_donor++) {
    InterCoeffPad* icd = &(m_InterCoeffData[i_donor]);
    m_NumReceivingCellsConcat += icd->m_NumRecCells;
  }

  // Build m_ReceivingCellIndiceesConcat
  m_ReceivingCellIndicesConcat = new size_t[m_NumReceivingCellsConcat];
  size_t count = 0;
  for (size_t i_donor = 0; i_donor < m_NumDonorPatches; i_donor++) {
    InterCoeffPad* icd = &(m_InterCoeffData[i_donor]);
    for (size_t i_rec=0; i_rec < icd->m_NumRecCells; i_rec++) {
      m_ReceivingCellIndicesConcat[count] = icd->m_RecCells[i_rec];
      count++;
    }
  }
#ifdef DEBUG
  if (count != m_NumReceivingCellsConcat) {
    BUG;
  }
#endif

  // Build unique index field m_ReceivingCellIndexUnique
  // Note: use data array stored in m_receive_cells at present.
  //       May also unify m_ReceivingCellIndexConcat .
  m_NumReceivingCellsUnique = m_receive_cells.size();
  m_ReceivingCellIndicesUnique = new size_t[m_NumReceivingCellsUnique];
  for (size_t i_rec_all = 0; i_rec_all < m_receive_cells.size(); i_rec_all++) {
    m_ReceivingCellIndicesUnique[i_rec_all] = m_receive_cells[i_rec_all];
  }

  // Build m_Donors
  m_Donors = new donor_t[m_NumDonorPatches];
  size_t count_rec_concat = 0;
  size_t count_donor_concat = 0;
  for (size_t i_donor = 0; i_donor < m_NumDonorPatches; i_donor++) {
    Patch* donor_patch = m_neighbours[i_donor].first;
    InterCoeffPad* icd = &(m_InterCoeffData[i_donor]);
    m_Donors[i_donor].variable_size = donor_patch->m_VariableSize;
    m_Donors[i_donor].data = donor_patch->m_Data;
    m_Donors[i_donor].num_receiver_cells = icd->m_NumRecCells;
    m_Donors[i_donor].stride = icd->m_StrideGivePerRec;
    m_Donors[i_donor].receiver_index_field_start = count_rec_concat;
    m_Donors[i_donor].donor_wi_field_start = count_donor_concat;
    count_rec_concat += icd->m_NumRecCells;
    count_donor_concat += icd->m_NumRecCells * icd->m_StrideGivePerRec;
  }

  m_NumDonorWIConcat = count_donor_concat;
#ifdef DEBUG
  if (count_rec_concat != m_NumReceivingCellsConcat) {
    BUG;
  }
#endif
  // Build concatenated donor index and weight fields
  m_DonorIndexConcat = new size_t[m_NumDonorWIConcat];
  m_DonorWeightConcat = new real[m_NumDonorWIConcat];
  size_t count_concat = 0;
  for (size_t i_donor = 0; i_donor < m_NumDonorPatches; i_donor++) {
    InterCoeffPad* icd = &(m_InterCoeffData[i_donor]);
    size_t count_in_donor = 0;
    for (size_t i_rec = 0; i_rec < icd->m_NumRecCells; i_rec++) {
      for (size_t i_s = 0; i_s < icd->m_StrideGivePerRec; i_s++) {
        m_DonorIndexConcat[count_concat] = icd->m_DonorCells[count_in_donor];
        m_DonorWeightConcat[count_concat] = icd->m_DonorWeights[count_in_donor];
        count_concat++;
        count_in_donor++;
      }
    }
  }
#ifdef DEBUG
  if (count_concat != m_NumDonorWIConcat) {
    BUG;
  }
#endif
};


void Patch::readSolverCodes(istringstream& iss_input)
{
  m_solvercodes.buildFrom(iss_input);
}


void Patch::setupTransformation(vec3_t xyzoref,
                                vec3_t base_i, vec3_t base_j)
{
  // Build up transformation matrix
  //.. build preliminary transformation this --> inertial
  //   transforms point or freevec given in coord-syst of "this" to
  //   parental coord-syst.
  /// @todo needs testing
  CoordTransform this2inertial;
  this2inertial.setVector(xyzoref);
  this2inertial.setMatrixFromBaseIJ(base_i, base_j);
  //.. build transformation matrix inertial --> this
  m_transformInertial2This.setAll(this2inertial.inverse());

  Transformation t;    /// @todo for compatibility only. Get rid if no longer needed
  t.setVector(xyzoref);
  setTransformation(t.inverse());
}

void Patch::insertNeighbour(Patch* neighbour_patch) {
  // Check, if boundary cells are yet extracted. If not, do so.
  if(!m_receiveCells_OK) {
    m_receive_cells.clear();
    //NEW_SEEK_EXCEPTION extractReceiveCells();
    buildRegions();  /// @todo need better method call handling, when to do what
    extractSeekCells();
    compactReceiveCellLists();
    m_receiveCells_OK = true;
  }
  // Insert donor patch and relative coord transformation into m_neighbours and get neighbourship index
  pair<Patch*, CoordTransformVV> dependency_h;
  dependency_h.first = neighbour_patch;
  CoordTransformVV ctvv_relative;
  ctvv_relative.setTransFromTo(m_transformInertial2This, neighbour_patch->m_transformInertial2This);
  dependency_h.second = ctvv_relative;
  m_neighbours.push_back(dependency_h);
  size_t i_neighbour = m_neighbours.size() - 1;  // checked: OK, since size() >= 1

  // Depending on calculation type, insert donor and relative coord transformation into InterCoeffWS lists
  if (m_InterpolateData) {
    m_InterCoeffData_WS.resize(i_neighbour + 1);
    //CHANGE_DONOR! m_InterCoeffData_WS[i_neighbour].setDonorPatch(neighbour_patch);
    m_InterCoeffData_WS[i_neighbour].setCoordTransform(ctvv_relative.extractReverse());
  }   /// @todo would be nice to have a link_vector<T>

  //  if (m_InterpolateGrad1N) {
  //    m_InterCoeffGrad1N_WS.resize(i_neighbour + 1);
  //    //CHANGE_DONOR! m_InterCoeffGrad1N_WS[i_neighbour].setDonorPatch(neighbour_patch);
  //    m_InterCoeffGrad1N_WS[i_neighbour].setCoordTransform(ctvv_relative.extractReverse());
  //  }

  bool found_dependency = computeDependencies(i_neighbour);

  // discard donor neighbourship relation, if not any was found
  if (!found_dependency) {
    m_neighbours.pop_back();
    //    if (m_InterpolateGrad1N) {
    //      m_InterCoeffGrad1N_WS.pop_back();
    //    }
    if (m_InterpolateData) {
      m_InterCoeffData_WS.pop_back();
    }
  }
  /** @todo Must apply InterCoeffxxxx_WS::adjust2Average AFTER all neighbouring donor patches have been inserted.
    *       This ensures, the correct number of hits (in case of multiple donor patches contributing to single cells)
    *       is known. */  // DONE!!!
}

void Patch::finalizeDependencies()
{
  // NOTE: All grad1N stuff ommited 2013_06_12
  //
  // Eliminate pretended receiving cells, if these received no hits at all,
  // m_receive_cell_data_hits[i] == 0 && m_receive_cell_grad1N_hits[i] == 0
  //
  // Do the following:
  //  0) Nothing to do, if no receive cells were found
  //  1) Build a scratch list v_shift with index shifts
  //  2) Build a scratch list index_new_from_old, so that an operation of type
  //     m_receive_cells..new[ll] = m_receive_cells..old[index_new_from_old[ll]]
  //     eliminates all zeros
  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  //  4) Eliminate non-hit cells by re-indexing in InterCoeffWS and reduce donor weights for
  //     receiving cells, being influenced by more than one donor patch. (More than one hit for
  //     same cell).

  //    ATTENTION: Actually indirect_receiveindex is needed only until InterCoeffWS::adjust2Average is done

  //  0) Nothing to do, if no receive cells were found
  if(m_receive_cells.size() == 0) {
    return;
  }

  //  1) Build a scratch list v_shift with index shifts
  size_t cumulated_shift = 0;
  vector<size_t> v_shift;
  v_shift.resize(m_receive_cells.size(), 0);
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    v_shift[i_rec] = cumulated_shift;
    bool any_hit = false;
    if(m_InterpolateData) {
      if(m_receive_cell_data_hits[i_rec] > 0) {
        any_hit = true;
      }
    }
    //    if(m_InterpolateGrad1N) {
    //      if(m_receive_cell_grad1N_hits[i_rec] > 0) {
    //        any_hit = true;
    //      }
    //    }
    if(!any_hit) { // not any donor contribution found
      cumulated_shift++;
    }
  }
  //  2) Build a scratch list index_new_from_old, so that an operation of type
  //     m_receive_cells..new[ll] = m_receive_cells..old[index_new_from_old[ll]]
  //     eliminates all zeros
  vector<size_t> index_new_from_old;
  index_new_from_old.resize(m_receive_cells.size(), 0);
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    index_new_from_old[i_rec] = i_rec - v_shift[i_rec];
  }
  //.. Total number of cells with neighbour contribs
  size_t num_cells_w_contrib = m_receive_cells.size() - cumulated_shift;

  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    m_receive_cells[index_new_from_old[i_rec]] = m_receive_cells[i_rec];
    if(m_InterpolateData) {
      m_receive_cell_data_hits[index_new_from_old[i_rec]] = m_receive_cell_data_hits[i_rec];
    }
    //    if(m_InterpolateGrad1N) {
    //      m_receive_cell_grad1N_hits[index_new_from_old[i_rec]] = m_receive_cell_grad1N_hits[i_rec];
    //    }
  }
  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  m_receive_cells.resize(num_cells_w_contrib);
  if(m_InterpolateData) {
    m_receive_cell_data_hits.resize(num_cells_w_contrib);
  }
  //  if(m_InterpolateGrad1N) {
  //    m_receive_cell_grad1N_hits.resize(num_cells_w_contrib);
  //  }

  //  4) Eliminate non-hit cells by re-indexing in InterCoeffWS and reduce donor weights for
  //     receiving cells, being influenced by more than one donor patch. (More than one hit for
  //     same cell).
  for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
    if(m_InterpolateData) {
      m_InterCoeffData_WS[i_donor].reindexIndReceive(index_new_from_old);
      m_InterCoeffData_WS[i_donor].adjust2Average(m_receive_cell_data_hits);
    }
    //    if(m_InterpolateGrad1N) {
    //      m_InterCoeffGrad1N_WS[i_donor].reindexIndReceive(index_new_from_old);
    //      m_InterCoeffGrad1N_WS[i_donor].adjust2Average(m_receive_cell_grad1N_hits);
    //    }
  }
  // Transfer data to padded data sets, if required.
  if(m_TransferPadded) {    /// @todo I dislike the if-trees. Too many options. Get experience and discard options.
    if(m_InterpolateData) {
      m_InterCoeffData.resize(m_InterCoeffData_WS.size());
      for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
        m_InterCoeffData[i_donor].build(m_InterCoeffData_WS[i_donor]);
        // postponed delete m_InterCoeffData_WS[i_donor];  // erase previous WS-version to save some memory on CPU.
      }

      // Transfer to direct lists
      /** @todo test for direct lists. Improve list build up and delete obsolete instances.
        *  Currently build up in sequence:
        *  m_InterCoeffData_WS -> m_InterCoeffData -> direct lists
        */
      buildDonorTransferData();

    }
    //    if(m_InterpolateGrad1N) {
    //      m_InterCoeffGrad1N.resize(m_InterCoeffGrad1N_WS.size());
    //      for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
    //        m_InterCoeffGrad1N[i_donor].build(m_InterCoeffGrad1N_WS[i_donor]);
    //        // postponed delete m_InterCoeffGrad1N_WS[i_donor];  // erase previous WS-version to save some memory on CPU.
    //      }
    //    }
  }
}

void Patch::compactReceiveCellLists()
{
  // Clean up m_receive_cells and hit counter lists
  // Sort cell indicees
  // Note: indirect indexing, but piecewise constant strides depending on patch type
  //       example structured: strides with (m_NumI*m_NumJ or m_NumJ*m_NumK or m_NumK*m_NumI)
  sort(m_receive_cells.begin(), m_receive_cells.end());
  // Remove duplicates
  vector<size_t>::iterator it;
  it = unique(m_receive_cells.begin(), m_receive_cells.end());
  m_receive_cells.resize(it - m_receive_cells.begin());
  // Create hit counter lists with same size
  if(m_InterpolateData) {
    m_receive_cell_data_hits.resize(m_receive_cells.size(), 0);
  }
  //  if(m_InterpolateGrad1N) {
  //    m_receive_cell_grad1N_hits.resize(m_receive_cells.size(), 0);
  //  }
}

/// @todo new mem-structure: must change access to borrowed pointers
void Patch::accessDonorData_WS(const size_t& field)
{
  vector<real*> donor_vars;
  vector<real*> this_vars;
  // assign variable pointers to work on
  for(size_t i_v=0; i_v<numVariables(); i_v++) {
    this_vars.push_back(getVariable(field, i_v));
  }
  donor_vars.resize(this_vars.size());

  // set all receiving data variables to 0, as donors will add their contributions onto
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
    //
    // NOTE: Following if-cond. is no longer needed to set variables to zero.
    // Reason: ommitting grad1n stuff ensures that m_receive_cells is a list of
    //         receiving nodes that refers only to data (not grad1n!!!) transfers.
    //
    //    if(m_receive_cell_data_hits[ll_rc] != 0) {
    //      size_t l_rc = m_receive_cells[ll_rc];  /// @todo indirect operation!!!
    //      for(size_t i_v=0; i_v<numVariables(); i_v++) {
    //        this_vars[i_v][l_rc] = 0.;
    //      }
    //    }
    // REPLACED BY BLOCK BELOW
#ifdef DEBUG
    if(m_receive_cell_data_hits[ll_rc] == 0) {
      BUG; // since 2013_06_12 (ommitting grad1n stuff) zero hits are no longer possible.
    }
#endif
    size_t l_rc = m_receive_cells[ll_rc];  /// @todo indirect operation!!!
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      this_vars[i_v][l_rc] = 0.;
    }
  }

  // loop through neighbouring donor patches
  for(size_t i_pd=0; i_pd<m_neighbours.size(); i_pd++) {
    Patch* donor = m_neighbours[i_pd].first;
    //.. assign foreign variable pointers (same sequence as in "this_vars"
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      donor_vars[i_v] = donor->getVariable(field, i_v);
    }
    //.. execute transfer contribution
    InterCoeffWS icws = m_InterCoeffData_WS[i_pd];
    icws.transFromTo(donor_vars, this_vars);    // NOTE: NO turning of vectorial variables
  }
}

void Patch::accessDonorDataPadded(const size_t& field)
{
  vector<real*> donor_vars;
  vector<real*> this_vars;
  // assign variable pointers to work on
  for(size_t i_v=0; i_v<numVariables(); i_v++) {
    this_vars.push_back(getVariable(field, i_v));
  }
  donor_vars.resize(this_vars.size());

  // set all receiving data variables to 0, as donors will add their contributions onto
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
#ifdef DEBUG
    if(m_receive_cell_data_hits[ll_rc] == 0) {
      BUG; // since 2013_06_12 (ommitting grad1n stuff) zero hits are no longer possible.
    }
#endif
    size_t l_rc = m_receive_cells[ll_rc];  /// @todo indirect operation!!!
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      this_vars[i_v][l_rc] = 0.;
    }
  }

  // loop through neighbouring donor patches
  for(size_t i_pd=0; i_pd<m_neighbours.size(); i_pd++) {
    Patch* donor = m_neighbours[i_pd].first;
    //.. assign foreign variable pointers (same sequence as in "this_vars"
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      donor_vars[i_v] = donor->getVariable(field, i_v);
    }
    //.. execute transfer contribution
    InterCoeffPad icp = m_InterCoeffData[i_pd];
    icp.transFromTo(donor_vars, this_vars);    // NOTE: NO turning of vectorial variables
  }
}

void Patch::accessDonorDataDirect(const size_t &field) {

  // assign variable pointers to work on
  /// @todo no need to set this_data all times, if fields are allways the same. Do at start up.
  vector<real*> this_vars;
  vector<real*> donor_vars;
  //.. assign this_vars
  for(size_t i_v=0; i_v<numVariables(); i_v++) {
    this_vars.push_back(getVariable(field, i_v));
  }
  //.. set size of donor_vars
  donor_vars.resize(this_vars.size());

  // Set all receiving data variables to 0, as donors will add their contributions onto
  #ifndef DEBUG
  #pragma omp parallel for
  #endif
  for(size_t ll_rc=0; ll_rc < m_NumReceivingCellsUnique; ll_rc++) {
    size_t l_rc = m_ReceivingCellIndicesUnique[ll_rc];
    for (size_t i_v=0; i_v<numVariables(); i_v++) {
      this_vars[i_v][l_rc] = 0.;
    }
  }

  // Get all donor contributions
  // Loop through neighbouring donor patches
  for (size_t i_pd=0; i_pd<m_NumDonorPatches; i_pd++) {
    donor_t& donor = m_Donors[i_pd];
    //.. assign foreign variable pointers (same sequence as in "this_vars")
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      donor_vars[i_v] = donor.data + i_v * donor.variable_size;
    }
    //.. loop for indirect receiving cells
    #ifndef DEBUG
    #pragma omp parallel for
    #endif
    for (size_t ll_rec = 0; ll_rec < donor.num_receiver_cells; ll_rec++) {
      //.... receiving cells index
      size_t i_rec = m_ReceivingCellIndicesConcat[donor.receiver_index_field_start + ll_rec];
      //.... start address in m_DonorCells/m_DonorWeights pattern
      size_t l_doner_cells_start = donor.donor_wi_field_start + ll_rec * donor.stride;
      //.... loop for contributing cells
      for(size_t i_contrib = 0; i_contrib < donor.stride; i_contrib++) {
        size_t l_wi = l_doner_cells_start + i_contrib;      // index of donor cell in concatenated lists
        size_t donor_cell_index = m_DonorIndexConcat[l_wi];
        real donor_cell_weight = m_DonorWeightConcat[l_wi];
        //...... loop for variables
        for(size_t i_v = 0; i_v < m_NumVariables; i_v++) {
          *(this_vars[i_v]+i_rec) += donor_vars[i_v][donor_cell_index] * donor_cell_weight;  // contribute to receiving cell
        }
      }
    }
  }
}


void Patch::accessTurnDonorData_WS(const size_t& field,
                                   const size_t& i_vx, const size_t& i_vy, const size_t& i_vz)
{
  /// @todo must check performance issues of these little vector allocations, assume OK
  vector<real*> donor_vars;
  vector<real*> this_vars;
  // sort variable sequences to get vectorial quantity at position 0,1,2 in sequence of donor_vars/this_vars
  vector<bool> blocked(numVariables(), false);
  vector<size_t> ind_iv;
  ind_iv.push_back(i_vx);
  ind_iv.push_back(i_vy);
  ind_iv.push_back(i_vz);
  blocked[i_vx] = true;
  blocked[i_vy] = true;
  blocked[i_vz] = true;
  for (size_t i_v = 0; i_v < numVariables(); i_v++) {
    if (!blocked[i_v]) {
      ind_iv.push_back(i_v);
    }
  }
  // assign variable pointers to work on
  for(size_t ii_v=0; ii_v<numVariables(); ii_v++) {
    size_t i_v = ind_iv[ii_v];
    this_vars.push_back(getVariable(field, i_v));
  }
  // set all receiving data variables to 0, as donors will add their contributions onto
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
    size_t l_rc = m_receive_cells[ll_rc];
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      this_vars[i_v][l_rc] = 0.;
    }
  }
  // loop through neighbouring donor patches
  for(size_t i_pd=0; i_pd<m_neighbours.size(); i_pd++) {
    Patch* donor = m_neighbours[i_pd].first;
    //.. assign foreign variable pointers (same sequence as in "this_vars"
    for(size_t ii_v=0; ii_v<numVariables(); ii_v++) {
      size_t i_v = ind_iv[ii_v];
      donor_vars.push_back(donor->getVariable(field, i_v));
    }
    //.. execute transfer contribution
    InterCoeffWS icws = m_InterCoeffData_WS[i_pd];
    icws.transTurnFromTo(donor_vars, this_vars);
  }
}

void Patch::allocateData()
{
  m_FieldSize = m_NumVariables * m_VariableSize;
  m_Data = new real [m_NumFields*m_FieldSize];
}

void Patch::deleteData()
{
  if(m_Data) {   /// @todo condition needed ??
    delete [] m_Data;
  }
}

void Patch::resize(size_t variable_size)
{
  m_VariableSize = variable_size;
  allocateData();
}

void Patch::computeVariableDifference(size_t i_field1, size_t i_var1, size_t i_field2, size_t i_var2, real &max_norm, real &l2_norm)
{
  RESTRICT real *var1 = getVariable(i_field1, i_var1);
  RESTRICT real *var2 = getVariable(i_field2, i_var2);
  max_norm = 0.0;
  l2_norm  = 0.0;
  for (size_t i = 0; i < m_VariableSize; ++i) {
    real diff = sqr(var2[i] - var1[i]);
    max_norm = max(max_norm, diff);
    l2_norm += diff;
  }
  max_norm = sqrt(max_norm);
  l2_norm  = sqrt(l2_norm);
}
