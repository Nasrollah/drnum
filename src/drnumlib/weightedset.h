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
#if !defined(WEIGHTEDSET_HH)
#define WEIGHTEDSET_HH

/**
be careful by using this class, when you allocate and delete thousands of objects of this type
**/
template<class T> class WeightedSet;

#include <algorithm>
#include <vector>
#include <iostream>
//#include <complex>
#include "drnum.h"

using namespace std;


/// @todo clean up this class.

/**
 * Class around vector<pair<size_t, T> > to keep pairs (index, a_T)
 *
 * Features simple arithmetics like adding, sorting, while concatenating a unique
 * set of indicees.
 *
 * operator+ and operator* must be defined for type T
 *
 * Example:
 *   WeightedSet<float>* a_ws;
 *   WeightedSet<float>* b_ws;
 *   a_ws = new WeightedSet<float>(2, 2.);        // creates a WS with one pair (2, 2.)
 *   vector<size_t> b_v; b_v.push_back(1); b_v.push_back(2); b_v.push_back(3);
 *   b_ws = new WeightedSet<float> b_ws(b_v, 1.); // creates a WS with the pairs (1, 1.), (2, 1.), (3, 1.)
 *   b_ws->pushBack(0, 100.);                     // appends pair (0, 100.) to the end of b_ws
 *   *a_ws += *b_ws;                              // creates a_ws = ((0, 100.), (1, 1.), (2, 3.), (3, 1.));
 *   *a_ws *= 3.;                                 // creates a_ws = ((0, 300.), (1, 3.), (2, 9.), (3, 3.));
 *
 */
template<class T>
class WeightedSet
{
public:     /// @todo need usefull protection rules
  vector<pair<size_t, T> > v;
  bool sure_sorted;
  T t_zero;
  T t_one;

public:

  /** Dummy construction
   */
  WeightedSet();

  /** Construction considering a single pair (j,w)
   * @param j the index for the single pair
   * @param t value of pair
   */
  WeightedSet(size_t j, T t);

  /** Construction considering a single index with weight "1"
   * @param j the index for the single pair
   */
  WeightedSet(size_t j);

  /** Construction considering equiweighted average for an index vector
   * @param vn the index vector
   */
  WeightedSet(vector<size_t> vn);

  /** Construction considering same weight for an index vector
   * @param vn the index vector
   * @param weight the weight of all entries
   */
  WeightedSet(vector<size_t> vn, T t);

  /** Construction considering equiweighted average for an int-index vector
   * @param vn the index vector
   */
  WeightedSet(vector<int> vn);

  /** Construction considering same weight for an int-index vector
   * @param vn the index vector
   * @param weight the weight of all entries
   */
  WeightedSet(vector<int> vn, T t);

  // Member functions
  void pushBack(const size_t& j, const T& t);            ///< appends a pair (j,s)
  void setup();                                          ///< initializes
  void sortWS();                                         ///< sort contents by index (first)
  void clearWS();                                        ///< Erase all contents
  void PrintNumSets();                                   ///< print to screen
  void MultScalar(const T& s);                           ///< multiply     v = s*v
  void Add(const WeightedSet<T>& w);                     ///< add:         v = v + w
  void MultiplyAdd(const WeightedSet<T>& w, const T& s); ///< mult_add:    v = v + s*w
  void Concatenate(const WeightedSet<T>& w);             ///< concatenate: v = [v , w]
  void Concatenate(const WeightedSet<T>& w, const T& s); ///< concatenate: v = [v , s*w]
  void Unify();    ///< adds up weights for multiple indicees and makes indicees unique
  size_t getSize() const;                  ///< return number of pairs.
  size_t getIndex(const size_t& ll) const; ///< return the index of the ll-th pair (index, weight)
  T getWeight(const size_t& ll) const;     ///< return the weight of the ll-th pair (index, weight)
  size_t HighNodeAddress() const;   ///< Find highest index stored v[..].first (compatibility)
  size_t HighestIndex() const;      ///< Find highest index stored v[..].first
  T WeightSum();                    ///< Compute sum of weights
  T WeightAbsMax();                 ///< Compute maximum absolute weight
  T RealValue(T* a);                ///< Compute sum(a[j[i]] * t[i]) for all i (compatibility)
  T computeValue(const T* a) const; ///< Compute sum(a[j[i]] * t[i]) for all i

  // Operators
  void operator=(const WeightedSet<T>& w);               ///< copy:        v = w
  void operator+=(const WeightedSet<T>& w);              ///< add:         v = v + w
  void operator*=(const T& s);                           ///< multiply     v = s*v

  /**
   * Eliminate all entries with absolute values of weight below eps.
   * @param eps the threshold value
   * @param relative if true, causes function to work with eps relative to largest absolute weight
   * @param keep_weight_sum if true, causes distribution of defects on remaining weights
   */
  void EliminateBelowEps(T eps, const bool& relative, const bool& keep_weight_sum);

  /**
   * Adjust sum of weights.
   * Shift-version: act by adding a constant to all contributors
   * @param aimed sum of the weights
   * @return true, if no error occured
   */
  bool adjustWeightSumShift(const T& weight_sum);

  /**
    * Transfer to fixed sized arrays for "indices" and "weights".
    * @param n_dimension fixed size of arrays
    * @param indices pointer on array for indices (will write onto)
    * @param weights pointer on array for weights (will write onto)
    * @return bool true indicatinf success
    */
  bool tranferToFixedArrays(size_t n_dimension,
                            size_t* indices, T* weights);

  /**
    * Transfer to fixed sized arrays for "indices" and "weights".
    * Overload version with ints
    * @param n_dimension fixed size of arrays
    * @param indices pointer on array for indices (will write onto)
    * @param weights pointer on array for weights (will write onto)
    * @return bool true indicatinf success
    */
  bool tranferToFixedArrays(size_t n_dimension,
                            int* indices, T* weights);


};

template<class T>
inline WeightedSet<T>::WeightedSet()
{
  setup();
  sure_sorted = true;
}

template<class T>
inline WeightedSet<T>::WeightedSet(size_t j, T t)
{
  setup();
  pair<size_t, T> new_pair;
  new_pair.first = j;
  new_pair.second = t;
  v.push_back(new_pair);
  sure_sorted = true;
}

template<class T>
inline WeightedSet<T>::WeightedSet(size_t j)
{
  setup();
  WeightedSet(j, t_one);
}

template<class T>
inline WeightedSet<T>::WeightedSet(vector<size_t> vn)
{
  setup();
  T t = t_one/vn.size();
  for(size_t i=0; i<vn.size(); i++) {
    pair<size_t, T> new_pair;
    new_pair.first = vn[i];
    new_pair.second = t;
    v.push_back(new_pair);
  }
  sortWS();
}

//  template<class T>
//  inline WeightedSet<T>::WeightedSet(size_t k_node) {
//    setup();
//    pair<size_t, T> new_pair;
//    new_pair.first = k_node;
//    new_pair.second = t_one;
//    v.push_back(new_pair);
//    sure_sorted = true;
//  }

template<class T>
inline WeightedSet<T>::WeightedSet(vector<size_t> vn, T t)
{
  setup();
  for(size_t i=0; i<vn.size(); i++) {
    pair<size_t, T> new_pair;
    new_pair.first = vn[i];
    new_pair.second = t;
    v.push_back(new_pair);
  }
  sortWS();
}

template<class T>
inline WeightedSet<T>::WeightedSet(vector<int> vn)
{
  setup();
  T t = t_one/vn.size();
  for(size_t i=0; i<vn.size(); i++) {
    pair<size_t, T> new_pair;
    new_pair.first = vn[i];
    new_pair.second = t;
    v.push_back(new_pair);
  }
  sortWS();
}

template<class T>
inline WeightedSet<T>::WeightedSet(vector<int> vn, T t)
{
  setup();
  for(size_t i=0; i<vn.size(); i++) {
    pair<size_t, T> new_pair;
    new_pair.first = vn[i];
    new_pair.second = t;
    v.push_back(new_pair);
  }
  sortWS();
}

template<class T>
inline void WeightedSet<T>::setup()
{
  t_zero = T(0);
  t_one = T(1);
  v.clear();
}

template<>
inline void WeightedSet<float>::setup()
{
  t_zero = 0.;
  t_one = 1.;
  v.clear();
}

template<>
inline void WeightedSet<double>::setup()
{
  t_zero = 0.;
  t_one = 1.;
  v.clear();
}

template<>
inline void WeightedSet<int>::setup()
{
  t_zero = 0;
  t_one = 1;
  v.clear();
}

template<>
inline void WeightedSet<size_t>::setup()
{
  t_zero = 0;
  t_one = 1;
  v.clear();
}

/// @todo Problem with ambigous real definition
//template<>
//inline void WeightedSet<complex<float> >::setup()
//{
//  t_zero = complex<float>(0.0, 0.0);
//  t_one = complex<float>(1.0, 0.0); //ATTENTION
//  v.clear();
//}

//template<>
//inline void WeightedSet<complex<double> >::setup()
//{
//  t_zero = complex<double>(0.0, 0.0);
//  t_one = complex<double>(1.0, 0.0); //ATTENTION
//  v.clear();
//}

template<class T>
inline void WeightedSet<T>::sortWS()
{
  sort(v.begin(), v.end());
  sure_sorted = true;
}

template<class T>
inline void WeightedSet<T>::clearWS()
{
  v.clear();
}

template<class T>
void WeightedSet<T>::pushBack(const size_t& j, const T& s) {
  pair<size_t, T> new_pair;
  new_pair.first = j;
  new_pair.second = s;
  v.push_back(new_pair);
  sure_sorted = false;
}

template<class T>
inline void WeightedSet<T>::PrintNumSets()
{
  T weight_sum = WeightSum();
  for(size_t i=0; i<v.size(); i++) {
    cout << "i = " << v[i].first << " ; weight = " << v[i].second  << " rel_weight = " << (v[i].second/weight_sum) << endl;
  }
  cout << "sum of weights = " << weight_sum << endl;
}

template<class T>
inline void WeightedSet<T>::operator=(const WeightedSet<T>& a_ws)
{
  v.clear();             /// @todo not sure, if this frees the mem
  v.reserve(a_ws.v.size());
  v.resize(a_ws.v.size());
  for(size_t i=0; i<a_ws.v.size(); i++) {
    v[i] = a_ws.v[i];
  }
  sure_sorted = a_ws.sure_sorted;
}

template<class T>
inline void WeightedSet<T>::operator+=(const WeightedSet<T>& a_ws) {
  Add(a_ws);
}

template<class T>
inline void WeightedSet<T>::MultScalar(const T& scalar)
{
  for(size_t i=0; i<v.size(); i++) {
    v[i].second *= scalar;
  }
}

template<class T>
inline void WeightedSet<T>::operator*=(const T& scalar) {
  MultScalar(scalar);
}

/// @todo modified 2012_07_30: test carefully.
//  fixed bug 2012_08_09

template<class T>
inline void WeightedSet<T>::Unify()
{
  sort(v.begin(), v.end());
  size_t it = 0;
  size_t hold_it = 0;
  bool initial = true;
  while(it < v.size()) {
    if(initial) {
      initial = false;
      v[hold_it] = v[it];
      //  it++;
    } else {
      if(v[it].first == v[hold_it].first) {
        // same address: add weights
        v[hold_it].second += v[it].second;
        //        it = erase(it);
      } else {
        // other address: copy
        hold_it++;
        v[hold_it] = v[it];
      }
    }
    it++;
  }
  v.resize(hold_it+1);

  // test explicit memory shrinking
  vector<pair<size_t, T> >(v).swap(v);

  sure_sorted = true;
}


template<class T>
inline void WeightedSet<T>::Concatenate(const WeightedSet<T>& w)
{
  for (size_t i = 0; i < w.v.size(); ++i) {
    push_back(w[i]);
  }
  sure_sorted = false;
}

template<class T>
inline void WeightedSet<T>::Concatenate(const WeightedSet<T>& w, const T& s)
{
  for (size_t i = 0; i < w.v.size(); ++i) {
    pair<size_t, T> pp;
    pp.first = w[i].first;
    pp.second = w[i].second * s;
    v.push_back(pp);
  }
  sure_sorted = false;
}

template<class T>
inline T WeightedSet<T>::WeightSum()
{
  T weight_sum = t_zero;
  for(size_t i=0; i < v.size(); i++) {
    weight_sum += v[i].second;
  }
  return weight_sum;
}

/// @todo new since 2012_07_30: test carefully.
template<class T>
inline T WeightedSet<T>::WeightAbsMax()
{
  T weight_absmax = t_zero;
  for(size_t i=0; i < v.size(); i++) {
    T test_w = v[i].second;
    if(test_w < t_zero) test_w = -test_w;  // no general abs-function for float, double, ...
    if(test_w > weight_absmax) weight_absmax = test_w;
  }
  return weight_absmax;
}

template<class T>
inline void WeightedSet<T>::Add(const WeightedSet<T>& addend_2) {
  Concatenate(addend_2);
  Unify();
}

template<class T>
inline void WeightedSet<T>::MultiplyAdd(const WeightedSet<T>& addend_2, const T& scalar) {
  Concatenate(addend_2, scalar);
  Unify();
}

// Include functions previously in WeightedSet.cc below

/** old heritage version.
* @todo this one was better than current Add, as it saved mem on large data sets.
* restore to new attribute-version, if needed.
*/
//template<class T>
//void WeightedSet<T>::Add(WeightedSet<T> addend_2)
//{
//  // check input for being sorted
//  //not in stl standard?  if(!is_sorted(begin(), end())) {
//  //  if(!sure_sorted) {
//  sort(begin(), end());
//  sure_sorted = true;
//  //  };
//  //not in stl standard?  if(!is_sorted(addend_2.begin(), addend_2.end())) {
//  //  if(!addend_2.sure_sorted) {
//  sort(addend_2.begin(), addend_2.end());
//  addend_2.sure_sorted = true;
//  //  };
//  //
//  WeightedSet<T> addend_1(*this);
//  //  addend_1 = *this;
//  clear();

//  // add go through sorted sequences of addends
//  pair<size_t, T> new_pair;
//  WeightedSet<T>::iterator addend_1_it = addend_1.begin();
//  WeightedSet<T>::iterator addend_2_it = addend_2.begin();
//  //  WeightedSet<T>::iterator this_it = begin();
//  while((addend_1_it != addend_1.end()) && (addend_2_it != addend_2.end())) {
//    if(addend_1_it->first < addend_2_it->first) {
//      // take the entry of addend_1
//      new_pair.first = addend_1_it->first;
//      new_pair.second = addend_1_it->second;
//      addend_1_it++;
//    }
//    else if(addend_1_it->first > addend_2_it->first) {
//      // take the entry of addend_2
//      new_pair.first = addend_2_it->first;
//      new_pair.second = addend_2_it->second;
//      addend_2_it++;
//    }
//    else {
//      // both addresses are the same: add weight
//      new_pair.first = addend_1_it->first;  // same as *addend_2->first
//      new_pair.second = addend_1_it->second + addend_2_it->second;
//      addend_1_it++;
//      addend_2_it++;
//    };
//    push_back(new_pair);
//  };
//  // append the rest of the addend that is still running
//  while(addend_1_it != addend_1.end()) {
//    new_pair.first = addend_1_it->first;
//    new_pair.second = addend_1_it->second;
//    addend_1_it++;
//    push_back(new_pair);
//  };
//  while(addend_2_it != addend_2.end()) {
//    new_pair.first = addend_2_it->first;
//    new_pair.second = addend_2_it->second;
//    addend_2_it++;
//    push_back(new_pair);
//  };
//}

template<class T>
inline size_t WeightedSet<T>::getSize() const
{
  return v.size();
}

template<class T>
inline size_t WeightedSet<T>::getIndex(const size_t& ll) const
{
  return v[ll].first;
}

template<class T>
inline T WeightedSet<T>::getWeight(const size_t& ll) const
{
  return v[ll].second;
}

template<class T>
inline size_t WeightedSet<T>::HighNodeAddress() const
{
  return HighestIndex();
}

template<class T>
inline size_t WeightedSet<T>::HighestIndex() const
{
  // get the maximum internal index of all entries
  size_t high_index = 0;
  for (size_t i=0; i<v.size(); i++) {
    if(v[i].first > high_index) {
      high_index = v[i].first;
    };
  };
  return high_index;
}

/// @todo new since 2012_07_30: test carefully.
template<class T>
inline void WeightedSet<T>::EliminateBelowEps(T eps, const bool& relative, const bool& keep_sum)
{
  // eliminate all entries with weights below eps
  Unify();   /// @todo Unify should check itself wether to do or not
  //  T weight_sum = WeightSum();
  if(relative) {
    eps *= WeightAbsMax();
  }
  T defect = t_zero;
  typename vector<pair<size_t, T> >::iterator it = v.begin();
  while(it != v.end()) {
    T test_w = (*it).second;
    if (test_w < t_zero) {
      test_w = -test_w;  // no general abs-function for float and double
    }
    if(test_w < eps) {
      defect += (*it).second;
      it = v.erase(it);
    } else {
      it++;
    }
  }
  /// @todo Correction below is potentially unsafe.
  if(keep_sum && (v.size()>0) ) {
    //unsafe T correction = weight_sum/(weight_sum-defect);
    //unsafe for(vector<pair<size_t, T> >::iterator i = begin(); i!=end(); i++) {
    //unsafe   (*i).second *= correction;
    // }
    T correction = defect / v.size();
    for(size_t i=0; i<v.size(); i++) {
      v[i].second += correction;
    }
  }
}

template<class T>
inline bool WeightedSet<T>::adjustWeightSumShift(const T& weight_sum)
{
  bool success = true;
  if(v.size() > 0) {
    T old_weight_sum = WeightSum();
    T correction = (weight_sum-old_weight_sum) / v.size();
    typename vector<pair<size_t, T> >::iterator it;
    for(it = v.begin(); it!=v.end(); it++) {
      (*it).second += correction;
    }
  } else {
    success = false;
  }
  return success;
}


template<class T>
inline bool WeightedSet<T>::tranferToFixedArrays(size_t n_dimension,
                                                 size_t* indices, T* weights)
{
  // Find the n_dimension highest weights in the set and store these on
  // fixed size arrays "indices" and "weights".
  // However, do not alter the sequence of the remaining entries in the set
  // Options:
  // 0) v.size() == 0 return as error
  // 1) matching size:
  //    * copy 1:1 onto fixes sized arrays
  // 2) fixed size array dimension "n_dimension" larger than v.size() :
  //    * apply padding: fill up indices with last index in v (v[v.size()-1].first)
  //                     fill up weights with t_zero
  // 3) fixed size array dimension "n_dimension" smaller than v.size() :
  //    * not implemented yet, issue a BUG

  // Option 0;
  if(v.size() == 0) {
    return false;
  }

  // Option 1: matching size
  if(n_dimension == v.size()) {
    for(size_t i=0; i<v.size(); i++) {
      indices[i] = v[i].first;
      weights[i] = v[i].second;
    }
  }

  // Option 2: fixed size array dimension "n_dimension" larger than v.size()
  else if(n_dimension > v.size()) {
    //.. copy
    for(size_t i=0; i<v.size(); i++) {
      indices[i] = v[i].first;
      weights[i] = v[i].second;
    }
    //.. apply padding
    for(size_t i=v.size(); i<n_dimension; i++) {
      indices[i] = v[v.size()].first;
      weights[i] = t_zero;
    }
  }

  // Option 3: fixed size array dimension "n_dimension" smaller than v.size()
  else {
    // not implemented yet, issue a BUG
    BUG;
  }

  return true;
}


template<class T>
inline bool WeightedSet<T>::tranferToFixedArrays(size_t n_dimension,
                                                 int* indices, T* weights)
{
  // Find the n_dimension highest weights in the set and store these on
  // fixed size arrays "indices" and "weights".
  // However, do not alter the sequence of the remaining entries in the set
  // Options:
  // 0) v.size() == 0 return as error
  // 1) matching size:
  //    * copy 1:1 onto fixes sized arrays
  // 2) fixed size array dimension "n_dimension" larger than v.size() :
  //    * apply padding: fill up indices with last index in v (v[v.size()-1].first)
  //                     fill up weights with t_zero
  // 3) fixed size array dimension "n_dimension" smaller than v.size() :
  //    * not implemented yet, issue a BUG

  // Option 0;
  if(v.size() == 0) {
    return false;
  }

  // Option 1: matching size
  if(n_dimension == v.size()) {
    for(size_t i=0; i<v.size(); i++) {
      indices[i] = v[i].first;
      weights[i] = v[i].second;
    }
  }

  // Option 2: fixed size array dimension "n_dimension" larger than v.size()
  else if(n_dimension > v.size()) {
    //.. copy
    for(size_t i=0; i<v.size(); i++) {
      indices[i] = v[i].first;  // auto-cast
      weights[i] = v[i].second;
    }
    //.. apply padding
    for(size_t i=v.size(); i<n_dimension; i++) {
      indices[i] = v[v.size()].first;  // auto-cast
      weights[i] = t_zero;
    }
  }

  // Option 3: fixed size array dimension "n_dimension" smaller than v.size()
  else {
    // not implemented yet, issue a BUG
    BUG;
  }

  return true;
}


template<class T>
inline T WeightedSet<T>::RealValue(T* a)
{
  T ret_val = t_zero;
  for(size_t i=0; i<v.size(); i++) {
    ret_val += a[v[i].first] * v[i].second;
  }
  return ret_val;
}

template<class T>
inline T WeightedSet<T>::computeValue(const T* a) const
{
  T ret_val = t_zero;
  for(size_t i=0; i<v.size(); i++) {
    ret_val += a[v[i].first] * v[i].second;
  }
  return ret_val;
}

#endif
