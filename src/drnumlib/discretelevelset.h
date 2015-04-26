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

#ifndef DISCRETELEVELSET_H
#define DISCRETELEVELSET_H

#include "drnum.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkPLYReader.h>

#include "patchgrid.h"
#include "postprocessingvariables.h"
#include "cartesianpatch.h"
#include "geometrytools.h"

#ifdef GPU
//#include "gpu_cartesianpatch.h"
#endif

class DiscreteLevelSet
{

private: // types

  typedef CGAL::Simple_cartesian<double> K;

  typedef K::FT         FT;
  typedef K::Ray_3      Ray;
  typedef K::Line_3     Line;
  typedef K::Point_3    Point;
  typedef K::Triangle_3 Triangle;

  typedef QVector<Triangle>::iterator                       TriangleIterator;
  typedef CGAL::AABB_triangle_primitive<K,TriangleIterator> TrianglePrimitive;
  typedef CGAL::AABB_traits<K, TrianglePrimitive>           TriangleTraits;
  typedef CGAL::AABB_tree<TriangleTraits>                   TriangleTree;
  typedef TriangleTree::Point_and_primitive_id              TrianglePointAndPrimitiveId;


private: // attributes

  PatchGrid*        m_PatchGrid;
  QVector<Triangle> m_Triangles;
  QVector<vec3_t>   m_TriNormals;
  QVector<vec3_t>   m_NodeNormals;
  QVector<vec3_t>   m_Nodes;
  QVector<ijk_t>    m_TriangleNodes;
  TriangleTree      m_TriangleTree;
  real              m_Tol;
  size_t            m_ExtraVarIndex;


protected: //

  void computeLevelSet(vtkPolyData* poly);
  real computePointLevelSet(vec3_t x);
  void levelSetPerCell(size_t i_patch);
  void recursiveLevelSet(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m);
  void interpolate(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m);


public:

  DiscreteLevelSet(PatchGrid *patch_grid, size_t i_extra_var = 0);

  void readGeometry(QString geometry_file_name);

};



class StoredLevelSet
{

protected: // attributes

  size_t m_IVar;


public: // methods

  CUDA_DH StoredLevelSet(size_t i_var) { m_IVar = i_var; }
  CUDA_DH StoredLevelSet() {}

  template <typename T_Patch>
  CUDA_DH real G(T_Patch& patch, size_t i, size_t j, size_t k, size_t i_field = 0)
  {
    return patch.f(i_field, m_IVar, i, j, k);
  }

};



class LevelSetPlotVars : public PostProcessingVariables
{
protected: // attributes

  size_t m_IVar;


public:

  LevelSetPlotVars(size_t i_var) { m_IVar = i_var; }

  virtual int numScalars() const { return 1; }
  virtual int numVectors() const { return 0; }

  virtual string getScalarName(int) const { return "G"; }
  virtual string getVectorName(int) const { BUG; }
  virtual real   getScalar(int, Patch* patch, int index, vec3_t) const { return patch->getExtraCPUVarset(m_IVar)[index]; }
  virtual vec3_t getVector(int, Patch*, int, vec3_t)             const { BUG; return vec3_t(0,0,0); }
};



#endif // DISCRETELEVELSET_H
