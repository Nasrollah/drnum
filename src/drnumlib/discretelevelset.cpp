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

#include "discretelevelset.h"

DiscreteLevelSet::DiscreteLevelSet(PatchGrid *patch_grid, size_t i_extra_var)
{
  m_PatchGrid = patch_grid;
  m_Tol = 0.55;
  m_ExtraVarIndex = i_extra_var;
}

void DiscreteLevelSet::readGeometry(QString geometry_file_name)
{
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  if (geometry_file_name.endsWith(".stl")) {
    cout << "Reading STL geometry" << endl;
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(qPrintable(geometry_file_name));
    reader->MergingOn();
    reader->Update();
    poly->DeepCopy(reader->GetOutput());
  } else if (geometry_file_name.endsWith(".ply")) {
    cout << "Reading PLY geometry" << endl;
    vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(qPrintable(geometry_file_name));
    reader->Update();
    poly->DeepCopy(reader->GetOutput());
  } else {
   BUG;
  }
  poly->BuildCells();
  computeLevelSet(poly);
}

void DiscreteLevelSet::computeLevelSet(vtkPolyData *poly)
{
  // build triangle tree
  {
    int num_faces = poly->GetNumberOfCells();
    int num_nodes = poly->GetNumberOfPoints();
    m_Triangles.clear();
    m_TriNormals.clear();
    m_Triangles.fill(Triangle(), num_faces);
    m_TriNormals.fill(vec3_t(0,0,0), num_faces);
    m_NodeNormals.clear();
    m_NodeNormals.fill(vec3_t(0,0,0), num_nodes);
    m_TriangleNodes.resize(num_faces);
    m_Nodes.resize(num_nodes);
    QVector<real> node_weights(num_nodes, 0.0);
    for (vtkIdType id_cell = 0; id_cell < num_faces; ++id_cell) {
      vtkIdType num_pts, *pts;
      poly->GetCellPoints(id_cell, num_pts, pts);
      if (num_pts != 3) {
        ERROR("only triangulated geometries are allowed");
      }
      dvec3_t a, b, c;
      poly->GetPoint(pts[0], a.data());
      poly->GetPoint(pts[1], b.data());
      poly->GetPoint(pts[2], c.data());
      m_Nodes[pts[0]] = vec3_t(a[0], a[1], a[2]);
      m_Nodes[pts[1]] = vec3_t(b[0], b[1], b[2]);
      m_Nodes[pts[2]] = vec3_t(c[0], c[1], c[2]);
      m_Triangles[id_cell] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      vec3_t u = b - a;
      vec3_t v = c - a;
      m_TriNormals[id_cell] = u.cross(v);
      m_TriNormals[id_cell].normalise();
      real w0 = fabs(GeometryTools::angle(c-a, b-a));
      real w1 = fabs(GeometryTools::angle(a-b, c-b));
      real w2 = fabs(GeometryTools::angle(b-c, a-c));
      node_weights[pts[0]] += w0;
      node_weights[pts[1]] += w1;
      node_weights[pts[2]] += w2;
      m_NodeNormals[pts[0]] += w0*m_TriNormals[id_cell];
      m_NodeNormals[pts[1]] += w1*m_TriNormals[id_cell];
      m_NodeNormals[pts[2]] += w2*m_TriNormals[id_cell];
      m_TriangleNodes[id_cell].i = pts[0];
      m_TriangleNodes[id_cell].j = pts[1];
      m_TriangleNodes[id_cell].k = pts[2];
    }
    m_TriangleTree.rebuild(m_Triangles.begin(), m_Triangles.end());
    m_TriangleTree.accelerate_distance_queries();
    for (int i = 0; i < m_NodeNormals.size(); ++i) {
      m_NodeNormals[i] *= 1.0/node_weights[i];
    }
  }

  for (size_t i_patch = 0; i_patch < m_PatchGrid->getNumPatches(); ++i_patch) {

    CartesianPatch* patch = dynamic_cast<CartesianPatch*>(m_PatchGrid->getPatch(i_patch));
    if (patch) {
      size_t i_m = patch->sizeI() - 1;
      size_t j_m = patch->sizeJ() - 1;
      size_t k_m = patch->sizeK() - 1;

      if ( i_m == 0 || j_m == 0 || k_m == 0) {
         cout << endl << endl << "Error! i-> " << i_m << " j-> " << j_m << " k-> " << k_m << endl;
         BUG;
      }
      recursiveLevelSet(patch, 0, 0, 0, i_m, j_m, k_m);
    }
    else {
      cout << "Computing levelset cell per cell for patch-> " << i_patch << "." << endl;
      levelSetPerCell(i_patch);
    }
  }
}

void DiscreteLevelSet::recursiveLevelSet(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m)
{
  size_t index[8];
  index[0] = patch->index(i_o, j_o, k_o);
  index[1] = patch->index(i_m, j_o, k_o);
  index[2] = patch->index(i_o, j_m, k_o);
  index[3] = patch->index(i_o, j_o, k_m);
  index[4] = patch->index(i_m, j_m, k_o);
  index[5] = patch->index(i_m, j_o, k_m);
  index[6] = patch->index(i_o, j_m, k_m);
  index[7] = patch->index(i_m, j_m, k_m);

  real g[8];
  for (int i = 0; i < 8; ++i) {
    vec3_t x = patch->xyzoCell(index[i]);
    g[i] = computePointLevelSet(x);
  }
  vec3_t x_ooo = patch->xyzoCell(index[0]);
  vec3_t x_ijk = patch->xyzoCell(index[7]);
  real d = m_Tol*(x_ooo - x_ijk).abs();

  if ( fabs(g[0]) < d || fabs(g[1]) < d || fabs(g[2]) < d || fabs(g[3]) < d ||
       fabs(g[4]) < d || fabs(g[5]) < d || fabs(g[6]) < d || fabs(g[7]) < d)
  {
    int trunc = 5;
    if ( (i_m-i_o) <= trunc && (j_m-j_o) <= trunc && (k_m-k_o) <= trunc) {
      for (size_t i = i_o; i <= i_m; ++i) {
        for (size_t j = j_o; j <= j_m; ++j) {
          for (size_t k = k_o; k <= k_m; ++k) {
            size_t i_cell = patch->index(i,j,k);
            vec3_t x = patch->xyzoCell(i_cell);
            real g = computePointLevelSet(x);
            patch->getExtraCPUVarset(m_ExtraVarIndex)[i_cell] = g;
          }
        }
      }
      return;
    } else {
      if ( (i_m-i_o) == 0 || (j_m-j_o) == 0 || (k_m-k_o) == 0) {
        cout << "WARNING at least one x_m-x_o == 0" << endl;
      }
      if ( (i_m - i_o) >= max(j_m-j_o, k_m-k_o)) {
        //cout << "half i" << endl;
        int i_new = (i_m - i_o)/2 + i_o;
        recursiveLevelSet(patch, i_o  , j_o, k_o, i_new, j_m, k_m );
        recursiveLevelSet(patch, i_new, j_o, k_o, i_m  , j_m, k_m );
        return;
      }
      if ( (j_m - j_o) >= max(i_m-i_o, k_m-k_o)) {
        //cout << "half j" << endl;
        int j_new = (j_m - j_o)/2 + j_o;
        recursiveLevelSet(patch, i_o, j_o  , k_o, i_m, j_new, k_m );
        recursiveLevelSet(patch, i_o, j_new, k_o, i_m, j_m, k_m );
        return;
      }
      if ( (k_m - k_o) >= max(i_m-i_o, j_m-j_o)) {
        //cout << "half k" << endl;
        int k_new = (k_m - k_o)/2 + k_o;
        recursiveLevelSet(patch, i_o, j_o, k_o  , i_m, j_m, k_new );
        recursiveLevelSet(patch, i_o, j_o, k_new, i_m, j_m, k_m   );
        return;
      }
    }
  }
  interpolate(patch, i_o, j_o, k_o, i_m, j_m, k_m);
}

void DiscreteLevelSet::interpolate(CartesianPatch* patch, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  real G_111 = computePointLevelSet(patch->xyzoCell(patch->index(i1, j1, k1)));
  real G_112 = computePointLevelSet(patch->xyzoCell(patch->index(i1, j1, k2)));
  real G_121 = computePointLevelSet(patch->xyzoCell(patch->index(i1, j2, k1)));
  real G_122 = computePointLevelSet(patch->xyzoCell(patch->index(i1, j2, k2)));
  real G_211 = computePointLevelSet(patch->xyzoCell(patch->index(i2, j1, k1)));
  real G_212 = computePointLevelSet(patch->xyzoCell(patch->index(i2, j1, k2)));
  real G_221 = computePointLevelSet(patch->xyzoCell(patch->index(i2, j2, k1)));
  real G_222 = computePointLevelSet(patch->xyzoCell(patch->index(i2, j2, k2)));

  for (size_t i = i1; i <= i2; ++i) {
    real di = real(i - i1)/real(i2 - i1);
    real G_i11 = (1-di)*G_111 + di*G_211;
    real G_i12 = (1-di)*G_112 + di*G_212;
    real G_i21 = (1-di)*G_121 + di*G_221;
    real G_i22 = (1-di)*G_122 + di*G_222;
    for (size_t j = j1; j <= j2; ++j) {
      real dj = real(j - j1)/real(j2 - j1);
      real G_ij1 = (1-dj)*G_i11 + dj*G_i21;
      real G_ij2 = (1-dj)*G_i12 + dj*G_i22;
      for (size_t k = k1; k <= k2; ++k) {
        real dk = real(k - k1)/real(k2 - k1);
        real G_ijk = (1-dk)*G_ij1 + dk*G_ij2;
        patch->getExtraCPUVarset(m_ExtraVarIndex)[patch->index(i,j,k)] = G_ijk;
      }
    }
  }
}

real DiscreteLevelSet::computePointLevelSet(vec3_t x)
{
  real G;
  try {
    Point p(x[0], x[1], x[2]);
    TrianglePointAndPrimitiveId result = m_TriangleTree.closest_point_and_primitive(p);
    Point cp = result.first;
    Triangle* T = result.second;
    int id = (T - m_Triangles.begin());
    vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
    vec3_t v = x - x_snap;
    G = v.abs();

    vec3_t a  = m_Nodes[m_TriangleNodes[id].i];
    vec3_t b  = m_Nodes[m_TriangleNodes[id].j];
    vec3_t c  = m_Nodes[m_TriangleNodes[id].k];
    vec3_t na = m_NodeNormals[m_TriangleNodes[id].i];
    vec3_t nb = m_NodeNormals[m_TriangleNodes[id].j];
    vec3_t nc = m_NodeNormals[m_TriangleNodes[id].k];

    mat3_t M;
    M[0] = b - a;
    M[1] = c - a;
    M[2] = M[0].cross(M[1]);
    M = M.transp();
    M = M.inverse();
    vec3_t r = x_snap - a;
    r = M*r;
    vec3_t n = na + r[0]*(nb - na) + r[1]*(nc - na);

    if (v*n < 0) {
      G *= -1;
    }
  } catch (...) {
    ERROR("cannot compute distance");
  }
  return G;
}

void DiscreteLevelSet::levelSetPerCell(size_t i_patch)
{
  for (size_t i_cell = 0; i_cell < m_PatchGrid->getPatch(i_patch)->variableSize(); ++i_cell) {
    vec3_t x = m_PatchGrid->getPatch(i_patch)->xyzoCell(i_cell);
    m_PatchGrid->getPatch(i_patch)->getExtraCPUVarset(m_ExtraVarIndex)[i_cell] = computePointLevelSet(x);
  }
}

