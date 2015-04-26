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
#ifndef COMPRESSIBLEVARIABLESANDG_H
#define COMPRESSIBLEVARIABLESANDG_H

#include "drnum.h"
#include "compressiblevariables.h"

#include <QString>


/** @todo
  * This is an intermediate test variable set to allow post processing on blockobject
  * design and simulations.
  */

template <typename TGas>
class CompressibleVariablesAndG : public CompressibleVariables<TGas>
{

  int m_ExtraVarIndex;


public:

  CompressibleVariablesAndG(int extra_var_index = 0);

  virtual int numScalars() const { return 6; }

  virtual string getScalarName(int i) const;
  virtual real   getScalar(int i, Patch* patch, int index, vec3_t x) const;

};

template <typename TGas>
CompressibleVariablesAndG<TGas>::CompressibleVariablesAndG(int extra_var_index)
{
  m_ExtraVarIndex = extra_var_index;
}

template <typename TGas>
string CompressibleVariablesAndG<TGas>::getScalarName(int i) const
{
  if (i <= 4) return CompressibleVariables<TGas>::getScalarName(i);
  if (i == 5) return "G";
  BUG;
  return "N/A";
}

template <typename TGas>
real CompressibleVariablesAndG<TGas>::getScalar(int i, Patch* patch, int index, vec3_t x) const
{
  if (i <= 4) return CompressibleVariables<TGas>::getScalar(i, patch, index, x);
  if (i == 5) return patch->getExtraCPUVarset(m_ExtraVarIndex)[index];

  BUG;
  return 0;
}

#endif // COMPRESSIBLEVARIABLESANDG_H
