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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "combiobjectand.h"

CombiObjectAnd::CombiObjectAnd(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


CombiObjectAnd::CombiObjectAnd(ObjectDefinition* object_a)
  : CombiObject(object_a)
{
}


bool CombiObjectAnd::evalBool()
{
  bool inside = true;
  for (size_t i_o = 0; i_o < m_Objects.size(); i_o++) {
    if(!m_Objects[i_o]->evalBool()) {
      inside = false;
      break;
    }
  }
  return inside;





//  bool a_inside = m_ObjectA->evalBool();
//  bool b_inside = m_ObjectB->evalBool();

//  return a_inside && b_inside;
}
