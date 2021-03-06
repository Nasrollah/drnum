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

template <class T>
struct mv_p
{
  static T apply(const T &a, const T &b) { return a+b; }
};

template <class T>
struct mv_m
{
  static T apply(const T &a, const T &b) { return a-b; }
};

template <class T>
struct mv_ml
{
  static T apply(const T &a, const T &b) { return a*b; }
};

template <class L, class O, class R>
struct ParseNode
{
  typedef typename R::value_type value_type;
  const L &l;
  const R &r;
  ParseNode(const L &a, const R &b) : l(a), r(b) {}
  value_type operator[](const unsigned int &i) const { return O::apply(l[i], r[i]); }
  unsigned int size() const { return r.size(); }
  value_type abs() const;
  value_type abs2() const;
};

template <class O, class R>
struct ParseNode<real, O, R>
{
  typedef typename R::value_type value_type;
  const real l;
  const R &r;
  ParseNode(const real a, const R &b) : l(a), r(b) {}
  value_type operator[](const unsigned int &i) const { return O::apply(l, r[i]); }
  unsigned int size() const { return r.size(); }
  value_type abs() const;
  value_type abs2() const;
};
