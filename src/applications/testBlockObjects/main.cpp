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
#ifndef GPU
#include "main.h"
#endif

extern "C" void GPU_main();

int main()
{
#ifdef GPU
  GPU_main();
#else

  cout << "hi" << endl;
#ifdef DEBUG
  cout << endl;
  cout << "DEBUG-MODE !!" << endl;
#else
  cout << "RELEASE-MODE !!" << endl;
#endif

#ifdef OPEN_MP
  int num_threads = 2;
  omp_set_num_threads(num_threads);
  cout << endl;
  cout << "*** NUMBER THREADS: " << num_threads << endl;
  cout << endl;
#endif

  run();
#endif
}
