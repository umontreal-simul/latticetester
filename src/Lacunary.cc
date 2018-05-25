// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <sstream>
#include <iostream>

#include "latticetester/Lacunary.h"


namespace LatticeTester
{

  /*=========================================================================*/

  std::string Lacunary::toString () const
  {
    std::ostringstream out;
    out << "dim = " << m_dim;
    out << "\nLac = {\n   ";
    for (int i = 0; i < m_dim; i++)
      out << m_lac[i] << "\n   ";
    out << "}\n";
    return out.str ();
  }


  /*=========================================================================*/


  bool Lacunary::calcIndicesStreams (int s, int w, int & minDim,
      int maxDim, int order)
  {
    if (w == 0) {
      if (minDim <= order)
        minDim = order + 1;
      return false;
    }

    if (m_dim < maxDim) {
      LatticeTester::DeleteVect (m_lac);
      LatticeTester::CreateVect (m_lac, maxDim);
      m_dim = maxDim;
    }

    BScal t1;
    LatticeTester::power2 (t1, (long) w);
    BScal t;
    t = 0;
    int i = 1;
    while (true) {
      for (int j = 0; j < s; j++) {
        if (i <= maxDim) {
          m_lac[i] = t + j;
          i++;
        } else
          return true;
      }
      t += t1;
    }
  }

} // namespace LatticeTester
