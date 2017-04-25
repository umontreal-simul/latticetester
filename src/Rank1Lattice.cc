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

#include "latticetester/Rank1Lattice.h"
#include "latticetester/Util.h"
#include <cassert>
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/matrix_proxy.hpp>

#ifdef WITH_NTL
#else
using namespace boost::numeric::ublas;
#endif

using namespace std;


namespace LatticeTester
{

Rank1Lattice::Rank1Lattice (const MScal & n, const MVect & a, int maxDim,
                           NormType norm):
                  IntLattice::IntLattice (n, maxDim, maxDim, norm)
{
   m_a.resize(1 + maxDim);
   m_a = a;
   init();
}


//=========================================================================

Rank1Lattice::~Rank1Lattice()
{
   m_a.clear ();
}


//=========================================================================

void Rank1Lattice::init()
{
   IntLattice::init();
   for (int r = 2; r <= getMaxDim(); r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}


//=========================================================================

Rank1Lattice & Rank1Lattice::operator= (const Rank1Lattice & lat)
{
   if (this == &lat)
      return * this;
   copy (lat);
   init ();
   int maxDim = lat.getMaxDim ();
   m_a.resize (1 + maxDim);
   m_a = lat.m_a;
   return *this;
}


//=========================================================================

Rank1Lattice::Rank1Lattice (const Rank1Lattice & lat):
      IntLattice::IntLattice (lat.m_m, lat.getOrder (),
                              lat.getMaxDim (), lat.getNorm ())
{
   // MyExit (1, "Rank1Lattice:: constructeur n'est pas terminé " );
   init ();
   int maxDim = lat.getMaxDim ();
   m_a.resize (1 + maxDim);
   m_a = lat.m_a;
}


//===========================================================================

void Rank1Lattice::buildBasis (int d)
{
   // assert(d <= getMaxDim());
   setDim (d);

   // conv(m_v[1][1], 1);

   for (int j = 1; j <= d; j++) {
      m_v (1, j) = m_a[j];
   }

   for (int i = 2; i <= d; i++) {
      for (int j = 1; j <= d; j++) {
         if (i == j) {
            m_v (i, j) = m_m;
         } else {
            m_v (i, j) = 0;
         }
      }
   }

   // if a[1] != 1, the basis must be triangularized
   if (m_v (1, 1) != 1) {
      Triangularization < Base > (m_v, m_w, d, d, m_m);
      dualize ();
   }
   CalcDual < Base > (m_v, m_w, d, m_m);
   m_v.setNegativeNorm (true);
   m_w.setNegativeNorm (true);
}

}
