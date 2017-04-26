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

#ifdef WITH_NTL
#include "NTL/ZZ.h"
#include "NTL/tools.h"
#endif

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "latticetester/Types.h"
#include "latticetester/Normalizer.h"


namespace LatticeTester
{

const int Normalizer::MAX_DIM;

/*-------------------------------------------------------------------------*/

Normalizer::Normalizer (const MScal & m0, int k0, int maxDim, std::string name,
                        NormType norm, double beta0) :
      m_name(name), m_norm(norm), m_m(m0), m_rank(k0), m_maxDim(maxDim),
      m_beta(beta0)
{
   m_cst = new double[maxDim];
}


/*-------------------------------------------------------------------------*/

void Normalizer::init (const MScal &m0, int k0, double beta0)
/*
 * Computes the vector Cst that corresponds to G, for a lattice of
 * density m^k for t >= k, and m^t for t < k.
 */
{
   double x, y;
   double logBeta;
   double logm;
   double k = k0;
   m_rank = k0;
   m_m = m0;
   m_beta = beta0;

   y = 1.0;
   logBeta = log (m_beta);
#ifdef WITH_NTL
   logm = log(NTL::to_ZZ(m_m));
#else
   logm = log(m_m);
#endif
   for (int j = 0; j < m_maxDim; j++) {
      if (j > k)
         y = k / j;
      x = log (getGamma(j)) + j * logBeta + y * logm;
      if (m_norm == L2NORM)
         x = x + x;
      m_cst[j] = exp (x);
   }
}


/*-------------------------------------------------------------------------*/

std::string Normalizer::ToString () const
{
   std::ostringstream os;
   os << "-----------------------------\n"
   << "Content of Normalizer object:\n\n Normalizer = " << m_name;
   os << "\n m = " << m_m;
   os << "\n rank = " << m_rank;
   os << "\n beta = " << std::setprecision (4) << m_beta << "\n\n";

   //   os.setf(std::ios::left);
   os << std::setprecision (13);
   for (int t = 0; t < m_maxDim; t++) {
      os << " Cst[" << std::setw(2) << std::right << t << "] = "
      << std::setw(14) << std::left << m_cst[t] << "\n";
   }
   os << "\n";
   return os.str ();
}


/*-------------------------------------------------------------------------*/

double Normalizer::getGamma (int) const
{
   return 1.0;
}


/*-------------------------------------------------------------------------*/

double & Normalizer::getCst (int j)
{
   assert (j >= 0 && j < m_maxDim); //fred
   return m_cst[j];
}

}
