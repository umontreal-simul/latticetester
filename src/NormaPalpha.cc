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

#include "latticetester/NormaPalpha.h"
#include "latticetester/IntFactor.h"
#include <stdexcept>


namespace LatticeTester
{

NormaPalpha::NormaPalpha (const MScal & m, int alpha, int s, NormType norm)
      : Normalizer (s, "Palpha", norm, 1.0)
{
   if (s > MAX_DIM)
      throw std::invalid_argument("NormaPalpha:   dimension > MAX_DIM");
   m_m = m;
   m_alpha = alpha;
}

/*=========================================================================*/

void NormaPalpha::init (int alpha)
/*
 * Computes the vector m_bounds that corresponds to the upper bound for a 
 * rank 1 lattice of density \f$m\f$ (prime number). The bound doesn't exit 
 * for dimension < 2.
 */
{
   m_alpha = alpha;
   for (int j = 2; j <= m_maxDim; j++)
      m_bounds[j] = calcBound (alpha, j);
}

/*=========================================================================*/

double NormaPalpha::calcBound (int alpha, int dim)
{
   double Res;

   const double eBasis = 2.71828182845904523536;
   double MM;
   conv (MM, m_m);

   if (dim <= 1) {
      std::cout << "NormaPalpha::calcBound:  dim < 2.   Returns -1" << std:: endl;
      return -1;
   }
   if (alpha <= 1) {
      std::cout << "NormaPalpha::calcBound:  alpha < 2.   Returns -1" << std:: endl;
      return -1;
   }

   int stat = IntFactor::isPrime (m_m,0);
   if (stat != PRIME) {
      std::cout << "NormaPalpha::calcBound:  m is not prime.   Returns -1" << std:: endl;
      return -1;
   }
   
   double Term1 = log (MM);
   if (Term1 <= alpha*dim /(alpha - 1)) {
      std::cout << "NormaPalpha::calcBound:" << std::endl;
      std::cout << "   m < exp(alpha*dim/(alpha - 1)) for dim = " << dim << std::endl;
      std::cout << "   Assumption required for existence of theoretical calcBound is not validated." << std::endl;
      std::cout << "   Returns -1" << std::endl;
      return -1;
   }

   Term1 = (2.0 * Term1 + dim) * eBasis / dim;   
   Res = alpha * dim * log(Term1) - alpha * log(MM);

   return exp(Res);
}

//========================================================================

}
