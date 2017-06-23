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

NormaPalpha::NormaPalpha (const MScal & n, int alp, int s, NormType norm)
      : Normalizer (n, s, "Palpha", norm, 1.0)
{
   if (s > MAX_DIM)
      throw std::invalid_argument("NormaPalpha:   dimension > MAX_DIM");
   for (int j = 2; j <= s; j++)
      m_cst[j] = calcBound (alp, j);
   alpha = alp;
}


/*=========================================================================*/

int NormaPalpha::getAlpha () const
{
   return alpha;
}


/*=========================================================================*/

double NormaPalpha::calcBound (int alpha, int dim)
{
   const double eBasis = 2.71828182845904523536;
   double MM;
   conv (MM, m_n);
   MScal mm;
   conv (mm, m_n);
   if (dim <= 1) {
      std::cout << "calcBound:  dim < 2.   Returns -1" << std:: endl;
      return -1;
   }
   if (alpha <= 1) {
      std::cout << "calcBound:  alpha < 2.   Returns -1" << std:: endl;
      return -1;
   }
   int stat = IntFactor::isPrime (mm,0);
   if (stat != PRIME) {
      std::cout << "calcBound:  m is not prime.   Returns -1" << std:: endl;
      return -1;
   }
   double Term1 = log (MM);
   if (Term1 <= alpha*dim /(alpha - 1)) {
      std::cout << "calcBound:  log (M) <= alpha*dim /(alpha - 1) for dim = "
         << dim << ".   Returns -1" << std:: endl;
      return -1;
   }

   Term1 = (2.0 * Term1 + dim) * eBasis / dim;
   int j;

   double Res = Term1;
   for (j = 1; j < dim; ++j)
      Res *= Term1;

   Res = Term1 = Res / MM;
   for (j = 1; j < alpha; ++j)
      Res *= Term1;

   return Res;
}

//========================================================================

}
