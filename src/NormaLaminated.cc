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

#include "latticetester/NormaLaminated.h"

namespace LatticeTester
{

/**
 * These laminated gamma constants are calculated as defined in 
 * Conway and Sloane book (Sphere packing, Lattices and groups) : 
 *    - equation (47) page 20 of chapter 1
 *    - table 6.1 page 158 of chapter 6
 */

const double NormaLaminated::m_gamma[] =
   {
      /* Gamma[0] = */    0.00000000000000,
      /* Gamma[1] = */    1.00000000000000,
      /* Gamma[2] = */    1.15470053837925,
      /* Gamma[3] = */    1.25992104989487,
      /* Gamma[4] = */    1.41421356237309,
      /* Gamma[5] = */    1.51571656651040,
      /* Gamma[6] = */    1.66536635531121,
      /* Gamma[7] = */    1.81144732852781,
      /* Gamma[8] = */    2.00000000000000,
      /* Gamma[9] = */    2.00000000000000,
      /* Gamma[10] = */   2.05837201792952,
      /* Gamma[11] = */   2.13008217887993,
      /* Gamma[12] = */   2.24492409661875,
      /* Gamma[13] = */   2.34692092000925,
      /* Gamma[14] = */   2.48864391982238,
      /* Gamma[15] = */   2.63901582154579,
      /* Gamma[16] = */   2.82842712474619,
      /* Gamma[17] = */   2.88668115405991,
      /* Gamma[18] = */   2.98682599936104,
      /* Gamma[19] = */   3.09851928453331,
      /* Gamma[20] = */   3.24900958542494,
      /* Gamma[21] = */   3.39145596751014,
      /* Gamma[22] = */   3.57278019514216,
      /* Gamma[23] = */   3.76602735259556,
      /* Gamma[24] = */   4.00000000000000,
      /* Gamma[25] = */   3.89061978964914,
      /* Gamma[26] = */   3.83450381188673,
      /* Gamma[27] = */   3.79980642833440,
      /* Gamma[28] = */   3.80678061204248,
      /* Gamma[29] = */   3.81328532378654,
      /* Gamma[30] = */   3.85616802789424,
      /* Gamma[31] = */   3.91155414534173,
      /* Gamma[32] = */   4.00000000000000,
      /* Gamma[33] = */   4.00000000000000,
      /* Gamma[34] = */   4.03398853947441,
      /* Gamma[35] = */   4.08000643768480,
      /* Gamma[36] = */   4.15703690412737,
      /* Gamma[37] = */   4.23124164832276,
      /* Gamma[38] = */   4.33546037046114,
      /* Gamma[39] = */   4.45012590438595,
      /* Gamma[40] = */   4.59479341998814,
      /* Gamma[41] = */   4.65735933059734,
      /* Gamma[42] = */   4.75016320908659,
      /* Gamma[43] = */   4.85364895344187,
      /* Gamma[44] = */   4.98703323713956,
      /* Gamma[45] = */   5.11791276058506,
      /* Gamma[46] = */   5.27922773469898,
      /* Gamma[47] = */   5.45208672393488,
      /* Gamma[48] = */   5.65685424949238
   };



/*=========================================================================*/


NormaLaminated::NormaLaminated (RScal & logDensity, int t, double beta)
      : Normalizer (logDensity, t, "Laminated", L2NORM, beta)
{
   if (t > MAX_DIM)
      throw std::invalid_argument("NormaLaminated:   dimension > MAX_DIM");
   Normalizer::init (logDensity, beta);
}


/*=========================================================================*/


inline double NormaLaminated::getGamma (int j) const
throw(std::out_of_range)
{
   if (j < 1 || j > MAX_DIM)
      throw std::out_of_range("NormaLaminated::getGamma");
   return m_gamma[j];
}

/*=======================================================================*/

}
