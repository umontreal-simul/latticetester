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

#include "latticetester/NormaRogers.h"

#include <cmath>

namespace LatticeTester
{

/**
 * These Rogers gamma constants are calculated as defined in 
 * Conway and Sloane book (Sphere packing, Lattices and groups) : 
 *    - equation (47) page 20 of chapter 1
 *    - table 1.2 page 15 of chapter 1
 */

const double NormaRogers::m_gamma0[ ] =
   {
      /* GamRogers[0] = */    0.0,
      /* GamRogers[1] = */    1.0,
      /* GamRogers[2] = */    1.15472,
      /* GamRogers[3] = */    1.2972925427178,
      /* GamRogers[4] = */    1.4492480809026,
      /* GamRogers[5] = */    1.5916002961306,
      /* GamRogers[6] = */    1.7315537290705,
      /* GamRogers[7] = */    1.8696088064088,
      /* GamRogers[8] = */    2.006052470232,
      /* GamRogers[9] = */    2.1411671718503,
      /* GamRogers[10] = */    2.2751349805586,
      /* GamRogers[11] = */    2.4081055004162,
      /* GamRogers[12] = */    2.5401903576369,
      /* GamRogers[13] = */    2.671499016465,
      /* GamRogers[14] = */    2.8020630856483,
      /* GamRogers[15] = */    2.9320505407083,
      /* GamRogers[16] = */    3.0614381882081,
      /* GamRogers[17] = */    3.1903070449466,
      /* GamRogers[18] = */    3.318714864331,
      /* GamRogers[19] = */    3.4466883426431,
      /* GamRogers[20] = */    3.5742655437525,
      /* GamRogers[21] = */    3.7014670196163,
      /* GamRogers[22] = */    3.8283274848644,
      /* GamRogers[23] = */    3.9548705630986,
      /* GamRogers[24] = */    4.0811157647776,
      /* GamRogers[25] = */    4.2071543016103,
      /* GamRogers[26] = */    4.3328598061492,
      /* GamRogers[27] = */    4.4583196677731,
      /* GamRogers[28] = */    4.583548484021,
      /* GamRogers[29] = */    4.7085595260287,
      /* GamRogers[30] = */    4.8333649016765,
      /* GamRogers[31] = */    4.9579756932973,
      /* GamRogers[32] = */    5.0824020747592,
      /* GamRogers[33] = */    5.2066534116689,
      /* GamRogers[34] = */    5.3307383476426,
      /* GamRogers[35] = */    5.454664878987,
      /* GamRogers[36] = */    5.5784404196715,
      /* GamRogers[37] = */    5.7020718581143,
      /* GamRogers[38] = */    5.8255656070255,
      /* GamRogers[39] = */    5.9489276473284,
      /* GamRogers[40] = */    6.0721635670068,
      /* GamRogers[41] = */    6.1952785955803,
      /* GamRogers[42] = */    6.3182776348,
      /* GamRogers[43] = */    6.4411652860615,
      /* GamRogers[44] = */    6.5639458749555,
      /* GamRogers[45] = */    6.6866234733141,
      /* GamRogers[46] = */    6.8092019190592,
      /* GamRogers[47] = */    6.9316848341156,
      /* GamRogers[48] = */    7.0540756406128
   };


/*=======================================================================*/


double NormaRogers::calcGamma (int dim)
{
   static const double pi = 3.1415926535897932384;
   static const double e = 2.7182818284590452353;
   static const double t = log2(e / sqrt(pi));
   static const double s = 4.0 * e * pi;
   const double dimr = dim;
   double r;

   r = 0.5 * dimr * log2 (dimr / s) + 1.5 * log2 (dimr) - t + 5.25 / (dimr + 2.5);
   r = 4 * exp2(2 * r / dimr);
   // PW_TODO : pourquoi *4 et pas *2 ? Ça colle pas avec le cours de Pierre page 289
   return r;
}


/*=========================================================================*/

NormaRogers::NormaRogers (const RScal & n, int t, double beta)
      : Normalizer (n, t, "Rogers", L2NORM, beta)
{
   m_gamma = new double[t + 1];
   int t0 = t;
/*   if (t0 > MAX_DIM)
      t0 = MAX_DIM;*/
   int i;
   for (i = 0; i <= t0; i++)
      m_gamma[i] = m_gamma0[i];
   for (i = t0 + 1; i <= t; i++)
      m_gamma[i] = calcGamma(i);

   Normalizer::init (n, beta);
}

/*=========================================================================*/

NormaRogers::~NormaRogers()
{
   delete[] m_gamma;
}

/*=========================================================================*/


inline double NormaRogers::getGamma (int j) const throw(std::out_of_range)
{
   if (j < 1 || j > m_maxDim) //fred
      throw std::out_of_range("NormaRogers::getGamma");
   return m_gamma[j];
}

/*=======================================================================*/

}
