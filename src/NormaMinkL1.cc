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

#include "latticetester/NormaMinkL1.h"
#include <cmath>

namespace LatticeTester
{

  const double NormaMinkL1::m_gamma0[] =
  {         //  GamMinkL1[t] = (t!)^{2/t}
    /* GamMinkL1[0] = */  0.00000000000000,
      /* GamMinkL1[1] = */  1.00000000000000,
      /* GamMinkL1[2] = */  2.00000000000000,
      /* GamMinkL1[3] = */  3.301927248894626,
      /* GamMinkL1[4] = */  4.898979485566356,
      /* GamMinkL1[5] = */  6.786916380543178,
      /* GamMinkL1[6] = */  8.962809493114328,
      /* GamMinkL1[7] = */  11.42450247602496,
      /* GamMinkL1[8] = */  14.17033543597957,
      /* GamMinkL1[9] = */  17.19898810749518,
      /* GamMinkL1[10] = */  20.50938353057181,
      /* GamMinkL1[11] = */  24.10062539497529,
      /* GamMinkL1[12] = */  27.97195541552144,
      /* GamMinkL1[13] = */  32.12272326936208,
      /* GamMinkL1[14] = */  36.55236475376784,
      /* GamMinkL1[15] = */  41.2603855160832,
      /* GamMinkL1[16] = */  46.24634867434057,
      /* GamMinkL1[17] = */  51.50986522412917,
      /* GamMinkL1[18] = */  57.05058648500341,
      /* GamMinkL1[19] = */  62.868198068697,
      /* GamMinkL1[20] = */  68.96241500217117,
      /* GamMinkL1[21] = */  75.33297774027004,
      /* GamMinkL1[22] = */  81.97964887293777,
      /* GamMinkL1[23] = */  88.90221038131129,
      /* GamMinkL1[24] = */  96.10046133233949,
      /* GamMinkL1[25] = */  103.5742159272685,
      /* GamMinkL1[26] = */  111.3233018382898,
      /* GamMinkL1[27] = */  119.3475587818154,
      /* GamMinkL1[28] = */  127.64683728756,
      /* GamMinkL1[29] = */  136.2209976308053,
      /* GamMinkL1[30] = */  145.069908901562,
      /* GamMinkL1[31] = */  154.1934481892748,
      /* GamMinkL1[32] = */  163.5914998656133,
      /* GamMinkL1[33] = */  173.2639549509683,
      /* GamMinkL1[34] = */  183.2107105527401,
      /* GamMinkL1[35] = */  193.4316693654913,
      /* GamMinkL1[36] = */  203.9267392246444,
      /* GamMinkL1[37] = */  214.6958327067107,
      /* GamMinkL1[38] = */  225.7388667701241,
      /* GamMinkL1[39] = */  237.0557624316304,
      /* GamMinkL1[40] = */  248.6464444739227,
      /* GamMinkL1[41] = */  260.5108411808306,
      /* GamMinkL1[42] = */  272.6488840968785,
      /* GamMinkL1[43] = */  285.0605078084658,
      /* GamMinkL1[44] = */  297.7456497442783,
      /* GamMinkL1[45] = */  310.7042499928673,
      /* GamMinkL1[46] = */  323.9362511355723,
      /* GamMinkL1[47] = */  337.44159809321,
      /* GamMinkL1[48] = */  351.220237985132,
  };

  /*=========================================================================*/

  double NormaMinkL1::calcGamma (int dim)
  {
    double gamma = 0.0;
    for (int i = 1; i <= dim; i++)
      gamma += log(i);
    gamma *= 2.0 / dim;

    return exp(gamma);
  }

  /*=========================================================================*/


  NormaMinkL1::NormaMinkL1 (RScal & logDensity, int t, double beta)
    : Normalizer (logDensity, t, "MinkL1", L1NORM, beta)
  {
    m_gamma = new double[t + 1];

    int t0 = t;
    if (t0 > MAX_DIM)
      t0 = MAX_DIM;
    for (int i = 0; i <= t0; i++)
      m_gamma[i] = m_gamma0[i];
    for (int i = t0 + 1; i <= t; i++)
      m_gamma[i] = calcGamma(i);

    Normalizer::init (logDensity, beta);
  }

  /*=========================================================================*/

  NormaMinkL1::~NormaMinkL1()
  {
    delete[] m_gamma;
  }
  /*=========================================================================*/


  inline double NormaMinkL1::getGamma (int j) const throw(std::out_of_range)
  {
    if (j < 1 || j > m_maxDim)
      throw std::out_of_range("NormaMinkL1::getGamma");
    return m_gamma[j];
  }

  /*=======================================================================*/

}
