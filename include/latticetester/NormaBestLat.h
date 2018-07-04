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

/* NormaBestLat.h for ISO C++ */
#ifndef LATTICETESTER__NORMABESTLAT_H
#define LATTICETESTER__NORMABESTLAT_H
#include "latticetester/Normalizer.h"
#include <stdexcept>


namespace LatticeTester {

  /**
   * This class implements the *best* theoretical bounds on the length of the
   * shortest nonzero vector in a lattice, based on the densest sphere packing
   * in lattices. The length of a vector is computed using the \f${\mathcal{L}}_2\f$
   * norm. The bounding lengths for a lattice containing \f$n\f$ points per 
   * unit volume in dimension \f$t\f$, are given by 
   * \f$\ell_t^* = \gamma_t^{1/2} n^{-1/t}\f$, where the \f$\gamma_t\f$ are 
   * the lattice constants for the *best* known lattices \cite mCON99a&thinsp;.
   * Note this class stores the log value of the density to handle larger values.
   */
  template<typename RedDbl>
    class NormaBestLat : public Normalizer<RedDbl> {
      public:

        /**
         * Constructor for the best bounds obtained for lattices. The lattices have
         * \f$Density\f$ points per unit volume, for all dimensions \f$\le t\f$. The bias 
         * factor `beta` \f$= \beta\f$ gives more weight to some of the dimensions. 
         * Restriction: \f$t \le48\f$.
         */
        NormaBestLat (RedDbl & logDensity, int t, double beta = 1);

        /**
         * Returns the value of the lattice constant \f$\gamma_j\f$ in
         * dimension \f$j\f$.
         */
        double getGamma (int j) const;
      private:

        /**
         * Lattice constants \f$\gamma_j\f$ for the most general lattices in each
         * dimension \f$j\f$.
         */
        static const double m_gamma[1 + Normalizer<RedDbl>::MAX_DIM];
    }; // End class NormaBestLat

  template<typename RedDbl>
    const double NormaBestLat<RedDbl>::m_gamma[] =
    {
      /* GamBestLat[0] = */   0.0,
      /* GamBestLat[1] = */   1.0,
      /* GamBestLat[2] = */   1.1547005383793,
      /* GamBestLat[3] = */   1.2599210498949,
      /* GamBestLat[4] = */   1.4142135623731,
      /* GamBestLat[5] = */   1.5157165665104,
      /* GamBestLat[6] = */   1.6653663553112,
      /* GamBestLat[7] = */   1.8114473285278,
      /* GamBestLat[8] = */   2.0,
      /* GamBestLat[9] = */   2.0,
      /* GamBestLat[10] = */  2.0583720179295,
      /* GamBestLat[11] = */  2.140198065871,
      /* GamBestLat[12] = */  2.3094010767585,
      /* GamBestLat[13] = */  2.3563484301065,
      /* GamBestLat[14] = */  2.4886439198224,
      /* GamBestLat[15] = */  2.6390158215458,
      /* GamBestLat[16] = */  2.8284271247462,
      /* GamBestLat[17] = */  2.8866811540599,
      /* GamBestLat[18] = */  2.986825999361,
      /* GamBestLat[19] = */  3.0985192845333,
      /* GamBestLat[20] = */  3.2490095854249,
      /* GamBestLat[21] = */  3.3914559675101,
      /* GamBestLat[22] = */  3.5727801951422,
      /* GamBestLat[23] = */  3.7660273525956,
      /* GamBestLat[24] = */  4.0,
      /* GamBestLat[25] = */  3.8906197896491,
      /* GamBestLat[26] = */  3.8345038118867,
      /* GamBestLat[27] = */  3.8405094116889,
      /* GamBestLat[28] = */  3.8858143186426,
      /* GamBestLat[29] = */  3.8513016372256,
      /* GamBestLat[30] = */  3.890079350856,
      /* GamBestLat[31] = */  4.0493929444608,
      /* GamBestLat[32] = */  4.2426406871193,
      /* GamBestLat[33] = */  4.1983166567599,
      /* GamBestLat[34] = */  4.1923458021689,
      /* GamBestLat[35] = */  4.2448520933335,
      /* GamBestLat[36] = */  4.3453285925836,
      /* GamBestLat[37] = */  4.2312416483228,
      /* GamBestLat[38] = */  4.4626316710462,
      /* GamBestLat[39] = */  4.5228010665648,
      /* GamBestLat[40] = */  4.6661029086385,
      /* GamBestLat[41] = */  4.8084724701927,
      /* GamBestLat[42] = */  4.9619948528877,
      /* GamBestLat[43] = */  5.1129393316586,
      /* GamBestLat[44] = */  5.2613041578794,
      /* GamBestLat[45] = */  5.4070956951517,
      /* GamBestLat[46] = */  5.5851474972462,
      /* GamBestLat[47] = */  5.7755698526865,
      /* GamBestLat[48] = */  6.0
    };



  /*=========================================================================*/

  template<typename RedDbl>
    NormaBestLat<RedDbl>::NormaBestLat (RedDbl & logDensity, int t, double beta)
    : Normalizer<RedDbl> (logDensity, t, "BestLat", L2NORM, beta)
    {
      if (t > this->MAX_DIM)
        throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
      Normalizer<RedDbl>::init (logDensity, beta);
    }


  /*=========================================================================*/

  template<typename RedDbl>
    inline double NormaBestLat<RedDbl>::getGamma (int j) const
    {
      if (j < 1 || j > this->MAX_DIM)
        throw std::out_of_range("NormaBestLat::getGamma");
      return m_gamma[j];
    }

} // End namespace LatticeTester
#endif
