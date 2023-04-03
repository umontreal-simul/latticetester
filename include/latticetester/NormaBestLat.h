// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
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

#ifndef LATTICETESTER_NORMABESTLAT_H
#define LATTICETESTER_NORMABESTLAT_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

  /**
   * This Normalizer class implements upper bounds on the length of the shortest nonzero
   * vector in a lattice. The Hermite constants \f$\gamma_s\f$ are approximated
   * by using the largest known density values for lattices \cite mCON99a.
   * These approximations are stored in a table accessible via the
   * `getGamma(int)` const method. In \cite mCON99a, Table I.1(a) gives the center
   * density \f$\delta_s\f$ of the densest known packings in dimension \f$s\f$.
   * We then have that
   * \f[
   *    \gamma_s = 4 \delta_s^{2/s}.
   * \f]
   * This class is to be used with the L2NORM (the Euclidian norm) exclusively.
   */

    class NormaBestLat : public Normalizer {
      public:

        /**
         * Constructs a `NormaBestLat` for up to `maxDim` dimensions, by assuming that the
         * log density is `logDensity` in all dimensions.
         * Restriction: `maxDim`\f$ \le 48\f$.
         */
        NormaBestLat (double logDensity, int maxDim);

    	/**
    	 * This constructor assumes that the primal lattice has scaling factor \f$m\f$
    	 * and order \f$k\f$, so its density is \f$m^k\f$ for \f$t\geq k\f$, and cannot
    	 * exceed  \f$m^s\f$ for projections in \f$s < k\f$ dimensions.
    	 */
    	NormaBestLat (double logm, int k, int maxDim);

    	/**
         * Constructs a `NormaBestLat` for up to `maxDim` dimensions, without computing the bounds.
         * Restriction: `maxDim`\f$ \le 48\f$.
         */
        NormaBestLat (int maxDim);

        /**
         * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
         * in dimension \f$j\f$.
         */
        double getGamma (int j) const;
        
      private:

        /**
         * Lattice constants \f$\gamma_j\f$ for the most general lattices in each
         * dimension \f$j\f$.
         */
        static const double m_gamma[1 + Normalizer::MAX_DIM];
    }; // End class NormaBestLat

    const double NormaBestLat::m_gamma[] =
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

    NormaBestLat::NormaBestLat (double logDensity, int maxDim)
    : Normalizer (maxDim, "BestLat", L2NORM)
    {
      if (maxDim > this->MAX_DIM)
        throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
      Normalizer::computeBounds (logDensity);
    }

    /*=========================================================================*/

      NormaBestLat::NormaBestLat (double logm, int k, int maxDim)
      : Normalizer (maxDim, "BestLat", L2NORM)
      {
        if (maxDim > this->MAX_DIM)
          throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
        Normalizer::computeBounds (logm, k);
      }

    /*=========================================================================*/

    NormaBestLat::NormaBestLat (int maxDim)
    : Normalizer (maxDim, "BestLat", L2NORM)
    {
      if (maxDim > this->MAX_DIM)
        throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
    }

  /*=========================================================================*/

    inline double NormaBestLat::getGamma (int j) const
    {
      if (j < 1 || j > this->MAX_DIM)
        throw std::out_of_range("NormaBestLat::getGamma");
      return m_gamma[j];
    }

} // End namespace LatticeTester

#endif

