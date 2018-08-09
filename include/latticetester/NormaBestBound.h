// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
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

/* NormaBestBound.h for ISO C++ */
#ifndef LATTICETESTER__NORMAROGERS_H
#define LATTICETESTER__NORMAROGERS_H
#include "latticetester/Normalizer.h"
#include <stdexcept>


namespace LatticeTester {

  /**
   * \Remark The current implementation does not implement the bounds for now
   *
   * This class implements upper bounds on the lenght of the shortest nonzero
   * vector in a lattice. To obtain these bounds, this class contains hard-coded
   * values of an approximation of the Hermite's constants \f$\gamma_s\f$
   * that are the best upper bounds that we are aware of. These Hermite's
   * constants are stored in a table accessible via the getGamma(int) const
   * method.
   *
   * For dimensions 0 through 8, the Hermite's constant are known exactly.
   * For dimensions 9 through 36, see https://doi.org/10.4007/annals.2003.157.689.
   * For dimensions 37 through 48 this uses Rogers's bound.
   *
   * This class is to be used with the L2NORM (the Euclidian norm) exclusively.
   * Note this class stores the log value of the density to handle larger values.
   */
  template<typename RedDbl>
    class NormaBestBound : public Normalizer<RedDbl> {
      public:

        /**
         * Constructor for this class. Suppose we want to use this normalizer
         * on a lattice with it's basis in the lines of \f$V\f$ of dimension
         * \f$t\f$. We can call this constructor as `NormaBestBound(abs(det(V)), t)`.
         * getPreComputedBound(t) will then return an upper bound on the
         * lenght of the shortest non-zero vector in dimension `t`. In the case
         * where the lattice also has the same density in lower dimensions than
         * `t`, pre-computed bounds will also be available.
         *
         * The bias factor `beta` gives more or less weight to some of the
         * dimensions (see Normalizer for details). It is recommended to keep it
         * at its default value because its usage is deprecated.
         * 
         * There is a restriction for `t` to be \f$\le48\f$.
         */
        NormaBestBound (RedDbl & logDensity, int t, double beta = 1);

        /**
         * Destructor.
         */
        ~NormaBestBound();

        /**
         * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
         * in dimension \f$j\f$.
         */
        double getGamma (int j) const;
      private:

        /**
         * Precomputed lattice constants \f$\gamma_j\f$ for the Rogers bounds
         * in each dimension \f$j \le48\f$.
         */
        static const double m_gamma[1 + Normalizer<RedDbl>::MAX_DIM];

    }; // End class NormaBestBound

  //===========================================================================

  /**
   * These Rogers gamma constants are calculated as defined in 
   * Conway and Sloane book (Sphere packing, Lattices and groups) : 
   *    - equation (47) page 20 of chapter 1
   *    - table 1.2 page 15 of chapter 1
   */
  template<typename RedDbl>
          const double NormaBestBound<RedDbl>::m_gamma[ ] =
  {
    /* GamBestBound[0] = */    1.,
    /* GamBestBound[1] = */    1.,
    /* GamBestBound[2] = */    1.,
    /* GamBestBound[3] = */    1.,
    /* GamBestBound[4] = */    1.,
    /* GamBestBound[5] = */    1.,
    /* GamBestBound[6] = */    1.,
    /* GamBestBound[7] = */    1.,
    /* GamBestBound[8] = */    1.,
    /* GamBestBound[9] = */    1.,
    /* GamBestBound[10] = */    1.,
    /* GamBestBound[11] = */    1.,
    /* GamBestBound[12] = */    1.,
    /* GamBestBound[13] = */    1.,
    /* GamBestBound[14] = */    1.,
    /* GamBestBound[15] = */    1.,
    /* GamBestBound[16] = */    1.,
    /* GamBestBound[17] = */    1.,
    /* GamBestBound[18] = */    1.,
    /* GamBestBound[19] = */    1.,
    /* GamBestBound[20] = */    1.,
    /* GamBestBound[21] = */    1.,
    /* GamBestBound[22] = */    1.,
    /* GamBestBound[23] = */    1.,
    /* GamBestBound[24] = */    1.,
    /* GamBestBound[25] = */    1.,
    /* GamBestBound[26] = */    1.,
    /* GamBestBound[27] = */    1.,
    /* GamBestBound[28] = */    1.,
    /* GamBestBound[29] = */    1.,
    /* GamBestBound[30] = */    1.,
    /* GamBestBound[31] = */    1.,
    /* GamBestBound[32] = */    1.,
    /* GamBestBound[33] = */    1.,
    /* GamBestBound[34] = */    1.,
    /* GamBestBound[35] = */    1.,
    /* GamBestBound[36] = */    1.,
    /* GamBestBound[37] = */    1.,
    /* GamBestBound[38] = */    1.,
    /* GamBestBound[39] = */    1.,
    /* GamBestBound[40] = */    1.,
    /* GamBestBound[41] = */    1.,
    /* GamBestBound[42] = */    1.,
    /* GamBestBound[43] = */    1.,
    /* GamBestBound[44] = */    1.,
    /* GamBestBound[45] = */    1.,
    /* GamBestBound[46] = */    1.,
    /* GamBestBound[47] = */    1.,
    /* GamBestBound[48] = */    1.
  };


  /*=======================================================================*/

  template<typename RedDbl>
    NormaBestBound<RedDbl>::NormaBestBound (RedDbl & logDensity, int t, double beta)
    : Normalizer<RedDbl> (logDensity, t, "Rogers", L2NORM, beta)
    {
      Normalizer<RedDbl>::init (logDensity, beta);
    }

  /*=========================================================================*/

  template<typename RedDbl>
    NormaBestBound<RedDbl>::~NormaBestBound()
    {
      delete[] m_gamma;
    }

  /*=========================================================================*/

  template<typename RedDbl>
    inline double NormaBestBound<RedDbl>::getGamma (int j) const
    {
      if (j < 1 || j > this->m_maxDim)
        throw std::out_of_range("NormaBestBound::getGamma");
      return m_gamma[j];
    }

}
#endif
