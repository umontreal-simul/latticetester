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

#ifndef LATTICETESTER_NORMALAMINATED_H
#define LATTICETESTER_NORMALAMINATED_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

  /**
   * This class implements upper bounds on the lenght of the shortest nonzero
   * vector in a lattice. To obtain these bounds, this class contains hard-coded
   * values of an approximation of the Hermite's constant \f$\gamma_s\f$ with
   * the density of the laminated lattice of dimension `s`. These Hermite's
   * constants are stored in a table accessible via the getGamma(int) const
   * method.
   *
   * A laminated lattice in dimension `s` is, roughly speaking, the densest
   * possible packing of a laminated lattice of dimension `s-1` with copies of
   * itself. These lattices are not always the densest lattices available, but
   * are a more instinctive construction of a dense lattice. In \cite mCON99a,
   * Table 6.1 gives the determinant \f$\lambda_s\f$ that can be used to recover
   * the center density \f$\delta_s = \lambda_s^{-1/2}\f$ of the densest
   * laminated lattice in dimension \f$s\f$.
   * From there, we get
   * \f[
   *    \gamma_s = 4 \delta_s^{2/s}.
   * \f]
   * The number \f$\gamma_n\f$ is the number returned by calling `getGamma(n)`.
   * The init() method of this class inherited from `Normalizer` computes the
   * upper bound on the shortest non-zero vector of the lattice as
   * \f[
   *    \gamma_s^{1/2} n^{-1/s}.
   * \f]
   * Here \f$n = \exp(\text{logDensity})\f$ is the density of the lattice to be
   * analyzed passed as an argument to the constructor of this object.
   *
   * This class is to be used with the L2NORM (the Euclidian norm) exclusively.
   * Note this class stores the log value of the density to handle larger values.
   */
  template<typename RealRed>
    class NormaLaminated : public Normalizer {
      public:

        /**
         * Constructor for this class. Suppose we want to use this normalizer
         * on a lattice with it's basis in the lines of \f$V\f$ of dimension
         * \f$t\f$. We can call this constructor as `NormaLaminated(abs(det(V)), t)`.
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
        NormaLaminated (RealRed & logDensity, int t, double beta = 1);

        /**
         * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
         * in dimension \f$j\f$.
         */
        double getGamma (int j) const;
      private:

        /**
         * Lattice constants \f$\gamma_j\f$ for the laminated lattices in each
         * dimension \f$j\f$.
         */
        static const double m_gamma[1 + Normalizer::MAX_DIM];
    }; // End class NormaLaminated

  //=============================================================================

  /**
   * These laminated gamma constants are calculated as defined in 
   * Conway and Sloane book (Sphere packing, Lattices and groups) : 
   *    - equation (47) page 20 of chapter 1
   *    - table 6.1 page 158 of chapter 6
   */
  template<typename RealRed>
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

  template<typename RealRed>
    NormaLaminated::NormaLaminated (RealRed & logDensity, int t,
        double beta):
      Normalizer (logDensity, t, "Laminated", L2NORM, beta)
    {
      if (t > this->MAX_DIM)
        throw std::invalid_argument("NormaLaminated:   dimension > this->MAX_DIM");
      Normalizer::init (logDensity, beta);
    }


  /*=========================================================================*/

  template<typename RealRed>
    inline double NormaLaminated::getGamma (int j) const
    {
      if (j < 1 || j > this->MAX_DIM)
        throw std::out_of_range("NormaLaminated::getGamma");
      return m_gamma[j];
    }

  extern template class NormaLaminated<double>;
  extern template class NormaLaminated<NTL::RR>;

}

#endif
