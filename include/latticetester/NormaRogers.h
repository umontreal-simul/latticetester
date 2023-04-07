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

#ifndef LATTICETESTER_NORMAROGERS_H
#define LATTICETESTER_NORMAROGERS_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

  /**
   * This class implements upper bounds on the lenght of the shortest nonzero
   * vector in a lattice. To obtain these bounds, this class contains hard-coded
   * values of an approximation of the Hermite's constants \f$\gamma_s\f$
   * calculated with Rogers's bound. These Hermite's constants are stored in a
   * table accessible via the getGamma(int) const method.
   *
   * Rogers bound has been introduced by Rogers in The Packing of Equal Spheres
   * in 1958 (citation to add in bibliography). This is a classical bound that
   * is implemented mainly for historical reasons. For thighter bounds, please
   * use NormaBestBound. Note that, since the value of \f$\gamma_n\f$ is known
   * exactly for \f$n \leq 8\f$, the Hermite's constant for these \f$n\f$ are
   * **not** upper bounds.
   * 
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
  template<typename RedDbl>
    class NormaRogers : public Normalizer<RedDbl> {
      public:

        /**
         * Constructor for this class. Suppose we want to use this normalizer
         * on a lattice with it's basis in the lines of \f$V\f$ of dimension
         * \f$t\f$. We can call this constructor as `NormaRogers(abs(det(V)), t)`.
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
        NormaRogers (RedDbl & logDensity, int t, double beta = 1);

        /**
         * Destructor.
         */
        ~NormaRogers();

        /**
         * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
         * in dimension \f$j\f$.
         */
        double getGamma (int j) const;
      private:

        /**
         * The lattice constants \f$\gamma_j\f$ are the Rogers bounds in each
         * dimension \f$j\f$.
         */
        double *m_gamma;

        /**
         * Precomputed lattice constants \f$\gamma_j\f$ for the Rogers bounds
         * in each dimension \f$j \le48\f$.
         */
        static const double m_gamma0[1 + Normalizer<RedDbl>::MAX_DIM];

        /**
         * Computes the Rogers bound in dimension \f$d\f$.
         */
        double calcGamma (int d);
    }; // End class NormaRogers

  //===========================================================================

  /**
   * This is not all the story, this has to be updated.
   *
   * These Rogers gamma constants are calculated as defined in 
   * Conway and Sloane book (Sphere packing, Lattices and groups) : 
   *    - equation (47) page 20 of chapter 1
   *    - table 1.2 page 15 of chapter 1
   */
  template<typename RedDbl>
          const double NormaRogers<RedDbl>::m_gamma0[ ] =
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

  template<typename RedDbl>
    double NormaRogers<RedDbl>::calcGamma (int dim)
    {
      static const double pi = 3.1415926535897932384;
      static const double e = 2.7182818284590452353;
      static const double t = log2(e / sqrt(pi));
      static const double s = 4.0 * e * pi;
      const double dimr = dim;
      double r;

      r = 0.5 * dimr * log2 (dimr / s) + 1.5 * log2 (dimr) - t + 5.25 / (dimr + 2.5);
      r = 4 * exp2(2 * r / dimr);

      /*
       * remark: 
       * Why "r = 4 * ..." and not "r = 2 * ..." as described in Pierre's 
       * course IFT6561 page 273 ?
       */

      return r;
    }


  /*=========================================================================*/

  template<typename RedDbl>
    NormaRogers<RedDbl>::NormaRogers (RedDbl & logDensity, int t, double beta)
    : Normalizer<RedDbl> (logDensity, t, "Rogers", L2NORM, beta)
    {
      m_gamma = new double[t + 1];

      int t0 = t;
      if (t0 > this->MAX_DIM)
        t0 = this->MAX_DIM;
      for (int i = 0; i <= t0; i++)
        m_gamma[i] = m_gamma0[i];
      for (int i = t0 + 1; i <= t; i++)
        m_gamma[i] = calcGamma(i);

      Normalizer<RedDbl>::init (logDensity, beta);
    }

  /*=========================================================================*/

  template<typename RedDbl>
    NormaRogers<RedDbl>::~NormaRogers()
    {
      delete[] m_gamma;
    }

  /*=========================================================================*/

  template<typename RedDbl>
    inline double NormaRogers<RedDbl>::getGamma (int j) const
    {
      if (j < 1 || j > this->m_maxDim)
        throw std::out_of_range("NormaRogers::getGamma");
      return m_gamma[j];
    }

  extern template class NormaRogers<double>;
  extern template class NormaRogers<NTL::RR>;

}

#endif
