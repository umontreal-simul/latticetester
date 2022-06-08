// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
// of Pierre L'Ecuyer at Université de Montréal.
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

#ifndef LATTICETESTER_NORMAMINKOWSKI_H
#define LATTICETESTER_NORMAMINKOWSKI_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

  /**
   * This class implements Minkowskiâ€™s theoretical **LOWER** bound on the length
   * of the shortest non-zero vector in a lattice. The length of a vector is
   * computed using the \f${\mathcal{L}}_2\f$ norm. The bounding lengths, for a
   * lattice containing \f$n\f$ points per unit volume in dimension \f$t\f$, 
   * are given by \f$\ell_t^* = \gamma_t^{1/2} n^{-1/t}\f$, where the 
   * \f$\gamma_t\f$ are the *Minkowski* lattice constants.
   */
  template<typename RealRed>
    class NormaMinkowski : public Normalizer<RealRed> {
      public:

        /**
         * Constructor for the bounds obtained for Minkowski lattices. The lattices
         * have \f$n\f$ points per unit volume, in all dimensions \f$\le t\f$. 
         * The bias factor `beta` \f$= \beta\f$ gives more weight to some of the 
         * dimensions. 
         * Note this class stores the log value of the density to handle larger values.
         * Restriction: \f$t \le48\f$.
         */
        NormaMinkowski (RealRed & logDensity, int t, double beta = 1);

        /**
         * Returns the value of the lattice constant \f$\gamma_j\f$ in
         * dimension \f$j\f$.
         */
        double getGamma (int j) const;
      private:

        /**
         * Lattice constants \f$\gamma_j\f$ for the Minkowski lattices in each
         * dimension \f$j\f$.
         */
        static const double m_gamma[1 + Normalizer<RealRed>::MAX_DIM];
    }; // End class NormaMinkowski

  //===========================================================================

  /*
   * This is (2/V_n)^(2/n) which seems wrong.
   * */
  template<typename RealRed>
    const double NormaMinkowski<RealRed>::m_gamma[ ] =
    {
      /* GamMinkowski[0] = */     0.00000000000000,
      /* GamMinkowski[1] = */     0.00000000000000,
      /* GamMinkowski[2] = */     1.04719756392815,
      /* GamMinkowski[3] = */     0.69062814882789,
      /* GamMinkowski[4] = */     0.66230588438641,
      /* GamMinkowski[5] = */     0.68895682802115,
      /* GamMinkowski[6] = */     0.73293650323622,
      /* GamMinkowski[7] = */     0.78407967437680,
      /* GamMinkowski[8] = */     0.83869147754813,
      /* GamMinkowski[9] = */     0.89517405737802,
      /* GamMinkowski[10] = */    0.95274948000853,
      /* GamMinkowski[11] = */    1.01100354118401,
      /* GamMinkowski[12] = */    1.06969923726369,
      /* GamMinkowski[13] = */    1.12869260847180,
      /* GamMinkowski[14] = */    1.18789173406828,
      /* GamMinkowski[15] = */    1.24723543612412,
      /* GamMinkowski[16] = */    1.30668158848970,
      /* GamMinkowski[17] = */    1.36620037299243,
      /* GamMinkowski[18] = */    1.42577021185918,
      /* GamMinkowski[19] = */    1.48537521407207,
      /* GamMinkowski[20] = */    1.54500351420012,
      /* GamMinkowski[21] = */    1.60464615768560,
      /* GamMinkowski[22] = */    1.66429633246555,
      /* GamMinkowski[23] = */    1.72394882699329,
      /* GamMinkowski[24] = */    1.78359964036430,
      /* GamMinkowski[25] = */    1.84324569710679,
      /* GamMinkowski[26] = */    1.90288463550686,
      /* GamMinkowski[27] = */    1.96251464853672,
      /* GamMinkowski[28] = */    2.02213436300702,
      /* GamMinkowski[29] = */    2.08174274687719,
      /* GamMinkowski[30] = */    2.14133903755970,
      /* GamMinkowski[31] = */    2.20092268604477,
      /* GamMinkowski[32] = */    2.26049331306133,
      /* GamMinkowski[33] = */    2.32005067447401,
      /* GamMinkowski[34] = */    2.37959463382297,
      /* GamMinkowski[35] = */    2.43912514042718,
      /* GamMinkowski[36] = */    2.49864221184901,
      /* GamMinkowski[37] = */    2.55814591979856,
      /* GamMinkowski[38] = */    2.61763637876553,
      /* GamMinkowski[39] = */    2.67711373682499,
      /* GamMinkowski[40] = */    2.73657816818394,
      /* GamMinkowski[41] = */    2.79602986712778,
      /* GamMinkowski[42] = */    2.85546904309690,
      /* GamMinkowski[43] = */    2.91489591667914,
      /* GamMinkowski[44] = */    2.97431071634673,
      /* GamMinkowski[45] = */    3.03371367580026,
      /* GamMinkowski[46] = */    3.09310503180911,
      /* GamMinkowski[47] = */    3.15248502245848,
      /* GamMinkowski[48] = */    3.21185388573071
    };



  /*=========================================================================*/

  template<typename RealRed>
    NormaMinkowski<RealRed>::NormaMinkowski (RealRed & logDensity, int t,
        double beta):
      Normalizer<RealRed> (logDensity, t, "Minkowski", L2NORM, beta)
    {
      if (t > this->MAX_DIM)
        throw std::invalid_argument("NormaMinkowski:   dimension > MAX_DIM");
      Normalizer<RealRed>::init (logDensity, beta);
    }


  /*=========================================================================*/

  template<typename RealRed>
    inline double NormaMinkowski<RealRed>::getGamma (int j) const
    {
      if (j < 1 || j > this->MAX_DIM)
        throw std::out_of_range("NormaMinkowski::getGamma");
      return m_gamma[j];
    }

  extern template class NormaMinkowski<double>;
  extern template class NormaMinkowski<NTL::RR>;

}

#endif
