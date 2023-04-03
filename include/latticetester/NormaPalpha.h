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

#ifndef LATTICETESTER_NORMAPALPHA_H
#define LATTICETESTER_NORMAPALPHA_H

#include "latticetester/Normalizer.h"
#include "latticetester/IntFactor.h"
//#include "NTL/ZZ.h"

namespace LatticeTester {

  /**
   * This class implements some theoretical bounds on the values of
   * \f$P_{\alpha}\f$ for a lattice (see class <tt>Palpha</tt>).
   */

    class NormaPalpha : public Normalizer {
      public:

        /**
         * Constructor for the bounds \f$B_{\alpha}(s)\f$ obtained for lattices, in
         * all dimensions \f$\le s\f$, where \f$\alpha= {}\f$<tt>alpha</tt>. The
         * lattices have rank \f$1\f$, with \f$m\f$ points per unit volume.
         * Restriction: \f$2 \le s \le48\f$, \f$\alpha\ge2\f$, and \f$m\f$ prime.
         */
        NormaPalpha (const std::int64_t m, int alpha, int s, NormType norm = L2NORM);

        /**
         * Computes and returns the bound \f$B_{\alpha}(s)\f$ given in
         * \cite mSLO94a&thinsp; (p. 83, Theorem 4.4). Given \f$s > 1\f$,
         * \f$\alpha> 1\f$, \f$m\f$ prime, and \f$m >
         * e^{\alpha s/(\alpha-1)}\f$, then there exists an integer vector
         * \f$\mathbf{a} \in\mathbb Z^s\f$ such that
         * \f[
         *   P_{\alpha}(s, \mathbf{a})  \le  B_{\alpha}(s)  = \frac{e}{s}^{\alpha s} \frac{(2\ln m + s)^{\alpha s}}{m^{\alpha}}.
         * \f]
         * If the conditions for the existence of the bound are not satisfied,
         * the function returns \f$-1\f$.
         */
        double calcBound (int alpha, int s);

        /**
         * Preventing the compiler from raising the warning: 
         * ''
         * latticetester/NormaPalpha.h:59:9: warning: 
         * 'LatticeTester::NormaPalpha::init' hides overloaded virtual function [-Woverloaded-virtual]
         * void init (int alpha);
         * ''
         */
        using Normalizer::computeBounds;

   
        /**
         * Initializes the bounds for the Palpha normalization.
         */
        void computeBounds (int alpha);

           /**
         * Initializes the bounds for the Palpha normalization.
         */
        void computeBounds ();



        /**
         * Returns the value of \f$\alpha\f$.
         */
        int getAlpha () const { return m_alpha; }

      private:

        /**
         * The value of \f$\m\f$.
         */
        std::int64_t m_m;

        /**
         * The value of \f$\alpha\f$.
         */
        int m_alpha;
    }; // End class NormaPalpha

  //===========================================================================

    NormaPalpha::NormaPalpha (const std::int64_t, int alpha, int maxDim,
        NormType norm):
      Normalizer (maxDim, "Palpha", norm)
    {
      if (maxDim > this->MAX_DIM)
        throw std::invalid_argument("NormaPalpha:   dimension > MAX_DIM");
      m_m = maxDim;
      m_alpha = alpha;
      computeBounds();
    }

  /*=========================================================================*/

        /*
     * Computes the vector m_bounds that corresponds to the upper bound for a 
     * rank 1 lattice of density \f$m\f$ (prime number). The bound does not exist
     * for dimension < 2.
     */
    void NormaPalpha::computeBounds ()
    {
      for (int j = 2; j <= this->m_maxDim; j++)
        this->m_bounds[j] = calcBound (m_alpha, j);
    }

  /*=========================================================================*/

    double NormaPalpha::calcBound (int alpha, int dim)
    {
      double Res;
      const double eBasis = 2.71828182845904523536;
      if (dim <= 1) {
        std::cout << "NormaPalpha::calcBound:  dim < 2.   Returns -1" << std:: endl;
        return -1;
      }
      if (alpha <= 1) {
        std::cout << "NormaPalpha::calcBound:  alpha < 2.   Returns -1" << std:: endl;
        return -1;
      }
      int stat = IntFactor<int64_t>::isPrime (m_m, 0);
      if (stat != PRIME) {
        std::cout << "NormaPalpha::calcBound:  m is not prime.   Returns -1" << std:: endl;
        return -1;
      }
      double Term1 = log (m_m);
      if (Term1 <= alpha*dim /(alpha - 1)) {
        std::cout << "NormaPalpha::calcBound:" << std::endl;
        std::cout << "   m < exp(alpha*dim/(alpha - 1)) for dim = " << dim 
          << std::endl;
        std::cout << "   Assumption required for existence of theoretical calcBound is not validated." << std::endl;
        std::cout << "   Returns -1" << std::endl;
        return -1;
      }
      Term1 = (2.0 * Term1 + dim) * eBasis / dim;   
      Res = alpha * dim * log(Term1) - alpha * log(m_m);
      return exp(Res);
    }

} // End namespace LatticeTester

#endif
