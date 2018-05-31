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

/* Normalizer.h for ISO C++ */
#ifndef LATTICETESTER__NORMALIZER_H
#define LATTICETESTER__NORMALIZER_H
#include "latticetester/Util.h"
#include "latticetester/Const.h"

#include <cassert>
#include <string>


namespace LatticeTester {

  /**
   * Classes which inherit from this base class are used in implementing bounds
   * on the length of the shortest nonzero vector in a lattice
   * \cite mCON99a&thinsp;. These bounds are used to normalize the length of
   * the shortest vectors. Tight lower bounds are available for all dimensions
   * for many important cases. In most cases, the \f${\mathcal{L}}_2\f$ norm is
   * used to compute the length of vectors.
   *
   * For some figures of merit, no useful bounds are known to normalize the
   * length of the shortest vector. In these cases, this base class will be
   * used as normalizer since it simply sets all normalization constants to 1.
   * This is necessary because the tests compare the normalized values of the
   * merit when searching for good lattices.
   *
   */
  template<typename RedDbl>
    class Normalizer {

      public:
        static const int MAX_DIM = 48;

        /**
         * Constructor for the bounds. Deals with lattices having
         * \f$n\f$ points per unit volume, in all dimensions \f$\le t\f$. `Name` is
         * the name of the Normalizer. The bias factor `beta` \f$= \beta\f$ gives
         * more weight to some of the dimensions: taking \f$\beta< 1\f$ inflates the
         * figure of merit by \f$(1/\beta)^t\f$, thus weakening the requirements for
         * large \f$t\f$ in a worst-case figure of merit. One normally uses
         * \f$\beta= 1\f$.
         * Note that the log value of the density is stored (instead of the density 
         * itself) so it is easier to manipulate really large values of density.
         *
         * \remark **Richard:** Je crois que ce facteur `beta` devrait
         * disparaître car des poids beaucoup plus généraux sont maintenant
         * implantés dans les classes `*Weights`.
         */

        Normalizer (RedDbl & logDensity, int t, std::string Name,
            NormType norm = L2NORM, double beta = 1);

        /**
         * Constructor only used by the NormaPalpha class. It doesn't take any
         * log density argument. This only works for rank1 lattices, having m points 
         * per unit of volume (m being a prime number), normalized with NormaPalpha.
         */
        Normalizer (
            int t, std::string Name, NormType norm = L2NORM, double beta = 1);

        /**
         * Destructor.
         */
        virtual ~Normalizer ()
        { delete [] m_bounds; }

        /**
         * Initializes the bounds on the length of the shortest vector. The
         * lattices have \f$Density\f$ points per unit volume and the bias factor 
         * is `beta` for all dimensions \f$j\le\f$ `maxDim`.
         */
        virtual void init (RedDbl & logDensity, double beta);

        /**
         * Returns this object as a string.
         */
        std::string ToString () const;

        /**
         * Returns the norm associated with this object.
         */
        NormType getNorm () const
        { return m_norm; }

        /**
         * Sets the log-density associated with this object to `logDensity`.
         */
        void setLogDensity (RedDbl logDensity)
        { m_logDensity = logDensity; }

        /**
         * Returns the `logDensity` associated with this object.
         */
        RedDbl getLogDensity () const
        { return m_logDensity; }

        /**
         * Sets the norm associated with this object to `norm`.
         */
        void setNorm (NormType norm)
        { m_norm = norm; }

        /**
         * Returns the maximal dimension for this object.
         */
        int getDim () const
        { return m_maxDim; }

        /**
         * Returns the bound on the length of the shortest nonzero vector in
         * dimension \f$j\f$ as computed in Normalizer::init.
         */
        double getPreComputedBound (int j) const;

        /**
         * Calculates and returns the bound on the length of the shortest nonzero vector in
         * dimension \f$j\f$.
         */
        double getBound (int j) const;

        /**
         * Returns the value of the lattice constant \f$\gamma_j\f$ in
         * dimension \f$j\f$. For this base class, always returns 1.
         */
        virtual double getGamma (int j) const;

      protected:
        /**
         * Name of the normalizer.
         */
        std::string m_name;

        /**
         * Norm associated with this object.
         */
        NormType m_norm;

        /**
         * log of the density, ie log of the number of points of the lattice 
         * per unit of volume.
         */
        RedDbl m_logDensity;

        /**
         * Only elements 1 to <tt>m_maxDim</tt> (inclusive) of arrays are
         * defined.
         */
        int m_maxDim;

        /**
         * Beta factor.
         */
        double m_beta;

        /**
         * Contains the bounds on the length of the shortest nonzero vector in
         * the lattice in each dimension.
         */
        double *m_bounds;

      private:
        /**
         * Use of the copy-constructor is forbidden.
         */
        Normalizer (const Normalizer<RedDbl> &);

        /**
         * Use of assigment is forbidden.
         */
        Normalizer<RedDbl> & operator= (const Normalizer<RedDbl> &);

    }; // End class Normalizer

  //===========================================================================

  // template<typename RedDbl>
  //   const int Normalizer<RedDbl>::MAX_DIM;

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    Normalizer<RedDbl>::Normalizer (RedDbl & logDensity0, int maxDim,
        std::string name, NormType norm, double beta0) :
      m_name(name), m_norm(norm), m_logDensity(logDensity0), m_maxDim(maxDim),
      m_beta(beta0)
  {
    m_bounds = new double[maxDim + 1];
  }

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    Normalizer<RedDbl>::Normalizer (int maxDim, std::string name,
        NormType norm, double beta0) :
      m_name(name), m_norm(norm), m_logDensity(0), m_maxDim(maxDim),
      m_beta(beta0)
  {
    m_bounds = new double[maxDim + 1];
  }

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    void Normalizer<RedDbl>::init (RedDbl &logDensity0, double beta0)
    /*
     * Computes the vector m_bounds that corresponds to the upper bound for a lattice of
     * log-density \f$logDensity_0\f$.
     */

    {
      double x, y;
      double logBeta;
      m_logDensity = logDensity0;
      m_beta = beta0;

      y = 1.0;
      logBeta = log (m_beta);

      for (int j = 1; j <= m_maxDim; j++) {
        y =  1. / j;
        if (typeid(RedDbl) != typeid(double)) {
          // We have to convert to double
          x = 0.5 * log (
              getGamma(j)) + j * logBeta - y * NTL::conv<double>(logDensity0);
        } else {
          // We suppose we already have NTL::ZZ as integer type
          x = 0.5 * log (getGamma(j)) + j * logBeta - y * logDensity0;
        }
        // #if NTL_TYPES_CODE == 3
        //         x = 0.5 * log (getGamma(j)) + j * logBeta - y * conv<double>(logDensity0);
        // #else 
        //         x = 0.5 * log (getGamma(j)) + j * logBeta - y * logDensity0;
        // #endif
        //log calculation to handle large values of n

        m_bounds[j] = exp(x); 
      }
    }

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    std::string Normalizer<RedDbl>::ToString () const
    {
      std::ostringstream os;
      os << "-----------------------------\n"
        << "Content of Normalizer object:\n\n Normalizer = " << m_name;
      os << "\n n = " << exp(m_logDensity);
      os << "(log(n) = " << m_logDensity << ")";
      os << "\n beta = " << std::setprecision (4) << m_beta << "\n\n";

      //   os.setf(std::ios::left);
      os << std::setprecision (13);
      for (int t = 1; t <= m_maxDim; t++) {
        os << " Bound[" << std::setw(2) << std::right << t << "] = "
          << std::setw(14) << std::left << m_bounds[t] << "\n";
      }
      os << "\n";
      return os.str ();
    }


  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    double Normalizer<RedDbl>::getGamma (int) const
    {
      return 1.0;
    }


  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    double Normalizer<RedDbl>::getPreComputedBound (int j) const
    {
      assert (j >= 1 && j <= m_maxDim);
      return m_bounds[j];

      /*
remark: 
in the init method, the bounds are pre-computed for the dimensions of
the projection, and are accessible throw this function. But in the code
a call to function getBound (below) is made. This means the pre-computed
bounds are not used and the bounds are calculated again at each step with 
the function below. Could be improved.
*/
    }

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    double Normalizer<RedDbl>::getBound (int j) const
    {
      /*
         assert (j >= 1 && j <= m_maxDim);
         if (j >= 1 && j <= Normalizer::MAX_DIM)
         return getPreComputedBound (j);
         else {
         */

      double x,y;
      double logBeta;
      y = 1./j;
      logBeta = log(m_beta);
      if (typeid(RedDbl) != typeid(double)) {
        // We have to convert to double
        x = 0.5 * log (
            getGamma(j)) + j * logBeta - y * NTL::conv<double>(m_logDensity);
      } else {
        // We suppose we already have NTL::ZZ as integer type
        x = 0.5 * log (getGamma(j)) + j * logBeta - y * m_logDensity;
      }
      //log calculation to handle large values of n
      // #if NTL_TYPES_CODE == 3
      //       x = 0.5 * log (getGamma(j)) + j * logBeta - y * conv<double>(m_logDensity);
      // #else
      //       x = 0.5 * log (getGamma(j)) + j * logBeta - y * m_logDensity;
      // #endif
      //log calculation to handle large values of n

      return exp(x);
      //}
    }


} // end namespace LatticeTester

#endif
