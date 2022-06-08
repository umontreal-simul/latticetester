// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
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

#ifndef LATTICETESTER_NORMALIZER_H
#define LATTICETESTER_NORMALIZER_H

#include "latticetester/Util.h"
#include "latticetester/Const.h"

#include <cassert>
#include <string>
#include <typeinfo>

namespace LatticeTester {

  /**
   * This is a base class for implementing normalization constants used in figures of merit,
   * to normalize the length of the shortest nonzero vector in either the primal or dual lattice.
   * These constants are based on upper bounds (or approximations) on the best possible length,
   * for given lattice density and dimension.  
   * The various subclasses of this class implement specific bounds.
   *
   * Given a lattice of density \f$ \eta \f$ in \f$t\f$ dimensions,
   * the Euclidean length \f$ d_t \f$ of a shortest nonzero lattice vector is upper-bounded as follows:
   * \f[
   *    d_t \le  d_t^*(\eta) \eqdef \gamma_t^{1/2} \eta^{-1/t},
   * \f]
   * where the constants  \f$ \gamma_t \f$ are known exactly only for  \f$ t\leq 8\f$.
   * For larger values of  \f$ t \f$, we use bounds or approximations of these constants,
   * which are defined in subclasses of `Normalizer`.
   * These bounds also hold for the  \f$ L^1\f$ lengths of the shortest vectors,
   * but with different constants  \f$ \gamma_t \f$.
   * In the dual lattice, the density is \f$1/\eta\f$ instead, and the Euclidean length
   * \f$ \ell_t \f$ of a shortest nonzero lattice vector then obeys:
   * \f[
   *    \ell_t \le  \ell_t^*(\eta) \eqdef \gamma_t^{1/2} \eta^{1/t}.
   * \f]
   * The bounds \f$ d_t^*(\eta)\f$ and \f$ \ell_t^*(\eta)\f$ are used to normalize
   * the values of \f$ d_t \f$ or \f$ \ell_t \f$ computed by the software for given lattices
   * to obtain standardized measures that lie between 0 and 1.
   * For the best lattices, these measures should be close to 1.
   *
   * Subclasses of this abstract class provide facilities to compute the bounds
   * \f$ d_t^*(\eta)\f$ or \f$ \ell_t^*(\eta)\f$ for a selected range of dimensions \f$ t \f$.
   * The first constructor will take this range (the maximum dimension) and the log of the density,
   * \f$ \log \eta \f$, as inputs. We work with the log density instead of the density itself
   * because the latter is sometimes extremely large or extremely small.
   * This constructor works fine when the density \f$ \eta\f$ is the same in all dimensions.
   *
   * But sometimes, the density may depend on the dimension.
   * This occurs for example for the lattice obtained from an MRG with modulus \f$ m\f$ and order \f$ k\f$,
   * whose density is typically \f$ m^k\f$ in \f$ s\ge k\f$ dimensions, and \f$ m^s\f$
   * in \f$ s < k\f$ dimensions.  The bounds can then be computed by taking this into account.
   * This is done by passing \f$ m\f$ and \f$ k\f$ as inputs to the constructor (in subclasses).
   *
   * In subclasses of this one, it is important to implement the
   * getGamma(int) const, getBound(int) const and init() methods and to call
   * the init() method in the constructor to pre-compute the bounds this class
   * will use for normalization.
   * The constructors here do not compute any bounds, they basically just reserve 
   * the space (an array) for the bounds.
   *
   * The prefered usage for this class is to declare a pointer to a Normalizer
   * and to instantiate a subclass with dynamically allocated memory:
   * \code{.cpp}
   * Normalizer<RScal>* norma;
   * norma = new NormaBestLat<RScal>(logDensity, t);
   * delete norma;
   * \endcode
   * 
   * Important: when making a search and examining millions of lattices, it is important 
   * NOT to construct a new Normalizer object and to recompute the constants for each lattice.
   */

  template<typename RealRed>
    class Normalizer {

      public:
        /**
         * The maximum dimension of the lattices for which this class can give
         * an upper bound.
         * */
        static const int MAX_DIM = 48;

        /**
         * This constructor creates a `Normalizer` by assuming that the density is
         * \f$\eta=\exp(\text{logDensity})\f$ in all dimensions \f$\leq t\f$.
         * The `name` parameter gives name that will be printed by the ToString() method.
         * The `norm` parameter is the `NormType` used by this object.
         */
        Normalizer (RealRed & logDensity, int t, std::string name,
            NormType norm = L2NORM);

        /**
         * This constructor assumes that the primal lattice has scaling factor \f$m\f$
         * and order \f$k\f$, so its density is \f$m^k\f$ for \f$t\geq k\f$, and cannot
         * exceed  \f$m^s\f$ for projections in \f$s < k\f$ dimensions.
         * The bounds \f$ d_t^*(\eta)\f$ (if `dualF=0`) or \f$ \ell_t^*(\eta)\f$ (if `dualF=1`)
         * will then be computed by assuming those densities.
         */
        Normalizer (Int & m, int k, bool dualF, int t, std::string name,
            NormType norm = L2NORM);

        /**
         * A constructor that does not take the density as an argument. The fields
         * are essentially the same as for `Normalizer(RealRed, int, std::string,
         * NormType, double)`.
         *
         * This is only used in the case of rank 1 lattices in the `NormaPalpha`
         * class with a prime density.   *** NEEEDED?  ***
         */
        Normalizer (int t, std::string Name, NormType norm = L2NORM);

        /**
         * Destructor.
         */
        virtual ~Normalizer ()
           { delete [] m_bounds; }

        /**
         * This is a method that will initialize the bounds this normalizer can
         * return. This will change the `logDensity` 
         * that are stored in this object. This will compute bounds for all
         * dimensions smaller than `t` (the parameter passed to the constructors)
         * that can be retrived with the getPreComputedBounds(int) method.
         */
        virtual void init (RealRed & logDensity);

        /**
         * Returns a string that describes this object.
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
        void setLogDensity (RealRed logDensity)
        { m_logDensity = logDensity; }

        /**
         * Returns the `logDensity` associated with this object.
         */
        RealRed getLogDensity () const
        { return m_logDensity; }

        /**
         * Sets the norm associated with this object to `norm`.
         */
        void setNorm (NormType norm)
        { m_norm = norm; }

        /**
         * Returns the maximal dimension for this object. This is the `t`
         * parameter of the constructors.
         */
        int getDim () const
        { return m_maxDim; }

        /**
         * Returns the bound for dimension `j` as computed in Normalizer::init().
         */
        double getPreComputedBound (int j) const;

        /**
         * Calculates and returns the bound on the length of the shortest nonzero vector in
         * dimension `j`.
         */
        virtual RealRed getBound (int j) const;

        /**
         * Returns the value of a lattice constant \f$\gamma\f$ in
         * dimension \f$j\f$. These constants can be used by subclasses to
         * implement the init() and the getBound(int) const methods. For this base
         * class, always returns 1.
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
        RealRed m_logDensity;

        /**
         * Only elements 1 to <tt>m_maxDim</tt> (inclusive) of m_bounds below
         * will be pre-computed. This stores the `t` parameter of the constructors.
         */
        int m_maxDim;

        /**
         * Contains the bounds on the length of the shortest nonzero vector in
         * the lattice in each dimension. This array is initialized by the init()
         * method, and it's values are returned with getPreComputedBound(int).
         */
        double *m_bounds;

      private:
        /**
         * Use of the copy-constructor is forbidden.
         */
        Normalizer (const Normalizer<RealRed> &);

        /**
         * Use of assigment is forbidden.
         */
        Normalizer<RealRed> & operator= (const Normalizer<RealRed> &);

    }; // End class Normalizer

  //===========================================================================

  template<typename RealRed>
    Normalizer<RealRed>::Normalizer (RealRed & logDensity0, int maxDim,
        std::string name, NormType norm) :
      m_name(name), m_norm(norm), m_logDensity(logDensity0), m_maxDim(maxDim),
  {
    m_bounds = new double[maxDim + 1];
  }

  /*-------------------------------------------------------------------------*/

  template<typename RealRed>
    Normalizer<RealRed>::Normalizer (int maxDim, std::string name,
        NormType norm) :
      m_name(name), m_norm(norm), m_logDensity(0), m_maxDim(maxDim),
  {
    m_bounds = new double[maxDim + 1];
  }

  /*-------------------------------------------------------------------------*/

  template<typename RealRed>
    void Normalizer<RealRed>::init (RealRed &logDensity0)
    /*
     * Computes the vector m_bounds that corresponds to the upper bounds on the 
     * best possible length of a shortest vector for a lattice of
     * log-density \f$logDensity_0\f$, in all dimensions up to maxDim.
     * This method assumes the same density in all dimensions. 
     */
    {
      double x;
      m_logDensity = logDensity0;
      for (int j = 1; j <= m_maxDim; j++) {
        x = 0.5*log (getGamma(j)) - (1.0/j) * NTL::conv<double>(logDensity0);
        m_bounds[j] = exp(x); 
      }
    }

  /*-------------------------------------------------------------------------*/

  template<typename RealRed>
    std::string Normalizer<RealRed>::ToString () const
    {
      std::ostringstream os;
      os << "-----------------------------\n"
        << "Content of Normalizer object:\n\n Normalizer = " << m_name;
      os << "\n n = " << exp(m_logDensity);
      os << "(log(n) = " << m_logDensity << ")";
      os << "\n\n";

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

  template<typename RealRed>
    double Normalizer<RealRed>::getGamma (int) const
    {
	  // In this abstract class, the gamma_t's are undefined.
      return -1.0;
    }


  /*-------------------------------------------------------------------------*/

  template<typename RealRed>
    double Normalizer<RealRed>::getPreComputedBound (int j) const
    {
      assert (j >= 1 && j <= m_maxDim);
      return m_bounds[j];

      /*
       * remark: 
       * in the init method, the bounds are pre-computed for the dimensions of
       * the projection, and are accessible through this function. But in the code
       * a call to function getBound (below) is made. This means the pre-computed
       * bounds are not used and the bounds are calculated again at each step with 
       * the function below. MUST be improved.
       */
    }

  /*-------------------------------------------------------------------------*/

  template<typename RealRed>
    RealRed Normalizer<RealRed>::getBound (int j) const
    {
      /*
         assert (j >= 1 && j <= m_maxDim);
         if (j >= 1 && j <= Normalizer::MAX_DIM)
         return getPreComputedBound (j);
         else {
         */
      if (getGamma(j) < 0) return RealRed(1.0);
      RealRed x,y;
      y = NTL::inv(RealRed(j));
      x = 0.5*NTL::log (getGamma(j)) - y*m_logDensity;
      return exp(x);
    }

  extern template class Normalizer<double>;
  extern template class Normalizer<NTL::RR>;

} // end namespace LatticeTester

#endif
