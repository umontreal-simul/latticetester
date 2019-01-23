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

#ifndef LATTICETESTER_NORMALIZER_H
#define LATTICETESTER_NORMALIZER_H

#include "latticetester/Util.h"
#include "latticetester/Const.h"

#include <cassert>
#include <string>
#include <typeinfo>

namespace LatticeTester {

  /**
   * Classes which inherit from this base class are used in implementing bounds
   * on the length of the shortest nonzero vector in a lattice.
   *
   * Given a lattice in dimension \f$t\f$, it is possible to center
   * non-interlapping spheres of radius \f$ d_t \f$ where \f$ d_t \f$ is the
   * lenght of the shortest vector in the lattice. It turns out that the
   * proportion of the space covered by these spheres can be used as a figure of
   * merit for both the dual and the primal lattice. Suppose \f$ V \f$ contains
   * the basis of the lattice in its lines. Then the density of the lattice is
   * defined by \f$ |\text{det}(V)| \f$ (this determinant is the volume/area contained
   * in the smallest parallelotope that we can fit between the points of the
   * lattice) it also is the number of points that the (not rescalled) lattice
   * will have in \f$[0,1)^t\f$. It is easy to see that the amount of space that
   * will be covered by the spheres is
   * \f[
   *    \frac{\text{Volume of an }n\text{-sphere of radius }d_n}{\text{Density of the lattice}}.
   * \f]
   * There are known upper bounds on the greatest proportion of the space that
   * can be covered in that way \cite mCON99a. If we divide the proportion of
   * the space covered in the current lattice by such an upper bound, the number
   * obtained is what we call **normalized**. This means that the number we will
   * get will be somewhere between 0 and 1 and that it can be used to compare the
   * distribution of the points of the lattice with the points of other lattices.
   * Having a value closer to 1 means a more evenly distributed lattice in some
   * way (see the section on the figures of merit to see which tests use this in
   * the intro). This class is intended as an interface to implement such
   * bounds. It is possible that subclasses of this one implement bounds for a
   * different interpretation of this problem, please look at their documentation
   * before using them. The bounds also sometimes use different norms.
   *
   * This base classe initializes the bounds at 1 and can be used if no
   * normalization can be done, or has to be done. This can be usefull in a few
   * implementations if there is a switch at runtime to instanciate a
   * `Normalizer` subclass because there will be no need to duplicate the code
   * for the case with no normalization.
   *
   * To create a subclass of this one, it is important to implement the
   * getGamma(int) const, getBound(int) const and init() methods and to call
   * the init() method in the constructor to pre-compute the bounds this class
   * will use for normalization.
   *
   * To instanciate this class or its subclasses, the usage of the copy
   * constructor or of the assignment = is prohibited by the fact that these
   * methods are private. The prefered usage for this class is to declare a
   * pointer to a Normalizer and to instanciate subclasses with dynamically 
   * allocated memory:
   * \code{.cpp}
   * Normalizer<RScal>* norma;
   * norma = new NormaBestLat<RScal>(logDensity, t);
   * delete norma;
   * \endcode
   */
  template<typename RedDbl>
    class Normalizer {

      public:
        /**
         * The maximum dimension of the lattices for which this class can give
         * an upper bound.
         * */
        static const int MAX_DIM = 48;

        /**
         * Complete constructor for a `Normalizer`. This will create a
         * normalizer for lattices with a density of \f$\exp(\text{logDensity})\f$
         * in all dimensions \f$\leq t\f$. Ususally, such a normalizer will
         * return bounds that take into a account the density. `Name` allows the
         * user to give a name to a normalizer object. This name will be printed
         * by the ToString() method. It serves no purpose implementation-wise,
         * but this can be usefull while debugging code. `norm` is the NormType
         * that will be used by this object. It cannot be changed. The usage of
         * `beta` is deprecated. This is a bias factor that can be usefull in
         * the case where a figure of merit with numerous projections is
         * computed. It can be used to give more weight to the first dimensions
         * by taking \f$\beta< 1\f$. It inflates the figures of merit by
         * \f$(1/\beta)^t\f$, thus weakening the requirements for good results
         * in large dimensions in a worst-case figure of merit. One normally
         * uses \f$\beta= 1\f$.
         *
         * Note that the log value of the density is stored (instead of the density 
         * itself) so it is easier to manipulate really large values of density.
         *
         * \remark **Richard:** Je crois que ce facteur `beta` devrait
         * disparaître car des poids beaucoup plus généraux sont maintenant
         * implantés dans les classes `*Weights`.
         * **Marc-Antoine:** Je crois que le design de cette classe est à
         * repenser. Il est probablement intéressant de considérer la logDensity
         * en tableaux parce que dans certaines applications la densité change
         * en fonction de la dimension.
         */
        Normalizer (RedDbl & logDensity, int t, std::string Name,
            NormType norm = L2NORM, double beta = 1);

        /**
         * Constructor that does not take the density as an argument. The fields
         * are essentially the same as for Normalizer(RedDbl, int, std::string,
         * NormType, double)
         *
         * This is only used in the case of rank 1 lattices in the NormaPalpha
         * class with a prime density.
         */
        Normalizer (int t, std::string Name, NormType norm = L2NORM,
            double beta = 1);

        /**
         * Destructor.
         */
        virtual ~Normalizer ()
        { delete [] m_bounds; }

        /**
         * This is a method that will initialize the bounds this normalizer can
         * return. This will change the `logDensity` and the `beta` variables
         * that are stored in this object. This will compute bounds for all
         * dimensions smaller than `t` (the parameter passed to the constructors)
         * that can be retrived with the getPreComputedBounds(int) method.
         */
        virtual void init (RedDbl & logDensity, double beta);

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
        virtual RedDbl getBound (int j) const;

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
        RedDbl m_logDensity;

        /**
         * Only elements 1 to <tt>m_maxDim</tt> (inclusive) of m_bounds bellow
         * will be pre-computed. This stores the `t` parameter of the constructors.
         */
        int m_maxDim;

        /**
         * Beta factor used to give more or less importance to some of the
         * dimensions.
         */
        double m_beta;

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
        Normalizer (const Normalizer<RedDbl> &);

        /**
         * Use of assigment is forbidden.
         */
        Normalizer<RedDbl> & operator= (const Normalizer<RedDbl> &);

    }; // End class Normalizer

  //===========================================================================

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
        x = 0.5*log (getGamma(j)) + j*logBeta - y*NTL::conv<double>(logDensity0);
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
       * remark: 
       * in the init method, the bounds are pre-computed for the dimensions of
       * the projection, and are accessible through this function. But in the code
       * a call to function getBound (below) is made. This means the pre-computed
       * bounds are not used and the bounds are calculated again at each step with 
       * the function below. Could be improved.
       */
    }

  /*-------------------------------------------------------------------------*/

  template<typename RedDbl>
    RedDbl Normalizer<RedDbl>::getBound (int j) const
    {
      /*
         assert (j >= 1 && j <= m_maxDim);
         if (j >= 1 && j <= Normalizer::MAX_DIM)
         return getPreComputedBound (j);
         else {
         */

      RedDbl x,y;
      RedDbl logBeta;
      y = NTL::inv(RedDbl(j));
      logBeta = NTL::log(m_beta);
      x = 0.5*NTL::log (getGamma(j)) + j*logBeta - y*m_logDensity;

      return exp(x);
    }

  extern template class Normalizer<double>;
  extern template class Normalizer<NTL::RR>;

} // end namespace LatticeTester

#endif
