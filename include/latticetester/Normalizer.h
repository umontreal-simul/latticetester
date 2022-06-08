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
   * Given a lattice in dimension \f$t\f$, it is possible to center
   * non-interlapping spheres of radius \f$ d_t \f$ where \f$ d_t \f$ is the
   * lenght of the shortest vector in the lattice. It turns out that the
   * proportion of the space covered by these spheres can be used as a figure of
   * merit for both the dual and the primal lattice. Suppose \f$ V \f$ contains
   * the basis of the lattice in its lines. Then the density of the lattice is
   * defined by \f$ \frac{1}{|\text{det}(V)|} \f$. It is the number of points
   * that the (not rescalled) lattice
   * will have in \f$[0,1)^t\f$. It is easy to see that the amount of space that
   * will be covered by the spheres is
   * \f[
   *    \text{Volume of an }n\text{-sphere of radius }d_n\times\text{Density of the lattice}}.
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
   * normalization can be done, or has to be done. This can be useful in a few
   * implementations if there is a switch at runtime to instanciate a
   * `Normalizer` subclass because there will be no need to duplicate the code
   * for the case with no normalization.
   *
   * To create a subclass of this one, it is important to implement the
   * getGamma(int) const, getBound(int) const and init() methods and to call
   * the init() method in the constructor to pre-compute the bounds this class
   * will use for normalization.
   * 
   * The constructors here do not compute any bounds, they basically just reserve 
   * the space (an array) for the bounds.
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
   * 
   * Important: when making a search and examining millions of lattices, it is important 
   * NOT to construct a new Normalizer object and to recompute the constants for each
   * lattice.  
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
         * Complete constructor for a `Normalizer`. This will create a
         * normalizer for lattices with a density of \f$\exp(\text{logDensity})\f$
         * in all dimensions \f$\leq t\f$. Ususally, such a normalizer will
         * return bounds that take into a account the density. `Name` allows the
         * user to give a name to a normalizer object. This name will be printed
         * by the ToString() method. It serves no purpose implementation-wise,
         * but this can be usefull while debugging code. `norm` is the NormType
         * that will be used by this object. It cannot be changed. 
         *
         * The log value of the density is stored (instead of the density 
         * itself) so it is easier to manipulate really large values of density.
         *
         * **Marc-Antoine:** Je crois que le design de cette classe est à
         * repenser. Il est probablement intéressant de considérer la logDensity
         * en tableaux parce que dans certaines applications la densité change
         * en fonction de la dimension.
         */
        Normalizer (RealRed & logDensity, int t, std::string Name,
            NormType norm = L2NORM);

        /**
         * Constructor that does not take the density as an argument. The fields
         * are essentially the same as for Normalizer(RealRed, int, std::string,
         * NormType, double)
         *
         * This is only used in the case of rank 1 lattices in the NormaPalpha
         * class with a prime density.
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
