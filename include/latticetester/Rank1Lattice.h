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

#ifndef LATTICETESTER_RANK1LATTICE_H
#define LATTICETESTER_RANK1LATTICE_H

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLattice.h"

namespace LatticeTester {

  /**
   * This subclass of `IntLattice` defines a general rank 1 lattice rule in \f$d\f$ dimensions,
   * whose points \f$\mathbb{u}_i\f$ are defined by
   * \f{equation}{
   *    \mathbf{u}_i = (i \mathbf{a} \mod m)/m,
   * \f}
   * where \f$\mathbf{a} \in \mathbb{Z}_m^d\f$ is called the generating vector.
   * The lattice is rescaled simply by removing the division by \f$m\f$.
   * Given \f$\mathbf{a}\f$, a basis of the rescaled (integer) lattice is given by
   * \f{align*}{
   *    \mathbf{v}_1 & = a \\
   *    \mathbf{v}_2 & = m \mathbf{e}_2 \\
   *    \vdots & \\
   *    \mathbf{v}_d & = m \mathbf{e}_v \\
   * \f}
   * where \f$\mathbf{e}_i\f$ is the \f$i^\text{th}\f$ unit vector.
   * 
   * A condition that is often required when building a rank 1 lattice is that
   * \f$\gcd(a_i, m) = 1,\ \forall 1\leq i \leq d\f$. When this condition is
   * verified, each lower-dimensional projection of the lattice contains the same number
   * of points (has the same density \f$m\f$) as the full lattice.
   * When searching for lattices that satisfy this condition, one may assume
   * without loss of generality generality that \f$a_1 = 1\f$.
   */
  template<typename Int, typename Real, typename RealRed>
    class Rank1Lattice: public IntLattice<Int, Real, RealRed> {

      private:
        typedef NTL::vector<Int>  IntVec;
        typedef NTL::matrix<Int>  IntMat;
        typedef NTL::vector<Real> RealVec;

      public:

        /**
         * This constructor takes as input the modulus `m`, the generating vector `a`,
         * the (maximal) dimension `maxDim`, and the norm used to measure the vector lengths.
         * The dimension of the vector `a` should be `maxDim`.
         * This constructor does not build the basis, to leave
         * more flexibility in the dimension when doing so.
         */
        Rank1Lattice (const Int & m, const IntVec & a, int maxDim,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Copy constructor.
         */
        Rank1Lattice (const Rank1Lattice<Int, Real, RealRed> & Lat);

        /**
         * Assigns `Lat` to this object.
        */
        Rank1Lattice & operator= (const Rank1Lattice<Int, Real, RealRed>
            & Lat);

        /**
         * Destructor.
         */
        ~Rank1Lattice();

        /**
         * Returns the first components of the generating vector \f$a\f$ as a string.
         * The number of components in the string will be the current dimension of the lattice.
         */
        std::string toStringCoef() const;

        /**
         * Builds a basis in `d` dimensions. This `d` must not exceed `this->maxDim()`.
         */
        void buildBasis (long d);

        /**
         * Dualizes the lattice by exchanging the primal and dual bases.
         */
        void dualize ();

        /**
         * Increases the current dimension by 1.
         */
        void incDim ();

      protected:

        /**
         * Initializes the rank 1 lattice. This just invokes `IntLattice::init()`.
         */
        void init();

        /**
         * The multipliers (generating vector) of the rank 1 lattice rule.
         * They are stored for up to `maxDim()` dimensions.
         * The first dimension has index 0.
         */
        IntVec m_a;
    };

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    Rank1Lattice<Int, Real, RealRed>::Rank1Lattice (
        const Int & m, const IntVec & a, int maxDim, NormType norm):
      IntLattice<Int, Real, RealRed> (m, maxDim, true, norm)
  {
    this->m_a = a;
    init();
  }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    Rank1Lattice<Int, Real, RealRed>::~Rank1Lattice()
    {
      this->m_a.clear ();
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    void Rank1Lattice<Int, Real, RealRed>::init()
    {
      IntLattice<Int, Real, RealRed>::init();
      // for (int r = 1; r < this->getDim(); r++)
      //   this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    Rank1Lattice<Int, Real, RealRed> &
    Rank1Lattice<Int, Real, RealRed>::operator= (
        const Rank1Lattice<Int, Real, RealRed> & lat)
    {
      if (this == &lat)
        return * this;
      this->copy (lat);
      init ();
      this->m_a = lat.m_a;
      return *this;
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    Rank1Lattice<Int, Real, RealRed>::Rank1Lattice (
        const Rank1Lattice<Int, Real, RealRed> & lat):
      IntLattice<Int, Real, RealRed> (
          lat.m_modulo, lat.getDim (), lat.getNormType ())
  {
    // MyExit (1, "Rank1Lattice:: constructor is incomplete" );
    init ();
    this->m_a = lat.m_a;
  }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    std::string Rank1Lattice<Int, Real, RealRed>::toStringCoef ()const
    {
      return toString (this->m_a, 0, this->getDim ());
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    void Rank1Lattice<Int, Real, RealRed>::incDim ()
    {
      assert(1 + this->getDim() <= this->m_maxDim);
      buildBasis (1 + this->getDim ());
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    void Rank1Lattice<Int, Real,RealRed>::buildBasis (long d)
    {
      assert(d <= this->m_maxDim);
      this->setDim (d);
      this->m_basis.SetDims(d,d);
      this->m_dualbasis.SetDims(d,d);

      // conv(m_v[1][1], 1);

      for (int j = 0; j < d; j++) {
        this->m_basis (0, j) = this->m_a[j];
      }

      for (int i = 1; i < d; i++) {
        for (int j = 0; j < d; j++) {
          if (i == j) {
            this->m_basis (i, j) = this->m_modulo;
          } else {
            this->m_basis (i, j) = 0;
          }
        }
      }

      // if a[0] != 1, the basis must be triangularized
      //BasisConstruction<Int> constr;
      if (this->m_basis (0, 0) != 1) {
        //constr.GCDConstruction(this->m_basis);
         Triangularization (
             this->m_basis, this->m_dualbasis, d, d, this->m_modulo);
         dualize ();
      }
      //constr.mDualTriangular(this->m_basis, this->m_dualbasis, this->m_modulo);
      calcDual (this->m_basis, this->m_dualbasis, d, this->m_modulo);
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    void Rank1Lattice<Int, Real, RealRed>::dualize ()
    {
      IntMat tmps(this->m_basis);
      this->m_basis = this->m_dualbasis;
      this->m_dualbasis = tmps;
    }

  extern template class Rank1Lattice<std::int64_t, double, double>;
  extern template class Rank1Lattice<NTL::ZZ, double, double>;
  extern template class Rank1Lattice<NTL::ZZ, NTL::RR, NTL::RR>;

} // End namespace LatticeTester

#endif
