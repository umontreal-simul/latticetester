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
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"

namespace LatticeTester {

  /**
   * This subclass of `IntLatticeExt` defines a general rank 1 lattice rule in \f$d\f$ dimensions,
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
template<typename Int, typename Real>
class Rank1Lattice: public IntLatticeExt<Int, Real> {

      private:
        typedef NTL::vector<Int>  IntVec;
        typedef NTL::matrix<Int>  IntMat;
        typedef NTL::vector<Real> RealVec;

      public:

        /**
         * This constructor takes as input the modulus `m`, the generating vector `aa`,
         * the (maximal) dimension `maxDim`, and the norm used to measure the vector lengths.
         * The length of the vector `aa` should be `maxDim`.
         * This constructor does not build the basis, to leave
         * more flexibility in the dimension when doing so.
         */
        Rank1Lattice (const Int & m, const IntVec & aa, int maxDim,
            // LatticeTester::NormType norm = LatticeTester::L2NORM);
            NormType norm = L2NORM);

        /**
         * Constructor for the special case of a Korobov lattice.
         * Here the generating vector has the form aa = (1, a, a^2 mod m, a^3 mod m, ...)
         * where a is an integer such that 1 < a < m.
         */
        Rank1Lattice (const Int & m, const Int & a, int maxDim,
            NormType norm = L2NORM);

        /**
         * Copy constructor.
         */
        Rank1Lattice (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Assigns `Lat` to this object.
        */
        Rank1Lattice & operator= (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Destructor.
         */
        ~Rank1Lattice();

        /**
         * Returns the first components of the generating vector \f$\ba\f$ as a string.
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
         * Initializes the rank 1 lattice. This just invokes `IntLatticeExt::init()`.
         */
        void init();

        /**
         * Vector of multipliers (generating vector) of the rank 1 lattice rule.
         * They are stored for up to `maxDim()` dimensions.
         * The first dimension has index 0.
         */
        IntVec m_a;
    };


//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
         const Int & m, const IntVec & aa, int maxDim, NormType norm):
         IntLatticeExt<Int, Real> (m, maxDim, true, norm) {
    this->m_a = aa;
    init();
  }

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
        const Int & m, const Int & a, int maxDim, NormType norm):
        IntLatticeExt<Int, Real> (m, maxDim, true, norm) {
    m_a.SetDim(maxDim);
	Int powa(1);  m_a[0] = powa;
    for (long i=1; i < maxDim; i++) {
    	powa = (a * powa) % m;
    	m_a[i] = powa;
    }
    init();
  }

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::~Rank1Lattice() {
      this->m_a.kill ();
    }

  //============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::init() {
    IntLatticeExt<Int, Real>::init();
      // for (int r = 1; r < this->getDim(); r++)
      //   this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
    }

  //============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real> & Rank1Lattice<Int, Real>::operator= (
         const Rank1Lattice<Int, Real> & lat) {
      if (this == &lat)
        return * this;
      this->copy (lat);
      init ();
      this->m_a = lat.m_a;
      return *this;
    }

  //============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice (
        const Rank1Lattice<Int, Real> & lat):
      IntLatticeExt<Int, Real> (
          lat.m_modulo, lat.getDim (), lat.getNormType ()) {
    // MyExit (1, "Rank1Lattice:: constructor is incomplete" );
    init ();
    this->m_a = lat.m_a;
  }

  //============================================================================

template<typename Int, typename Real>
std::string Rank1Lattice<Int, Real>::toStringCoef ()const {
      return toString (this->m_a, 0, this->getDim ());
    }

  //============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDim () {
      assert(1 + this->getDim() <= this->m_maxDim);
      buildBasis (1 + this->getDim ());
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }

  //============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildBasis (long d) {
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

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::dualize () {
      IntMat tmps(this->m_basis);
      this->m_basis = this->m_dualbasis;
      this->m_dualbasis = tmps;
    }

//============================================================================

template class Rank1Lattice<std::int64_t, double>;
template class Rank1Lattice<NTL::ZZ, double>;
template class Rank1Lattice<NTL::ZZ, NTL::RR>;

} // End namespace LatticeTester

#endif
