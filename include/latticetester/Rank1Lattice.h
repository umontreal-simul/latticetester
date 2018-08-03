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

#ifndef RANK1LATTICE_H
#define RANK1LATTICE_H
#include "latticetester/Const.h"
#include "latticetester/IntLattice.h"


namespace LatticeTester {

  /**
   * This class implements a general rank 1 lattice basis. For the values
   * \f$a_1, a_2, …, a_d\f$ given, the \f$d\f$-dimensional lattice basis is
   * formed as:
   * \f[
   * \mathbf{b_1} = (a_1, a_2, …, a_d),\quad\mathbf{b_2} = (0, n, 0, …, 0),\quad…, \quad\mathbf{b_d} = (0, …, 0, n)
   * \f]
   * Without loss of generality, one may choose \f$a_1 = 1\f$.
   *
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    class Rank1Lattice: public IntLattice<Int, BasInt, Dbl, RedDbl> {
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::vector<BasInt> BasIntVec;
        typedef NTL::matrix<BasInt> BasIntMat;
        typedef NTL::vector<Dbl> DblVec;
      public:

        /**
         * Constructor. \f$d\f$ represents the number of multipliers in the array
         * `a`.
         */
        Rank1Lattice (const Int & n, const IntVec & a, int d,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Copy constructor.
         */
        Rank1Lattice (const Rank1Lattice<Int, BasInt, Dbl, RedDbl> & Lat);

        /**
         * Assigns `Lat` to this object.
         */
        Rank1Lattice & operator= (const Rank1Lattice<Int, BasInt, Dbl, RedDbl>
            & Lat);

        /**
         * Destructor.
         */
        ~Rank1Lattice();

        /**
         * Returns the vector of multipliers \f$a\f$ as a string.
         */
        std::string toStringCoef() const;

        /**
         * Builds the basis in dimension \f$d\f$.
         */
        void buildBasis (int d);

        /**
         * Dualize the matrix. The matrix entered need to have
         * the particular shape describe ERWAN
         */
        void dualize ();

        /**
         * Increases the dimension by 1.
         */
        void incDim ();
      protected:

        /**
         * Initializes the rank 1 lattice.
         */
        void init();

        /**
         * The multipliers of the rank 1 lattice rule.
         */
        IntVec m_a;
    };

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    Rank1Lattice<Int, BasInt, Dbl, RedDbl>::Rank1Lattice (
        const Int & n, const IntVec & a, int maxDim, NormType norm):
      IntLattice<Int, BasInt, Dbl, RedDbl> (n, 1, maxDim, norm)
  {
    this->m_a = a;
    init();
  }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    Rank1Lattice<Int, BasInt, Dbl, RedDbl>::~Rank1Lattice()
    {
      this->m_a.clear ();
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void Rank1Lattice<Int, BasInt, Dbl, RedDbl>::init()
    {
      IntLattice<Int, BasInt, Dbl, RedDbl>::init();
      for (int r = 1; r < this->getDim(); r++)
        this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    Rank1Lattice<Int, BasInt, Dbl, RedDbl> &
    Rank1Lattice<Int, BasInt, Dbl, RedDbl>::operator= (
        const Rank1Lattice<Int, BasInt, Dbl, RedDbl> & lat)
    {
      if (this == &lat)
        return * this;
      copy (lat);
      init ();
      this->m_a = lat.m_a;
      return *this;
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    Rank1Lattice<Int, BasInt, Dbl, RedDbl>::Rank1Lattice (
        const Rank1Lattice<Int, BasInt, Dbl, RedDbl> & lat): 
      IntLattice<Int, BasInt, Dbl, RedDbl> (
          lat.m_modulo, lat.getOrder (), lat.getDim (), lat.getNorm ())
  {
    // MyExit (1, "Rank1Lattice:: constructeur n'est pas terminé " );
    init ();
    this->m_a = lat.m_a;
  }


  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    std::string Rank1Lattice<Int, BasInt, Dbl, RedDbl>::toStringCoef ()const
    {
      return toString (this->m_a, 0, this->getDim ());
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void Rank1Lattice<Int, BasInt, Dbl, RedDbl>::incDim ()
    {
      // kill();
      buildBasis (1 + this->getDim ());
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void Rank1Lattice<Int, BasInt, Dbl,RedDbl>::buildBasis (int d)
    {
      // assert(d <= getMaxDim());
      this->setDim (d);

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
      if (this->m_basis (0, 0) != 1) {
        Triangularization < BasIntMat > (
            this->m_basis, this->m_dualbasis, d, d, this->m_modulo);
        dualize ();
      }
      CalcDual < BasIntMat > (
          this->m_basis, this->m_dualbasis, d, this->m_modulo);
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void Rank1Lattice<Int, BasInt, Dbl, RedDbl>::dualize ()
    {
      BasIntMat tmps(this->m_basis);
      this->m_basis = this->m_dualbasis;
      this->m_dualbasis = tmps;
    }

} // End namespace LatticeTester
#endif
