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


#ifndef LATTICETESTER__INTLATTICEBASIS_H
#define LATTICETESTER__INTLATTICEBASIS_H
#include "latticetester/Const.h"
#include "latticetester/Util.h"
#include "latticetester/ntlwrap.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

namespace LatticeTester {

  /**
   * \class IntLatticeBasis
   *
   * \brief This class represents a Lattice and its basis
   *
   * This class offers tools to manipulate lattice bases. Each lattice is
   * at least represented by a basis \f$V\f$, a dimension and a Norm.
   * Users can beside precise the dual lattice \f$W\f$ and the modulo.
   * In that case, the flag \f$m_withDual\f$ is set to \f$true\f$.
   */

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    class IntLatticeBasis {
      public:

        /**
         * Constructor. The primal basis is initialize with identity,
         * the dimension of the lattice with dim and the Norm used for
         * reduction with norm.
         */
        IntLatticeBasis (const int dim, NormType norm = L2NORM);

        /**
         * Constructor. The primal basis is initialize with \f$basis\f$,
         * the dimension of the lattice with dim and the Norm used for
         * reduction with norm.
         */
        IntLatticeBasis (const BasIntMat basis, const int dim,
            NormType norm = L2NORM);

        /**
         * Constructor. The primal basis is initialize with \f$primalbasis\f$,
         * the dual basis is initualize with \f$dualbasis\f$, the
         * dimension of the lattice with \f$dim\f$ and the Norm used for
         * reduction with norm.
         */
        IntLatticeBasis (const BasIntMat primalbasis, const BasIntMat dualbasis,
            const Int modulo, const int dim, NormType norm = L2NORM);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        IntLatticeBasis (const IntLatticeBasis<Int, BasInt, BasIntVec,
            BasIntMat, Dbl, DblVec, RedDbl> & Lat);

        /**
         * Destructor
         */
        ~IntLatticeBasis ();

        /**
         * Cleans and releases all the memory allocated for this lattice.
         */
        void kill ();

        /**
         * Copy the lattice
         */
        void copyBasis (const IntLatticeBasis<Int, BasInt, BasIntVec,
            BasIntMat, Dbl, DblVec, RedDbl> & lat);

        /**
         * Copy the n first elements of the lattice lat
         */
        void copyBasis (const IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat,
            Dbl, DblVec, RedDbl> & lat, int n);

        /**
         * Set all the norm of the matrix to -1
         */
        void initVecNorm ();

        /**
         * Return the basis, a BMat-type matrix
         */
        BasIntMat & getBasis () { return m_basis; }

        /**
         * Return the dual basis, a BMat-type Matrix
         */
        BasIntMat & getDualBasis () { return m_dualbasis; }

        /**
         * Return the dimension of the lattice, a int type
         */
        int getDim () const { return m_dim; }

        /**
         * Return the Norm type used by the lattice
         */
        NormType getNorm () const { return m_norm; }

        /**
         * Return the norm of the i-th line of the basis
         */
        Dbl getVecNorm (const int & i) { return m_vecNorm[i]; }

        /**
         * Return the norm of each line of the basis
         */
        DblVec getVecNorm () const { return m_vecNorm; }

        /**
         * Return the norm of the i-th line of the basis
         */
        Dbl getDualVecNorm (const int & i) { return m_dualvecNorm[i]; }

        /**
         * Return the norm of each line of the basis
         */
        DblVec getDualVecNorm () const { return m_dualvecNorm; }

        /**
         * Return the modulo used for the Dual
         */
        Int getModulo () const { return m_modulo; }

        /**
         * Set the dimension of the basis
         */
        void setDim (const int & d) { if(d>0) m_dim = d;}

        /**
         * Set the norm used by the lattice
         */
        void setNorm (const NormType & norm) { m_norm = norm; }

        /**
         * Set the norm i egal to the value of NScal
         */
        void setVecNorm ( const Dbl & value, const int & i)
        {
          m_vecNorm[i] = value;
        }

        /**
         * Set the dual norm i egal to the value of NScal
         */
        void setDualVecNorm ( const Dbl & value, const int & i)
        {
          m_dualvecNorm[i] = value;
        }

        /**
         * Return True if we use Dual.
         */
        bool withDual() { return m_withDual; }

        /**
         * Set the Dual Flag.
         */
        void setDualFlag(bool flag) { m_withDual = flag; }

        /**
         * Get and det m_xx but I don't now what means m_xx
         */
        bool getXX (int i) const { return m_xx[i]; } //???
        void setXX (bool val, int i) { m_xx[i] = val; } //??

        /**
         * Set the norm of all vectors to -1
         */
        void setNegativeNorm ();

        /**
         * Set the norm of the i-eme vector to -1
         */
        void setNegativeNorm (const int & i){ m_vecNorm[i] = -1; }

        /**
         * Set the norm of all vectors of the dual basis to -1
         */
        void setDualNegativeNorm ();

        /**
         * Set the norm of the i-eme vector in the dual Basis to -1
         */
        void setDualNegativeNorm (const int & i){ m_dualvecNorm[i] = -1; }

        /**
         * Recalculates the norm of each vector in the basis of
         * the lattice
         */
        void updateVecNorm ();

        /**
         * Recalculates the norm of each vector in the basis of
         * the lattice from d to dim
         */
        void updateVecNorm (const int & d);

        /**
         * Recalculates the norm of each vector in the basis of
         * the lattice
         */
        void updateDualVecNorm ();

        /**
         * Recalculates the norm of each vector in the dual basis of
         * the lattice from d to dim
         */
        void updateDualVecNorm (const int & d);

        /**
         * Updates the norm of vector at dimension `d` using the `L2NORM`.
         */
        void updateScalL2Norm (const int i);

        /**
         * Updates the norm of all basis vectors from dimensions `d1` to `d2`
         * (exclusive) using the `L2NORM`.
         */
        void updateScalL2Norm (const int k1, const int k2);

        /**
         * Updates the norm of vector in the Dual Basis at dimension `d`
         * using the `L2NORM`.
         */
        void updateDualScalL2Norm (const int i);

        /**
         * Updates the norm of all dual basis vectors from dimensions `d1`
         * to `d2` (exclusive) using the `L2NORM`.
         */
        void updateDualScalL2Norm (const int k1, const int k2);

        /**
         * Exchanges vectors \f$i\f$ and \f$j\f$ in the basis.
         */
        void permute (int i, int j);

        /**
         * Check duality
         */
        bool checkDuality();

        /**
         * Sorts the basis vectors with indices from \f$d\f$ to the dimension
         * of the basis by increasing length. The dual vectors are permuted
         * accordingly. Assumes that the lengths of the corresponding basis
         * vectors are up to date.
         */
        void sort (int d);

        /**
         * Return a string with the primal basis and its norms
         */
        std::string toStringBasis() const;

        /**
         * Return a string with the dual basis and its norms
         */
        std::string toStringDualBasis() const;

        /**
         * Writes the lattice and the parameters on standard output
         */
        void write () const;

      protected:

        /**
         * Represent the basis of the lattice. BMat is defined in Types.h
         * This is a matrix where the type of the number can change.
         */
        BasIntMat m_basis;

        /**
         * Represent the dual basis of the lattice. BMat is defined in Types.h
         * This is a matrix where the type of the number can change.
         */
        BasIntMat m_dualbasis;


        /**
         * The dimension of the space.
         */
        int m_dim;

        /**
         * The norm used in the reduction and for this lattice
         */
        NormType m_norm;

        /**
         * The norm of each vector in the basis.
         */
        DblVec m_vecNorm;

        /**
         * The norm of each vector in the dual basis.
         */
        DblVec m_dualvecNorm;

        /**
         * The modulo linked to the m-dual.
         */
        Int m_modulo;

        /**
         * If m_withDual is true, we can use Dual Basis.
         */
        bool m_withDual;

        /**
         * This table is used in the Minkowski reduction.
         */
        bool *m_xx;


    }; // class IntLatticeBasis

  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::IntLatticeBasis (const int dim, NormType norm):
      m_dim (dim),
      m_norm (norm),
      m_modulo(0),
      m_withDual(false),
      m_xx(0)

  {
    this->m_basis.resize(dim,dim);
    this->m_vecNorm.resize (dim);
    initVecNorm();
  }

  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::IntLatticeBasis (const BasIntMat basis, const int dim, 
        NormType norm):
      m_basis (basis),
      m_dim (dim),
      m_norm (norm),
      m_modulo(0),
      m_withDual(false),
      m_xx(0)
  {
    this->m_vecNorm.resize (dim);
    initVecNorm();
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::IntLatticeBasis (
        const BasIntMat primalbasis,
        const BasIntMat dualbasis,
        const Int modulo,
        const int dim,
        NormType norm):
      IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl>(
          primalbasis, dim, norm)
  {
    this->m_dualbasis = dualbasis;
    this->m_withDual = true;
    this->m_dualvecNorm.resize (dim);
    setDualNegativeNorm();
    this->m_modulo = modulo;
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::IntLatticeBasis (const IntLatticeBasis<Int, BasInt, BasIntVec,
        BasIntMat, Dbl, DblVec, RedDbl> & lat):
      m_dim (lat.getDim ()),
      m_norm (lat.getNorm ()),
      m_xx(0)
  {
    copyBasis (lat);
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::~IntLatticeBasis ()
    {
      kill();
      this->m_basis.BasIntMat::clear ();
      this->m_dualbasis.BasIntMat::clear ();
      this->m_vecNorm.clear ();
      this->m_dualvecNorm.clear ();
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::kill ()
    {
      delete [] this->m_xx;
      this->m_xx = 0;
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::copyBasis (const IntLatticeBasis<Int, BasInt, BasIntVec,
        BasIntMat, Dbl, DblVec, RedDbl> & lat)
    {
      //if(m_dim == lat.m_dim)

      this->m_basis = lat.m_basis;
      this->m_dualbasis = lat.m_dualbasis;
      this->m_vecNorm = lat.m_vecNorm;
      this->m_dualvecNorm = lat.m_dualvecNorm;
      this->m_withDual = lat.m_withDual;
      this->m_modulo = lat.m_modulo;
      this->m_xx = new bool[m_dim];
      for (int i = 0; i < this->m_dim; i++)
        this->m_xx[i] = lat.getXX(i);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::copyBasis (const IntLatticeBasis<Int, BasInt, BasIntVec, 
        BasIntMat, Dbl, DblVec, RedDbl> & lat, int n)
    {
      if(this->m_dim == n) {
        CopyMatr(this->m_basis, lat.m_basis, n);
        CopyVect(this->m_vecNorm, lat.m_vecNorm, n);
        this->m_withDual = lat.m_withDual;
        if(this->m_withDual){
          this->m_dualbasis.resize(this->m_basis.size1(),
              this->m_basis.size1());
          this->m_dualvecNorm.resize(this->m_basis.size1());
          CopyMatr(this->m_dualbasis, lat.m_dualbasis, n);
          CopyVect(this->m_dualvecNorm, lat.m_dualvecNorm, n);
        }
        this->m_modulo = lat.m_modulo;
        this->m_xx = new bool[n];
        for (int i = 0; i < n; i++)
          this->m_xx[i] = lat.getXX(i);
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::initVecNorm ()
    {
      this->m_xx = new bool[this->m_dim];
      for(int i = 0; i < this->m_dim; i++){
        this->m_vecNorm[i] = -1;
        this->m_xx[i] = true;
      }
    }


  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::setNegativeNorm ()
    {
      for (int i = 0; i<this->m_dim; i++){
        this->m_vecNorm[i] = -1;
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::setDualNegativeNorm ()
    {
      for(int i = 0; i < this->m_dim; i++){
        this->m_dualvecNorm[i] = -1;
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::updateVecNorm ()
    {
      updateVecNorm (0);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::updateVecNorm (const int & d)
    {
      assert (d >= 0);

      for (int i = d; i < this->m_dim; i++) {
        if (this->m_vecNorm[i] < 0) {
          NTL::matrix_row<BasIntMat> row(this->m_basis, i);
          if (this->m_norm == L2NORM) {
            ProdScal<Int> (row, row, this->m_dim, this->m_vecNorm[i]);
          } else {
            CalcNorm <BasIntVec, Dbl> (row, this->m_dim, this->m_vecNorm[i],
                this->m_norm);
          }
        }
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::updateDualVecNorm ()
    {
      updateDualVecNorm (0);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::updateDualVecNorm (const int & d)
    {
      assert (d >= 0);

      for (int i = d; i < this->m_dim; i++) {
        if (this->m_dualvecNorm[i] < 0) {
          NTL::matrix_row<BasIntMat> row(this->m_dualbasis, i);
          if (this->m_norm == L2NORM) {
            ProdScal<Int> (row, row, this->m_dim, this->m_dualvecNorm[i]);
          } else {
            CalcNorm <BasIntVec, Dbl> (row, this->m_dim, this->m_dualvecNorm[i],
                this->m_norm);
          }
        }
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::updateScalL2Norm (const int i)
    {
      NTL::matrix_row<BasIntMat> row(this->m_basis, i);
      ProdScal<Int> (row, row, this->m_dim, this->m_vecNorm[i]);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::updateScalL2Norm (const int k1, const int k2)
    {
      for (int i = k1; i < k2; i++) {
        updateScalL2Norm(i);
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::updateDualScalL2Norm (const int i)
    {
      NTL::matrix_row<BasIntMat> row(this->m_dualbasis, i);
      ProdScal<Int> (row, row, this->m_dim, this->m_dualvecNorm[i]);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::updateDualScalL2Norm (const int k1,
        const int k2)
    {
      for (int i = k1; i < k2; i++) {
        updateDualScalL2Norm(i);
      }
    }


  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::permute (int i, int j)
    {
      if (i == j)
        return ;
      for (int k = 0; k < this->m_dim; k++){
        swap9 (this->m_basis(j,k), this->m_basis(i,k));
        if(this->m_withDual){
          swap9(this->m_dualbasis(j,k), this->m_dualbasis(i,k));
        }
      }
      swap9 (this->m_vecNorm[i], this->m_vecNorm[j]);
      if(this->m_withDual){
        swap9 (this->m_dualvecNorm[i], this->m_dualvecNorm[j]);
      }
      bool b = this->m_xx[j];
      this->m_xx[j] = this->m_xx[i];
      this->m_xx[i] = b;
    }


  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    bool IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, 
    RedDbl>::checkDuality ()
    {
      if(!this->m_withDual) {
        std::cout << "DO NOT USE IntLatticeBasis::checkDuality without dual"
          << std::endl;
        return false;
      }
      BasInt S;
      int dim = getDim ();

      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          NTL::matrix_row<const BasIntMat> row1(this->m_basis, i);
          NTL::matrix_row<const BasIntMat> row2(this->m_dualbasis, j);
          ProdScal<Int> (row1, row2, dim, S);
          if (j != i) {
            if (S != 0) {
              std::cout << "******  checkDuality failed for V[" << i <<
                "] and W[" << j << "]" << std::endl;
              return false;
            }
          } else if (S != this->m_modulo) {
            std::cout << "******  checkDuality failed for i, j = " << i << " , " <<
              j << std::endl;
            return false;
          }
        }
      }
      return true;


    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::sort (int d)
    /*
     * We assume that the square lengths are already updated.
     * This gives flexibility to the user to put something else than
     * the square Euclidean length in V.vecNorm, W.vecNorm, etc.
     */
    {
      int dim = getDim ();
      for (int i = 0; i < dim; i++){
        if (getVecNorm(i) < 0) {
          std::cout << "\n***** ERROR: sort   Negative norm for i = " << i <<
            ",  dim = " << dim << std::endl;
        }
      }

      for (int i = d; i < dim; i++){
        int k = i;
        for (int j = i + 1; j < dim; j++) {
          if (getVecNorm (j) < getVecNorm (k))
            k = j;
        }
        if (i != k)
          permute (i, k);
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::write () const
    {
      std::cout << "Dim = " << this->m_dim << " \n \n";
      for (int i = 0; i < this->m_dim; i++) {
        std::cout << "  |     ";
        for (int j = 0; j < this->m_dim; j++) {
          std::cout << std::setprecision (15) << this->m_basis(i,j) << "\t";
        }
        std::cout << "    |";
        if(this->m_withDual){
          std::cout <<"  |     ";
          for (int j = 0; j < this->m_dim; j++) {
            std::cout << std::setprecision (15) <<
              this->m_dualbasis(i,j) << "\t";
          }
          std::cout << "    |";
        }
        std::cout << "\n";

      }
      std::cout << "\n";
      std::cout << "Norm used : " << toStringNorm(this->m_norm) << "\n"
        << std::endl;
      std::cout << "Norm of each Basis vector : \n";
      std::cout << " Primal     ";
      if(this->m_withDual)
        std::cout << "\t Dual \n";
      std::cout << "\n";

      for (int i = 0; i < this->m_dim; i++) {
        std::cout << "   ";
        if (this->m_vecNorm[i] < 0) {
          std::cout << "NaN OR Not computed";
        } else {
          std::cout << this->m_vecNorm[i];
        }
        if(this->m_withDual){
          std::cout << "\t \t \t ";
          if (this->m_dualvecNorm[i] < 0)
            std::cout << "NaN OR Not computed";
          else
            std::cout << this->m_dualvecNorm[i];
        }
        std::cout << "\n";
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    std::string IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
           RedDbl>::toStringBasis () const
    {
      std::ostringstream os;
      os << "Primal Basis:\n";
      os << "  Dim = " << this->m_dim << " \n";
      for (int i = 0; i < this->m_dim; i++) {
        os << "    [";
        for (int j = 0; j < this->m_dim; j++)
          os << " " <<  std::setprecision (15) << this->m_basis(i,j);
        os << " ]\n";
      }

      os << "  Norm:\n";
      for (int i = 0; i < this->m_dim; i++) {
        os << "    ";
        if (this->m_vecNorm[i] < 0) {
          os << "-1" << std::endl;
        } else {
          os << this->m_vecNorm[i] << std::endl;
        }
      }
      os << std::endl;
      return os.str ();
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename BasIntVec, 
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
    std::string IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec,
    RedDbl>::toStringDualBasis () const
    {
      std::ostringstream os;
      os << "Dual Basis:\n";
      os << "  Dim = " << this->m_dim << " \n";
      for (int i = 0; i < this->m_dim; i++) {
        os << "    [";
        for (int j = 0; j < this->m_dim; j++)
          os << " " <<  std::setprecision (15) << this->m_dualbasis(i,j);
        os << " ]\n";
      }

      os << "  Norm:\n";
      for (int i = 0; i < this->m_dim; i++) {
        os << "    ";
        if (this->m_dualvecNorm[i] < 0) {
          os << "-1" << std::endl;
        } else {
          os << this->m_dualvecNorm[i] << std::endl;
        }
      }
      os << std::endl;
      return os.str ();
    }

} // namespace LatticeTester

#endif


