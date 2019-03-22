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


#ifndef LATTICETESTER_INTLATTICEBASIS_H
#define LATTICETESTER_INTLATTICEBASIS_H

#include "latticetester/Const.h"
#include "latticetester/Util.h"
#include "latticetester/ntlwrap.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

namespace LatticeTester {

  /**
   * This class represents a lattice and its basis and offers tools to do basic 
   * manipulations on lattice bases. Lattices are always stored rescaled. That
   * is, we only consider lattices with rational coordinates such that `m` is a
   * lcm for all the denominators in basis coordinates. It is then possible to
   * represent the lattice in the integers instead of the real numbers by
   * multiplying all it's vectors by `m`. We call this lattice a rescaled
   * lattice. In practice, this allows an exact representation of the arithmetic
   * on the basis. This is usefull in all the use cases of this software.
   *
   * There are numerous ways to represent a lattice. Depending on the
   * calculations that need to be done, it is possible to only provide the basis
   * vectors, the norm and the dimension. There are also fields to provide a
   * dual basis and an `m` (the variables named `modulo`). Both these fields are
   * needed for most applications using the dual basis. The dimension \f$d\f$
   * specifies that the basis contains \f$d\f$ vectors in \f$\mathbb{Z}^d\f$.
   * The norm is one of the norm of NormType.
   *
   * This class also has methods and attributes that can be used to store and
   * compute the norms of the basis and dual basis vectors. It can also permute
   * vectors in the basis or sort them by length and do the corresponding 
   * changes in the dual. Finally, it can also check if the dual really is a
   * dual to the primal basis.
   *
   * This class is made to be built upon. It only offers the bare basics of what
   * a user could want to do with a Lattice. It is most likely that the
   * utilities presented here do not suffice common needs. For a more extensive
   * representation of a lattice look for the `IntLattice` and the `Rank1Lattice`
   * classes.
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    class IntLatticeBasis {

      private:
        // Forward definition of types to be used in this class.
        typedef NTL::vector<BasInt> BasIntVec;
        typedef NTL::matrix<BasInt> BasIntMat;
        typedef NTL::vector<Dbl> DblVec;

      public:

        /**
         * Constructor initializing the primal basis with the identity matrix.
         * The dimension of the lattice is set to `dim` and the norm used for
         * reduction to `norm`.
         */
        IntLatticeBasis (const int dim, NormType norm = L2NORM);

        /**
         * Constructor taking all three needed component of an IntLatticeBasis.
         * The primal basis is initialized with `basis`, the dimension of the
         * lattice with `dim` and the norm used for reduction with `norm`.
         */
        IntLatticeBasis (const BasIntMat basis, const int dim,
            NormType norm = L2NORM);

        /**
         * Complete constructor. The primal basis is initialized with `primalbasis`,
         * the dual basis with `dualbasis`, the \f$m\f$ used for rescaling with 
         * `modulo`, the dimension of the lattice with `dim` and the norm used 
         * for reduction with `norm`.
         */
        IntLatticeBasis (const BasIntMat primalbasis, const BasIntMat dualbasis,
            const Int modulo, const int dim, NormType norm = L2NORM);

        /**
         * Copy constructor. This will copy the entirety of Lat into `*this`.
         */
        IntLatticeBasis (const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & Lat);

        /**
         * Destructor.
         */
        ~IntLatticeBasis ();

        /**
         * Cleans and releases all the memory allocated to this lattice.
         */
        void kill ();

        /**
         * Copy the lattice `lat`, except it's NormType and dimension, into this
         * object. This does not check if the dimensions match.
         */
        void copyBasis (const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat);

        /**
         * Copy the n first elements of the basis of the lattice `lat` into this
         * object. The object into which `lat` is copied has to be of dimension n.
         */
        void copyBasis (const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat,
            long n);

        /**
         * Initializes a vector containing the norms of the basis vectors to -1
         * at all components.
         */
        void initVecNorm ();

        /**
         * Returns the basis represented in a matrix.
         */
        BasIntMat & getBasis () { return m_basis; }

        /**
         * Returns the dual basis represented in a matrix.
         */
        BasIntMat & getDualBasis () { return m_dualbasis; }

        /**
         * Returns the dimension of the lattice. (Both the number of vectors in
         * the basis and the number of coordinates of those vectors.
         */
        int getDim () const { return m_dim; }

        /**
         * Returns the NormType used by this lattice.
         */
        NormType getNorm () const { return m_norm; }

        /**
         * Returns the norm of the i-th vector of the basis.
         */
        Dbl getVecNorm (const int & i) { return m_vecNorm[i]; }

        /**
         * Returns the norm of each vector of the basis in a vector.
         */
        DblVec getVecNorm () const { return m_vecNorm; }

        /**
         * Returns the norm of the i-th vector of the dual basis.
         */
        Dbl getDualVecNorm (const int & i) { return m_dualvecNorm[i]; }

        /**
         * Returns the norm of each vector of the dual basis in a vector.
         */
        DblVec getDualVecNorm () const { return m_dualvecNorm; }

        /**
         * Returns the `m` used for rescaling if it has been defined. Returns `0`
         * otherwise.
         */
        Int getModulo () const { return m_modulo; }

        /**
         * Sets the dimension of the basis to `d`. This won't change any of the
         * vectors of the basis by itself. This method should not be called
         * directly on an object of this class except in a function specifically
         * changing the dimension of this object.
         */
        void setDim (const int & d) { if(d>0) m_dim = d;}

        /**
         * Sets the NormType used by this lattice to `norm`.
         */
        void setNorm (const NormType & norm) { m_norm = norm; }

        /**
         * Sets the norm of the `i`-th component of the basis to `value`. The
         * usage of `updateVecNorm(const int&)` is recommended over this function.
         */
        void setVecNorm ( const Dbl & value, const int & i)
        {
          m_vecNorm[i] = value;
        }

        /**
         * Sets the norm of the `i`-th component of the dual basis to `value`.
         * The usage of `updateDualVecNorm(const int&)` is recommended over this function.
         */
        void setDualVecNorm ( const Dbl & value, const int & i)
        {
          m_dualvecNorm[i] = value;
        }

        /**
         * Returns `true` if a dual has been defined and `false` otherwise.
         */
        bool withDual() { return m_withDual; }

        /**
         * Sets the `withDual` flag to `flag`. This flag indicates whether or 
         * not this IntLatticeBasis contains a dual basis. It is the flag
         * returned by `withDual()`.
         */
        void setDualFlag(bool flag) { m_withDual = flag; }

        /**
         * Sets all the values in the array containing the norms of the basis 
         * vectors to -1.
         */
        void setNegativeNorm ();

        /**
         * Sets the value of the `i`-th component in the array containing the 
         * norms of the basis vectors to -1.
         */
        void setNegativeNorm (const int & i){ m_vecNorm[i] = -1; }

        /**
         * Sets all the values in the array containing the norms of the dual basis 
         * vectors to -1.
         */
        void setDualNegativeNorm ();

        /**
         * Sets the value of the `i`-th component in the array containing the 
         * norms of the dual basis vectors to -1.
         */
        void setDualNegativeNorm (const int & i){ m_dualvecNorm[i] = -1; }

        /**
         * Updates the array containing the basis vectors norms by recomputing 
         * them.
         * */
        void updateVecNorm ();

        /**
         * Updates the array containing the basis vectors norms from the `d`-th 
         * component to the last by recomputing them.
         * */
        void updateVecNorm (const int & d);

        /**
         * Updates the array containing the dual basis vectors norms by recomputing 
         * them.
         * */
        void updateDualVecNorm ();

        /**
         * Updates the array containing the dual basis vectors norms from the `d`-th 
         * component to the last by recomputing them.
         * */
        void updateDualVecNorm (const int & d);

        /**
         * Updates the `i`-th value of the array containing the norms of the 
         * basis vectors by recomputing it using the `L2NORM`.
         */
        void updateScalL2Norm (const int i);

        /**
         * Updates the `k1`-th to the `k2-1`-th values of the array containing 
         * the norms of the basis vectors by recomputing them using the `L2NORM`.
         */
        void updateScalL2Norm (const int k1, const int k2);

        /**
         * Updates the `i`-th value of the array containing the norms of the 
         * dual basis vectors by recomputing it using the `L2NORM`.
         */
        void updateDualScalL2Norm (const int i);

        /**
         * Updates the `k1`-th to the `k2-1`-th values of the array containing 
         * the norms of the dual basis vectors by recomputing them using the `L2NORM`.
         */
        void updateDualScalL2Norm (const int k1, const int k2);

        /**
         * Exchanges vectors `i` and `j` in the basis. This also changes the 
         * dual basis vectors and the arrays containing secondary information
         * about the two basis (like the norms) accordingly.
         */
        void permute (int i, int j);

        /**
         * Returns `true` if the dual basis contained in the object really is 
         * the dual of the basis, and `false` otherwise. This also returns false
         * if no dual has been specified.
         */
        bool checkDuality();

        /**
         * Sorts the basis vectors with indices greater of equal to \f$d\f$ by 
         * increasing length. The dual vectors are permuted accordingly. Assumes 
         * that the lengths of the corresponding basis vectors are up to date.
         */
        void sort (int d);

        /**
         * Returns a string with the primal basis and its norms.
         */
        std::string toStringBasis() const;

        /**
         * Returns a string with the dual basis and its norms.
         */
        std::string toStringDualBasis() const;

        /**
         * Writes the lattice and its parameters on standard output. This prints
         * the dimension, the norm used, the basis and dual basis vectors and
         * the basis and dual basis vector norms.
         */
        void write () const;

      protected:

        /**
         * Each row of this matrix represents a vector in the basis of the
         * lattice.
         */
        BasIntMat m_basis;

        /**
         * Each row of this matrix represents a vector in the dual basis of the
         * lattice.
         */
        BasIntMat m_dualbasis;


        /**
         * The dimension of the lattice. The dimension is both the number of 
         * vectors in the basis, and the dimension \f$d\f$ of \f$\mathbb{Z}^d\f$
         * containing the lattice.
         */
        int m_dim;

        /**
         * The NormType used in the reduction and for this lattice.
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
         * The `m` used for rescaling the lattice.
         */
        BasInt m_modulo;

        /**
         * If m_withDual is `true` a dual basis has been specified, otherwise it
         * is `false`.
         */
        bool m_withDual;

        /**
         * This table is used in the Minkowski reduction, but it's usage is quite
         * obscure.
         */
        //bool *m_xx;

    }; // class IntLatticeBasis

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::IntLatticeBasis (const int dim,
        NormType norm):
      m_dim (dim),
      m_norm (norm),
      m_modulo(0),
      m_withDual(false)
      //m_xx(0)

  {
    this->m_basis.resize(dim,dim);
    this->m_vecNorm.resize (dim);
    initVecNorm();
  }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::IntLatticeBasis (
        const BasIntMat basis, const int dim, NormType norm):
      m_basis (basis),
      m_dim (dim),
      m_norm (norm),
      m_modulo(0),
      m_withDual(false)
      //m_xx(0)
  {
    this->m_vecNorm.resize (dim);
    initVecNorm();
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::IntLatticeBasis (
        const BasIntMat primalbasis,
        const BasIntMat dualbasis,
        const Int modulo,
        const int dim,
        NormType norm):
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>(
          primalbasis, dim, norm)
  {
    this->m_dualbasis = BasIntMat(dualbasis);
    this->m_withDual = true;
    this->m_dualvecNorm.resize (dim);
    setDualNegativeNorm();
    this->m_modulo = modulo;
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::IntLatticeBasis (
        const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat):
      m_dim (lat.getDim ()),
      m_norm (lat.getNorm ())
      //m_xx(0)
  {
    copyBasis (lat);
  }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::~IntLatticeBasis ()
    {
      kill();
      this->m_basis.BasIntMat::clear ();
      this->m_dualbasis.BasIntMat::clear ();
      this->m_vecNorm.clear ();
      this->m_dualvecNorm.clear ();
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::kill ()
    {
      //delete [] this->m_xx;
      //this->m_xx = 0;
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::copyBasis (
        const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat)
    {
      //if(m_dim == lat.m_dim)

      this->m_basis = BasIntMat(lat.m_basis);
      this->m_dualbasis = BasIntMat(lat.m_dualbasis);
      this->m_vecNorm = DblVec(lat.m_vecNorm);
      this->m_dualvecNorm = DblVec(lat.m_dualvecNorm);
      this->m_withDual = lat.m_withDual;
      this->m_modulo = lat.m_modulo;
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::copyBasis (
        const IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat, long n)
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
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::initVecNorm ()
    {
      //this->m_xx = new bool[this->m_dim];
      for(int i = 0; i < this->m_dim; i++){
        this->m_vecNorm[i] = -1;
        //this->m_xx[i] = true;
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::setNegativeNorm ()
    {
      for (int i = 0; i<this->m_dim; i++){
        this->m_vecNorm[i] = -1;
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::setDualNegativeNorm ()
    {
      for(int i = 0; i < this->m_dim; i++){
        this->m_dualvecNorm[i] = -1;
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateVecNorm ()
    {
      updateVecNorm (0);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateVecNorm (const int & d)
    {
      assert (d >= 0);

      for (int i = d; i < this->m_dim; i++) {
        NTL::matrix_row<BasIntMat> row(this->m_basis, i);
        if (this->m_norm == L2NORM) {
          ProdScal<Int> (row, row, this->m_dim, this->m_vecNorm[i]);
        } else {
          CalcNorm <BasIntVec, Dbl> (row, this->m_dim, this->m_vecNorm[i],
              this->m_norm);
        }
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateDualVecNorm ()
    {
      updateDualVecNorm (0);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateDualVecNorm (const int & d)
    {
      assert (d >= 0);

      for (int i = d; i < this->m_dim; i++) {
        NTL::matrix_row<BasIntMat> row(this->m_dualbasis, i);
        if (this->m_norm == L2NORM) {
          ProdScal<Int> (row, row, this->m_dim, this->m_dualvecNorm[i]);
        } else {
          CalcNorm <BasIntVec, Dbl> (row, this->m_dim, this->m_dualvecNorm[i],
              this->m_norm);
        }
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateScalL2Norm (const int i)
    {
      NTL::matrix_row<BasIntMat> row(this->m_basis, i);
      ProdScal<Int> (row, row, this->m_dim, this->m_vecNorm[i]);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateScalL2Norm (
        const int k1, const int k2)
    {
      for (int i = k1; i < k2; i++) {
        updateScalL2Norm(i);
      }
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateDualScalL2Norm (
        const int i)
    {
      NTL::matrix_row<BasIntMat> row(this->m_dualbasis, i);
      ProdScal<Int> (row, row, this->m_dim, this->m_dualvecNorm[i]);
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::updateDualScalL2Norm (
        const int k1,
        const int k2)
    {
      for (int i = k1; i < k2; i++) {
        updateDualScalL2Norm(i);
      }
    }


  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::permute (int i, int j)
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
      //bool b = this->m_xx[j];
      //this->m_xx[j] = this->m_xx[i];
      //this->m_xx[i] = b;
    }


  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    bool IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::checkDuality ()
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

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::sort (int d)
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

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    void IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::write () const
    {
      std::cout << "Dim = " << this->m_dim << " \n \n";
      std::cout << std::setprecision(10) << "Primal basis vectors:\n";
      for (int i = 0; i < this->m_dim; i++) {
        std::cout << this->m_basis[i];
        //for (int j = 0; j < this->m_dim; j++) {
        //  std::cout <<  this->m_basis(i,j);
        //}
        std::cout << "\n";
      }
      std::cout << "\nDual basis vectors:\n";
      for (int i = 0; i < this->m_dim; i++) {
        if(this->m_withDual){
          std::cout << this->m_dualbasis[i];
          //for (int j = 0; j < this->m_dim; j++) {
          //  std::cout << this->m_dualbasis(i,j);
          //}
          std::cout << "\n";
        }
      }
      std::cout << "\n";
      std::cout << "Norm used: " << toStringNorm(this->m_norm) << "\n"
        << std::endl;
      std::cout << "Norm of each Basis vector: \n";
      std::cout << "Primal";
      if(this->m_withDual)
        std::cout << "\t\tDual\n";
      std::cout << "\n";

      for (int i = 0; i < this->m_dim; i++) {
        if (this->m_vecNorm[i] < 0) {
          std::cout << "NaN OR Not computed";
        } else {
          if (this->m_norm == L2NORM) {
            std::cout << NTL::sqrt(this->m_vecNorm[i]);
          } else {
            std::cout << this->m_vecNorm[i];
          }
        }
        std::cout << "\t";
        if(this->m_withDual){
          if (this->m_dualvecNorm[i] < 0)
            std::cout << "NaN OR Not computed";
          else{
            if (this->m_norm == L2NORM) {
              std::cout << NTL::sqrt(this->m_dualvecNorm[i]);
            } else {
              std::cout << this->m_dualvecNorm[i];
            }
          }
        }
        std::cout << "\n";
      }
      std::cout << std::endl;
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    std::string IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::toStringBasis () const
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

      os << "  Norms:\n";
      os << "    [";
      for (int i = 0; i < this->m_dim; i++) {
        if (this->m_vecNorm[i] < 0) {
          os << "-1" << " ";
        } else {
          os << this->m_vecNorm[i] << " ";
        }
      }
      os << "]" << std::endl;
      return os.str ();
    }

  /*=========================================================================*/

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
    std::string IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::toStringDualBasis 
    () const
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

      os << "  Norms:\n";
      os << "    [";
      for (int i = 0; i < this->m_dim; i++) {
        if (this->m_dualvecNorm[i] < 0) {
          os << "-1" << " ";
        } else {
          os << this->m_dualvecNorm[i] << " ";
        }
      }
      os << "]" << std::endl;
      return os.str ();
    }

  extern template class IntLatticeBasis<std::int64_t, std::int64_t, double, double>;
  extern template class IntLatticeBasis<NTL::ZZ, NTL::ZZ, double, double>;
  extern template class IntLatticeBasis<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>;

} // namespace LatticeTester

#endif
