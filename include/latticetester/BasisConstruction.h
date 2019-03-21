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

#ifndef LATTICETESTER_BASISCONSTRUCTION_H
#define LATTICETESTER_BASISCONSTRUCTION_H

#include "NTL/LLL.h"

#include "latticetester/IntLatticeBasis.h"
#include "latticetester/ntlwrap.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"

namespace LatticeTester {

  template <typename BasIntMat>
    struct LLLConstr {
      void LLLConstruction(BasIntMat& matrix);
    };

/**
 * This class implements general methods to perform a lattice basis
 * construction from a set of vectors, as well as general methods to obtain
 * the dual lattice basis depending on the current lattice basis.
 *
 * This module only works with NTL matrices because those are objects
 * aware of their shape.
 *
 * What is done here seems like very simple linear algebra BUT it is not. The
 * fact that we operate only on the ring of integers makes it so that these
 * simple algorithms can become way too slow very rapidly because the numbers
 * they manipulate grow very fast. This also means that in many cases the
 * usage of standard `long` type integers will overflow.
 *
 * ### A note on basis construction
 * Although basis construction is mathematically simple (and also
 * computationnally simple when in the field of real numbers), things can
 * get messy. In higher dimensions (~30), the coefficients in the matrix start
 * getting really big (\f$\gg 2^{64}\f$) and needing a lot of memory
 * (run the BasisConstruction example and see this for yourself). Therefore, we
 * give users a few tips about the usage of this module.
 * - Prefer the usage of NTL types when using this module. The functions don't
 *   have any kind of overflow detection.
 * - Reduce the basis before doing a triangularization. Reducing a basis with
 *   LLL is much faster than GCDConstruction and seems to make this operation
 *   easier to perform.
 * - Use specialized methods. With a more in depth knowledge of your problem, it
 *   is possible that there are much more efficient ways to build a basis and its
 *   dual (and/or that those matrices will already be triangular).
 * */
template<typename BasInt> class BasisConstruction{

  private:
    typedef NTL::vector<BasInt> BasIntVec;
    typedef NTL::matrix<BasInt> BasIntMat;
    struct LLLConstr<BasIntMat> spec;

  public:
    /**
     * This functions takes a set of generating vectors of a vector space and
     * finds a basis for this space whilst applying LLL reduction. This is
     * much faster than applying GCDConstruction, but it doesn't help building
     * the dual.
     * */
    void LLLConstruction(BasIntMat& matrix);

    /**
     * This function does some kind of Gaussian elimination on
     * \f$\text{span}(\text{matrix}^t)\f$ as a vector space on \f$\mathbb{Z}\f$ to
     * obtain a basis of this vector space (lattice).
     * This takes a matrix `matrix` where each row `v[0], ..., v[n]` is
     * considered as a vector. This basically performs, for each column,
     * Euclid's algorithm on the rows under the diagonal to set all
     * coefficients to zero except the one in the diagonal. Finding the GCD
     * of all the coefficients from the diagonal and under allows to perform
     * Gaussian elimination using only the 3 allowed operations that do not
     * change the span of the vectors:
     *   - Multiply a line by \f$-1\f$,
     *   - Add an integer multiple of line \f$i\f$ to line \f$j\f$ for \f$i \neq j\f$,
     *   - Swap line \f$i\f$ and line \f$j\f$.
     * After constructing this basis, the algorithm eliminates negative
     * coefficients in the matrix.
     *
     * WATCH OUT. This function (and building mecanism) are very memory heavy
     * has the numbers below the diagonal can grow very big.
     * */
    void GCDConstruction(BasIntMat& matrix);

    /**
     * This method builds the dual of `matrix` in `dualMatrix` and multiplies
     * modulo by the rescaling factor used. This function uses the clever method
     * developed in \cite rCOU96a that notes that for `B` to be a `m`-dual to
     * `A`, we have to have that \f$AB^t = mI\f$. It is quite easy to show that,
     * knowing `A` is upper triangular, `B` will be a lower triangular matrix
     * with `A(i,i)*B(i,i) = m` for all `i` and \f$ A_i \cdot B_j = 0\f$ for
     * \f$i\neq j\f$. And that to get the second condition, we simply have to
     * recursively take for each line
     * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
     *
     * This does the computation much faster than doing a traditional solving
     * of a linear system.
     * */
    void DualConstruction(BasIntMat& matrix, BasIntMat& dualMatrix,
        BasInt& modulo);

    /**
     * This does the same thing as DualConstruction(), but is much slower. This
     * is here simply for the sake of comparison and should not be used in
     * practice.
     *
     * Suppose the basis matrix contains basis vectors on its lines and is
     * \f$p\times q\f$ where \f$q \geq p\f$. We can compute the \f$m\f$-dual
     * as follows.
     * Let's note the basis matrix `V` and the dual matrix `W` and have lines
     * of `W` also contain dual basis vectors in its lines. We know that
     * \f$VW^t = mI_{p\times p}\f$ where \f$I\f$ is the identity matrix. Now, in
     * the case of a \f$m\f$-dual computation, we can assume that all the
     * arithmetic is done modulo \f$m\f$ and that \f$m\f$ is prime??
     * */
    void DualSlow(BasIntMat& matrix, BasIntMat& dualMatrix,
        BasInt& modulo);

    /**
     * This is a method that does the general construction of the projection
     * `proj` of the basis `in` and puts it in `out`. This will completely
     * overwrite the lattice basis in `out`, changing the dimension. If one
     * wants to compute the dual for this projection, it has to be done
     * afterwards with the DualConstruction() method.
     * */
    template<typename Int, typename Dbl, typename RedDbl>
      void ProjectionConstruction(IntLatticeBasis<Int, BasInt, Dbl, RedDbl>& in,
          IntLatticeBasis<Int, BasInt, Dbl, RedDbl>& out, Coordinates& proj);
  };

  //============================================================================
  // Implementation

  template<typename BasInt>
    void BasisConstruction<BasInt>::LLLConstruction(BasIntMat& matrix){
      spec.LLLConstruction(matrix);
    }

  template<typename BasInt>
    void BasisConstruction<BasInt>::GCDConstruction(BasIntMat& matrix)
  {
    // It is important to note that the lines of matrix are the basis vectors
    long rows = matrix.NumRows();
    long cols = matrix.NumCols();
    long max_rank = rows < cols ? rows : cols;
    long rank = 0;
    // The basis will have at most max_rank vectors.
    BasInt q;
    for (long i = 0; i < max_rank; i++) {
      // We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
      // algorithm and applying transformations to the matrix
      for (long j = i+1; j < rows; j++) {
        while (matrix[j][i] != 0) {
          matrix[i].swap(matrix[j]);
          q = matrix[j][i]/matrix[i][i];
          matrix[j] -= q * matrix[i];
        }
      }
      // We make sure that the coefficients are positive.
      // This is because the algorithms we use work for positive vectors.
      // if (matrix[i][i] < 0) matrix[i] *= -1;
      // for (long j = i-1; j >= 0; j--) {
      //   if (matrix[j][i] < 0) {
      //     if (-matrix[j][i] > matrix[i][i]){
      //       q = -matrix[j][i]/matrix[i][i] + 1;
      //       matrix[j] += q * matrix[i];
      //     } else {
      //       matrix[j] += matrix[i];
      //     }
      //   }
      // }
      if (matrix[i][i] != 0) rank++;
    }
    // We remove zero vectors from the basis.
    matrix.SetDims(rank, cols);
  }

  //============================================================================

  /**
   * This algorithm calculates the dual as well as the `m` used for rescalling.
   * It checks if this `m` divides the modulo given to the algorithm.
   * Right now, this assumes the basis is triangular, might need to change it.
   * */
  template<typename BasInt>
    void BasisConstruction<BasInt>::DualSlow(BasIntMat& matrix,
        BasIntMat& dualMatrix, BasInt& modulo)
  {
    // We need to have a triangular basis matrix
    if (! CheckTriangular(matrix, matrix.NumRows(), modulo))
      GCDConstruction(matrix);
    long dim = matrix.NumRows();
    if (dim != matrix.NumCols()) {
      std::cout << "matrix has to be square, but dimensions do not fit.\n";
      return;
    }
    BasInt m;
    m = 1;
    NTL::ident(dualMatrix, dim);
    BasInt gcd;
    for (long i = dim-1; i>=0; i--) {
      for (long j = i+1; j < dim; j++) {
        dualMatrix[i] -= matrix[i][j] * dualMatrix[j];
        matrix[i] -= matrix[i][j] * matrix[j];
      }
      gcd = matrix[i][i];
      for (long j = i; j < dim; j++) {
        gcd = NTL::GCD(gcd, dualMatrix[i][j]);
      }
      gcd = matrix[i][i] / gcd;
      if (gcd != 1) {
        dualMatrix *= gcd;
        m *= gcd;
      }
      for (long j = i; j < dim; j++) {
        dualMatrix[i][j] /= matrix[i][i];
      }
      matrix[i][i] = 1;
    }
  }

  template<typename BasInt>
    void BasisConstruction<BasInt>::DualConstruction(BasIntMat& matrix,
        BasIntMat& dualMatrix, BasInt& modulo)
  {
    // We need to have a triangular basis matrix
    if (! CheckTriangular(matrix, matrix.NumRows(), BasInt(0)))
      GCDConstruction(matrix);
    long dim = matrix.NumRows();
    if (dim != matrix.NumCols()) {
      std::cout << "matrix has to be square, but dimensions do not fit.\n";
      return;
    }
    if (modulo < 1) {
      std::cerr << "modulo has to be a positive integer.\n";
      exit(1);
      return;
    }
    dualMatrix.SetDims(dim, dim);
    for (int i = 0; i < dim; i++) {
      for (int j = i + 1; j < dim; j++)
        NTL::clear (dualMatrix(i,j));

      if (!NTL::IsZero(matrix(i,i))) {
        BasInt gcd = NTL::GCD(modulo, matrix(i,i));
        modulo *= matrix(i,i) / gcd;
        dualMatrix *= matrix(i,i) / gcd;
      }

      DivideRound (modulo, matrix(i,i), dualMatrix(i,i));

      for (int j = i - 1; j >= 0; j--) {

        NTL::clear (dualMatrix(i,j));

        for (int k = j + 1; k <= i; k++)
          dualMatrix(i,j) += matrix(j,k) * dualMatrix(i,k);

        if (dualMatrix(i,j) != 0)
          dualMatrix(i,j) = -dualMatrix(i,j);


        if (!NTL::IsZero(dualMatrix(i,j) % matrix(j,j))) {
          BasInt gcd = NTL::GCD(dualMatrix(i,j), matrix(j,j));
          modulo *= matrix(j,j) / gcd;
          dualMatrix *= matrix(j,j) / gcd;
        }

        DivideRound (dualMatrix(i,j), matrix(j,j), dualMatrix(i,j));
      }

    }
  }

  //============================================================================

    template<typename BasInt>
    template<typename Int, typename Dbl, typename RedDbl>
      void BasisConstruction<BasInt>::ProjectionConstruction(
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>& in,
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>& out,
      Coordinates& proj) {
        std::size_t dim = proj.size();
        if (dim > in.getDim())
            MyExit(1, "Coordinates do not match 'in' dimension.");
        BasIntMat new_basis, tmp(in.getBasis());
        new_basis.SetDims(dim, tmp.NumRows());
        tmp = NTL::transpose(tmp);
        auto it = proj.cbegin();
        for (std::size_t i = 0; i < dim; i++) {
          if (*it <= in.getDim())
            new_basis[i] = tmp[*it];
          else
            MyExit(1, "Coordinates do not match 'in' dimension.");
        }
        new_basis = NTL::transpose(new_basis);
        LLLConstruction(new_basis);
        out = IntLatticeBasis<Int, BasInt, Dbl, RedDbl>(new_basis, dim, in.getNorm());
      }

  extern template class BasisConstruction<std::int64_t>;
  extern template class BasisConstruction<NTL::ZZ>;

} // end namespace LatticeTester

#endif
