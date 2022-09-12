// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER_BASISCONSTRUCTION_H
#define LATTICETESTER_BASISCONSTRUCTION_H

#include "NTL/LLL.h"

#include "latticetester/IntLatticeBase.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"

namespace LatticeTester {

template<typename IntMat>
struct LLLConstr {
	void LLLConstruction(IntMat &matrix);
};

/**
 * This class offers methods to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 * The implementation relies on NTL and uses NTL matrices.
 *
 * NTL already offers a very efficient method to construct an LLL-reduced basis from a set
 * of generating vectors.  This is the most effective way of constructing a basis
 * and it is encapsulated in the `LLLConstruction` method given below.
 * We also offer an alternative that constructs a triangular basis, in `GCDTriangularBasis`.
 * To compute the $m$-dual of a given basis, we have a general method implemented in
 * `mDualComputation`, and a faster method in `mDualTriangular` that works only when
 * the basis is upper-triangular.
 * The methods `Util::Triangularization` and `Util::CalcDual` do essentially the same
 * things; however, the methods given here perform more verifications.
 *  ***  We should compare the speeds
 *
 * A few tips about the usage of this class:
 * - Prefer the usage of NTL types when using this module. The methods here do not
 *   have any kind of overflow detection.
 * - Reduce the basis before doing a triangularization. Reducing a basis with
 *   LLL is much faster than the GCDConstruction and seems to make this operation
 *   easier to perform.    ***  To be tested again.
 * - Use specialized methods. With a more in depth knowledge of your problem, it
 *   is possible that there are much more efficient ways to build a basis and its
 *   dual (and/or those matrices may already be triangular).
 */

template<typename Int> class BasisConstruction {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	struct LLLConstr<IntMat> spec;   // Why this?

public:

	/**
	 * This functions takes a set of generating vectors of a vector space and
	 * finds a basis for this space by applying LLL reduction, using the NTL implementation.
	 * This is much faster than applying GCDConstruction, but it does not provide a triangular basis.
	 * It is a good idea to use it before computing a triangular basis.
	 */
	void LLLConstruction(IntMat &matrix);

	/**
	 * This function does essentially the same thing as `Util::Triangularization`.
	 * It uses a form of Gaussian elimination to obtain an upper triangular basis
	 * for the smallest lattice that contains the generating vectors which are the
	 * rows of the given matrix `matrix`.
	 * In each column, it applies Euclid's algorithm to the elements under the
	 * diagonal and change the corresponding rows to set all these elements to
	 * zero except the one in the diagonal. The allowed row operations are only
	 *   - Multiply a row by \f$-1\f$,
	 *   - Add an integer multiple of row \f$i\f$ to row \f$j\f$ for \f$i \neq j\f$,
	 *   - Swap row \f$i\f$ with row \f$j\f$.
	 * All these operations can be performed modulo the scaling factor `m`.
	 * After constructing this basis, the algorithm eliminates negative
	 * coefficients in the matrix.
	 * Warning: In this implementation, the numbers below the diagonal can grow
	 * very large, so the method may require a lot of memory.
	 *  ***  But everything should be done modulo m ???
	 */
	void GCDTriangularBasis(IntMat &matrix);

	/**
	 * This function does essentially the same thing as `Util::CalcDual`.
	 * It assumes that `matrix` contains a triangular basis of the primal lattice
	 * scaled by the factor `m`.  It computes and returns the `m`-dual basis
	 * in `dualMatrix`.  This function uses the method described in \cite rCOU96a.
	 * Since `A=matrix` is upper triangular, `B=mdualMatrix` will be lower triangular
	 * with `A(i,i)*B(i,i) = m` for all `i` and \f$ A_i \cdot B_j = 0\f$ for
	 * \f$i\neq j\f$. To get the second condition, we simply have to
	 * recursively take for each line
	 * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
	 * This is much faster than a traditional solving of a linear system.
	 */
	void mDualTriangular(IntMat &matrix, IntMat &dualMatrix, Int m);

	/**
	 * This does the same thing as DualConstruction(), but is much slower. This
	 * is here simply for the sake of comparison and should not be used in practice.
	 * Suppose the basis matrix contains basis vectors on its lines and is
	 * \f$p\times q\f$ where \f$q \geq p\f$. We can compute the \f$m\f$-dual
	 * as follows.
	 * Let's note the basis matrix `V` and the dual matrix `W` and have lines
	 * of `W` also contain dual basis vectors in its lines. We know that
	 * \f$VW^t = mI_{p\times p}\f$ where \f$I\f$ is the identity matrix. Now, in
	 * the case of a \f$m\f$-dual computation, we can assume that all the
	 * arithmetic is done modulo \f$m\f$.

	 void DualSlow(IntMat& matrix, IntMat& dualMatrix, Int& m);
	 */

	/**
	 * This function assumes that `matrix` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular.
	 * It computes and returns the `m`-dual basis in `dualMatrix`.
	 * ***  TO BE IMPLEMENTED **
	 */
	void mDualComputation(IntMat &matrix, IntMat &dualMatrix, Int m);

	/**
	 * This method constructs a basis for the projection `proj` of the basis `in`.
	 * using the `LLLConstruction` method, and puts it in `out`. This will
	 * overwrite the lattice basis in `out`, changing also the dimension.
	 * The returned basis is not triangular in general.
	 */
	template<typename Int, typename Real, typename RealRed>
	void ProjectionConstruction(IntLatticeBase<Int, Real, RealRed> &in,
			IntLatticeBase<Int, Real, RealRed> &out, const Coordinates &proj);
};

//============================================================================
// Implementation

template<>
void LLLConstr<NTL::matrix<std::int64_t>>::LLLConstruction(
		NTL::matrix<std::int64_t> &matrix);

template<>
void LLLConstr<NTL::matrix<NTL::ZZ>>::LLLConstruction(
		NTL::matrix<NTL::ZZ> &matrix);

template<typename Int>
void BasisConstruction<Int>::LLLConstruction(IntMat &matrix) {
	spec.LLLConstruction(matrix);
}

template<typename Int>
void BasisConstruction<Int>::GCDTriangularBasis(IntMat &matrix) {
	// It is important to note that the lines of matrix are the basis vectors
	long rows = matrix.NumRows();
	long cols = matrix.NumCols();
	long max_rank = rows < cols ? rows : cols;
	long rank = 0;
	// The basis will have at most max_rank vectors.
	Int q;
	for (long i = 0; i < max_rank; i++) {
		// We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
		// algorithm and applying transformations to the matrix
		for (long j = i + 1; j < rows; j++) {
			while (matrix[j][i] != 0) {
				matrix[i].swap(matrix[j]);
				q = matrix[j][i] / matrix[i][i];
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
		if (matrix[i][i] != 0) {
			rank++;
			if (matrix[i][i] < 0)
				matrix[i] *= Int(-1);
		}
	}
	// We remove zero vectors from the basis.
	matrix.SetDims(rank, cols);
}

//============================================================================

/**
 * This algorithm calculates the dual as well as the `m` used for rescaling.
 * It checks if this `m` divides the modulo given to the algorithm.
 * Right now, this assumes the basis is triangular, might need to change it.

 template<typename Int>
 void BasisConstruction<Int>::DualSlow(IntMat& matrix,
 IntMat& dualMatrix, Int& modulo)
 {
 // We need to have a triangular basis matrix
 if (! CheckTriangular(matrix, matrix.NumRows(), modulo))
 GCDConstruction(matrix);
 long dim = matrix.NumRows();
 if (dim != matrix.NumCols()) {
 std::cout << "matrix has to be square, but dimensions do not fit.\n";
 return;
 }
 Int m(1);
 NTL::ident(dualMatrix, dim);
 Int gcd;
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
 */

template<typename Int>
void BasisConstruction<Int>::mDualTriangular(IntMat &matrix, IntMat &dualMatrix,
		Int m) {
	// We need to have a triangular basis matrix
	if (!CheckTriangular(matrix, matrix.NumRows(), Int(0)))
		GCDConstruction(matrix);
	long dim = matrix.NumRows();
	if (dim != matrix.NumCols()) {
		std::cout << "matrix has to be square, but dimensions do not fit.\n";
		return;
	}
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);
		return;
	}
	dualMatrix.SetDims(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(dualMatrix(i, j));
		if (!NTL::IsZero(matrix(i, i))) {
			Int gcd = NTL::GCD(m, matrix(i, i));
			m *= matrix(i, i) / gcd;
			dualMatrix *= matrix(i, i) / gcd;
		}

		DivideRound(m, matrix(i, i), dualMatrix(i, i));
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(dualMatrix(i, j));
			for (int k = j + 1; k <= i; k++)
				dualMatrix(i, j) += matrix(j, k) * dualMatrix(i, k);
			if (dualMatrix(i, j) != 0)
				dualMatrix(i, j) = -dualMatrix(i, j);
			if (!NTL::IsZero(dualMatrix(i, j) % matrix(j, j))) {
				Int gcd = NTL::GCD(dualMatrix(i, j), matrix(j, j));
				m *= matrix(j, j) / gcd;
				dualMatrix *= matrix(j, j) / gcd;
			}
			DivideRound(dualMatrix(i, j), matrix(j, j), dualMatrix(i, j));
		}
	}
}

//============================================================================

template<typename Int>
void BasisConstruction<Int>::mDualComputation(IntMat &matrix,
		IntMat &dualMatrix, Int m) {
	// **  TO BE IMPLEMENTED  **
}

//============================================================================

template<typename Int>
template<typename Int, typename Real, typename RealRed>
void BasisConstruction<Int>::ProjectionConstruction(
		IntLatticeBase<Int, Real, RealRed> &in,
		IntLatticeBase<Int, Real, RealRed> &out, const Coordinates &proj) {
	std::size_t dim = proj.size();
	unsigned int lat_dim = in.getDim();
	if (dim > lat_dim)
		MyExit(1, "Coordinates do not match the dimensions of `in`.");
	IntMat new_basis, tmp(NTL::transpose(in.getBasis()));
	new_basis.SetDims(dim, tmp.NumRows());
	tmp = NTL::transpose(tmp);
	auto it = proj.cbegin();
	for (std::size_t i = 0; i < dim; i++) {
		if (*it <= lat_dim)
			new_basis[i] = tmp[*it];
		else
			MyExit(1, "Coordinates do not match the dimensions of `in`.");
		it++;
	}
	new_basis = NTL::transpose(new_basis);
	LLLConstruction(new_basis);
	out = IntLatticeBase<Int, Real, RealRed>(new_basis, dim, in.getNorm());
}

extern template class BasisConstruction<std::int64_t> ;
extern template class BasisConstruction<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
