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
#include <NTL/mat_GF2.h>
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"
#include "NTL/LLL.h"
#include "NTL/HNF.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/IntLattice.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <type_traits>

namespace LatticeTester {

/**
 * This class offers methods to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 * The implementation relies on NTL and uses NTL matrices.
 * When the basis turns out to have fewer rows than columns, some of the methods
 * add implicitly the rescaled unit vectors to the set of generating vectors.
 * Then the basis matrix is always square.
 *
 * NTL already offers a very efficient method to construct an LLL-reduced basis from a set
 * of generating vectors.  This is the most effective way of constructing a basis
 * and it is encapsulated in the `LLLConstruction0` method given below.
 * We also offer an alternative that constructs a triangular basis, in `GCDTriangularBasis`.
 * To compute the $m$-dual of a given basis, we have a general method implemented in
 * `mDualComputation`, and a faster method in `mDualTriangular` that works only when
 * the basis is upper-triangular.
 * The methods `Util::Triangularization` and `Util::calcDual` do essentially the same
 * things; however, the methods given here perform more verifications.
 *
 *  ***  In the end, we want to remove these methods from Util and move them here!
 *
 *  UPDATE THIS AFTER WE ARE DONE WITH THE TESTS:
 *
 * A few tips about the usage of this class:
 * - Prefer the usage of NTL types when using this module. The methods here do not
 *   have any kind of overflow detection.
 * - Reduce the basis before doing a triangularization. Reducing a basis with
 *   LLL is much faster than the GCDTriangularBasis and seems to make this operation
 *   easier to perform.    ***  To be tested again.
 * - Use specialized methods. With a more in depth knowledge of your problem, it
 *   is possible that there are much more efficient ways to build a basis and its
 *   dual (and/or those matrices may already be triangular).
 */

template<typename Int> class BasisConstruction {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;

public:

	/**
	 * This function takes a set of generating vectors of a lattice in matrix `gen` and
	 * finds a lattice basis by applying LLL reduction with the given value of `delta`,
	 * using the NTL implementation specified by `prec`. See the class `EnumTypes` and
	 * the documentation of LLL in NTL for the meaning and choices for `prec`.
	 * This function is implemented only for the \texttt{NTL::ZZ} type, because it uses NTL.
	 * It *does not* assume that all vectors m e_i belong to the lattice, so
	 * it may return a basis matrix that has fewer rows than columns!
	 * To make sure that these vectors belong to the lattice, we can add them
	 * explicitly beforehand to the set of generating vectors, or call the next method.
	 */
	void LLLConstruction0(IntMat &gen, double delta = 0.999999,
			PrecisionType precision = DOUBLE);

   /**
    * Similar to `LLLConstruction`, except that in case the set of generating
    * vectors do not generate a full-dimensional lattice, it adds the vectors
    * m e_i to the generating set, so it always returns a square matrix.
    */
	void LLLBasisConstruction(IntMat &gen, Int &m, double delta = 0.999999,
			PrecisionType precision = DOUBLE);

	/**
	 * Takes a set of generating vectors in the matrix `gen` and iteratively
	 * transforms it into a lower triangular lattice basis into the matrix `basis`.
	 * `gen` and `basis` must have the same number of rows and the same number of columns.
	 * All the computations are done modulo the scaling factor `m`.
	 * After the execution, `gen` will contain irrelevant information (garbage)
	 * and `basis` will contain an upper triangular basis.
	 * Perhaps with zero rows at the end, in general, unless we assume implicitly that
	 * all vectors of the form m e_i are in the generating set.
	 * The algorithm is explained in the \lattester{} guide.
	 */
	void lowerTriangularBasis(IntMat &gen, IntMat &basis, Int &m);

	/**
	 * Similar to `lowerTriangularBasis`, except that the returned basis is upper triangular.
	 */
	void upperTriangularBasis(IntMat &gen, IntMat &basis, Int &m);

	/**
	 * This is an old implementationSame that uses a form of Gaussian elimination to
	 * obtain an upper triangular basis for the smallest lattice that contains the generating
	 * vectors which are the rows of the given `matrix`. It returns the basis in the same `matrix`.
	 * This function *does not* assume that all vectors `m e_i` belong to the lattice, and
	 * it may return a basis matrix that has fewer rows than columns!
	 */
	void GCDTriangularBasis(IntMat &gen, Int &m);

	/**
	 * Takes an upper triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
	 * The method assumes that each coefficient on the diagonal of `basis` is nonzero and divides `m`.
	 * That is, the basis matrix must be square and invertible.
	 * The algorithm is described in the Lattice Tester guide \cite iLEC22l.
	 * Since the basis is upper triangular, its m-dual will be lower triangular.
	 */
	void mDualUpperTriangular(const IntMat &basis, IntMat &basisDual,
			const Int &m);

	/**
	 * This function does essentially the same thing as `mDualUpperTriangular`, but it is
	 * slightly different and slower. It uses the method described in \cite rCOU96a.
	 */
	void mDualUpperTriangular96(IntMat &basis, IntMat &basisDual, Int &m);

	/**
	 * This function assumes that `basis` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular, and it returns in `basisDual`
	 * the m-dual basis.
	 */
	void mDualBasis(IntMat &basis, IntMat &basisDual, Int &m);

	/**
	 * Constructs a basis for the projection `proj` of the lattice `in`,
	 * using LLLConstruction, and puts it in `out`. The basis is not triangular.
	 * This will overwrite the lattice basis in `out` and change the dimension.
	 * It does not update the dual.
	 */
	template<typename Real>
	void projectionConstructionLLL(IntLattice<Int, Real> &in,
			IntLattice<Int, Real> &out, const Coordinates &proj);
	

	class ZZWatcher {
	public:
	   NTL::ZZ& watched;
	   explicit
	   ZZWatcher(NTL::ZZ& _watched) : watched(_watched) {}

	   ~ZZWatcher() { watched.KillBig(); }
	};

	#define NTL_ZZNewRegister(x) NTL_TLS_LOCAL(NTL::ZZ, x); ZZWatcher _WATCHER__ ## x(x)
	
};

//============================================================================
// Implementation

template<typename Int>
void BasisConstruction<Int>::LLLConstruction0(IntMat &gen, double delta,
		PrecisionType prec) {
	std::cerr << "LLLConstruction0 can only be done with NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

template<>
void BasisConstruction<NTL::ZZ>::LLLConstruction0(NTL::matrix<NTL::ZZ> &gen,
		double delta, PrecisionType prec) {
	long rank;
	switch (prec) {
	case DOUBLE:
		rank = NTL::LLL_FP(gen, delta);
		break;
	case QUADRUPLE:
		rank = NTL::LLL_QP(gen, delta);
		break;
	case XDOUBLE:
		rank = NTL::LLL_XD(gen, delta);
		break;
	case RR:
		rank = NTL::LLL_RR(gen, delta);
		break;
	default:
		std::cerr << "LLLConstruction0: unknown precision type.\n";
	}
	long num = gen.NumRows();
	for (long i = 0; i < rank; i++) {
		NTL::swap(gen[i], gen[num - rank + i]);
	}
	gen.SetDims(rank, gen.NumCols());
}

//============================================================================
// Implementation

template<typename Int>
void BasisConstruction<Int>::LLLBasisConstruction(IntMat &gen, Int &m, double delta,
		PrecisionType prec) {
	std::cerr << "LLLBasisConstruction can only be done with NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

template<>
void BasisConstruction<NTL::ZZ>::LLLBasisConstruction(NTL::matrix<NTL::ZZ> &gen,
		NTL::ZZ &m, double delta, PrecisionType prec) {
	LLLConstruction0 (gen, delta, prec);
	int rank = gen.NumRows();
	int dim = gen.NumCols();
	if (rank == dim)
		return;
    gen.SetDims(rank + dim, dim);
    // We now add the m e_i vectors, and we redo the LLL.
    int i, j;
    for (i=rank; i < rank+dim; i++) {
        for (j=0; j < dim; j++) {
        	if (i==j) gen[i][j] = m; else gen[i][j] = 0;
        }
    }
	LLLConstruction0 (gen, delta, prec);
	std::cerr << "LLLBasisConstruction: we had to add some rows!\n";
}

//===================================================================

template<typename Int>
void BasisConstruction<Int>::upperTriangularBasis(IntMat &gen, IntMat &basis,
		Int &m) {
	std::cerr << "LLLBasisConstruction can only be done with NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

template<>
void BasisConstruction<NTL::ZZ>::upperTriangularBasis(NTL::matrix<NTL::ZZ> &gen, NTL::matrix<NTL::ZZ> &basis,
		NTL::ZZ &n) {
	
	//The following variables are used in the following
	IntVec coeff_gcd, coeff_xi, xi;
	NTL::ZZ gcd_tower, C, D;
	int dim1 = gen.NumRows();
	int dim2 = gen.NumCols();
	int k;
		
	//If necessary, the entries of the input matrix can be reduced. However, this costs some run time.
	//for (int k = 0; k < dim1; k++) {
	//	for (int j = 0; j < dim2; j++) {
	//		gen(k,j) = gen(k, j) % m; 
	//	}
	//}
	
	//Define dimensions of vectors
	coeff_gcd.SetLength(dim2);
	coeff_xi.SetLength(dim2);
	xi.SetLength(dim2);
	
	NTL_ZZNewRegister(gcd); 
	NTL_ZZNewRegister(m);
	
	m = n;
	
	for (int i = 0; i < dim2; i++) {
		
		//Delete values set during the last execution of the loop
		NTL::clear(xi);
		NTL::clear(coeff_gcd);
		
		//Search for the first non-zero element in the row		
		for (k = 0; k < dim1; k++)
		{
			if (gen[k][i] != 0)
				break;
		}
		
		//Reduce the other generators as they are used often in the following
		for (int j = k; j < dim1; j++) {
				gen[j][i] = gen[j][i] % m; 
		}
		
		//The else-case adds m e_i to the basis matrix
		if (k < dim1) {
			gcd = m; //GCD(m,gen[k][i]);
			coeff_gcd[k] = 1;
			gcd_tower = gcd;  			
			
			//Find the other coefficients by applying the Euclidean algorithm multiple times
			for (int j = k; j < dim1; j++) {	
				if (gen[j][i] == 0) 
					coeff_gcd[j] = 0;
				else {
					XGCD(gcd, C, D, gcd_tower, gen[j][i]);
					coeff_gcd[j] = D;
					for (int l = 0; l < j; l++) {
						mul(coeff_gcd[l], coeff_gcd[l], C);
					}
					gcd_tower = gcd;
				}		
			}		
			
			//Reduce the coefficients found during the Euclidean algorithm
			for (int j = 0; j < dim1; j++) {
				coeff_gcd[j] = coeff_gcd[j] % m;
			}
			
			//We have now found all the coefficients and can compute the vector x_i		
			for (int k = 0; k < dim1; k++) {
					if (coeff_gcd[k] != 0) {							
						for (int j = i; j < dim2; j++) {
							MulAddTo(xi[j], gen[k][j], coeff_gcd[k]);
						}
					}
			}
			
			//Next we calculate the new vectors v_i
			//For that purpose, we at first calculate the coefficients with which x_i needs to be multiplied
			for (int j = 0; j < dim1; j++) {
				div(coeff_xi[j], gen[j][i], gcd);
				coeff_xi[j] =  coeff_xi[j] % m;
				xi[j] = xi[j] % m;
			}
			
			//Update the v_i
			for (int k = 0; k < dim1; k++) {
				if (coeff_xi[k] != 0) {
					for (int j = i; j < dim2; j++) {
						MulSubFrom(gen[k][j], coeff_xi[k], xi[j]);
					}
				}
			}		
			
			//Set the base vectors	
			basis[i] = xi;
		} else {
			for (int j = 0; j < dim2; j++) {
				if (j != i)
					basis[i][j] = 0;
				else
					basis[i][j] = m;
			}
		}
	}
	//std::cout << basis;
}


//==============================================================================

template<typename Int>
void BasisConstruction<Int>::lowerTriangularBasis(IntMat &gen, IntMat &basis,
		Int &m) {
	IntVec coeff, vl, v2;
	Int C, D, val, gcd;
	int pc, pl, k;
	int dim1 = gen.NumRows();
	int dim2 = gen.NumCols();
	pl = dim1 - 1;
	pc = dim2 - 1;
	while (pl >= 0 && pc >= 0) {
		for (int i = 0; i < dim1; i++)
			Modulo(gen(i, pc), m, gen(i, pc));
		coeff.SetLength(dim2);
		k = 0;
		while (k < dim1 && gen(k, pc) == 0) {
			coeff[k] = 0;
			k++;
		}
		if (k < dim1) {
			gcd = gen(k, pc);
			coeff[k] = 1;
			val = gcd;
			for (int i = k + 1; i < dim1; i++) {
				if (gen(i, pc) == 0) {
					coeff[i] = 0;
					continue;
				}
				Euclide(val, gen(i, pc), C, D, gcd);
				coeff[i] = D;
				for (int j = 0; j < i; j++)
					coeff[j] *= C;
				val = gcd;
			}
			int coeffN[dim2];
			int nb = 0;
			for (int a = 0; a < dim1; a++) {
				if (coeff[a] != 0) {
					coeffN[nb] = a;
					nb++;
				}
			}
			vl.SetLength(dim2);
			int ind = 0;
			for (int j = 0; j < dim2; j++) {
				for (int i = 0; i < nb; i++) {
					ind = coeffN[i];
					vl[j] = vl[j] + coeff[ind] * gen(ind, j);
				}
				Modulo(vl[j], m, vl[j]);
			}
			for (int i = 0; i < dim1; i++) {
				if (gen(i, pc) != 0) {
					v2 = (gen(i, pc) / gcd) * vl;
					for (int j = 0; j < dim2; j++)
						Modulo(v2[j], m, v2[j]);
					for (int j = 0; j < dim2; j++) {
						gen(i, j) = gen(i, j) - v2[j];
						Modulo(gen(i, j), m, gen(i, j));
					}
				}
			}
			basis[pl] = vl;
		} else {
			for (int j1 = 0; j1 < dim2; j1++) {
				if (j1 != pl)
					NTL::clear(basis(pl, j1));
				else
					basis(pl, j1) = m;
			}
		}
		coeff.clear();
		vl.clear();
		pl--;
		pc--;
	}
}

//======================================================

template<typename Int>
void BasisConstruction<Int>::GCDTriangularBasis(IntMat &gen, Int &m) {
	// On exit, the rows of matrix are the basis vectors.
	long rows = gen.NumRows();
	long cols = gen.NumCols();
	long max_rank = rows < cols ? rows : cols;
	long rank = 0;
	// The basis will have at most max_rank vectors.
	Int q;
	for (long i = 0; i < max_rank; i++) {
		// We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
		// algorithm and applying transformations to the matrix
		for (long j = i + 1; j < rows; j++) {
			while (gen[j][i] != 0) {
				gen[i].swap(gen[j]);
				q = gen[j][i] / gen[i][i];
				gen[j] -= q * gen[i];
				for (int k = 0; k < max_rank; k++)
					Modulo(gen[j][k], m, gen[j][k]);
			}
		}
		if (gen[i][i] != 0) {
			rank++;
			if (gen[i][i] < 0)
				gen[i] *= Int(-1);
		}
	}
	// We remove zero vectors from the basis.
	gen.SetDims(rank, cols);
}

//======================================================

// This is the old version from Couture and L'Ecuyer (1996).
template<typename Int>
void BasisConstruction<Int>::mDualUpperTriangular96(IntMat &basis,
		IntMat &basisDual, Int &m) {
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);  return;
	}
	// We must have a triangular basis matrix in the first place.
	if (!CheckTriangular(basis, basis.NumRows(), Int(0)))
		GCDTriangularBasis(basis, m);
	long dim = basis.NumRows();
	if (dim != basis.NumCols()) {
		std::cout << ":mDualUpperTriangular96: basis matrix must be square.\n";
		return;
	}
	basisDual.SetDims(dim, dim);
	Int mm = m;            // Local copy of m that can be changed.
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(basisDual(i, j));
		if (!NTL::IsZero(basisDual(i, i))) {
			Int gcd = NTL::GCD(mm, basis(i, i));
			mm *= basis(i, i) / gcd;
			basisDual *= basis(i, i) / gcd;
		}

		DivideRound(mm, basis(i, i), basisDual(i, i));
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(basisDual(i, j));
			for (int k = j + 1; k <= i; k++)
				basisDual(i, j) += basis(j, k) * basisDual(i, k);
			if (basisDual(i, j) != 0)
				basisDual(i, j) = -basisDual(i, j);
			if (!NTL::IsZero(basisDual(i, j) % basis(j, j))) {
				Int gcd = NTL::GCD(basisDual(i, j), basis(j, j));
				mm *= basis(j, j) / gcd;
				basisDual *= basis(j, j) / gcd;
			}
			DivideRound(basisDual(i, j), basis(j, j), basisDual(i, j));
		}
	}
}

template<>
void BasisConstruction<NTL::ZZ>::mDualUpperTriangular96( NTL::matrix<NTL::ZZ>  &basis,
		NTL::matrix<NTL::ZZ>  &basisDual, NTL::ZZ &m) {
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);  return;
	}
	
	NTL::ZZ gcd, fac;
	NTL::ZZ mm = m;            // Local copy of m that can be changed.
	
	int dim;
	
	// We must have a triangular basis matrix in the first place.
	//if (!CheckTriangular(basis, basis.NumRows(), 0))
	//	GCDTriangularBasis(basis, m);
	//long dim = basis.NumRows();
	//if (dim != basis.NumCols()) {
	//	std::cout << ":mDualUpperTriangular96: basis matrix must be square.\n";
	//	return;
	//}
	
	dim = basis.NumRows();
	
	basisDual.SetDims(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(basisDual[i][j]);
		if (!NTL::IsZero(basisDual[i][i])) {
			gcd = NTL::GCD(mm, basis[i][i]);
			div(fac, basis[i][i], gcd); 
			mul(mm, mm, fac); 
			mul(basisDual[i][i], basisDual[i][i], fac);
		}

		DivideRound(mm, basis[i][i], basisDual[i][i]);
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(basisDual(i, j));
			for (int k = j + 1; k <= i; k++)
				basisDual[i][j] += basis[j][k] * basisDual[i][k];
			//if (basisDual(i, j) != 0)
				negate(basisDual[i][j], basisDual[i][j]);
			if (!NTL::IsZero(basisDual[i][j] % basis[j][j])) {
				gcd = NTL::GCD(basisDual(i, j), basis[j][j]);
				div(fac, basis[j][j], gcd); 
				mul(mm, mm, fac); 
				mul(basisDual[i][j], basisDual[i][j], fac);
			}
			div(basisDual[i][j], basisDual[i][j], basis[j][j]);
		}
	}
	
	//std::cout << basisDual;
}


//===================================================

/**
 * For `B` to be `m`-dual to `A`, we have to have that \f$AB^t = mI\f$.
 * Since `A` is upper triangular, `B` will be a lower triangular matrix
 * with `A(i,i)*B(i,i) = m` for all `i` and
 * \f$ A_i \cdot B_j = 0\f$ for \f$i\neq j\f$. To get the second condition,
 * we simply have to recursively take for each line
 * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
 */
template<typename Int>
void BasisConstruction<Int>::mDualUpperTriangular(const IntMat &A, IntMat &B,
		const Int &m) {
	int dim = A.NumRows();
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(B(i, j));
		DivideRound(m, A(i, i), B(i, i));
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(B(i, j));
			for (int k = j + 1; k <= i; k++)
				B(i, j) += A(j, k) * B(i, k);
			if (B(i, j) < 0)
				B(i, j) = -B(i, j);
			DivideRound(B(i, j), A(j, j), B(i, j));
		}
	}
}

template<>
void BasisConstruction<NTL::ZZ>::mDualUpperTriangular(const NTL::matrix<NTL::ZZ> &A, NTL::matrix<NTL::ZZ> &B,
		const NTL::ZZ &m) {
	int dim = A.NumRows();	
	
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(B[i][j]);
		div(B[i][i], m, A[i][i]); 
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(B[i][j]);
			for (int k = j + 1; k <= i; k++)
				MulSubFrom(B[i][j], A[j][k], B[i][k]);
			div(B[i][j], B[i][j], A[j][j]);		
		}
	}
	
	//std::cout << B;
}

template<typename Int>
void BasisConstruction<Int>::mDualBasis(IntMat &basis, IntMat &basisDual,
		Int &m) {
	std::cerr << "mDualBasis is implemented only for NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

// The specialization for the case where `Int = ZZ`.
template<>
void BasisConstruction<NTL::ZZ>::mDualBasis(
		NTL::matrix<NTL::ZZ> &basis, NTL::matrix<NTL::ZZ> &basisDual, NTL::ZZ &m) {
	NTL::ZZ d;
	NTL::Mat<NTL::ZZ> C;
	int dim = basis.NumRows();
	if (dim != basis.NumCols()) {
		std::cerr << "mDualBasis: the given basis matrix must be square.\n";
		exit(1);
	}
	C.SetDims(dim, dim);
	inv(d, basisDual, basis);
	NTL::ZZ m2 = m / d;
	transpose(C, basis);
	for (int i = 1; i < dim; i++) {
		for (int j = 1; j < dim; j++)
			basis(i, j) = (m2 * C(i, j));
	}
}

//=================================================================================

template<typename Int>
template<typename Real>
void BasisConstruction<Int>::projectionConstructionLLL(
		IntLattice<Int, Real> &in, IntLattice<Int, Real> &out,
		const Coordinates &proj) {
	std::size_t size = proj.size();
	unsigned int lat_dim = in.getDim();
	if (size > lat_dim)
		MyExit(1, "More projection coordinates than the dimension of `in`.");
	IntMat new_basis, tmp(NTL::transpose(in.getBasis()));
	new_basis.SetDims(size, tmp.NumRows());
	tmp = NTL::transpose(tmp);
	auto it = proj.cbegin();
	for (std::size_t i = 0; i < size; i++) {
		if (*it <= lat_dim)
			new_basis[i] = tmp[*it];
		else
			MyExit(1, "Coordinates do not match the dimensions of `in`.");
		it++;
	}
	new_basis = NTL::transpose(new_basis);
	LLLConstruction(new_basis);
	out = IntLattice<Int, Real>(new_basis, size, in.getNormType());
}

template class BasisConstruction<std::int64_t>;
template class BasisConstruction<NTL::ZZ>;
// extern template class BasisConstruction<std::int64_t> ;
// extern template class BasisConstruction<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
