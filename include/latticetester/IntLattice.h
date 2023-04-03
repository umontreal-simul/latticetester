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

#ifndef LATTICETESTER_INTLATTICE_H
#define LATTICETESTER_INTLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NTLWrap.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

namespace LatticeTester {

/**
 * An `IntLattice` object is an integral lattice, with its basis and `m`-dual basis
 * (the latter is optional). There are tools to perform simple manipulations on those lattice bases.
 * The value of `m` is always chosen in a way that all coordinates of the basis and
 * of its `m`-dual are integers, so they can be represented exactly.
 * The basis or its dual is rescaled by `m`, which is typically the smallest integer with this property.
 *
 * The dimension $t$ of the lattice is the number of independent vectors that form a basis.
 * Usually, these vectors also have $t$ coordinates, but in general they may have more.
 * A norm is also chosen in `NormType` to measure the vector lengths; by default it is the
 * Euclidean norm.
 * Methods and attributes are offered to compute and store the norms of the basis and dual basis vectors,
 * to permute basis vectors, sort them by length and do the corresponding changes in the dual, etc.
 * An `IntLattice` object contains several protected variables to store all these quantities.
 * For better efficiency, we should avoid creating too many of these objects, for example when
 * making searches for good lattices.
 *
 * The class `IntLatticeExt` extends this class and contains virtual methods that must
 * be defined in its subclasses.
 */
template<typename Int, typename Real>
class IntLattice {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;
	typedef NTL::matrix<Real> RealMat;

public:

	/**
	 * Constructs a lattice whose basis is the identity, in `dim` dimensions,
	 * with the specified norm type, and the scaling factor `m` and dual basis undefined.
	 */
	IntLattice(const int dim, NormType norm = L2NORM);

	/**
	 * Constructs a lattice with the given basis, in `dim` dimensions,
	 * and with the specified norm type. The dual basis and `m` are not initialized.
	 * The `basis` matrix must be a `dim` by `dim` square integer matrix.
	 */
	IntLattice(const IntMat basis, const int dim, NormType norm = L2NORM);

	/**
	 * Constructs a lattice with the given basis and given m-dual basis for the given `m`,
	 * in `dim` dimensions, and with the specified norm type.
	 */
	IntLattice(const IntMat primalbasis, const IntMat dualbasis,
			const Int m, const int dim, NormType norm = L2NORM);

	/**
	 * Copy constructor. Makes a deep copy of `lat` into `*this`.
	 */
	IntLattice(const IntLattice<Int, Real> &lat);

	/**
	 * Destructor.
	 */
	~IntLattice();

	/**
	 * Cleans and releases all the memory allocated to this lattice.
	 * **Deprecated: REMOVE!**
	 */
	// void kill();

	/**
	 * Makes a full deep copy of the lattice `lat` into this object.
	 * New matrix and vector objects are constructed to store the bases and norms.
	 */
	void copyLattice(const IntLattice<Int, Real> &lat);

	/*
	 * Previously named `copyLattice`.
	 * Overwrites the first `n` elements of the basis of the lattice `lat` over the elements
	 * of the basis of the current object. The latter must have dimension `n` already,
	 * otherwise an error message is printed and nothing also is done!
	 * The vector norms and the dual basis (if available) are also overwritten.
	 * The difference with `copyLattice` is that here, no new matrix or vector is constructed;
	 * the previous ones are re-used.
	 */
	 void overwriteLattice(const IntLattice<Int, Real> &lat, long n);
	 
	/**
	 * Initializes a vector containing the norms of the basis vectors to -1
	 * for all components.  It means the norms are no longer up to date.
	 */
	void initVecNorm();

	/**
	 * Returns the basis represented in a matrix.
	 */
	IntMat& getBasis() {
		return m_basis;
	}

	/**
	 * Returns the m-dual basis represented in a matrix.
	 */
	IntMat& getDualBasis() {
		return m_dualbasis;
	}

	/**
	 * Returns the dimension of the lattice, which is the dimension of the basis vectors,
	 * and also usually the number of independent vectors in the basis.
	 */
	int getDim() const {
		return m_dim;
	}

	/**
	 * Returns the `NormType` used by this lattice.
	 */
	NormType getNormType() const {
		return m_norm;
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of the i-th vector of the basis,
	 * with the index i starting at 0.
	 */
	Real getVecNorm(const int &i) {
		return m_vecNorm[i];
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of each basis vector, in a vector.
	 */
	RealVec getVecNorm() const {
		return m_vecNorm;
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of the i-th vector of the m-dual basis.
	 */
	Real getDualVecNorm(const int &i) {
		return m_dualvecNorm[i];
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of each vector of the m-dual basis, in a vector.
	 */
	RealVec getDualVecNorm() const {
		return m_dualvecNorm;
	}

	/**
	 * Returns the scaling factor `m`, which is 1 by default.
	 */
	Int getModulo() const {
		return m_modulo;
	}

	/**
	 * Sets the dimension of the basis to `dim`. This does not change any of the
	 * basis vectors, but only the dimension variable.
	 * Warning: After calling this method, the size of the basis matrix may no longer agree with the dimension.
	 */
	void setDim(const int &dim) {
		if (dim > 0)
			m_dim = dim;
	}

	/**
	 * Sets the `NormType` used by this lattice to `norm`.
	 */
	void setNormType(const NormType &norm) {
		m_norm = norm;
	}

	/**
	 * Sets the norm of the `i`-th component of the basis to `value`, which is assumed to
	 * be the correct value.  To recompute the norm, use `updateVecNorm(const int&)` instead.
	 */
	void setVecNorm(const Real &value, const int &i) {
		m_vecNorm[i] = value;
	}

	/**
	 * Sets the norm of the `i`-th component of the m-dual basis to `value`,  which is assumed to
	 * be the correct value.  To recompute the norm, use `updateDualVecNorm(const int&)` instead.
	 */
	void setDualVecNorm(const Real &value, const int &i) {
		m_dualvecNorm[i] = value;
	}

	/**
	 * Returns `true` iff an m-dual basis is available.
	 */
	bool withDual() const {
		return m_withDual;
	}

	/**
	 * Sets the `withDual` flag to `flag`. This flag indicates whether or
	 * not this IntLattice contains an up-to-date m-dual basis. It is the flag
	 * returned by `withDual()`.
	 */
	void setDualFlag(bool flag) {
		m_withDual = flag;
	}

	/**
	 * Sets all the values in the array containing the norms of the basis vectors to -1.
	 * This means that these norms are no longer up to date.
	 */
	void setNegativeNorm();

	/**
	 * Sets the value of the `i`-th component in the array containing the
	 * norms of the basis vectors to -1.
	 */
	void setNegativeNorm(const int &i) {
		m_vecNorm[i] = -1;
	}

	/**
	 * Sets all the values in the array containing the norms of the dual basis
	 * vectors to -1, to indicate that these norms are no longer up to date.
	 */
	void setDualNegativeNorm();

	/**
	 * Sets the value of the `i`-th component in the array containing the
	 * norms of the m-dual basis vectors to -1.
	 */
	void setDualNegativeNorm(const int &i) {
		m_dualvecNorm[i] = -1;
	}

	/**
	 * Updates the array containing the basis vectors norms by recomputing them.
	 */
	void updateVecNorm();

	/**
	 * Updates the array containing the basis vectors norms from the `d`-th
	 * component to the last, by recomputing them.
	 * Putting `d=0` recomputes all the norms.
	 */
	void updateVecNorm(const int &d);

	/**
	 * Updates the array containing the m-dual basis vectors norms by recomputing them.
	 * Assumes that the dual basis is available.
	 */
	void updateDualVecNorm();

	/**
	 * Updates the array containing the m-dual basis vectors norms from the `d`-th
	 * component to the last by recomputing them.
	 * */
	void updateDualVecNorm(const int &d);

	/**
	 * Updates the `i`-th value of the array containing the square norms of the
	 * basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateScalL2Norm(const int i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the square norms of the basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateScalL2Norm(const int k1, const int k2);

	/**
	 * Updates the `i`-th value of the array containing the square norms of the
	 * m-dual basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the square norms of the m-dual basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int k1, const int k2);

	/**
	 * Exchanges vectors `i` and `j` in the basis. This also changes the
	 * m-dual basis vectors and the arrays containing secondary information
	 * about the two basis (like the norms) accordingly.
	 */
	void permute(int i, int j);

	/**
	 * Exchanges vectors `i` and `j` in the basis without changing the m-dual.
	 * See `permute()`.
	 */
	void permuteNoDual(int i, int j);

    /**
     * Exchange the primal and m-dual bases.
     * If the dual is not defined, exits with an error message.
     */
    void dualize ();

    /**
	 * Returns `true` iff the m-dual basis contained in the object really is
	 * the m-dual of the basis.  This also returns false
	 * if no dual has been specified.
	 */
	bool checkDuality();

	/**
	 * Sorts the basis vectors with indices greater of equal to `d` by
	 * increasing length. The m-dual vectors are permuted accordingly. Assumes
	 * that the lengths (norms) of the corresponding basis vectors are up to date.
	 */
	void sort(int d);

	/**
	 * Sorts the basis vectors with indices greater of equal to `d` by
	 * increasing length. The m-dual vectors are **not** permuted. See `sort()`.
	 */
	void sortNoDual(int d);

	/**
	 * Returns a string that contains the primal basis vectors and their norms.
	 */
	std::string toStringBasis() const;

	/**
	 * Returns a string with the m-dual basis vectors and their norms.
	 */
	std::string toStringDualBasis() const;

	/**
	 * Returns a string that represents the lattice and its parameters.
	 * It contains the dimension, the norm used, the basis and m-dual basis vectors and
	 * the basis and dual basis vector norms.
	 */
	std::string toString() const;

	/**
	 * Writes on standard output the string returned by `toString`.
	 */
	void write() const;

protected:

	/**
	 * The rows of this matrix are the primal basis vectors.
	 */
	IntMat m_basis;

	/**
	 * The rows of this matrix are the m-dual basis vectors.  May not be initialized.
	 * When m_withDual = true, it must be initialized.
	 */
	IntMat m_dualbasis;

	/**
	 * The dimension of the lattice, which is the number of (independent) vectors
	 * in the basis. It cannot exceed the number of coordinates in those vectors.
	 */
	int m_dim;

	/**
	 * The NormType used to measure the vector lengths for this lattice.
	 * It is used for the basis reduction and compute a shortest vector, for example.
	 */
	NormType m_norm;

	/**
	 * A vector that stores the norm of each basis vector.
	 * In case of the L_2 norm, it contains the square norm instead.
	 * A value of -1 means that the norm is not up to date.
	 */
	RealVec m_vecNorm;

	/**
	 * Similar to vecNorm, but for the m-dual basis.
	 */
	RealVec m_dualvecNorm;

    /**
	 * The scaling factor `m` used for rescaling the lattice. It is 0 when undefined.
	 */
	Int m_modulo;

	/**
	 * This `m_withDual` variable is `true` iff an m-dual basis is available.
	 */
	bool m_withDual;
};

// class IntLattice

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const int dim, NormType norm) :
		m_dim(dim), m_norm(norm), m_modulo(0), m_withDual(false) {
	this->m_basis.resize(dim, dim);
	this->m_vecNorm.resize(dim);
	initVecNorm();
}

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntMat basis,
		const int dim, NormType norm) :
		m_basis(basis), m_dim(dim), m_norm(norm),  m_modulo(0), m_withDual(false) {
	this->m_vecNorm.resize(dim);
	initVecNorm();
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntMat primalbasis,
		const IntMat dualbasis, const Int m, const int dim, NormType norm) :
		IntLattice<Int, Real>(primalbasis, dim, norm) {
	this->m_dualbasis = IntMat(dualbasis);
	this->m_dualvecNorm.resize(dim);
	this->m_modulo = m;
	this->m_withDual = true;
	setDualNegativeNorm();
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(
		const IntLattice<Int, Real> &lat) {
	copyLattice(lat);
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::~IntLattice() {
	// kill();
	this->m_basis.IntMat::clear();
	this->m_dualbasis.IntMat::clear();
	this->m_vecNorm.clear();
	this->m_dualvecNorm.clear();
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::copyLattice(
		const IntLattice<Int, Real> &lat) {
	this->m_dim = lat.m_dim;
	this->m_basis = IntMat(lat.m_basis);
	this->m_dualbasis = IntMat(lat.m_dualbasis);
	this->m_norm = lat.m_norm;
	this->m_vecNorm = RealVec(lat.m_vecNorm);
	this->m_dualvecNorm = RealVec(lat.m_dualvecNorm);
	this->m_modulo = lat.m_modulo;
	this->m_withDual = lat.m_withDual;
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::overwriteLattice(
		const IntLattice<Int, Real> &lat, long n) {
	if (this->m_dim == n) {
		CopyMatr(this->m_basis, lat.m_basis, n);
		CopyVect(this->m_vecNorm, lat.m_vecNorm, n);
		this->m_withDual = lat.m_withDual;
		if (this->m_withDual) {
			this->m_dualbasis.resize(this->m_basis.size1(),
					this->m_basis.size1());
			this->m_dualvecNorm.resize(this->m_basis.size1());
			CopyMatr(this->m_dualbasis, lat.m_dualbasis, n);
			CopyVect(this->m_dualvecNorm, lat.m_dualvecNorm, n);
		}
		this->m_modulo = lat.m_modulo;
	}
	else
		std::cout << "Calling IntLattice::overwriteLattice with wrong dimension"
				<< std::endl;
	}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::initVecNorm() {
	for (int i = 0; i < this->m_dim; i++) {
		this->m_vecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::setNegativeNorm() {
	for (int i = 0; i < this->m_dim; i++) {
		this->m_vecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::setDualNegativeNorm() {
	for (int i = 0; i < this->m_dim; i++) {
		this->m_dualvecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateVecNorm() {
	updateVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateVecNorm(const int &d) {
	assert(d >= 0);
	for (int i = d; i < this->m_dim; i++) {
		NTL::matrix_row<IntMat> row(this->m_basis, i);
		if (this->m_norm == L2NORM) {
			ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
		} else {
			CalcNorm<IntVec, Real>(row, this->m_dim, this->m_vecNorm[i],
					this->m_norm);
		}
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm() {
	updateDualVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm(const int &d) {
	assert(d >= 0);
	for (int i = d; i < this->m_dim; i++) {
		NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
		if (this->m_norm == L2NORM) {
			ProdScal<Int>(row, row, this->m_dim, this->m_dualvecNorm[i]);
		} else {
			CalcNorm<IntVec, Real>(row, this->m_dim, this->m_dualvecNorm[i],
					this->m_norm);
		}
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int i) {
	NTL::matrix_row<IntMat> row(this->m_basis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int k1,
		const int k2) {
	for (int i = k1; i < k2; i++) {
		updateScalL2Norm(i);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int i) {
	NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_dualvecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int k1,
		const int k2) {
	for (int i = k1; i < k2; i++) {
		updateDualScalL2Norm(i);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permute(int i, int j) {
	if (i == j)
		return;
	for (int k = 0; k < this->m_dim; k++) {
		swap9(this->m_basis(j, k), this->m_basis(i, k));
		if (this->m_withDual) {
			swap9(this->m_dualbasis(j, k), this->m_dualbasis(i, k));
		}
	}
	swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
	if (this->m_withDual) {
		swap9(this->m_dualvecNorm[i], this->m_dualvecNorm[j]);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permuteNoDual(int i, int j) {
	if (i == j)
		return;
	for (int k = 0; k < this->m_dim; k++) {
		swap9(this->m_basis(j, k), this->m_basis(i, k));
	}
	swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
}

//===========================================================================

template<typename Int, typename Real>
    void IntLattice<Int, Real>::dualize () {
    if(!(this->m_withDual)) {
	    std::cout << "\n***** ERROR: calling dualize while dual basis is not defined."
          << std::endl;
	    return;
		}
    std::swap(this->m_basis, this->m_dualbasis);
    this->setNegativeNorm ();
    this->setDualNegativeNorm ();
  }

/*=========================================================================*/

template<typename Int, typename Real>
bool IntLattice<Int, Real>::checkDuality() {
	if (!this->m_withDual) {
		std::cout << "Calling IntLattice::checkDuality with undefined m-dual"
				<< std::endl;
		return false;
	}
	Int S;
	int dim = getDim();

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			NTL::matrix_row<const IntMat> row1(this->m_basis, i);
			NTL::matrix_row<const IntMat> row2(this->m_dualbasis, j);
			ProdScal<Int>(row1, row2, dim, S);
			if (j != i) {
				if (S != 0) {
					std::cout << "******  checkDuality failed for V[" << i
							<< "] and W[" << j << "]" << std::endl;
					return false;
				}
			} else if (S != this->m_modulo) {
				std::cout << "******  checkDuality failed for i, j = " << i
						<< " , " << j << std::endl;
				return false;
			}
		}
	}
	return true;

}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::sort(int d)
/*
 * We assume that the (square) lengths are already updated.
 * This gives flexibility to the user to put something else than
 * the square Euclidean length in vecNorm.
 */
{
	int dim = getDim();
	for (int i = 0; i < dim; i++) {
		if (getVecNorm(i) < 0) {
			std::cout << "\n***** ERROR: in sort, Negative norm for i = " << i
					<< ",  dim = " << dim << std::endl;
		}
	}

	for (int i = d; i < dim; i++) {
		int k = i;
		for (int j = i + 1; j < dim; j++) {
			if (getVecNorm(j) < getVecNorm(k))
				k = j;
		}
		if (i != k)
			permute(i, k);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::sortNoDual(int d)
/*
 * We assume that the (square) lengths are already updated.
 * This gives flexibility to the user to use something else than
 * the square Euclidean length.
 */
{
	int dim = getDim();
	for (int i = 0; i < dim; i++) {
		if (getVecNorm(i) < 0) {
			std::cout << "\n***** ERROR: sort   Negative norm for i = " << i
					<< ",  dim = " << dim << std::endl;
		}
	}
	for (int i = d; i < dim; i++) {
		int k = i;
		for (int j = i + 1; j < dim; j++) {
			if (getVecNorm(j) < getVecNorm(k))
				k = j;
		}
		if (i != k)
			permuteNoDual(i, k);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toString() const {
	std::ostringstream os;
	os << "Dim = " << this->m_dim << " \n \n";
	os << std::setprecision(10) << "Primal basis vectors:\n";
	for (int i = 0; i < this->m_dim; i++) {
		os << this->m_basis[i];
		//for (int j = 0; j < this->m_dim; j++) {
		//  os <<  this->m_basis(i,j);
		//}
		os << "\n";
	}
	os << "\nm-Dual basis vectors:\n";
	for (int i = 0; i < this->m_dim; i++) {
		if (this->m_withDual) {
			os << this->m_dualbasis[i];
			//for (int j = 0; j < this->m_dim; j++) {
			//  os << this->m_dualbasis(i,j);
			//}
			os << "\n";
		}
	}
	os << "\n";
	os << "Norm used: " << toStringNorm(this->m_norm) << "\n"
			<< std::endl;
	os << "Norm of each Basis vector: \n";
	os << "Primal";
	if (this->m_withDual)
		os << "\t\tDual\n";
	os << "\n";

	for (int i = 0; i < this->m_dim; i++) {
		if (this->m_vecNorm[i] < 0) {
			os << "NaN OR Not computed";
		} else {
			if (this->m_norm == L2NORM) {
				os << NTL::sqrt(this->m_vecNorm[i]);
			} else {
				os << this->m_vecNorm[i];
			}
		}
		os << "\t";
		if (this->m_withDual) {
			if (this->m_dualvecNorm[i] < 0)
				os << "NaN OR Not computed";
			else {
				if (this->m_norm == L2NORM) {
					os << NTL::sqrt(this->m_dualvecNorm[i]);
				} else {
					os << this->m_dualvecNorm[i];
				}
			}
		}
		os << "\n";
	}
	os << std::endl;
	return os.str();
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::write() const {
	std::cout << this->toString() << "\n";
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toStringBasis() const {
	std::ostringstream os;
	os << "Primal Basis:\n";
	os << "  Dim = " << this->m_dim << " \n";
	for (int i = 0; i < this->m_dim; i++) {
		os << "    [";
		for (int j = 0; j < this->m_dim; j++)
			os << " " << std::setprecision(15) << this->m_basis(i, j);
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
	return os.str();
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toStringDualBasis() const {
	std::ostringstream os;
	os << "m-Dual Basis:\n";
	os << "  Dim = " << this->m_dim << " \n";
	for (int i = 0; i < this->m_dim; i++) {
		os << "    [";
		for (int j = 0; j < this->m_dim; j++)
			os << " " << std::setprecision(15) << this->m_dualbasis(i, j);
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
	return os.str();
}

extern template class IntLattice<std::int64_t, double> ;
extern template class IntLattice<NTL::ZZ, double> ;
extern template class IntLattice<NTL::ZZ, NTL::RR> ;

} // namespace LatticeTester

#endif
