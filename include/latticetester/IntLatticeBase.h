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
#include "latticetester/NTLWrap.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

namespace LatticeTester {

/**
 * This class acts as a base class for `IntLattice`. It contains a subset of the methods
 * and offers constructors to build a lattice from an arbitrary basis.
 * The class `IntLattice` extends this class and contains virtual methods that must
 * be defined in its subclasses.
 * An object of this class is an integral lattice, with its basis and `m`-dual basis.
 * There are tools to perform simple manipulations on those lattice bases.
 * The value of `m` is always chosen in a way that all coordinates of the basis and
 * of its `m`-dual are integers, so they can be represented exactly.
 * The basis or its dual is rescaled by `m`.
 * Typically, `m` is the smallest integer with this property.
 *
 * The dimension $t$ of the lattice is the number of independent vectors that form a basis.
 * Usually, these vectors also have $t$ coordinates.
 * A norm is also chosen in `NormType` to measure the vector lengths; by default it is the
 * Euclidean norm.
 * Methods and attributes are offered to compute and store the norms of the basis and dual basis vectors,
 * to permute basis vectors, sort them by length and do the corresponding changes in the dual, etc.
 *
 * Note that in this class, the indices of basis vectors and coordinates start at 1.
 */
template<typename Int, typename Real, typename RealRed>
class IntLatticeBase {

private:
	// Forward definition of types to be used in this class.
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;

public:

	/**
	 * Constructs a lattice whose basis is the identity, in `dim` dimensions,
	 * with `m=1`, and with the specified norm type.
	 */
	IntLatticeBase(const int dim, NormType norm = L2NORM);

	/**
	 * Constructs a lattice with the given basis, in `dim` dimensions,
	 * and with the specified norm type. The dual basis and `m` are not initialized.
	 * The `basis` matrix must be a `dim` by `dim` square matrix.
	 */
	IntLatticeBase(const IntMat basis, const int dim, NormType norm = L2NORM);

	/**
	 * Constructs a lattice with the given basis and given m-dual basis for the given `m`,
	 * in `dim` dimensions, and with the specified norm type.
	 */
	IntLatticeBase(const IntMat primalbasis, const IntMat dualbasis,
			const Int m, const int dim, NormType norm = L2NORM);

	/**
	 * Copy constructor. Makes a deep copy of Lat into `*this`.
	 */
	IntLatticeBase(const IntLatticeBase<Int, Real, RealRed> &Lat);

	/**
	 * Destructor.
	 */
	~IntLatticeBase();

	/**
	 * Cleans and releases all the memory allocated to this lattice.
	 */
	void kill();

	/**
	 * Copy the lattice `lat`, except it's NormType and dimension, into this
	 * object. This does not check if the dimensions match.
	 */
	void copyLattice(const IntLatticeBase<Int, Real, RealRed> &lat);

	/**
	 * Copy the n first elements of the basis of the lattice `lat` into this
	 * object. The object into which `lat` is copied has to be of dimension n.
	 */
	void copyLattice(const IntLatticeBase<Int, Real, RealRed> &lat, long n);

	/**
	 * Initializes a vector containing the norms of the basis vectors to -1
	 * at all components.
	 */
	void initVecNorm();

	/**
	 * Returns the basis represented in a matrix.
	 */
	IntMat& getBasis() {
		return m_basis;
	}

	/**
	 * Returns the dual basis represented in a matrix.
	 */
	IntMat& getDualBasis() {
		return m_dualbasis;
	}

	/**
	 * Returns the dimension of the lattice. (Both the number of vectors in
	 * the basis and the number of coordinates of those vectors.
	 */
	int getDim() const {
		return m_dim;
	}

	/**
	 * Returns the NormType used by this lattice.
	 */
	NormType getNorm() const {
		return m_norm;
	}

	/**
	 * Returns the norm of the i-th vector of the basis.
	 */
	Real getVecNorm(const int &i) {
		return m_vecNorm[i];
	}

	/**
	 * Returns the norm of each vector of the basis in a vector.
	 */
	DblVec getVecNorm() const {
		return m_vecNorm;
	}

	/**
	 * Returns the norm of the i-th vector of the dual basis.
	 */
	Real getDualVecNorm(const int &i) {
		return m_dualvecNorm[i];
	}

	/**
	 * Returns the norm of each vector of the dual basis in a vector.
	 */
	DblVec getDualVecNorm() const {
		return m_dualvecNorm;
	}

	/**
	 * Returns the `m` used for rescaling if it has been defined. Returns `0`
	 * otherwise.
	 */
	Int getModulo() const {
		return m_modulo;
	}

	/**
	 * Sets the dimension of the basis to `d`. This won't change any of the
	 * vectors of the basis by itself. This method should not be called
	 * directly on an object of this class except in a function specifically
	 * changing the dimension of this object.
	 */
	void setDim(const int &d) {
		if (d > 0)
			m_dim = d;
	}

	/**
	 * Sets the NormType used by this lattice to `norm`.
	 */
	void setNorm(const NormType &norm) {
		m_norm = norm;
	}

	/**
	 * Sets the norm of the `i`-th component of the basis to `value`. The
	 * usage of `updateVecNorm(const int&)` is recommended over this function.
	 */
	void setVecNorm(const Real &value, const int &i) {
		m_vecNorm[i] = value;
	}

	/**
	 * Sets the norm of the `i`-th component of the dual basis to `value`.
	 * The usage of `updateDualVecNorm(const int&)` is recommended over this function.
	 */
	void setDualVecNorm(const Real &value, const int &i) {
		m_dualvecNorm[i] = value;
	}

	/**
	 * Returns `true` if a dual has been defined and `false` otherwise.
	 */
	bool withDual() const {
		return m_withDual;
	}

	/**
	 * Sets the `withDual` flag to `flag`. This flag indicates whether or
	 * not this IntLatticeBase contains a dual basis. It is the flag
	 * returned by `withDual()`.
	 */
	void setDualFlag(bool flag) {
		m_withDual = flag;
	}

	/**
	 * Sets all the values in the array containing the norms of the basis
	 * vectors to -1.
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
	 * vectors to -1.
	 */
	void setDualNegativeNorm();

	/**
	 * Sets the value of the `i`-th component in the array containing the
	 * norms of the dual basis vectors to -1.
	 */
	void setDualNegativeNorm(const int &i) {
		m_dualvecNorm[i] = -1;
	}

	/**
	 * Updates the array containing the basis vectors norms by recomputing
	 * them.
	 * */
	void updateVecNorm();

	/**
	 * Updates the array containing the basis vectors norms from the `d`-th
	 * component to the last by recomputing them.
	 * */
	void updateVecNorm(const int &d);

	/**
	 * Updates the array containing the dual basis vectors norms by recomputing
	 * them.
	 * */
	void updateDualVecNorm();

	/**
	 * Updates the array containing the dual basis vectors norms from the `d`-th
	 * component to the last by recomputing them.
	 * */
	void updateDualVecNorm(const int &d);

	/**
	 * Updates the `i`-th value of the array containing the norms of the
	 * basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateScalL2Norm(const int i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the norms of the basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateScalL2Norm(const int k1, const int k2);

	/**
	 * Updates the `i`-th value of the array containing the norms of the
	 * dual basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the norms of the dual basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int k1, const int k2);

	/**
	 * Exchanges vectors `i` and `j` in the basis. This also changes the
	 * dual basis vectors and the arrays containing secondary information
	 * about the two basis (like the norms) accordingly.
	 */
	void permute(int i, int j);

	/**
	 * Exchanges vectors `i` and `j` in the basis without changing the dual.
	 * See `permute()`.
	 */
	void permuteNoDual(int i, int j);

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
	void sort(int d);

	/**
	 * Sorts the basis vectors with indices greater of equal to \f$d\f$ by
	 * increasing length. The dual vectors are **not** permuted. See `sort()`.
	 */
	void sortNoDual(int d);

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
	void write() const;

protected:

	/**
	 * Each row of this matrix represents a vector in the basis of the
	 * lattice.
	 */
	IntMat m_basis;

	/**
	 * Each row of this matrix represents a vector in the dual basis of the
	 * lattice.
	 */
	IntMat m_dualbasis;

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
	Int m_modulo;

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
};
// class IntLatticeBase

//===========================================================================

template<typename Int, typename Real, typename RealRed>
IntLatticeBase<Int, Real, RealRed>::IntLatticeBase(const int dim, NormType norm) :
		m_dim(dim), m_norm(norm), m_modulo(0), m_withDual(false)
//m_xx(0)

{
	this->m_basis.resize(dim, dim);
	this->m_vecNorm.resize(dim);
	initVecNorm();
}

//===========================================================================

template<typename Int, typename Real, typename RealRed>
IntLatticeBase<Int, Real, RealRed>::IntLatticeBase(const IntMat basis,
		const int dim, NormType norm) :
		m_basis(basis), m_dim(dim), m_norm(norm), m_modulo(0), m_withDual(false)
//m_xx(0)
{
	this->m_vecNorm.resize(dim);
	initVecNorm();
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
IntLatticeBase<Int, Real, RealRed>::IntLatticeBase(const IntMat primalbasis,
		const IntMat dualbasis, const Int m, const int dim, NormType norm) :
		IntLatticeBase<Int, Real, RealRed>(primalbasis, dim, norm) {
	this->m_dualbasis = IntMat(dualbasis);
	this->m_withDual = true;
	this->m_dualvecNorm.resize(dim);
	setDualNegativeNorm();
	this->m_modulo = m;
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
IntLatticeBase<Int, Real, RealRed>::IntLatticeBase(
		const IntLatticeBase<Int, Real, RealRed> &lat) :
		m_dim(lat.getDim()), m_norm(lat.getNorm())
//m_xx(0)
{
	copyLattice(lat);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
IntLatticeBase<Int, Real, RealRed>::~IntLatticeBase() {
	kill();
	this->m_basis.IntMat::clear();
	this->m_dualbasis.IntMat::clear();
	this->m_vecNorm.clear();
	this->m_dualvecNorm.clear();
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::kill() {
	//delete [] this->m_xx;
	//this->m_xx = 0;
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::copyLattice(
		const IntLatticeBase<Int, Real, RealRed> &lat) {
	//if(m_dim == lat.m_dim)

	this->m_basis = IntMat(lat.m_basis);
	this->m_dualbasis = IntMat(lat.m_dualbasis);
	this->m_vecNorm = DblVec(lat.m_vecNorm);
	this->m_dualvecNorm = DblVec(lat.m_dualvecNorm);
	this->m_withDual = lat.m_withDual;
	this->m_modulo = lat.m_modulo;
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::copyLattice(
		const IntLatticeBase<Int, Real, RealRed> &lat, long n) {
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
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::initVecNorm() {
	//this->m_xx = new bool[this->m_dim];
	for (int i = 0; i < this->m_dim; i++) {
		this->m_vecNorm[i] = -1;
		//this->m_xx[i] = true;
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::setNegativeNorm() {
	for (int i = 0; i < this->m_dim; i++) {
		this->m_vecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::setDualNegativeNorm() {
	for (int i = 0; i < this->m_dim; i++) {
		this->m_dualvecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateVecNorm() {
	updateVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateVecNorm(const int &d) {
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

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateDualVecNorm() {
	updateDualVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateDualVecNorm(const int &d) {
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

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateScalL2Norm(const int i) {
	NTL::matrix_row<IntMat> row(this->m_basis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateScalL2Norm(const int k1,
		const int k2) {
	for (int i = k1; i < k2; i++) {
		updateScalL2Norm(i);
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateDualScalL2Norm(const int i) {
	NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_dualvecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::updateDualScalL2Norm(const int k1,
		const int k2) {
	for (int i = k1; i < k2; i++) {
		updateDualScalL2Norm(i);
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::permute(int i, int j) {
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
	//bool b = this->m_xx[j];
	//this->m_xx[j] = this->m_xx[i];
	//this->m_xx[i] = b;
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::permuteNoDual(int i, int j) {
	if (i == j)
		return;
	for (int k = 0; k < this->m_dim; k++) {
		swap9(this->m_basis(j, k), this->m_basis(i, k));
	}
	swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
bool IntLatticeBase<Int, Real, RealRed>::checkDuality() {
	if (!this->m_withDual) {
		std::cout << "DO NOT USE IntLatticeBase::checkDuality without dual"
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

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::sort(int d)
/*
 * We assume that the square lengths are already updated.
 * This gives flexibility to the user to put something else than
 * the square Euclidean length in V.vecNorm, W.vecNorm, etc.
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
			permute(i, k);
	}
}

/*=========================================================================*/

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::sortNoDual(int d)
/*
 * We assume that the square lengths are already updated.
 * This gives flexibility to the user to put something else than
 * the square Euclidean length in V.vecNorm, W.vecNorm, etc.
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

template<typename Int, typename Real, typename RealRed>
void IntLatticeBase<Int, Real, RealRed>::write() const {
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
		if (this->m_withDual) {
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
	if (this->m_withDual)
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
		if (this->m_withDual) {
			if (this->m_dualvecNorm[i] < 0)
				std::cout << "NaN OR Not computed";
			else {
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

template<typename Int, typename Real, typename RealRed>
std::string IntLatticeBase<Int, Real, RealRed>::toStringBasis() const {
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

template<typename Int, typename Real, typename RealRed>
std::string IntLatticeBase<Int, Real, RealRed>::toStringDualBasis() const {
	std::ostringstream os;
	os << "Dual Basis:\n";
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

extern template class IntLatticeBase<std::int64_t, std::int64_t, double, double> ;
extern template class IntLatticeBase<NTL::ZZ, NTL::ZZ, double, double> ;
extern template class IntLatticeBase<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR> ;

} // namespace LatticeTester

#endif
