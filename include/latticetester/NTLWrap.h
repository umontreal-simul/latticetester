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

#ifndef LATTICETESTER__NTLWRAP_H
#define LATTICETESTER__NTLWRAP_H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include <cstdint>

/**
 * This module extends the `Vec` and `Mat` classes of NTL. It was previously
 * necessary because NTL and boost (an old dependency) did not use the same 
 * function names and indices.
 * Most of the methods below only define alias names for some NTL methods.
 * This was done at the time to have the same names for methods that are in both
 * boost and NTL, allowing LatticeTester to work with either the boost or NTL library,
 * depending on pre-processing statements. But these alias could probably be removed
 * because we no longer use boost in LatticeTester.
 * However, they may still be used in LatNet Builder.                 **************
 *
 * New functions have also been implemented in this module as a way to overload a
 * few operators and methods of NTL (especially on matrix and vector types) to
 * the usage of `NTL::Mat<std::uint64_t>` because some basic utilities that we need
 * for those integers are not offered in NTL.
 */

namespace NTL {

/**
 * A subclass of the `NTL::Vec<T>` class. It extends its parent with a some additional
 * methods and overloads (changes) a few others for compatibility with boost.
 * Suggestion: Avoid using these redefined methods, because this can be very confusing
 * and dangerous. We should probably remove them!                            *********
 */
template<typename T> class vector: public Vec<T> {
public:

	typedef std::int64_t size_type;

	/**
	 * Empty constructor. Creates a vector of size 0.
	 */
	vector<T>() {
	}
	;
	/**
	 * Allocation constructor. Allocates the space needed by a vector of size
	 * `size`.
	 * @param size The size of the vector
	 */
	vector<T>(size_type size) :
			Vec<T>(INIT_SIZE, size) {
	}
	;
	/**
	 * Copy constructor. Makes a copy of `v`.
	 * @param v Vector to be copied.
	 */
	vector<T>(const Vec<T> &v) :
			Vec<T>(v) {
	}
	;

	/**
	 * Destructor. Frees memory of all T and destroys the vector.
	 * */
	~vector() {
	}
	;

	/**
	 * Set the vector lenght to size.
	 * New objects are initialized using the default contructor for T
	 * This uses `%NTL::%Vec<T>::%SetLength(size)`
	 * @param size The new size wanted
	 */
	void resize(size_type size) {
		this->SetLength(size);
	}

	/**
	 * Releases space and sets length to 0.
	 * This uses `%NTL::%Vec<T>::%kill()`.
	 * BAD RENAMING!   clear() should not be the same as kill().    **********
	 */
	void clear() {
		this->kill();
	}

	/*
	 * The function in this comment adds nothing since NTL::Vec<T> now
	 * implements a member function: `void swap(Vec<T>& other);`. It is kept
	 * here for history.
	 * Thanks to Ayman for pointing this out.
	 *
	 * inline void swap (matrix<T> &m) { NTL::swap(*this, m); }
	 */

	/**
	 * Returns the current lenght of the vector.
	 * This uses `%NTL::%Vec<T>::%length()`.
	 */
	size_type size() const {
		return this->length();
	}

	/**
	 * Returns the number of allocated and initialized elements in the vector.
	 * It is possible for the vector length to differ from this one if the max
	 * length has been preemptively expanded or if the vector has been
	 * shrinked because NTL does not free components until explicitly asked.
	 * This uses `%NTL::%Vec<T>::%MaxLength()`.
	 */
	size_type max_size() const {
		return this->MaxLength();
	}

	/**
	 * a copy from NTL::Vec<T>::append
	 */
	void push_back(const T &element) {
		this->append(element);
	}

	/**
	 * return the last element, this one is for STL compatibility.   ?????
	 */
	T back() {
		return (*(this->end() - 1));
	}

	/**
	 * Returns `true` if the vector has length 0 and `false` otherwise.
	 */
	bool empty() const {
		return size() == 0;
	}

	/**
	 * Change in the indexation reference for () operator to start from 0.
	 * In `NTL::Vec<T>` the () operator starts from 1 which is not compatible
	 * with boost.
	 */
	const T& operator()(size_type i) const {
		return (*this)[i];
	}
	T& operator()(size_type i) {
		return (*this)[i];
	}

};

/**
 * A subclass of the `NTL::Mat<T>` class. It extends its parent with a few
 * methods and overloads a few others with more compatible defaults.
 * */
template<typename T> class matrix: public Mat<T> {

public:

	typedef std::int64_t size_type;

	/**
	 * Empty constructor.
	 */
	matrix<T>() {
	}
	/**
	 * Copy constructor. Creates a new matrix that is a copy of a.
	 * @param a Matrix to be copied.
	 * */
	matrix<T>(const Mat<T> &a) :
			Mat<T>(a) {
	}
	/**
	 * Allocation constructor.
	 * Creates and allocates a `size1`\f$\times\f$`size2` matrix, initializing
	 * the elements T with their default constructor.
	 * @param size1 Height of the matrix
	 * @param size2 Width of the matrix
	 * */
	matrix<T>(size_type size1, size_type size2) :
			Mat<T>(INIT_SIZE, size1, size2) {
	}

	/**
	 * Set the matrix dimensions to `size1`\f$\times\f$`size2`.
	 * This is just an alis to `%NTL::%Mat<T>::%SetDims(size1, size2)`.
	 * @param size1 New number of rows in the matrix.
	 * @param size2 New number of columns in the matrix.
	 */
	void resize(size_type size1, size_type size2) {
		this->SetDims(size1, size2);
	}

	/**
	 * Releases space and sets the matrix this size \f$0\times 0\f$.
	 * This uses `NTL::Mat<T>::kill()`.
	 *
	 * NOT COMPATIBLE WITH clear() in NTL and ELSEWHERE !!!   *********
	 * It should initialize to 0 but not destroy the object !!
	 */
	void clear() {
		this->kill();
	}

	/**
	 * Returns the number of rows of the matrix.
	 */
	size_type size1() const {
		return this->NumRows();
	}

	/**
	 * Returns the number of columns of the matrix.
	 */
	size_type size2() const {
		return this->NumCols();
	}

	/**
	 * Overload to change the indexation reference for (i,j) operator to start
	 * from 0.
	 * In NTL::Vec<T> the (i,j) operator starts from 1 which is not compatible
	 * with boost.
	 */
	T& operator()(size_type i, size_type j) {
		return (*this)[i][j];
	}
	const T& operator()(size_type i, size_type j) const {
		return (*this)[i][j];
	}
};

/**
 * An extension of `NTL::vector<T>` implemented in this module to be used as
 * a matrix row.
 */
template<class M>
class matrix_row: public vector<typename M::value_type> {
public:
	/**
	 * Copy constructor? Creates an object that is a copy if `data`.
	 * */
	inline matrix_row(M &data, typename M::size_type i) {
		this->_vec__rep = (typename M::value_type*&) data[i]._vec__rep;
	}
	/**
	 * Empty constructor.
	 * */
	inline ~matrix_row() {
		this->_vec__rep = 0; /* avoid destruction in parent class */
	}
};

//============================================================================

/// Renaming this supported type for convenience.
typedef Mat<std::int64_t> Mat_64;
/// Renaming this supported type for convenience.
typedef Vec<std::int64_t> Vec_64;

//============================================================================

/**
 * Transforms `mat` into the identity matrix of dimensions
 * \f$\text{dim}\times\text{dim}\f$.
 * */
void ident(Mat_64 &mat, long dim);

/**
 * Transposes `A` into `X`.
 * This is a template overload of the transpose function of NTL. This
 * does basically the same thing as the implementation NTL uses. It might
 * be necessary to implement the swap function for the type `T` for this to work.
 * */
template<typename T>
void transpose(NTL::Mat<T> &X, const NTL::Mat<T> &A) {
	long n = A.NumRows();
	long m = A.NumCols();

	long i, j;

	// If both matrices have the same address, we need to transpose in place
	if (&X == &A) {
		if (n == m)
			for (i = 0; i < n; i++)
				for (j = 0; j < m; j++)
					std::swap(X[i][j], X[j][i]);
		else {
			NTL::Mat<T> tmp;
			tmp.SetDims(m, n);
			for (i = 0; i < n; i++)
				for (j = 0; j < m; j++)
					tmp[j][i] = A[i][j];
			X.kill();
			X = tmp;
		}
	} else {
		X.SetDims(m, n);
		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				X[j][i] = A[i][j];
	}
}

/**
 * Another implementation of the `transpose` function.  Returns the ranspose of `a`.
 * */
template<typename T>
inline NTL::Mat<T> transpose(const NTL::Mat<T> &a) {
	NTL::Mat<T> x;
	transpose(x, a);
	return x;
}

//============================================================================

/**
 * \name Operator overloads
 * @{
 * These are operator overloads for Mat_64 and Vec_64 types. Only the
 * overloads we currently use are defined.
 * */
Vec_64 operator*(const Vec_64 &vec, std::int64_t a);
Vec_64 operator*(std::int64_t a, const Vec_64 &vec);
std::int64_t operator*(const Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator+=(Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator-=(Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator*=(Vec_64 &vec, std::int64_t a);
Mat_64& operator*=(Mat_64 &mat, std::int64_t a);
Mat_64 operator*(const Mat_64 &mat1, const Mat_64 &mat2);
/**
 * @}
 * */

double determinant(const NTL::matrix<std::int64_t> &mat);

} // End namespace NTL

#endif // LATTICETESTER__NTLWRAP_H
