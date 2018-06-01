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

#ifndef __NTLWRAP_H__
#define __NTLWRAP_H__

#include <NTL/vector.h>
#include <NTL/matrix.h>

#include <cstdint>

// The next line has been removed because it makes us use the namespaces std and NTL
//NTL_CLIENT

/**
 * This module contains extensions of certain classes in NTL. It was previously
 * necessary because NTL and boost (an old dependency) did not use the same 
 * function names and indices.
 *
 * This name conversion was meant to have the same function names in boost and NTL
 * and allows us to have LatticeTester work with either boost library or NTL library
 * depending on preprocessing statements.
 */

namespace NTL
{
  /**
   * A subclass of the `NTL::Vec<T>` class. It extends its parent with a few 
   * methods and overloads a few others with more compatible defaults.
   * */
  template <typename T> class vector : public Vec<T>{
    public:

      /** \todo Check if this type is valid or should be reverted to `long`*/
      typedef std::int64_t size_type;

      /**
       * Empty constructor. Creates a vector of size 0.
       */
      vector<T>() {};
      /**
       * Allocation constructor. Allocates the space needed by a vector of size
       * `size`.
       * @param size The size of the vector
       */
      vector<T>(size_type size) : Vec<T>(INIT_SIZE, size) {};
      /**
       * Copy constructor. Makes a copy of `v`.
       * @param v Vector to be copied.
       */
      vector<T>(const Vec<T> &v) : Vec<T>(v) {};

      /**
       * Destructor. Frees memory of all T and destroys the vector.
       * */
      ~vector () {};

      /**
       * Set the vector lenght to size.
       * New objects are initialized using the default contructor for T
       * This uses `%NTL::%Vec<T>::%SetLength(size)`
       * @param size The new size wanted
       */
      void resize(size_type size) { this->SetLength(size); }

      /**
       * Releases space and sets length to 0.
       * This uses `%NTL::%Vec<T>::%kill()`.
       */
      void clear() { this->kill(); }

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
      size_type size() const { return this->length(); }

      /**
       * Returns the number of allocated and initialized elements in the vector.
       * It is possible for the vector length to differ from this one if the max
       * length has been preemptively expanded or if the vector has been 
       * shrinked because NTL does not free components until explicitly asked.
       * This uses `%NTL::%Vec<T>::%MaxLength()`.
       */
      size_type max_size() const { return this->MaxLength(); }

      /**
       * Returns `true` if the vector has length 0 and `false` otherwise.
       */
      bool empty() const { return size() == 0; }

      /**
       * Change in the indexation reference for () operator to start from 0.
       * In `NTL::Vec<T>` the () operator starts from 1 which is not compatible
       * with boost.
       */
      const T& operator()(size_type i) const { return (*this)[i]; }
      T&  operator()(size_type i) { return (*this)[i]; }

  };

  /**
   * A subclass of the `NTL::Mat<T>` class. It extends its parent with a few 
   * methods and overloads a few others with more compatible defaults.
   * */
  template <typename T> class matrix : public Mat<T>{

    public:

      /** \todo Check if this type is valid or should be reverted to `long`*/
      typedef std::int64_t size_type;

      /**
       * Empty constructor.
       */
      matrix<T>() {}
      /**
       * Copy constructor. Creates a new matrix that is a copy of a.
       * @param a Matrix to be copied.
       * */
      matrix<T>(const Mat<T>& a) : Mat<T>(a) {}
      /**
       * Allocation constructor.
       * Creates and allocates a `size1`\f$\times\f$`size2` matrix, initializing
       * the elements T with their default constructor.
       * @param size1 Height of the matrix
       * @param size2 Width of the matrix
       * */
      matrix<T>(size_type size1, size_type size2) :
        Mat<T>(INIT_SIZE, size1, size2) {}

      /**
       * Set the matrix dimensions to `size1`\f$\times\f$`size2`.
       * This uses `%NTL::%Mat<T>::%SetDims(size1, size2)`.
       * @param size1 New height of the matrix.
       * @param size2 New width of the matrix.
       */
      void resize(size_type size1, size_type size2)
      {
        this->SetDims(size1, size2);
      }

      /**
       * Releases space and sets the matrix this size \f$0\times 0\f$.
       * This uses `%NTL::%Mat<T>::%kill()`.
       */
      void clear() { this->kill(); }

      /*
       * The function in this comment adds nothing since NTL::Vec<T> now
       * implements a member function: `void swap(Mat<T>& other);`. It is kept
       * here for history.
       * Thanks to Ayman for pointing this out.
       *
       * inline void swap (matrix<T> &m) { NTL::swap(*this, m); }
       * */

      /**
       * Returns the number of rows of the matrix.
       */
      size_type size1() const { return this->NumRows(); }

      /**
       * Returns the number of columns of the matrix.
       */
      size_type size2() const { return this->NumCols(); }

      /**
       * Overload to change the indexation reference for (i,j) operator to start 
       * from 0.
       * In NTL::Vec<T> the (i,j) operator starts from 1 which is not compatible
       * with boost.
       */
      T& operator()(size_type i, size_type j) { return (*this)[i][j]; }
      const T& operator()(size_type i, size_type j) const {
        return (*this)[i][j];
      }

  };




  /**
   * An extension of `NTL::vector<T>` implemented in this module to be used as
   * a matrix row.
   */
  template <class M>
    class matrix_row : public vector<typename M::value_type> {
      public:
        inline matrix_row(M& data, typename M::size_type i) {
          this->_vec__rep = (typename M::value_type*&) data[i]._vec__rep;
        }
        inline ~matrix_row() {
          this->_vec__rep = 0; /* avoid destruction in parent class */
        }
    };
} // End namespace NTL

#endif // __NTLWRAP_H__
