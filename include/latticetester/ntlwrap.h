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
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>


NTL_CLIENT

typedef NTL::Vec<double> vec_double;
typedef NTL::Mat<double> mat_double;
typedef NTL::Vec<long>   vec_zz;
typedef NTL::Mat<long>   mat_zz;

// define a specialization of vector<T> as a wrapper around NTL vec_T
#define WRAP_NTL_VECTOR(S,T) \
   template <> \
   class vector<T> : public vec_##S \
   { \
   public: \
      typedef T value_type; \
      typedef T* pointer; \
      typedef T& reference; \
      typedef const T* const_pointer; \
      typedef const T& const_reference; \
      typedef long size_type; \
    \
      inline vector<T>() {} \
      inline vector<T>(size_type size) : vec_##S(INIT_SIZE, size) {} \
      inline vector<T>(const vec_##S &v) : vec_##S(v) {} \
    \
      inline void resize(size_type size) { SetLength(size); } \
      inline void clear() { kill(); } \
      /*    inline void swap (vector<T> &v) { NTL::swap(*this, v); } */ \
      inline size_type size() const { return length(); } \
      inline size_type max_size() const { return MaxLength(); } \
      inline bool empty() const { return size() == 0; } \
    \
      inline const_reference operator()(size_type i) const { return (*this)[i]; } \
      inline reference operator()(size_type i) { return (*this)[i]; } \
   }

// define a specialization of matrix<T> as a wrapper around NTL mat_T
#define WRAP_NTL_MATRIX(S,T) \
   template<> \
   class matrix<T> : public mat_##S \
   { \
   public: \
      typedef T value_type; \
      typedef T* pointer; \
      typedef T& reference; \
      typedef const T* const_pointer; \
      typedef const T& const_reference; \
      typedef long size_type; \
    \
      inline matrix<T>() {} \
      inline matrix<T>(const mat_##S& a) : mat_##S(a) {} \
      inline matrix<T>(size_type size1, size_type size2) : mat_##S(INIT_SIZE, size1, size2) {} \
    \
      inline void resize(size_type size1, size_type size2) { SetDims(size1, size2); } \
      inline void clear() { kill(); } \
      /*   inline void swap (matrix<T> &m) { NTL::swap(*this, m); }*/   \
      inline size_type size1() const { return NumRows(); } \
      inline size_type size2() const { return NumCols(); } \
    \
      inline reference operator()(size_type i, size_type j) { return (*this)[i][j]; } \
      inline const_reference operator()(size_type i, size_type j) const { return (*this)[i][j]; } \
   }

namespace NTL
{
   template <typename T> class vector { /* empty generic template: must be specialized */ };
   template <typename T> class matrix { /* empty generic template: must be specialized */ };

   // declare specializations
   WRAP_NTL_VECTOR(ZZ,ZZ);
   WRAP_NTL_MATRIX(ZZ,ZZ);
   WRAP_NTL_VECTOR(RR,RR);
   WRAP_NTL_MATRIX(RR,RR);
   WRAP_NTL_VECTOR(double,double);
   WRAP_NTL_MATRIX(double,double);
   WRAP_NTL_VECTOR(zz,long);
   WRAP_NTL_MATRIX(zz,long);

   // matrix proxy
   template <class M>
   class matrix_row : public vector<typename M::value_type> {
   public:
      inline matrix_row(M& data, typename M::size_type i) { this->_vec__rep = (typename M::value_type*&) data[i]._vec__rep; }
      inline ~matrix_row() { this->_vec__rep = 0; /* avoid destruction in parent class */ }
   };
}

#endif // __NTLWRAP_H__
