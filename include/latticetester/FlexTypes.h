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


#ifndef LATTICETESTER__FLEXTYPES_H
#define LATTICETESTER__FLEXTYPES_H

  /**
   * Here we describe the flexible types used by LatticeTester to large integers
   * and high-precision real numbers, as well as vectors and arrays of such numbers.
   * Those flexible types are used in templates when instantiating objects.
   *
   * The flexible type `Int` represents an integer. It can be either `std::int64_t`
   * or the NTL integer type `NTL::ZZ`. The types `IntVec` and `IntMat` represent
   * vectors and matrices of such flexible integers, respectively.
   * The flexible type `Real` represents a real number.
   * It can be either `double` or the NTL real type `NTL::RR`.
   * The types `RealVec` and `RealMat` represent vectors
   * and matrices of such real numbers.  The `Real` type is used for the Cholesky
   * decomposition (which can be very sensitive to accuracy), to compute bounds
   * in the BB algorithm, and for vector norms.
   */
   

     #define IntVec   NTL::vector<Int> 
     #define IntMat   NTL::matrix<Int> 
     #define RealVec  NTL::vector<Real> 
     #define RealMat  NTL::matrix<Real> 
     
  

#endif
