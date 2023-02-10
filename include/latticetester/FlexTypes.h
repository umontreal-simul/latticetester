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

#include <cstdint>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/RR.h>

#include "latticetester/NTLWrap.h"

  /**
   * Here we describe the flexible types used by LatticeTester to large integers
   * and high-precision real numbers, as well as vectors and arrays of such numbers.
   * Those flexible types are used in templates when instantiating objects.
   *
   * The type `Int` represents a flexible integer. It can be either `std::int64_t`
   * or the NTL integer type `NTL::ZZ`. The types `IntVec` and `IntMat` represent
   * vectors and matrices of such flexible integers, respectively.
   *
   * For the real numbers, we use the two flexible representations `Real` and `RealRed`.
   * Each one can be either `double` or the NTL real type `NTL::RR`.
   * The types `RealVec`, `RealRedVec`, `RealMat`, and `RealRedMat` represent vectors
   * and matrices of such real numbers.  The `RealRed` type is used for the Cholesky
   * decomposition (which can be very sensitive to accuracy) and to compute bounds
   * in the BB algorithm, whereas `Real` is used for other real numbers.
   */

namespace LatticeTester {

    typedef NTL::ZZ       Int;
    typedef std::int64_t  Int;
    typedef NTL::RR       Real;
    typedef double        Real;
    typedef double        RealRed;
    typedef NTL::RR       RealRed;


   typedef NTL::vector<Int> IntVec;
   typedef NTL::matrix<Int> IntMat;
   typedef NTL::vector<Real> RealVec;
   typedef NTL::matrix<Real> RealMat;
   typedef NTL::vector<RealRed> RealRedVec;
   typedef NTL::matrix<RealRed> RealRedMat;
}

#endif
