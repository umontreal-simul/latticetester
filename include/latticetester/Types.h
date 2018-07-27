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

/**
 * \file latticetester/Types.h
 * See the \ref Types "Types" module and \ref DIY
 * */
/**
 * \defgroup Types Types
 * @{
 *
 * Sets standard <tt>typedef</tt>’s for the types that can be used in
 * LatticeTester. Although it is not directly used anymore, it is provided to 
 * external users as a way to easily use supported types, as well as for
 * historical reasons. The types depend on how `NTL_TYPES_CODE` is defined: 
 * types used can be primitives, like `std::int64_t`, `double`, etc., or large
 * number types defined in the NTL dependency. The `NTL_TYPES_CODE` variable 
 * can be defined when compiling (for example, with `-DNTL_TYPES_CODE=x` passed
 * to a compiler like g++) or in the executables.
 *
 * Originally this file was created as a way to define the types used by the 
 * program by passing a definition of `NTL_TYPES_CODE` to the compiler. This 
 * trick was used in the first version of this software in Modula 2 to switch 
 * from hardware supported types to software implemented types to have access
 * to a much needed arbitrary precision in some calculations. With modern 
 * programming languages, such as C++, we have access to tools, such as
 * overloading and templates, that makes it much easier to avoid having to 
 * recompile every time we want to change types. Types are now declared when
 * intanciating a class rather than globally with a definition passed to a 
 * compiler, making it possible to create an executable that supports any 
 * combination of types (virtually as we only test and support for the types
 * combinations of this file as of now).
 *
 * Some of the `typedef`s define the same type and this is intentional. It comes 
 * from the old functionnalities of the program. There are different type 
 * prefixes: `M`, `B`, `N` and `R`. Originally, these prefixes and their 
 * respective `typedef`s where used in different modules to do different kind 
 * of calculations.
 * - `M` stands for modulo and served as a type for general computations in the 
 *   lattice (as most of the calculations on integer are made modulo m). Note
 *   that for an unknown reason the actual modulo types are not used anymore.
 *   These types must be some kind of integer types because we build the 
 *   lattices in \f$\mathbb{Z}_m\f$.
 * - `B` stands for basis. These types where used in basis representation,
 *    calculation and reductions. These types must be an integer type for the 
 *    same reason as `M`
 * - `N` stands for norm. These types are used when there is a norm or a scalar
 *   product to compute. Because these kind of calculations can, in some cases,
 *   give floating point numbers (especially in the case of a norm) these types 
 *   must be floating point numers.
 * - `R` stands for reduced. These types where used in calculations requiring 
 *   floating point numbers when reducing a basis (such as Cholesky 
 *   decomposition) as well as figures of merit. These types must be floating 
 *   point numbers.
 *
 * In the current version of the software, the types where replaced by template
 * typenames (that feel slightly less criptic to me) as follows:
 * - `MScal -> Int`
 * - `MVect -> IntVec`
 * - `MMat  -> IntMat`
 * - `M*P   -> nothing (never used)`
 * - `BScal -> BasInt`
 * - `BVect -> BasIntVec`
 * - `BMat  -> BasIntMat`
 * - `NScal -> Dbl`
 * - `NVect -> DblVec`
 * - `NMat  -> DblMat`
 * - `RScal -> RedDbl`
 * - `RVect -> RedDblVec`
 * - `RMat  -> RedDblMat`
 *
 * The typedefs are defined as follow for the different values of the
 * NTL_TYPES_CODE constant.
 * \code{.cpp}
 * #if NTL_TYPES_CODE == 1
 * // the case  "LLDD"
 *
 * #include "NTL/lzz_p.h"
 * #include "NTL/vec_lzz_p.h"
 * #include "NTL/mat_lzz_p.h"
 * #include "NTL/lzz_pX.h"
 * #include "NTL/lzz_pE.h"
 * #include "NTL/lzz_pEX.h"
 *
 * typedef std::int64_t              MScal;
 * typedef NTL::vector<std::int64_t> MVect;
 * typedef NTL::matrix<std::int64_t> MMat;
 * typedef zz_p     MScalP;
 * typedef NTL::vec_zz_p MVectP;
 * typedef NTL::mat_zz_p MMatP;
 * typedef std::int64_t                BScal;
 * typedef NTL::vector<std::int64_t>   BVect;
 * typedef NTL::matrix<std::int64_t>   BMat;
 * typedef double              NScal;
 * typedef NTL::vector<double> NVect;
 * typedef NTL::matrix<double> NMat;
 * typedef double              RScal;
 * typedef NTL::vector<double> RVect;
 * typedef NTL::matrix<double> RMat;
 * typedef zz_pX    PolX;
 * typedef zz_pE    PolE;
 *
 * #elif NTL_TYPES_CODE == 2
 * // the case  "ZZDD"
 * 
 * #include "NTL/ZZ.h"
 * #include "NTL/vec_ZZ.h"
 * #include "NTL/mat_ZZ.h"
 * #include "NTL/ZZ_p.h"
 * #include "NTL/vec_ZZ_p.h"
 * #include "NTL/mat_ZZ_p.h"
 * #include "NTL/ZZ_pE.h"
 * #include "NTL/ZZ_pX.h"
 * #include "NTL/ZZ_pEX.h"
 * typedef NTL::ZZ              MScal;
 * typedef NTL::vector<NTL::ZZ> MVect;
 * typedef NTL::matrix<NTL::ZZ> MMat;
 * typedef NTL::ZZ_p     MScalP;
 * typedef NTL::vec_ZZ_p MVectP;
 * typedef NTL::mat_ZZ_p MMatP;
 * typedef NTL::ZZ              BScal;
 * typedef NTL::vector<NTL::ZZ> BVect;
 * typedef NTL::matrix<NTL::ZZ> BMat;
 * typedef double              NScal;
 * typedef NTL::vector<double> NVect;
 * typedef NTL::matrix<double> NMat;
 * typedef double              RScal;
 * typedef NTL::vector<double> RVect;
 * typedef NTL::matrix<double> RMat;
 * typedef NTL::ZZ_pX    PolX;
 * typedef NTL::ZZ_pE    PolE;
 * #elif NTL_TYPES_CODE == 3
 * // the case  "ZZRR"
 *
 * #include "NTL/ZZ.h"
 * #include "NTL/vec_ZZ.h"
 * #include "NTL/mat_ZZ.h"
 * #include "NTL/RR.h"
 * #include "NTL/vec_RR.h"
 * #include "NTL/mat_RR.h"
 * #include "NTL/ZZ_p.h"
 * #include "NTL/vec_ZZ_p.h"
 * #include "NTL/mat_ZZ_p.h"
 * #include "NTL/ZZ_pE.h"
 * #include "NTL/ZZ_pX.h"
 * #include "NTL/ZZ_pEX.h"
 *
 * typedef NTL::ZZ              MScal;
 * typedef NTL::vector<NTL::ZZ> MVect;
 * typedef NTL::matrix<NTL::ZZ> MMat;
 * typedef NTL::ZZ_p     MScalP;
 * typedef NTL::vec_ZZ_p MVectP;
 * typedef NTL::mat_ZZ_p MMatP;
 * typedef NTL::ZZ              BScal;
 * typedef NTL::vector<NTL::ZZ> BVect;
 * typedef NTL::matrix<NTL::ZZ> BMat;
 * typedef NTL::RR              NScal;
 * typedef NTL::vector<NTL::RR> NVect;
 * typedef NTL::matrix<NTL::RR> NMat;
 * typedef NTL::RR              RScal;
 * typedef NTL::vector<NTL::RR> RVect;
 * typedef NTL::matrix<NTL::RR> RMat;
 * typedef NTL::ZZ_pX    PolX;
 * typedef NTL::ZZ_pE    PolE;
 * #endif
 * \endcode
 *
 * **Marc-Antoine** : I removed those because the using keyword
 * at file scope imports everything in the namespace NTL for every file 
 * including the header Types.h. I feel that, although LatticeTester does not 
 * use Types.h anymore, it could be a good idea to distribute it as a list of 
 * supported/default types combinations. In that use case it is an awful idea
 * to use the NTL namespace.
 * 
 *      using NTL::matrix_row;
 *      using namespace NTL;
 * @}
 * */

 ///\cond
#ifndef LATTICETESTER__TYPES_H
#define LATTICETESTER__TYPES_H

#include <cstdint>

#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "ntlwrap.h"
#if NTL_TYPES_CODE == 1
// the case  "LLDD"

#include "NTL/lzz_p.h"
#include "NTL/vec_lzz_p.h"
#include "NTL/mat_lzz_p.h"
#include "NTL/lzz_pX.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pEX.h"

typedef std::int64_t              MScal;
typedef zz_p     MScalP; // This appears nowhere
typedef NTL::vec_zz_p MVectP; // This appears nowhere
typedef NTL::mat_zz_p MMatP; // This appears only once
typedef std::int64_t                BScal;
typedef double              NScal;
typedef double              RScal;
typedef zz_pX    PolX; // This appears nowhere
typedef zz_pE    PolE; // This appears nowhere

#elif NTL_TYPES_CODE == 2
// the case  "ZZDD"

#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/vec_ZZ_p.h"
#include "NTL/mat_ZZ_p.h"
#include "NTL/ZZ_pE.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZ_pEX.h"

typedef NTL::ZZ              MScal;
typedef NTL::ZZ_p     MScalP;
typedef NTL::vec_ZZ_p MVectP;
typedef NTL::mat_ZZ_p MMatP;
typedef NTL::ZZ              BScal;
typedef double              NScal;
typedef double              RScal;
typedef NTL::ZZ_pX    PolX;
typedef NTL::ZZ_pE    PolE;

#elif NTL_TYPES_CODE == 3
// the case  "ZZRR"

#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/RR.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include "NTL/ZZ_p.h"
#include "NTL/vec_ZZ_p.h"
#include "NTL/mat_ZZ_p.h"
#include "NTL/ZZ_pE.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZ_pEX.h"

typedef NTL::ZZ              MScal;
typedef NTL::ZZ_p     MScalP;
typedef NTL::vec_ZZ_p MVectP;
typedef NTL::mat_ZZ_p MMatP;
typedef NTL::ZZ              BScal;
typedef NTL::RR              NScal;
typedef NTL::RR              RScal;
typedef NTL::ZZ_pX    PolX;
typedef NTL::ZZ_pE    PolE;
#endif

#if NTL_TYPES_CODE

typedef NTL::vector<MScal> MVect;
typedef NTL::matrix<MScal> MMat;
typedef NTL::vector<BScal> BVect;
typedef NTL::matrix<BScal> BMat;
typedef NTL::vector<NScal> NVect;
typedef NTL::matrix<NScal> NMat; // This appears nowhere
typedef NTL::vector<RScal> RVect;
typedef NTL::matrix<RScal> RMat;

#endif


namespace LatticeTester {
  typedef void ProcII (int, int);
}

#endif
///\endcond
