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

/**
 * \file latticetester/Types.h
 *
 * Sets the <tt>typedef</tt>’s for all the types used in LatticeTester. Depending on
 * how `NTL_TYPES_CODE` is defined, all types used will be primitives, like `long`,
 * `double`, etc., or the large number types defined in NTL will be used. The
 * `NTL_TYPES_CODE` variable is defined in the Makefile. \anchor REF__Types_mod_Types
 *
 */

#ifndef LATTICETESTER__TYPES_H
#define LATTICETESTER__TYPES_H

//#define WITH_NTL
//#define NTL_TYPES_CODE 2

#ifdef WITH_NTL

   #include <NTL/vector.h>
   #include <NTL/matrix.h>
   #include "ntlwrap.h"
   using NTL::matrix_row;
   using namespace NTL;

#if NTL_TYPES_CODE == 1
   // the case  "LLDD"

   #include "NTL/lzz_p.h"
   #include "NTL/vec_lzz_p.h"
   #include "NTL/mat_lzz_p.h"
   #include "NTL/lzz_pX.h"
   #include "NTL/lzz_pE.h"
   #include "NTL/lzz_pEX.h"

   typedef long              MScal;
   typedef NTL::vector<long> MVect;
   typedef NTL::matrix<long> MMat;
   typedef zz_p     MScalP;
   typedef vec_zz_p MVectP;
   typedef mat_zz_p MMatP;
   typedef long                BScal;
   typedef NTL::vector<long>   BVect;
   typedef NTL::matrix<long>   BMat;
   typedef double              NScal;
   typedef NTL::vector<double> NVect;
   typedef NTL::matrix<double> NMat;
   typedef double              RScal;
   typedef NTL::vector<double> RVect;
   typedef NTL::matrix<double> RMat;
   typedef zz_pX    PolX;
   typedef zz_pE    PolE;
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

   typedef ZZ              MScal;
   typedef NTL::vector<ZZ> MVect;
   typedef NTL::matrix<ZZ> MMat;
   typedef ZZ_p     MScalP;
   typedef vec_ZZ_p MVectP;
   typedef mat_ZZ_p MMatP;
   typedef ZZ              BScal;
   typedef NTL::vector<ZZ> BVect;
   typedef NTL::matrix<ZZ> BMat;
   typedef double              NScal;
   typedef NTL::vector<double> NVect;
   typedef NTL::matrix<double> NMat;
   typedef double              RScal;
   typedef NTL::vector<double> RVect;
   typedef NTL::matrix<double> RMat;
   typedef ZZ_pX    PolX;
   typedef ZZ_pE    PolE;
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

   typedef ZZ              MScal;
   typedef NTL::vector<ZZ> MVect;
   typedef NTL::matrix<ZZ> MMat;
   typedef ZZ_p     MScalP;
   typedef vec_ZZ_p MVectP;
   typedef mat_ZZ_p MMatP;
   typedef ZZ              BScal;
   typedef NTL::vector<ZZ> BVect;
   typedef NTL::matrix<ZZ> BMat;
   typedef RR              NScal;
   typedef NTL::vector<RR> NVect;
   typedef NTL::matrix<RR> NMat;
   typedef RR              RScal;
   typedef NTL::vector<RR> RVect;
   typedef NTL::matrix<RR> RMat;
   typedef ZZ_pX    PolX;
   typedef ZZ_pE    PolE;
#else
#error "NTL_TYPES_CODE is not defined"
#endif



   namespace LatticeTester {
      typedef void ProcII (int, int);
   }

#else

   // use Boost ublas
   #include <boost/numeric/ublas/vector.hpp>
   #include <boost/numeric/ublas/vector_proxy.hpp>
   #include <boost/numeric/ublas/matrix.hpp>
   #include <boost/numeric/ublas/matrix_proxy.hpp>


   typedef long                                 MScal;
   typedef boost::numeric::ublas::vector<long>  MVect;
   typedef boost::numeric::ublas::matrix<long>  MMat;

   typedef long                                 BScal;
   typedef boost::numeric::ublas::vector<long>  BVect;
   typedef boost::numeric::ublas::matrix<long>  BMat;

   typedef double                                 NScal;
   typedef boost::numeric::ublas::vector<double>  NVect;
   typedef boost::numeric::ublas::matrix<double>  NMat;

   typedef double                                 RScal;
   typedef boost::numeric::ublas::vector<double>  RVect;
   typedef boost::numeric::ublas::matrix<double>  RMat;

#endif
#endif


