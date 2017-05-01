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
 * \file latticetester/PolTypes.h
 *
 * Sets the <tt>typedef</tt>’s for all the types used in LatMRG. Depending on
 * how `NTL_TYPES_CODE` is defined, all types used will be primitives, like `long`,
 * `double`, etc., or the large number types defined in NTL will be used. The
 * `NTL_TYPES_CODE` variable is defined in the Makefile. \anchor REF__Types_mod_Types
 *
 */

#ifndef LATTICETESTER__TYPES_H
#define LATTICETESTER__TYPES_H


/** Ayman
    * PolBasis USES NTL 
    * for the moment we only consider the case NTL_TYPES_CODE = LLDD
    */


#include <NTL/vector.h>
#include <NTL/matrix.h>
#include "ntlwrap.h"
using NTL::matrix_row;
#include <NTL/GF2E.h>
using namespace NTL;


typedef NTL::matrix<GF2E>  EMat;
typedef NTL::vector<GF2E>   EVect;

// the case  "ZZDD"
typedef long MScal;
typedef NTL::vector<long> MVect;
typedef NTL::matrix<long> MMat;

typedef long BScal;
typedef NTL::vector<long>   BVect;
typedef NTL::matrix<long>   BMat;
typedef double              NScal;
typedef NTL::vector<double> NVect;
typedef NTL::matrix<double> NMat;
typedef double              RScal;
typedef NTL::vector<double> RVect;
typedef NTL::matrix<double> RMat;


   namespace LatMRG {
      typedef void ProcII (int, int);
   }


#endif
