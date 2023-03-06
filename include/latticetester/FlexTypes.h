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
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/NTLWrap.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#define   LD  1
#define   ZD  2
#define   ZR  3
#define   LR  4

     
     #if  TYPES_CODE  ==  LD 

	   typedef std::int64_t  Int;  
        typedef double  Real; 


     #elif  TYPES_CODE  == ZD 
         
	     typedef NTL::ZZ Int;
          typedef double Real;
    
     #elif   TYPES_CODE ==  ZR
   
	     typedef NTL::ZZ Int;
          typedef NTL::RR Real;

     
     #elif  TYPES_CODE  ==  LR 
  
	 	typedef std::int64_t  Int;  
          typedef NTL::RR Real;
	 
     #endif

     #ifdef TYPES_CODE
       typedef NTL::vector<Int> IntVec;
       typedef NTL::matrix<Int> IntMat;
       typedef NTL::vector<Real> RealVec;
       typedef NTL::matrix<Real> RealMat;
     #endif

#endif
