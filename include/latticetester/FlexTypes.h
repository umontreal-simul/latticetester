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

     
     #if  TYPES_CODE  ==   LD 

       /*
        #include "NTL/lzz_p.h"
        #include "NTL/vec_lzz_p.h"
        #include "NTL/mat_lzz_p.h"
        #include "NTL/lzz_pX.h"
        #include "NTL/lzz_pE.h"
        #include "NTL/lzz_pEX.h"
       */ 

	   typedef std::int64_t  Int;  
        typedef double  Real; 


     #elif  TYPES_CODE  == ZD 

         
          #include "NTL/ZZ.h"
          #include "NTL/vec_ZZ.h"
          #include "NTL/mat_ZZ.h"
          #include "NTL/ZZ_p.h"
          #include "NTL/vec_ZZ_p.h"
          #include "NTL/mat_ZZ_p.h"
          #include "NTL/ZZ_pE.h"
          #include "NTL/ZZ_pX.h"
          #include "NTL/ZZ_pEX.h"
         

	     typedef NTL::ZZ Int;
          typedef double Real;
    
     #elif   TYPES_CODE == ZR
         /*
          #include "NTL/ZZ.h"
          #include "NTL/vec_ZZ.h"
          #include "NTL/mat_ZZ.h"
          #include "NTL/ZZ_p.h"
          #include "NTL/vec_ZZ_p.h"
          #include "NTL/mat_ZZ_p.h"
          #include "NTL/ZZ_pE.h"
          #include "NTL/ZZ_pX.h"
          #include "NTL/ZZ_pEX.h"
          */

	     typedef NTL::ZZ Int;
          typedef NTL::RR Real;

     
     #elif  TYPES_CODE  ==  LR 

         /*
          #include "NTL/ZZ.h"
          #include "NTL/vec_ZZ.h"
          #include "NTL/mat_ZZ.h"
          #include "NTL/ZZ_p.h"
          #include "NTL/vec_ZZ_p.h"
          #include "NTL/mat_ZZ_p.h"
          #include "NTL/ZZ_pE.h"
          #include "NTL/ZZ_pX.h"
          #include "NTL/ZZ_pEX.h"

          */   

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
