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

#include "latticetester/Const.h"
#include <string>


using namespace std;

namespace LatticeTester
{

//===========================================================================

string toStringNorm (NormType norm)
{
   switch (norm) {
   case SUPNORM:
      return "SUPNORM";
   case L1NORM:
      return "L1NORM";
   case L2NORM:
      return "L2NORM";
   case ZAREMBANORM:
      return "ZAREMBANORM";
   default:
      return "***** NormType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringGen (GenType gen)
{
   switch (gen) {
   case MWC:
      return "MWC";
   case LCG:
      return "LCG";
   case MRG:
      return "MRG";
   case KOROBOV:
      return "KOROBOV";
   case RANK1:
      return "RANK1";
   default:
      return "***** GenType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringLattice (LatticeType lat)
{
   switch (lat) {
   case FULL:
      return "FULL";
   case RECURRENT:
      return "RECURRENT";
   case ORBIT:
      return "ORBIT";
   case PRIMEPOWER:
      return "PRIMEPOWER";
   default:
      return "***** LatticeType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringPrime (PrimeType stat)
{
   switch (stat) {
   case PRIME:
      return "PRIME";
   case PROB_PRIME:
      return "PROB_PRIME";
   case COMPOSITE:
      return "COMPOSITE";
   default:
      return "UNKNOWN";
   }
}


//===========================================================================

string toStringDecomp (DecompType deco)
{
   switch (deco) {
   case DECOMP:
      return "DECOMP";
   case DECOMP_WRITE:
      return "DECOMP_WRITE";
   case DECOMP_READ:
      return "DECOMP_READ";
   case DECOMP_PRIME:
      return "DECOMP_PRIME";
   default:
      return "***** DecompType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringOutput (OutputType sort)
{
   switch (sort) {
   case TERMINAL:
      return "TERMINAL";
   case RES:
      return "RES";
   case TEX:
      return "TEX";
   case GEN:
      return "GEN";
   default:
      return "***** OutputType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringCriterion (CriterionType criter)
{
   switch (criter) {
   case SPECTRAL:
      return "SPECTRAL";
   case BEYER:
      return "BEYER";
   case PALPHA:
      return "PALPHA";
   default:
      return "***** CriterionType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringNorma (NormaType norma)
{
   switch (norma) {
   case BESTLAT:
      return "BESTLAT";
   case LAMINATED:
      return "LAMINATED";
   case ROGERS:
      return "ROGERS";
   case MINKL1:
      return "MINKL1";
   case MINKOWSKI:
      return "MINKOWSKI";
   case NORMA_GENERIC:
      return "NORMA_GENERIC";
   default:
      return "***** NormaType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringCalc (CalcType calc)
{
   switch (calc) {
   case PAL:
      return "PAL";
   case NORMPAL:
      return "NORMPAL";
   case BAL:
      return "BAL";
   case SEEKPAL:
      return "SEEKPAL";
   default:
      return "***** CalcType: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringImplemCond (ImplemCond cond)
{
   switch (cond) {
   case NO_COND:
      return "NO_COND";
   case APP_FACT:
      return "APP_FACT";
   case EQUAL_COEF:
      return "EQUAL_COEF";
   case ZERO_COEF:
      return "ZERO_COEF";
   case POWER_TWO:
      return "POWER_TWO";
   default:
      return "***** ImplemCond: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

string toStringSearchMethod (SearchMethod method)
{
   switch (method) {
   case EXHAUST:
      return "EXHAUST";
   case RANDOM:
      return "RANDOM";
   default:
      return "***** SearchMethod: IMPOSSIBLE CASE ";
   }
}


//===========================================================================

}
