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

  string toStringCriterion (CriterionType criter)
  {
    switch (criter) {
      case LENGTH:
        return "LENGTH";
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

  string toStringProblem (ProblemType prob)
  {
    switch (prob) {
      case BASIS:
        return "BASIS";
      case DUAL:
        return "DUAL";
      case REDUCTION:
        return "REDUCTION";
      case SHORTEST:
        return "SHORTEST";
      case MERIT:
        return "MERIT";
      default:
        return "***** ProblemType: IMPOSSIBLE CASE ";
    }
  }


  //===========================================================================

  string toStringNorma (NormaType norma)
  {
    switch (norma) {
      case BESTLAT:
        return "BESTLAT";
      case BESTBOUND:
        return "BESTBOUND";
      case LAMINATED:
        return "LAMINATED";
      case ROGERS:
        return "ROGERS";
      case MINKL1:
        return "MINKL1";
      case MINK:
        return "MINK";
      case L1:
        return "L1";
      case L2:
        return "L2";
      case NONE:
        return "NONE";
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

  string toStringPreRed (PreReductionType prered)
  {
    switch (prered) {
      case FULL:
        return "FULL";
      case BKZ:
        return "BKZ";
      case DIETER:
        return "DIETER";
      case LLL:
        return "LLL";
      case NOPRERED:
        return "NOPRERED";
      default:
        return "***** PreReductionType: IMPOSSIBLE CASE ";
    }
  }


  //===========================================================================

  string toStringPrecision (PrecisionType precision)
  {
    switch (precision) {
      case DOUBLE:
        return "DOUBLE";
      case QUADRUPLE:
        return "QUADRUPLE";
      case EXPONENT:
        return "EXPONENT";
      case ARBITRARY:
        return "ARBITRARY";
      case EXACT:
        return "EXACT";
      default:
        return "***** PrecisionType: IMPOSSIBLE CASE ";
    }
  }

  //===========================================================================


  string toStringOutput (OutputType sort)
  {
    switch (sort) {
      case TERM:
        return "TERM";
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

  const std::array<unsigned int, NB_PRIMES> PRIMES_ARRAY = {{
    #include "../data/primes.dat"
  }};


}
