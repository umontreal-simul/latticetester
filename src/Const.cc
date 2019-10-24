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
#include <iostream>


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

  //===========================================================================

  int toNormString       (NormType& norm, const std::string& str) {
    if (str == "SUPNORM") {
      norm = SUPNORM;
    } else if (str == "L1NORM") {
      norm = L1NORM;
    } else if (str == "L2NORM") {
      norm = L2NORM;
    } else if (str == "ZAREMBANORM") {
      norm = ZAREMBANORM;
    } else {
      std::cerr << str << " is not a NormType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toPrimeString      (PrimeType& prime, const std::string& str) {
    if (str == "UNKNOWN") {
      prime = UNKNOWN;
    } else if (str == "PRIME") {
      prime = PRIME;
    } else if (str == "PROB_PRIME") {
      prime = PROB_PRIME;
    } else if (str == "COMPOSITE") {
      prime = COMPOSITE;
    } else {
      std::cerr << str << " is not a PrimeType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toCriterionString  (CriterionType& criterion, const std::string& str) {
    if (str == "LENGTH") {
      criterion = LENGTH;
    } else if (str == "SPECTRAL") {
      criterion = SPECTRAL;
    } else if (str == "BEYER") {
      criterion = BEYER;
    } else if (str == "PALPHA") {
      criterion = PALPHA;
    } else if (str == "BOUND_JS") {
      criterion = BOUND_JS;
    } else {
      std::cerr << str << " is not a CriterionType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toProblemString    (ProblemType& problem, const std::string& str) {
    if (str == "MERIT") {
      problem = MERIT;
    } else if (str == "BASIS") {
      problem = BASIS;
    } else if (str == "DUAL") {
      problem = DUAL;
    } else if (str == "REDUCTION") {
      problem = REDUCTION;
    } else if (str == "SHORTEST") {
      problem = SHORTEST;
    } else {
      std::cerr << str << " is not a ProblemType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toNormaString      (NormaType& norma, const std::string& str) {
    if (str == "BESTLAT") {
      norma = BESTLAT;
    } else if (str == "BESTBOUND") {
      norma = BESTBOUND;
    } else if (str == "LAMINATED") {
      norma = LAMINATED;
    } else if (str == "ROGERS") {
      norma = ROGERS;
    } else if (str == "MINK") {
      norma = MINK;
    } else if (str == "MINKL1") {
      norma = MINKL1;
    } else if (str == "L1") {
      norma = L1;
    } else if (str == "L2") {
      norma = L2;
    } else if (str == "NONE") {
      norma = NONE;
    } else {
      std::cerr << str << " is not a NormaType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toCalcString       (CalcType& calc, const std::string& str) {
    if (str == "PAL") {
      calc = PAL;
    } else if (str == "NORMPAL") {
      calc = NORMPAL;
    } else if (str == "BAL") {
      calc = BAL;
    } else if (str == "SEEKPAL") {
      calc = SEEKPAL;
    } else {
      std::cerr << str << " is not a CalcType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toPreRedString     (PreReductionType& prered, const std::string& str) {
    if (str == "FULL") {
      prered = FULL;
    } else if (str == "LLL") {
      prered = LLL;
    } else if (str == "BKZ") {
      prered = BKZ;
    } else if (str == "NOPRERED") {
      prered = NOPRERED;
    } else if (str == "DIETER") {
      prered = DIETER;
    } else {
      std::cerr << str << " is not a PreReductionType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toOutputString     (OutputType& output, const std::string& str) {
    if (str == "TERM") {
      output = TERM;
    } else if (str == "RES") {
      output = RES;
    } else if (str == "TEX") {
      output = TEX;
    } else if (str == "GEN") {
      output = GEN;
    } else {
      std::cerr << str << " is not a OutputType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toPrecisionString  (PrecisionType& precision, const std::string& str) {
    if (str == "DOUBLE") {
      precision = DOUBLE;
    } else if (str == "QUADRUPLE") {
      precision = QUADRUPLE;
    } else if (str == "EXPONENT") {
      precision = EXPONENT;
    } else if (str == "ARBITRARY") {
      precision = ARBITRARY;
    } else if (str == "EXACT") {
      precision = EXACT;
    } else {
      std::cerr << str << " is not a PrecisionType.\n";
      return 1;
    }
    return 0;
  }

}
