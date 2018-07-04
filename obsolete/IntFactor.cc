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

#include "latticetester/IntFactor.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>

#include <unistd.h>
#include <stdio.h>

#include "NTL/ZZ.h"


namespace LatticeTester
{


  std::string IntFactor::toString () const
  {
    char c;

    switch (m_status) {
      case PRIME:
        c = 'P';
        break;
      case PROB_PRIME:
        c = 'Q';
        break;
      case COMPOSITE:
        c = 'C';
        break;
      default:
        c = 'U';
        break;
    }

    std::ostringstream sortie;
    sortie << m_factor << std::setw(10) << m_multiplicity << std::setw(10) << c;
    return sortie.str ();
  }


  //===========================================================================

  inline
    std::string IntFactor::toString (PrimeType stat)
    {
      return toStringPrime (stat);
    }


  //===========================================================================

  PrimeType IntFactor::isPrime (const MScal & y, long k)
  {
    static constexpr unsigned int NB_PRIMES = 6543;
    // NbPrem has to be instanciated if we use NTL types
    MScal NbPrem;
    NbPrem = 2;
    NTL::ZZ LIM;
    LIM = NTL::conv<ZZ>("4295098369");
    //    unsigned int c[200];
    static const std::array<unsigned int, NB_PRIMES> primes =
    {{
#include "../data/prime.dat"
     }};
    // std::ifstream in ("../share/latticetester/data/prime.dat"); // contains primes < 2^16
    //    if (!(in.is_open())) {
    //       std::cerr << "Error:   cannot open file   prime.dat\n";
    //       exit(8);
    //    }
    //    in.ignore (100, '\n');  // ignore first line = number of primes in file

    MScal ys = NTL::SqrRoot (y);
    //    if (ys > 65536)
    //       ys = 65536;
    //    // The condition NbPrem < 65536 is necessary because there seems to
    //    // be a bug in ZZ input: does not detect end-of-file and gives an error.
    unsigned int i = 1;
    while (i < NB_PRIMES && (NbPrem <= ys)) {
      if (y % NbPrem == 0)
        return COMPOSITE;
      NbPrem = primes[i];
      i++;
    }
    if (y <= LIM)
      return PRIME;

    /* A n'est divisible par aucun des nombres premiers inf. a 2^16. */
    return isProbPrime (y, k);
  }

  PrimeType IntFactor::isPrime (long k)
  {
    return isPrime (m_factor, k);
  }


  //===========================================================================

  PrimeType IntFactor::isProbPrime (const MScal & y, long k)
  {
    PrimeType stat;
    long res = NTL::ProbPrime (y, k);

    switch (res) {
      case 0:
        stat = COMPOSITE;
        break;
      case 1:
        stat = PROB_PRIME;
        break;
      default:
        stat = UNKNOWN;
    }
    return stat;
  }

  //===========================================================================

}                                 // namespace LatticeTester
