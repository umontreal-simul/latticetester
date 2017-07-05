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

#ifdef WITH_NTL
#include "NTL/ZZ.h"
#endif


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


#ifdef WITH_NTL

PrimeType IntFactor::isPrime (const MScal & y, long k)
{
   MScal NbPrem;
   NTL::ZZ LIM;
   LIM = NTL::conv<ZZ>("4294967295");
   std::ifstream in ("/u/simardr/latmrg/prime.dat"); // contains primes < 2^16
   if (!(in.is_open())) {
      std::cerr << "Error:   cannot open file   prime.dat\n";
      exit(8);
   }
   in.ignore (100, '\n');  // ignore first line = number of primes in file

   MScal ys = NTL::SqrRoot (y);
   if (ys > 65536)
      ys = 65536;
   // The condition NbPrem < 65536 is necessary because there seems to
   // be a bug in ZZ input: does not detect end-of-file and gives an error.
   while ((in >> NbPrem) && (NbPrem <= ys)) {
      if (y % NbPrem == 0)
         return COMPOSITE;
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

#else

PrimeType IntFactor::isPrime (const MScal& n, long)
{
   if (n == 2)
      return PRIME;
   if ((n & 1) == 0) //means that n is even
      return COMPOSITE;
   long factmax = static_cast<long>(std::sqrt (n));
   for (long i = 3; i <= factmax; i += 2)
      if (n % i == 0)
         return COMPOSITE;
   return PRIME;
}

#endif


//===========================================================================

#ifdef WITH_NTL
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
#endif

//===========================================================================

}                                 // namespace LatMRG
