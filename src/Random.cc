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

#include "latticetester/Random.h"
#include <climits>


#define GERME   (ULONG_MAX / 54321)
#define NORM53  1.11022302462515654e-16  /* 1/2^53 */
#define NORM32  2.3283064365386963e-10   /* 1/2^32 */
#define MASK53  0x1fffffffffffffUL       /* 2^53 - 1 */

#if ULONG_MAX > 4294967295UL
#define ULONG_64_OUI
#endif

namespace {

  std::uint64_t etat1=GERME, etat2=GERME, etat3=GERME, etat4=GERME, etat5=GERME;

    std::uint64_t randValue ()
    {
      std::uint64_t b;

#ifdef ULONG_64_OUI
      // Générateur LFSR258 de L'Ecuyer
      b = ((etat1 << 1) ^ etat1) >> 53;
      etat1 = ((etat1 & 18446744073709551614UL) << 10) ^ b;
      b = ((etat2 << 24) ^ etat2) >> 50;
      etat2 = ((etat2 & 18446744073709551104UL) << 5) ^ b;
      b = ((etat3 << 3) ^ etat3) >> 23;
      etat3 = ((etat3 & 18446744073709547520UL) << 29) ^ b;
      b = ((etat4 << 5) ^ etat4) >> 24;
      etat4 = ((etat4 & 18446744073709420544UL) << 23) ^ b;
      b = ((etat5 << 3) ^ etat5) >> 33;
      etat5 = ((etat5 & 18446744073701163008UL) << 8) ^ b;
      return (etat1 ^ etat2 ^ etat3 ^ etat4 ^ etat5);

#else
      // Générateur LFSR113 de L'Ecuyer
      b  = ((etat1 << 6) ^ etat1) >> 13;
      etat1 = ((etat1 & 4294967294UL) << 18) ^ b;
      b  = ((etat2 << 2) ^ etat2) >> 27;
      etat2 = ((etat2 & 4294967288UL) << 2) ^ b;
      b  = ((etat3 << 13) ^ etat3) >> 21;
      etat3 = ((etat3 & 4294967280UL) << 7) ^ b;
      b  = ((etat4 << 3) ^ etat4) >> 12;
      etat4 = ((etat4 & 4294967168UL) << 13) ^ b;
      return (etat1 ^ etat2 ^ etat3 ^ etat4);

#endif
    }
}

namespace LatticeTester {
  namespace Random {

    std::int64_t randInt (std::int64_t i, std::int64_t j)
    {
#ifdef ULONG_64_OUI
      std::uint64_t d = j - i + 1;
      std::uint64_t q = 0x4000000000000000L / d;
      std::uint64_t r = 0x4000000000000000L % d;
      std::uint64_t v = 0x4000000000000000L - r;
      std::uint64_t res;

      do {
        res = randValue() >> 2;
      } while (res >= v);

      return i + (std::int64_t) (res / q);

#else
      // moins précise que ci-dessus
      double res = randU01();
      int d = j - i + 1;
      return i + static_cast<long>(d * res);
#endif
    }

    //=========================================================================

    NTL::ZZ randInt (NTL::ZZ i, NTL::ZZ j)
    {
      NTL::ZZ d = j - i + 1;
      long numbits = NTL::NumBits(d);
      NTL::ZZ maxvalue = NTL::power2_ZZ(numbits);
      NTL::ZZ q = maxvalue / d;
      NTL::ZZ r = maxvalue % d;
      NTL::ZZ v = maxvalue - r;
      NTL::ZZ res;
      do {
        long leftbits = numbits;
#ifdef ULONG_64_OUI
        // We have to go 63 bits at a time because the conversion is between long
        // and ZZ, not unsigned long and ZZ
        while (leftbits > 63) {
          res = (res << 63) + NTL::ZZ(randBits(63));
          leftbits -= 63;
        }
#else
        while (leftbits > 31) {
          res = res << 31 + NTL::ZZ(randBits(31));
          leftbits -= 31;
        }
#endif
        res = (res << leftbits) + NTL::ZZ(randBits(leftbits));
      } while (res >= v);

      return i + (res / q);
    }

    //============================================================================

    double randU01 ()
    {
#ifdef ULONG_64_OUI
      //   return NORM53 * (randValue () >> 11);
      return NORM53 * (randValue () & MASK53);
#else
      return NORM32 * randValue ();
#endif
    }


    //=========================================================================

    uint64_t randBits (int s)
    {
#ifdef ULONG_64_OUI
      return randValue () >> (64 - s);
#else
      return randValue () >> (32 - s);
#endif
    }


    //=========================================================================

    void setSeed (uint64_t seed)
    {
      // Set one high bit = 1 to make sure initial state is valid
      etat1 = seed | 0x40000000;
      etat2 = seed | 0x40000000;
      etat3 = GERME;
      etat4 = GERME;
      etat5 = GERME;
    }


    //========================================================================

  }
}
