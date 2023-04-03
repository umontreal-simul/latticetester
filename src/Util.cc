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

#include "latticetester/Util.h"
#include <cstdlib>
#include <cmath>
#include <climits>

#include "NTL/ZZ.h"
using namespace NTL;

namespace
{

#define GERME  9876543219876543L
#define NORM53  1.11022302462515654e-16   // 1 / 2^53
#define NORM52  2.22044604925031308e-16   // 1 / 2^52


  std::uint64_t etat1 = GERME,
                etat2 = GERME, etat3 = GERME, etat4 = GERME, etat5 = GERME;


  // This is the integer version of LFSR258 available at
  // http://simul.iro.umontreal.ca/ 
  std::uint64_t RandValue ()
  {
    std::uint64_t b;

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
  }

}  // namespace

//============================================================================


namespace LatticeTester
{

  const std::int64_t TWO_EXP[64] = {
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
    65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216,
    33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648,
    4294967296, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472,
    274877906944, 549755813888, 1099511627776, 2199023255552, 4398046511104,
    8796093022208, 17592186044416, 35184372088832, 70368744177664,
    140737488355328, 281474976710656, 562949953421312, 1125899906842624,
    2251799813685248, 4503599627370496, 9007199254740992, 18014398509481984,
    36028797018963968, 72057594037927936, 144115188075855872, 288230376151711744,
    576460752303423488, 1152921504606846976, 2305843009213693952,
    (std::int64_t)4611686018427387904U, (std::int64_t)9223372036854775808U };



  //============================================================================
  // binary GCD from wikipedia

  std::int64_t gcd (std::int64_t u, std::int64_t v)
  {
    int shift;
    if (u < 0) u = -u;
    if (v < 0) v = -v;

    // GCD(0,x) := x
    if (u == 0 || v == 0)
      return u | v;

    /* (u | v) & 1 is 1 if u or v is odd. This divides u and v by two while
     * both are even and stores the number of divisions in shift.*/
    for (shift = 0; ((u | v) & 1) == 0; ++shift) {
      u >>= 1;
      v >>= 1;
    }

    // Since at most one of u and v is even we can divide u by two if it is
    while ((u & 1) == 0)
      u >>= 1;

    // From now on, u is always odd.
    do {
      // We make v odd too since u is odd and gcd(u,v) = gcd(u,v/2) if v is even
      while ((v & 1) == 0)     // Loop X
        v >>= 1;

      /* Now u and v are both odd, so diff(u, v) is even. Let u = min(u,
         v), v = diff(u, v)/2. */
      if (u < v) {
        v -= u;
      } else {
        long diff = u - v;
        u = v;
        v = diff;
      }
      v >>= 1;
    } while (v != 0);

    return u << shift;
  }

  //===========================================================================

  void MyExit (int status, std::string msg)
  {
    std::cout << "\n***** Error " << msg << std::endl;
    exit (status);

  }


  //===========================================================================

  std::int64_t Factorial(int t)
  {
    if (t == 0 || t == 1)
      return 1;
    long fact = 1;
    for (int i = 2;i <= t;i++) {
      fact *= i;
    }
    return fact;
  }


  //===========================================================================
  // Générateur LFSR258 de L'Ecuyer

  std::int64_t RandInt (std::int64_t i, std::int64_t j)
  {
    // The range covered by the interval
    std::uint64_t d = j-i+1;
    // Gap between numbers of the output. That is, if we generate a number
    // in (0,q) with RandValue(), we will get i and in [q,2q) we will get i+1
    // and so on.
    std::uint64_t q = 0x4000000000000000L / d;
    // These two numbers make a correction to avoid a bias. Basically if we want
    // to generate an integers with the inverse density, we need to select each
    // integer with the same probaility q. Because of the rounding in
    // integer division d*q is not always 1. What we have instead is that q*v=d.
    // Because of that, we do not want to consider RandValue() bigger than v.
    std::uint64_t r = 0x4000000000000000L % d;
    std::uint64_t v = 0x4000000000000000L - r;
    std::uint64_t res;

    do {
      res = RandValue() >> 2;
    } while (res >= v);

    return i + (std::uint64_t) (res / q);
  }

  // Cuts the last bits of LFSR258 64 bits int and converts it to a double
  // precision number
  double RandU01()
  {
    return NORM53 * (RandValue() >> 11);
  }


  std::uint64_t RandBits (int s)
  {
    return RandValue () >> (64 - s);
  }

  void SetSeed (std::uint64_t seed)  {
    // Choose one bit = 1 to make sure initial state is valid
    /*
     * etat1 = seed | 0x40000000;
     * etat3 = seed | 0x40000000;
     */
    // Making sure the seed is valid
    etat1 = seed | 2;
    etat2 = seed | 512;
    etat3 = seed | 4096;
    etat4 = seed | 131072;
    etat5 = seed | 8388608;
  }


  //===========================================================================

}        // namespace LatticeTester
