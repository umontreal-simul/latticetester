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

#ifndef LATTICETESTER__RANDOM_H
#define LATTICETESTER__RANDOM_H

#include <cstdint>
#include <climits>
#define GERME   (ULONG_MAX / 54321)

#include <NTL/ZZ.h>

namespace LatticeTester {

  /**
   * This class generates random numbers (in fact pseudo-random numbers).
   * The generator used is the 64-bits generator `LFSR258`
   * from L'Ecuyer \cite rLEC99a  with period length near \f$2^{258}\f$
   * for 64-bits machines, and the 32-bits generator `LFSR113`
   * from L'Ecuyer \cite rLEC99a  with period length near \f$2^{113}\f$
   * on  32-bits machines. Thus the random numbers generated will be
   * different on 32-bits and 64-bits machines.
   */
  namespace Random {
    /**
     * Returns a random number in \f$[0, 1)\f$. The number has 53 random bits
     * of resolution on 64-bits machines, and 32 random bits
     * on 32-bits machines.
     */
    double randU01();

    /**
     * Return a random integer in \f$[i, j]\f$. The numbers \f$i\f$ and \f$j\f$ can occur.
     * Restriction: \f$i < j\f$. The generated number is limited to 62 bits.
     */
    std::int64_t randInt (std::int64_t i, std::int64_t j);

    /**
     * This is an overload of the previous method to create an NTL::ZZ integer
     * in the range \f$[i,j]\f$. NTL has a method to generate pseudo random
     * numbers in a given range, but it is not deterministic. This method
     * concatenates random bits until it has a big enough number and tries to
     * assign the correct integer to it.
     * */
    NTL::ZZ randInt(NTL::ZZ i, NTL::ZZ j);

    /**
     * Returns random blocks of \f$s\f$ bits (\f$s\f$-bit integers).
     */
    std::uint64_t randBits (int s);

    /**
     * Sets the seed of the generator. If not called, a default seed is used.
     */
    void setSeed (std::uint64_t seed = GERME);
  }

}

#endif

