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


namespace LatticeTester {
  /**
   *
   *
   * This class generates random numbers (in fact pseudo-random numbers).
   * The generator used is the 64-bits generator \c LFSR258
   * from L'Ecuyer \cite rLEC99a  with period length near \f$2^{258}\f$
   * for 64-bits machines, and the 32-bits generator \c LFSR113
   * from L'Ecuyer \cite rLEC99a  with period length near \f$2^{113}\f$
   * on  32-bits machines. Thus the random numbers generated will be
   * different on 32-bits and 64-bits machines.
   *
   *
   */
  class Random {
    public:

      /**
       * Constructor using a default seed for the generator.
       * One may reset the seed by calling \c setSeed.
       */
      Random();

      /**
       * Destructor.
       */
      ~Random()  {}

      /**
       * Returns a random number in \f$[0, 1)\f$. The number has 53 random bits
       * of resolution on 64-bits machines, and 32 random bits
       * on 32-bits machines.
       */
      double randU01();

      /**
       * Return a random integer in \f$[i, j]\f$. The numbers \f$i\f$ and \f$j\f$ can occur.
       * Restriction: \f$i < j\f$.
       */
      int randInt (int i, int j);

      /**
       * Returns random blocks of \f$s\f$ bits (\f$s\f$-bit integers).
       */
      std::uint64_t randBits (int s);

      /**
       * Sets the seed of the generator. If not called, a default seed is used.
       */
      void setSeed (std::uint64_t seed);


    private:

      std::uint64_t randValue();

      std::uint64_t etat1, etat2, etat3, etat4, etat5;
  };

}

#endif

