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

#ifndef LATTICETESTER_INTFACTOR_H
#define LATTICETESTER_INTFACTOR_H

#include <string>
#include <iomanip>
#include <cstdint>
#include <sstream>

#include "latticetester/Const.h"

#include <NTL/ZZ.h>

namespace LatticeTester {

  /**
   * The objects of this class are the "prime" factors in the decomposition of
   * a positive integer.
   * The class also contains functions to determine whether a number is prime,
   * probably prime or composite. These methods can be used externally to test
   * the primality of an integer, or to test if this factor is prime.
   */
  template<typename Int>
    class IntFactor {
      public:

        /**
         * Constructor for a factor \f$x\f$ of multiplicity `mult` that has a 
         * PrimeType `stat`
         */
        IntFactor (const Int & x, int mult = 1,
            PrimeType stat = UNKNOWN):
          m_factor (x), m_multiplicity (mult), m_status (stat) {   }

        /**
         * Returns the numeric value of this factor.
         */
        Int getFactor () const { return m_factor; }

        /**
         * Sets the value of this factor to \f$x\f$.
         */
        void setFactor (const Int & x) { m_factor = x; }

        /**
         * Returns the multiplicity of this object.
         */
        int getMultiplicity () const { return m_multiplicity; }

        /**
         * Sets the multiplicity of this object to \f$m\f$.
         */
        void setMultiplicity (int m) { m_multiplicity = m; }

        /**
         * Returns the PrimeType of this object.
         */
        PrimeType getStatus () const { return m_status; }

        /**
         * Sets the PrimeType of this object to \f$s\f$.
         */
        void setStatus (PrimeType s) { m_status = s; }

        /**
         * Tests whether \f$y\f$ is prime. First tests whether \f$y\f$ is
         * divisible by all small primes \f$p\f$ (\f$p < 2^{16}\f$) that are
         * kept in file `prime.dat`. Then applies the Miller-Rabin probability
         * test with \f$k\f$ trials.
         */
        static PrimeType isPrime (const Int & y, std::int64_t k);

        /**
         * Tests whether this factor is prime. Similar to `isPrime` above.
         */
        LatticeTester::PrimeType isPrime (std::int64_t k);

        /**
         * Transforms status `stat` in an easily readable string and returns
         * it.
         */
        static std::string toString (LatticeTester::PrimeType stat);

        /**
         * Returns this object as a string.
         */
        std::string toString () const;


  

      private:

        /**
         * The factor.
         */
        Int m_factor;

        /**
         * The multiplicity of this factor.
         */
        int m_multiplicity;

        /**
         * The status of this factor, i.e. whether it is prime, composite, ...
         */
        LatticeTester::PrimeType m_status;

        /**
         * Applies the Miller-Rabin probability test with \f$k\f$ trials to
         * \f$y\f$.
         */
        static LatticeTester::PrimeType isProbPrime (const Int & y, std::int64_t k);
    };    // class IntFactor

  //===========================================================================


  template<typename Int>
    std::string IntFactor<Int>::toString () const
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
      sortie << m_factor << std::setw(10) << m_multiplicity << std::setw(10) 
        << c;
      return sortie.str ();
    }


  //===========================================================================

  template<typename Int>
    inline
    std::string IntFactor<Int>::toString (PrimeType stat)
    {
      return toStringPrime (stat);
    }

  //===========================================================================

  template<typename Int>
    PrimeType IntFactor<Int>::isPrime (const Int & y, std::int64_t k)
    {
      // NbPrem has to be instantiated if we use NTL types
      Int NbPrem;
      NbPrem = 2;
      NTL::ZZ LIM;
      LIM = NTL::conv<NTL::ZZ>("4295098369");

      Int ys = NTL::SqrRoot (y);
      unsigned int i = 1;
      while (i < NB_PRIMES && (NbPrem <= ys)) {
        if (y % NbPrem == 0)
          return COMPOSITE;
        NbPrem = PRIMES_ARRAY[i];
        i++;
      }
      if (y <= LIM)
        return PRIME;

      /* A is divisible by none of the prime numbers smaller than 2^16. */
      return isProbPrime (y, k);
    }

  //===========================================================================

  template<typename Int>
    PrimeType IntFactor<Int>::isPrime (std::int64_t k)
    {
      return isPrime (m_factor, k);
    }


  //===========================================================================

  template<typename Int>
    PrimeType IntFactor<Int>::isProbPrime (const Int & y, std::int64_t k)
    {
      PrimeType stat;
      std::int64_t res = NTL::ProbPrime (y, k);

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

  extern template class IntFactor<NTL::ZZ>;
  extern template class IntFactor<std::int64_t>;

}     // namespace LatticeTester

#endif
