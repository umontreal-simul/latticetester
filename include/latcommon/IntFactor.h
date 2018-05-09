// This file is part of LatCommon.
//
// LatCommon
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

#ifndef LATCOMMON__INTFACTOR_H
#define LATCOMMON__INTFACTOR_H
#include <string>
#include "latcommon/Types.h"
#include "latcommon/Const.h"
 

namespace LatCommon {

/**
 * The objects of this class are the "prime" factors in the decomposition of
 * a positive integer. The class contains functions to determine whether a
 * number is prime, probably prime or composite. The related class
 * `IntFactorization` contains the list of "prime" factors. (This class
 * should be a private nested class inside `IntFactorization`.)
 *
 * <div class="LatSoft-bigskip"></div>
 */
class IntFactor {
public:

   /**
    * Constructor for a factor \f$x\f$, with multiplicity `mult` and status
    * `stat`.
    */
   IntFactor (const MScal & x, int mult = 1,
                 LatCommon::PrimeType stat = LatCommon::UNKNOWN):
                 m_factor (x), m_multiplicity (mult), m_status (stat) {   }

   /**
    * Returns the value of `factor`.
    */
   MScal getFactor () const { return m_factor; }

   /**
    * Sets the value of `factor` to \f$x\f$.
    */
   void setFactor (const MScal & x) { m_factor = x; }

   /**
    * Returns the multiplicity of this object.
    */
   int getMultiplicity () const { return m_multiplicity; }

   /**
    * Sets the multiplicity of this object to \f$m\f$.
    */
   void setMultiplicity (int m) { m_multiplicity = m; }

   /**
    * Returns the status of this object.
    */
   LatCommon::PrimeType getStatus () const { return m_status; }

   /**
    * Sets the status of this object to \f$s\f$.
    */
   void setStatus (LatCommon::PrimeType s) { m_status = s; }

   /**
    * Tests whether \f$y\f$ is prime. First tests whether \f$y\f$ is
    * divisible by all small primes \f$p\f$ (\f$p < 2^{16}\f$) that are
    * kept in file `prime.dat`. Then applies the Miller-Rabin probability
    * test with \f$k\f$ trials.
    */
   static LatCommon::PrimeType isPrime (const MScal & y, long k);

   /**
    * Tests whether this factor is prime. Similar to `isPrime` above.
    */
   LatCommon::PrimeType isPrime (long k);

   /**
    * Transforms status `stat` in an easily readable string and returns
    * it.
    */
   static std::string toString (LatCommon::PrimeType stat);

   /**
    * Returns this object as a string.
    */
   std::string toString () const;

private:

   /**
    * The factor.
    */
   MScal m_factor;

   /**
    * The multiplicity of this factor.
    */
   int m_multiplicity;

   /**
    * The status of this factor, i.e. whether it is prime, composite, ...
    */
   LatCommon::PrimeType m_status;

#ifdef WITH_NTL
   /**
    * Applies the Miller-Rabin probability test with \f$k\f$ trials to
    * \f$y\f$.
    */
   static LatCommon::PrimeType isProbPrime (const MScal & y, long k);
#endif
};    // class IntFactor

}     // namespace LatMRG
#endif
