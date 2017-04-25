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

#ifndef LATCOMMON__COORDINATES_H
#define LATCOMMON__COORDINATES_H

#include <iterator>
#include <set>
#include <map>
#include <iostream>
#include "latcommon/Util.h"


namespace LatCommon {

/**
 * Set of coordinates.
 *
 * \remark Coordinates indices start at 0.  To make the input/output
 * representation of coordinates start at 1 instead of 0, see #humanize.  For
 * code written with coordinate indices starting at 1, set #humanize to \c false.
 */
class Coordinates : public std::set<size_t> {
public:
   /**
    * &ldquo;Humanize&rdquo; the formatting of coordinate values.
    *
    * If set to \c true, the external representation of coordinate values if
    * shifted by one with respect to its internal representation.  More
    * precisely, an internal coordinate value \f$j\f$ is mapped to the external
    * representation \f$j+1\f$ during output, and an external coordinate value
    * \f$j+1\f$ is mapped to internal representation \f$j\f$ during input.
    *
    * Defaults to \c true;
    *
    * \sa #asOutput() #asInput()
    */
   static bool humanize;

   /**
    * Maps the internal representation of a coordinate value to its external
    * representation.
    */
   inline static value_type asOutput(value_type i)
   { return i + (humanize ? 1 : 0); }

   /**
    * Maps the external representation of a coordinate value to its internal
    * representation.
    */
   inline static value_type asInput(value_type i)
   { return i - (humanize ? 1 : 0); }

   /**
    * Constructs an empty coordinate set.
    */
   Coordinates():
      std::set<value_type>()
   { }

   /**
    * Copy-constructor.
    */
   Coordinates(const Coordinates& other):
        std::set<value_type>(other)
   { }

   /**
    * Constructs a coordinate set populated with the values from \c first
    * (inclusively) to \c last (exclusively).
    */
   template<typename InputIterator>
     Coordinates(InputIterator first, InputIterator last):
        std::set<value_type>(first, last)
   { }
};


/**
 * \relates Coordinates
 * Formats the coordinate set \c coords and outputs it to \c os.
 */
std::ostream& operator<< (std::ostream& os, const Coordinates& coords);

/**
 * \relates Coordinates
 * Reads a formatted coordinate set from \c is.
 *
 * The input must consist of positive integers separated by whitespace and/or by
 * commas, and optionally enclosed in braces.  The ordering is not important.
 * Repeated values are ignored.
 * For example, the following strings are valid input that would produce
 * equivalent Coordinates objects:
 * - <tt>1 2 5</tt>
 * - <tt>1, 2, 5</tt>
 * - <tt>{1 2 5}</tt>
 * - <tt>{1,2,5}</tt>
 * - <tt>{1, 2, 5}</tt>
 * - <tt>2 5 1</tt>
 * - <tt>2 1 5 1</tt>
 */
std::istream& operator>> (std::istream& is, Coordinates& coords);

}
#endif
