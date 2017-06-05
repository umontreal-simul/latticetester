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

#ifndef LATTICETESTER__COORDINATES_H
#define LATTICETESTER__COORDINATES_H

#include <iterator>
#include <set>
#include <map>
#include <iostream>
#include "latticetester/Util.h"


namespace LatticeTester {

/**
 * Set of coordinates.
 */

class Coordinates : public std::set<size_t> {
    public:
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

} //end namespace

#endif
