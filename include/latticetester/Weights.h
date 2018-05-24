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

#ifndef LATTICETESTER__WEIGHTS_H
#define LATTICETESTER__WEIGHTS_H

#include <string>
#include "latticetester/CoordinateSets.h"


namespace LatticeTester {

/**
 * Scalar weight type.
 *
 * \note We could have used \c Weight, but it might be wise to leave this \c
 * typedef in case we decide to use <tt>long Weight</tt> at some point.
 */
typedef double Weight;

/**
 * Abstract weights class.
 *
 * This abstract class is the basis for different kinds of weights used to
 * accentuate the importance given to some projections when computing
 * figures of merit for lattices or point sets.
 */
class Weights {
public:

   /**
    * Destructor.
    */
   virtual ~Weights()
   { }

   /**
    * Returns the weight of the projection specified by `projection`.
    */
   virtual Weight getWeight (const Coordinates & projection) const = 0;

   virtual std::string name() const = 0;

protected:
   /**
    * Identifies the type of weights, formats them and outputs them on \c os.
    *
    * \remark Deriving classes should identify themselves in the output.
    */
   virtual void format(std::ostream& os) const = 0;

   friend std::ostream & operator<< (std::ostream & out, const Weights & o);
};


/**
 * \relates Weights
 * Identifies the type of weights, formats them and outputs them on \c os.
 */
inline std::ostream & operator<< (std::ostream & os, const Weights & o)
{ o.format(os); return os; }

}

#endif
