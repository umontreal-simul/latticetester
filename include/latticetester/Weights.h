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

#ifndef LATTICETESTER__WEIGHTS_H
#define LATTICETESTER__WEIGHTS_H

#include <string>

#include "latticetester/Coordinates.h"

namespace LatticeTester {

  /**
   * \relates Weights
   * Scalar weight type.
   *
   * \note We could have used \c Weight, but it might be wise to leave this \c
   * typedef in case we decide to use <tt>std::int64_t Weight</tt> at some point.
   */
  typedef double Weight;

  /**
   * Abstract class representing Weights for figures of merit. Typically, if one
   * wants to analyze a lattice, there are multiple figures of merit available.
   * This class presents an interface to give weights to all those figures of
   * merit so that a global figure of merit can be calculated as a weighted mean
   * of many figures of merit.
   *
   * This abstract class is the basis for different kinds of weights used to
   * accentuate the importance of some projections when computing
   * figures of merit for lattices or point sets. For more insight on what can
   * be done with weights, look what the subclasses implement. These classes are
   * not used directly by LatticeTester since LatticeTester does not implement
   * calculations of figure of merit for different projections of a point set.
   * For example usages, one should look at the LatNetBuilder and the LatMRG
   * software.
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

      /**
       * Returns the interlacing factor of the weights. This is used in LatNet Builder
       * to parameterize the figures of merit for interlaced digital nets.
       */ 
      virtual unsigned int interlacingFactor() const
      {
        return 1;
      }

    protected:
      /**
       * Identifies the type of weights, formats them and outputs them on \c os.
       *
       * Subclasses that implement Weights should identify themselves in the output.
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
