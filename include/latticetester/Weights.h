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

#ifndef LATTICETESTER__WEIGHTS_H
#define LATTICETESTER__WEIGHTS_H

#include <string>

#include "latticetester/Coordinates.h"

namespace LatticeTester {

  /**
   * A scalar weight type.
   */
  typedef double Weight;

  /**
   * Abstract class that defines an interface to specify Weights given to projections
   * in figures of merit. This class and its subclasses are not used directly by
   * LatticeTester since LatticeTester does not implement figures of merit.
   * For examples of usages, one may look at LatNetBuilder and LatMRG.
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
