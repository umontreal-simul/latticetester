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

#ifndef LATTICETESTER__UNIFORM_WEIGHTS_H
#define LATTICETESTER__UNIFORM_WEIGHTS_H

#include "latticetester/Weights.h"

#include <sstream>
#include <vector>

namespace LatticeTester {
  /**
   * This class is used to implement the same weight for all projections. It
   * represents the trivial case of no weight. The weights can all be chosen as
   * 1.
   */
  class UniformWeights : public Weights {
    protected:

      Weight m_weight;

    public:

      /**
       * Constructs uniform weights.
       *
       * \param weight     Weight for all projections.
       */
      explicit UniformWeights (Weight weight)
      { m_weight = weight; }

      /**
       * Destructor.
       */
      virtual ~UniformWeights()
      { }

      /**
       * Returns the same weight regardless of the specified indices.
       */
      virtual Weight getWeight (const Coordinates &) const  { return m_weight; }

    protected:
      /// \copydoc LatticeTester::Weights::format()
      virtual void format(std::ostream& os) const {
        os << "UniformWeights(" << m_weight << ")";
      }
  };

}

#endif

