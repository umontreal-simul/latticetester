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

#ifndef LATTICETESTER__UNIFORM_WEIGHTS_H
#define LATTICETESTER__UNIFORM_WEIGHTS_H

#include "latticetester/Weights.h"

#include <sstream>
#include <vector>

namespace LatticeTester {
  /**
   * Specifies projection weights that are the same (usually 1) for all projections.
   */
  class WeightsUniform : public Weights {
    protected:

      Weight m_weight;

    public:

      /**
       * Constructs uniform weights.
       *
       * \param weight     Weight given to all projections.
       */
      explicit WeightsUniform (Weight weight)
      { m_weight = weight; }

      /**
       * Destructor.
       */
      virtual ~WeightsUniform()
      { }

      /**
       * Returns the same weight regardless of the specified indices.
       */
      virtual Weight getWeight (const Coordinates &) const  { return m_weight; }

    protected:
      /// \copydoc LatticeTester::Weights::format()
      virtual void format(std::ostream& os) const {
        os << "WeightsUniform(" << m_weight << ")";
      }
  };

}

#endif

