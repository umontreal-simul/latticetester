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

#ifndef LATTICETESTER__POD_WEIGHTS_H
#define LATTICETESTER__POD_WEIGHTS_H

#include "latticetester/Weights.h"
#include "latticetester/WeightsOrderDependent.h"
#include "latticetester/WeightsProduct.h"
#include <vector>

namespace LatticeTester {
  /**
   * Defines product and order-dependent (POD) weights.
   * The weight of a projection is the sum of a product weight and an order-dependent weight.
   */
  class WeightsPOD : public Weights {
    protected:

      WeightsOrderDependent m_orderDependentWeights;
      WeightsProduct m_productWeights;

    public:

      /**
       * Constructs POD weights with default weight.
       */
      WeightsPOD();

      /**
       * Destructor.
       */
      virtual ~WeightsPOD()
      { }
      
      /**
       * Returns the weight of the projection specified by `projection`.
       */
      virtual Weight getWeight (const Coordinates & projection) const;

      /**
       * Returns the order-dependent part of the weights.
       */
      WeightsOrderDependent& getOrderDependentWeights()
      { return m_orderDependentWeights; }

      const WeightsOrderDependent& getOrderDependentWeights() const
      { return m_orderDependentWeights; }

      /**
       * Returns the product part of the weights.
       */
      WeightsProduct& getProductWeights()
      { return m_productWeights; }

      const WeightsProduct& getProductWeights() const
      { return m_productWeights; }

      double getWeightForOrder (Coordinates::size_type order) const
    {
        return getOrderDependentWeights().getWeightForOrder(order);
    }

    double getWeightForCoordinate (Coordinates::size_type coordinate) const
    {
        return getProductWeights().getWeightForCoordinate(coordinate);
    }

// #ifdef WITH_XML
//       /**
//        * Static factory method; create a \c WeightsPOD object by
//        * parsing XML data.
//        */
//       static WeightsPOD* createFromXML (const pugi::xml_node & node);
// #endif

    protected:
      /// \copydoc LatticeTester::Weights::format()
      virtual void format(std::ostream& os) const;
  };

}

#endif
