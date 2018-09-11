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

#ifndef LATTICETESTER__PRODUCT_WEIGHTS_H
#define LATTICETESTER__PRODUCT_WEIGHTS_H

#include "latticetester/Weights.h"
#include <vector>
// #ifdef WITH_XML
// #include <pugixml.hpp>
// #endif


namespace LatticeTester {
  /**
   * Product weights.
   *
   * The weight of a projection is equal to the product of the individual weights
   * of each coordinate involved.
   */
  class ProductWeights : public Weights {
    protected:

      Weight m_defaultWeight;
      std::vector<Weight> m_weights;

    public:

      /**
       * Constructs projection-dependent weights with a default weight.
       *
       * \param defaultWeight    Default weight.
       */
      explicit ProductWeights (Weight defaultWeight = 0.0);

      /**
       * Destructor.
       */
      virtual ~ProductWeights()  {}

      /**
       * Returns the weight of the projection specified by \c projection.
       */
      virtual Weight getWeight (const Coordinates & projection) const;

      /**
       * Returns the weight associated to the given coordinate.
       */
      virtual Weight getWeightForCoordinate (
          Coordinates::size_type coordinate) const
      {
        return coordinate < m_weights.size() ?
          m_weights[coordinate] : m_defaultWeight;
      }

      /**
       * Sets the weight for the coordinate specified by \c coordinate.
       */
      virtual void setWeightForCoordinate (
          Coordinates::size_type coordinate, Weight weight);

      /**
       * Sets the default weight of all coordinates for which a weight
       * has not been set explicitly set using #setWeightForCoordinate().
       */
      virtual void setDefaultWeight (Weight weight)
      { m_defaultWeight = weight; }

      virtual const std::vector<Weight>& getWeights() const
      {
        return m_weights;
      }

      virtual Weight getDefaultWeight() const
      {
        return m_defaultWeight;
      }

    /**
    * Multiply each weight by a constant.
    */
   void multiplyWeights(Weight c){
       for (auto& w: m_weights){
           w *= c;
       }
       m_defaultWeight *= c;
   }

// #ifdef WITH_XML
//       /**
//        * Static factory method; create a \c ProductWeights object by
//        * parsing XML data.
//        */
//       static ProductWeights* createFromXML (const pugi::xml_node & node);
// #endif

    protected:
      /// \copydoc LatticeTester::Weights::format()
      virtual void format(std::ostream& os) const;
  };

}

#endif
