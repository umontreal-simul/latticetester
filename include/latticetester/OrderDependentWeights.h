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

#ifndef LATTICETESTER__ORDER_DEPENDENT_WEIGHTS_H
#define LATTICETESTER__ORDER_DEPENDENT_WEIGHTS_H

#include "latticetester/Weights.h"
#include <vector>
#ifdef WITH_XML
#include <pugixml.hpp>
#endif


namespace LatticeTester {
  /**
   * Order-dependent weights.
   *
   * This class implements order dependent  weights.
   * The weight of a projection depends only on the order of the projection,
   * i.e. the dimension of the projection.
   */
  class OrderDependentWeights : public Weights {
    protected:

      Weight m_defaultWeight;
      std::vector<Weight> m_weights;

    public:

      /**
       * Constructs order-dependent weights with default weight.
       *
       * \param defaultWeight   Default weight.
       */
      explicit OrderDependentWeights (Weight defaultWeight = 0.0);

      /**
       * Destructor.
       */
      virtual ~OrderDependentWeights()
      { }

      virtual std::string name() const { return "order-dependent"; }

      /**
       * Returns the weight of the projection specified by \c projection.
       */
      virtual Weight getWeight (const Coordinates & projection) const;

      /**
       * Returns the weight associated to the given order.
       */
      virtual Weight getWeightForOrder (Coordinates::size_type order) const
      { return order < m_weights.size() ? m_weights[order] : m_defaultWeight; }

      /**
       * Sets the weight for the order specified by \c order.
       */
      virtual void setWeightForOrder (Coordinates::size_type order, Weight weight);

      /**
       * Sets the default weight of all orders for which a weight
       * has not been set explicitly set using #setWeightForOrder().
       */
      virtual void setDefaultWeight (Weight weight)
      { m_defaultWeight = weight; }

      virtual Weight getDefaultWeight () const 
      { return m_defaultWeight; } 

      virtual unsigned int getSize () const 
      { return (unsigned int) m_weights.size(); } 

#ifdef WITH_XML
      /**
       * Static factory method; create a \c OrderDependentWeights object by
       * parsing XML data.
       */
      static OrderDependentWeights* createFromXML (const pugi::xml_node & node);
#endif

    protected:
      /// \copydoc LatticeTester::Weights::format()
      virtual void format(std::ostream& os) const;
  };

}

#endif
