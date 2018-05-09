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

#ifndef LATCOMMON__POD_WEIGHTS_H
#define LATCOMMON__POD_WEIGHTS_H

#include "latcommon/Weights.h"
#include "latcommon/OrderDependentWeights.h"
#include "latcommon/ProductWeights.h"
#include <vector>
#ifdef WITH_XML
#include <pugixml.hpp>
#endif


namespace LatCommon {
/**
 * Product and order-dependent (POD) weights.
 *
 * This class implements POD weights.
 * The weight of a projection depends on the order of the projection and on the
 * coordinates that are in the projection.
 */
class PODWeights : public Weights {
protected:

   OrderDependentWeights m_orderDependentWeights;
   ProductWeights m_productWeights;

public:

   /**
    * Constructs order-dependent weights with default weight.
    */
   PODWeights();

   /**
    * Destructor.
    */
   virtual ~PODWeights()
   { } 

   /**
    * Returns the weight of the projection specified by \c projection.
    */
   virtual Weight getWeight (const Coordinates & projection) const;

   virtual std::string name() const { return "POD"; }

   /**
    * Returns the order-dependent part of the weights.
    */
   OrderDependentWeights& getOrderDependentWeights()
   { return m_orderDependentWeights; }

   const OrderDependentWeights& getOrderDependentWeights() const
   { return m_orderDependentWeights; }

   /**
    * Returns the product part of the weights.
    */
   ProductWeights& getProductWeights()
   { return m_productWeights; }

   const ProductWeights& getProductWeights() const
   { return m_productWeights; }

#ifdef WITH_XML
   /**
    * Static factory method; create a \c PODWeights object by
    * parsing XML data.
    */
   static PODWeights* createFromXML (const pugi::xml_node & node);
#endif

protected:
   /// \copydoc LatCommon::Weights::format()
   virtual void format(std::ostream& os) const;
};

}

#endif
