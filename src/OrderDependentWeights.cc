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

#include "latticetester/OrderDependentWeights.h"
#include "latticetester/Weights.h"
#include "latticetester/Util.h"

#include <sstream>


namespace LatticeTester
{

  //===========================================================================

  OrderDependentWeights::OrderDependentWeights (Weight defaultWeight)
    : m_defaultWeight(defaultWeight)
  {
  }

  //===========================================================================

  void OrderDependentWeights::setWeightForOrder (Coordinates::size_type order, Weight weight)
  {
    if (order >= m_weights.size())
      m_weights.resize(order + 1);
    m_weights[order] = weight;
  }

  //===========================================================================

  Weight OrderDependentWeights::getWeight (const Coordinates& projection) const
  {
    return getWeightForOrder(projection.size());
  }

  //===========================================================================

  void OrderDependentWeights::format(std::ostream& os) const
  {
    using LatticeTester::operator<<;
    os << "OrderDependentWeights(" << m_weights << ", default=" << m_defaultWeight << ")";
  }

} // namespace LatticeTester

//===========================================================================

// #ifdef WITH_XML
// #include "xmlerror.hpp"
// 
// namespace LatticeTester
// {
// 
//   OrderDependentWeights* OrderDependentWeights::createFromXML (const pugi::xml_node& root)
//   {
//     OrderDependentWeights* o = new OrderDependentWeights();
// 
//     pugi::xml_node node;
// 
//     try {
//       // default weight
//       node = root.child("default").child("weight");
//       if (node)
//         o->setDefaultWeight(lexical_cast<Weight>(node.child_value()));
// 
//       // per-dimension weights
//       for (pugi::xml_node dnode = root.child("dimension"); dnode; dnode = dnode.next_sibling("dimension")) {
//         // weight
//         node = dnode.child("weight");
//         if (!node)
//           throw pugi::xml_error(dnode, "missing <weight> element");
//         Weight weight = lexical_cast<Weight>(node.child_value());
//         // dimension (projection order)
//         node = dnode.child("dimension");
//         if (!node)
//           throw pugi::xml_error(dnode, "missing <dimension> element");
//         int dimension = lexical_cast<int>(node.child_value());
// 
//         // store weight
//         o->setWeightForOrder(dimension, weight);
//       }
// 
//       return o;
//     }
//     catch (bad_lexical_cast& e) {
//       delete o;
//       throw pugi::xml_error(node, e.what());
//     }
//   }
// 
// } // namespace LatticeTester
// #endif
