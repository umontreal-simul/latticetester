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

#include "latticetester/ProductWeights.h"
#include "latticetester/Util.h"
#include <sstream>


namespace LatticeTester
{

  //===========================================================================

  ProductWeights::ProductWeights (Weight defaultWeight)
    : m_defaultWeight(defaultWeight)
  {
  }

  //===========================================================================

  void ProductWeights::setWeightForCoordinate (Coordinates::size_type coordinate, Weight weight)
  {
    if (coordinate >= m_weights.size())
      m_weights.resize(coordinate + 1);
    m_weights[coordinate] = weight;
  }

  //===========================================================================

  Weight ProductWeights::getWeight (const Coordinates & projection) const
  {
    if (projection.empty())
      return m_defaultWeight;

    Weight w = 1.0;
    Coordinates::const_iterator it = projection.begin();
    while (it != projection.end()) {
      w *= getWeightForCoordinate(*it);
      ++it;
    }
    return w;
  }

  //===========================================================================

  void ProductWeights::format(std::ostream& os) const
  {
    using LatticeTester::operator<<;
    os << "ProductWeights(" << m_weights << ", default=" << m_defaultWeight << ")";
  }

} // namespace LatticeTester

//===========================================================================

// #ifdef WITH_XML
// #include "xmlerror.hpp"
// 
// namespace LatticeTester
// {
// 
//   ProductWeights* ProductWeights::createFromXML (const pugi::xml_node & root)
//   {
//     ProductWeights* o = new ProductWeights();
// 
//     pugi::xml_node node;
// 
//     try {
//       // default weight (optional)
//       node = root.child("default").child("weight");
//       if (node)
//         o->setDefaultWeight(lexical_cast<Weight>(node.child_value()));
// 
//       // per-coordinate weights
//       for (pugi::xml_node cnode = root.child("coordinate"); cnode; cnode = cnode.next_sibling("coordinate")) {
//         // weight
//         node = cnode.child("weight");
//         if (!node)
//           throw pugi::xml_error(cnode, "missing <weight> element");
//         Weight weight = lexical_cast<Weight>(node.child_value());
//         // coordinate index
//         node = cnode.child("coordinate");
//         if (!node)
//           throw pugi::xml_error(cnode, "missing <coordinate> element");
//         int coordinate = lexical_cast<Weight>(node.child_value());
//         // store weight
//         o->setWeightForCoordinate(coordinate, weight);
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
