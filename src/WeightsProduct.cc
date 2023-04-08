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

#include "latticetester/WeightsProduct.h"
#include "latticetester/Util.h"
#include <sstream>


namespace LatticeTester
{

  //===========================================================================

  WeightsProduct::WeightsProduct (Weight defaultWeight)
    : m_defaultWeight(defaultWeight)
  {
  }

  //===========================================================================

  void WeightsProduct::setWeightForCoordinate (Coordinates::size_type coordinate, Weight weight)
  {
    if (coordinate >= m_weights.size())
      m_weights.resize(coordinate + 1);
    m_weights[coordinate] = weight;
  }

  //===========================================================================

  Weight WeightsProduct::getWeight (const Coordinates & projection) const
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

  void WeightsProduct::format(std::ostream& os) const
  {
    using LatticeTester::operator<<;
    os << "WeightsProduct(" << m_weights << ", default=" << m_defaultWeight << ")";
  }

} // namespace LatticeTester

//===========================================================================

// #ifdef WITH_XML
// #include "xmlerror.hpp"
// 
// namespace LatticeTester
// {
// 
//   WeightsProduct* WeightsProduct::createFromXML (const pugi::xml_node & root)
//   {
//     WeightsProduct* o = new WeightsProduct();
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
//         int64_t coordinate = lexical_cast<Weight>(node.child_value());
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
