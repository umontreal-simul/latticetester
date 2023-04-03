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

#include "latticetester/WeightsPOD.h"
#include <sstream>

namespace LatticeTester
{

  //===========================================================================

  WeightsPOD::WeightsPOD ()
  {
  }

  //===========================================================================

  Weight WeightsPOD::getWeight (const Coordinates& projection) const
  {
    return m_orderDependentWeights.getWeight(projection) * m_productWeights.getWeight(projection);
  }

  //===========================================================================

  void WeightsPOD::format(std::ostream& os) const
  {
    using LatticeTester::operator<<;
    os << "WeightsPOD(" << m_orderDependentWeights << ", " << m_productWeights << ")";
  }

} // namespace LatticeTester

//===========================================================================

// #ifdef WITH_XML
// #include "xmlerror.hpp"
// 
// namespace LatticeTester
// {
// 
//   WeightsPOD* WeightsPOD::createFromXML (const pugi::xml_node& root)
//   {
//     WeightsPOD* o = new WeightsPOD();
// 
//     pugi::xml_node node;
// 
//     try {
//       // order-dependent weights
//       node = root.child("order-dependent");
//       if (!node)
//         throw pugi::xml_error(root, "missing <order-dependent> element");
//       WeightsOrderDependent* odw = WeightsOrderDependent::createFromXML(*node);
//       o->m_orderDependentWeights = *odw;
//       delete odw;
// 
//       // product weights
//       node = root.child("product");
//       if (!node)
//         throw pugi::xml_error(root, "missing <product> element");
//       WeightsProduct* pw = WeightsProduct::createFromXML(*node);
//       o->m_productWeights = *pw;
//       delete pw;
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
