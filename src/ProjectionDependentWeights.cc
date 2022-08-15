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

#include "latticetester/Util.h"
#include "latticetester/ProjectionDependentWeights.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

namespace LatticeTester
{

  const ProjectionDependentWeights::WeightsMap ProjectionDependentWeights::m_emptyWeights;

  //===========================================================================

  ProjectionDependentWeights::ProjectionDependentWeights():
    Weights()
  {}

  //===========================================================================

  ProjectionDependentWeights::ProjectionDependentWeights (const ProjectionDependentWeights & o):
    Weights(),
    m_weights (o.m_weights)
  {}

  //===========================================================================

  void ProjectionDependentWeights::setWeight (const Coordinates & projection, Weight weight)
  {
    // set weight only if it is non-zero
    if (weight) {
      Coordinates::value_type largestIndex = *projection.rbegin();

      if (largestIndex >= m_weights.size())
        m_weights.resize(largestIndex + 1);

      m_weights[largestIndex][projection] = weight;
    }
  }

  //===========================================================================

  Weight ProjectionDependentWeights::getWeight (const Coordinates& projection) const
  {
    if (projection.size() == 0){
      return 0.0;
    }
    Coordinates::value_type largestIndex = *projection.rbegin();

    if (largestIndex < m_weights.size()) {
      // try explicit weights first
      WeightsMap::const_iterator it = m_weights[largestIndex].find(projection);
      if (it != m_weights[largestIndex].end())
        return it->second;
    }

    // fall back to zero
    return 0.0;
  }

  //===========================================================================

  const ProjectionDependentWeights::WeightsMap& ProjectionDependentWeights::getWeightsForLargestIndex (Coordinates::value_type largestIndex) const
  {
    if (largestIndex < m_weights.size())
      return m_weights[largestIndex];

    return m_emptyWeights;
  }

  //===========================================================================

  void ProjectionDependentWeights::format(ostream& os) const
  {
    using LatticeTester::operator<<;
    os << "ProjectionDependentWeights(" << m_weights << ")";
#if 0
    // line by line
    for (WeightsMap::const_iterator it = map.begin (); it != map.end (); it++)
      os << setw (30) << toString (it->first) << ": " << it->second
        << endl;
#endif
  }

  //===========================================================================
  // non-member functions
  //===========================================================================

  /**
   * Helper function to check the next characters from an input stream.
   *
   * Returns \c true if the next characters in \c is are \c token.  Upon a match, the characters
   * are removed from \c is; otherwise, they are left in \c is.
   */
  static bool check_next_chars(istream& is, const string& token)
  {
    string::const_iterator it = token.begin();
    while (it != token.end()) {
      int c = is.get();
      if (c != *it) {
        is.putback(c);
        break;
      }
      ++it;
    }
    if (it == token.end())
      return true;
    while (it != token.begin())
      is.putback(*(--it));
    return false;
  }

  /**
   * Helper function to skip all characters of a given class.
   */
  static void skip_any(istream& is, const string& characters)
  {
    while (characters.find(is.peek()) != string::npos)
      is.get();
  }

  //===========================================================================

  istream& operator>> (istream& is, ProjectionDependentWeights& weights)
  {
    string whitespace = " \t\n\r";
    string sep = whitespace + ",";

    // skip leading whitespace
    skip_any(is, whitespace);

    // check if coordinate set is enclosed in braces
    bool with_braces = false;
    if (is.peek() == '{') {
      is.get();
      with_braces = true;
      skip_any(is, whitespace);
    }

    while (is.good()) {

      // exit on closing brace
      if (with_braces && is.peek() == '}') {
        is.get();
        break;
      }

      // skip comments
      if (is.peek() == '#') {
        while (is.get() != '\n');
        skip_any(is, whitespace);
        continue;
      }

      // read coordinates
      Coordinates coords;
      is >> coords;

      // check for pair separator
      skip_any(is, whitespace);
      if (!(check_next_chars(is, ":") ||
            check_next_chars(is, "->") ||
            check_next_chars(is, "=>")))
        throw invalid_argument("invalid characters while parsing "
            "ProjectionDependentWeights from input stream");
      skip_any(is, whitespace);

      // read weight
      Weight weight;
      is >> weight;

      weights.setWeight(coords, weight);

      // skip across-pairs separator
      skip_any(is, sep);
    }
    return is;
  }

} // namespace

//===========================================================================

// #ifdef WITH_XML
// #include "xmlerror.hpp"
// 
// namespace LatticeTester
// {
// 
//   ProjectionDependentWeights* ProjectionDependentWeights::createFromXML (const pugi::xml_node& root)
//   {
//     ProjectionDependentWeights* o = new ProjectionDependentWeights();
// 
//     pugi::xml_node node;
// 
//     try {
//       // per-projection weights
//       for (pugi::xml_node pnode = root.child("projection"); pnode; pnode = pnode.next_sibling("projection")) {
// 
//         // weight
//         node = pnode.child("weight");
//         if (!node)
//           throw pugi::xml_error(pnode, "missing <weight> element");
//         Weight weight;
//         weight = lexical_cast<Weight>(node.child_value());
// 
//         // coordinates (projection)
//         node = pnode.child("coordinates");
//         if (!node)
//           throw pugi::xml_error(pnode, "missing <coordinates> element");
// 
//         // parse coordinate list as projection
//         Coordinates projection;
//         for (node = node.child("coordinate"); node; node = node.next_sibling("coordinate"))
//           projection.insert(lexical_cast<int>(node.child_value()));
// 
//         // store projection weight
//         o->setWeight(projection, weight);
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
// } // namespace
// #endif
