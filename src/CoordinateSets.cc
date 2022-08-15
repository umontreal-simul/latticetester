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
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

  namespace CoordinateSets {

    //==========================================================================
    // FromRanges
    //==========================================================================

    FromRanges::FromRanges (
        Coordinates::size_type minOrder, Coordinates::size_type maxOrder,
        Coordinates::value_type minCoord, Coordinates::value_type maxCoord)
    {
      for (Coordinates::size_type order = minOrder; order <= maxOrder; order++)
        includeOrder(order, minCoord, maxCoord);
    }

    //==========================================================================

    FromRanges::FromRanges ()
    {}

    //==========================================================================

    void FromRanges::includeOrder (Coordinates::size_type order,
        Coordinates::value_type minCoord,
        Coordinates::value_type maxCoord)
    {
      if (order > 0 && maxCoord < minCoord + order - 1)
        throw std::invalid_argument
          ("CoordinateSets::FromRanges::includeOrder(): maxCoord < minCoord + order - 1");

      m_ranges[order] = Range (minCoord, maxCoord);
    }

    //==========================================================================

    void FromRanges::excludeOrder (Coordinates::size_type order)
    {
      RangeMap::const_iterator it = m_ranges.find (order);
      if (it != m_ranges.end ())
        m_ranges.erase (it);
    }

    //==========================================================================

    void FromRanges::const_iterator::resetToOrder (
        const RangeMap::const_iterator& itRange)
    {
      // empty coordinate set
      m_value.clear ();

      // stop if all orders exhausted
      if (itRange == m_seq->ranges().end ()) {
        m_atEnd = true;
        return;
      }
      // new order
      const Coordinates::size_type order = itRange->first;

      // new range
      const Coordinates::value_type minCoord = itRange->second.first;
      const Coordinates::value_type maxCoord = itRange->second.second;

      // double check range
      while (order > 0 && maxCoord < minCoord + order - 1)
        // this should never happen
        throw std::logic_error
          ("CoordinateSets::FromRanges::resetToOrder(): maxCoord < minCoord + order - 1");

      // add initial coordinates to the set
      Coordinates::value_type i = minCoord;
      while (m_value.size () < order && i <= maxCoord)
        m_value.insert(i++);
    }

    //==========================================================================

    FromRanges::const_iterator& FromRanges::const_iterator::operator++()
    {
      if (m_atEnd)
        return *this;

      // current order
      const Coordinates::size_type order = m_value.size();

      // iterator on ranges pointing to current order
      RangeMap::const_iterator itRange = m_seq->ranges().find (order);

      if (itRange == m_seq->ranges().end ()) {
        // current order has been removed during iteration
        m_atEnd = true;
        return *this;
      }

      const Coordinates::value_type maxCoord = itRange->second.second;

      // maximum index value
      Coordinates::value_type maxval = maxCoord;

      // iterator before the min element
      Coordinates::reverse_iterator rend = m_value.rend();

      // iterator on the max element
      Coordinates::reverse_iterator rbegin = m_value.rbegin();
      Coordinates::reverse_iterator rit = rbegin;

      // find the maximum index that has not reached its maximum value yet
      while (rit != rend && *rit >= maxval) {
        ++rit;
        --maxval;
      }

      if (rit == rend) {
        // all indices have reached their maximum value
        // go to next order
        resetToOrder (++itRange);
        return *this;
      }
      // Increment the projection indices as if they were a sequence of
      // digits.

      // make note of the highest index that is below its maximum value (to be
      // incremented)
      Coordinates::value_type high = *rit;

      // erase all indices that are at their maximum value plus the one to be
      // incremented
      m_value.erase(--rit.base(), rbegin.base());

      // insert new values
      while (m_value.size () < order)
        m_value.insert(++high);

      return *this;
    }

    //==========================================================================

    FromRanges::const_iterator FromRanges::const_iterator::operator++(int)
    {
      FromRanges::const_iterator tmp(*this);
      ++(*this);
        return tmp;
    }
    //==========================================================================
    // CoordinateSubsets
    //==========================================================================

    Subsets::Subsets (const Coordinates & coords,
        Coordinates::size_type minOrder,
        Coordinates::size_type maxOrder)
      : m_coords (coords), m_minOrder (minOrder), m_maxOrder (maxOrder)
    {}

    //==========================================================================

    void Subsets::const_iterator::resetToOrder (Coordinates::size_type order)
    {
      m_value.clear();

      // check if any projections at all are available
      if (order > m_seq->maxOrder() || m_seq->coords().size () < order) {
        m_atEnd = true;
        return;
      }

      Coordinates::const_iterator itp = m_seq->coords().begin ();
      while (m_value.size () < order && itp != m_seq->coords().end ())
        m_value.insert(*itp++);
    }

    //==========================================================================

    Subsets::const_iterator& Subsets::const_iterator::operator++()
    {
      if (m_atEnd) return *this;

      // current order
      Coordinates::size_type order = m_value.size();

      // maximum index value
      Coordinates::reverse_iterator ritMaxVal = m_seq->coords().rbegin();

      // iterator before the min element
      Coordinates::reverse_iterator rend = m_value.rend();

      // iterator on the max element
      Coordinates::reverse_iterator rbegin = m_value.rbegin();
      Coordinates::reverse_iterator rit = rbegin;

      // find the maximum index that has not reached its maximum value yet
      while (rit != rend && *rit >= *ritMaxVal) {
        ++rit;
        ++ritMaxVal;
      }

      if (rit == rend) {
        // all indices have reached their maximum value
        // go to next order
        resetToOrder (order + 1);
        return *this;
      }
      // Increment the projection indices as if they were a sequence of
      // digits.

      // make note of the highest index that is below its maximum value (to be
      // incremented)
      Coordinates::const_iterator itNewVal = m_seq->coords().upper_bound (*rit);
      //if (itNewVal == m_seq->coords().end()) {
      //  m_atEnd = true;
      //  return *this;
      //}
      // erase all indices that are at their maximum value plus the one to be
      // incremented
      m_value.erase(--rit.base(), rbegin.base());

      // insert new values
      while (m_value.size() < order)
        m_value.insert(*itNewVal++);
      return *this;
    }

    //==========================================================================

    Subsets::const_iterator Subsets::const_iterator::operator++(int)
    {
      Subsets::const_iterator tmp(*this);
      ++(*this);
      return tmp;
    }

  } // End namespace CoordinatesSets

} // End namespace LatticeTester
