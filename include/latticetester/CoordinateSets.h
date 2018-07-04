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

#ifndef LATTICETESTER__COORDINATE_SETS_H
#define LATTICETESTER__COORDINATE_SETS_H

#include <set>
#include <map>
#include <iostream>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/config.hpp>

#include "latticetester/Coordinates.h"

namespace LatticeTester {

  /**
   * A classes containing multiple sets of coordinates.
   *
   * These can be used to specify a set of projections of lattices or point sets.
   * The classes in this namespace are iterable but are not containers, so they
   * require very little storage.  They virtually contain objects of type
   * LatticeTester::Coordinates.
   */
  namespace CoordinateSets {

    /**
     * A CoordinateSets for coordinates within a given range.
     * This contains ways to build all subsets of coordinates of a given size
     * in an interval \f$\{\mathtt{minCoord}, \dots, \mathtt{maxCoord}\}\f$.
     * The intended usage of this class is to generate subsets of different 
     * orders for a single interval.
     * 
     * This class keeps a map of intervals indiced with the order that will be
     * used to generate subsets from them. This means that this class only keeps
     * one interval for each order.
     * 
     * When you iterate through this class, it generates the next subset of 
     * coordinates in the interval associated with the current order. If there 
     * is none left, it moves on to the next order and starts generating subsets
     * from the interval associated with this order. It is not possible to 
     * generate coordinate sets of the same order for different intervals with 
     * the same instance of this class.
     */
    class FromRanges {

      public:
        /**
         * Constructs a set of all subsets of \f$\{\mathtt{minCoord}, \dots,
         * \mathtt{maxCoord}\}\f$ with minimum and maximum cardinality
         * specified by \c minOrder and \c maxOrder.
         * For example, to select  all 1, 2, and 3-tuples over coordinates 2, 3, 4,
         * one may use the declaration <tt>FromRanges range(1, 3, 2, 4)</tt>;
         * this gives the sets
         * <tt>range = {{2}, {3}, {4}, {2, 3}, {2, 4}, {3, 4}, {2, 3, 4}}</tt>.
         */
        FromRanges (Coordinates::size_type  minOrder, Coordinates::size_type  maxOrder,
            Coordinates::value_type minCoord, Coordinates::value_type maxCoord);

        /**
         * Constructs an empty set of coordinate sets.
         *
         */
        FromRanges();

        /**
         * Include all subsets \f$\mathfrak u\f$ of
         * \f$\{\mathtt{minCoord}, \dots, \mathtt{maxCoord}\}\f$ of `order`
         * \f$|\mathfrak u| = \mathtt{order}\f$.
         * For example, calling \c includeOrder(3, 1, 5) causes all 3-tuples over
         * coordinates \f$1, \dots, 5\f$ to be included.
         * If `order = 0` (corresponding to the empty set),
         * `minCoord` and `maxCoord` are ignored.
         * 
         * Except for the case where `order = 0`, an exception is thrown if
         * \f$\mathtt{maxCoord} < \mathtt{minCoord} + \mathtt{order} - 1\f$.
         *
         * If the object already contains an interval for size `order`,
         * this interval is overwritten by the new one.
         */
        void includeOrder (Coordinates::size_type order,
            Coordinates::value_type minCoord,
            Coordinates::value_type maxCoord);

        /**
         * Removes the element associated with `order` from the list of 
         * of intervals to generate from. This means that no coordinate set
         * \f$\mathfrak u\f$ of order \f$|\mathfrak u| = \mathtt{order}\f$ will
         * be returned by the iterator until a new interval is given to the 
         * object for this order.
         */
        void excludeOrder (Coordinates::size_type order);

      private:
        typedef std::pair<Coordinates::value_type, Coordinates::value_type> Range;
        typedef std::map<Coordinates::size_type, Range> RangeMap;

        // Why is there a private getter?
        const RangeMap& ranges() const { return m_ranges; }

      public:
        /**
         * An iterator class used internaly by the `FromRange` class. Given an
         * object of this class, it is possible to cycle through the element it 
         * contains with the increment (`++`) or decrement (`--`) operators.
         * */
        class const_iterator : public boost::iterators::iterator_facade<const_iterator,
        const Coordinates, // value
        boost::iterators::forward_traversal_tag> // traversal
      {
        public:
          struct end_tag {};

          /**
           * Constructor for an iterator at the begining of the list of sets that
           * `seq` contains.
           * */
          explicit const_iterator(const FromRanges& seq):
            m_seq(&seq), m_atEnd(false)
        { resetToOrder (m_seq->ranges().begin()); }

          /**
           * Constructor for an iterator at the end of the list of sets that
           * `seq` contains.
           * */
          const_iterator(const FromRanges& seq, end_tag):
            m_seq(&seq), m_atEnd(true) // end
        {}

          /**
           * Empty FromRange::const_iterator constructor
           * */
          const_iterator():
            m_seq(nullptr), m_atEnd(false)
        { }

        private:
          friend class boost::iterators::iterator_core_access;

          /**
           * The iterator advances of one value. It generates a new subset of
           * coordinates and stores it in m_value that can be accessed with 
           * deference().
           * */
          void increment();

          /**
           * Compares this instance with `other`, returning true if they are 
           * associated with the same FromRange object and if 
           * */
          bool equal(const const_iterator& other) const
          { 
            return m_seq == other.m_seq 
              && (other.m_atEnd ? m_atEnd : m_value == other.m_value);
          }

          /**
           * Returns the coordinates set at which the iterator is.
           * */
          const Coordinates& dereference() const
          { return m_value; }

          /**
           * Resets the iterator at the beginning of order.
           * */
          void resetToOrder (const RangeMap::const_iterator& it);

        private:
          const FromRanges* m_seq;
          bool m_atEnd;
          Coordinates m_value;

      }; // end FromRanges::const_iterator class

        /**
         * Returns a const_iterator pointing to the first element in the sequence
         * of coordinates sets that the object contains. It can then be used to
         * cycle through all the sets with `++`.
         */
        const_iterator begin() const
        { return const_iterator(*this); }

        /**
         * Returns a const_iterator pointing past the last element in the seq. 
         * It can then be used to cycle through all the sets with `--`.
         */
        const_iterator end() const
        { return const_iterator(*this, typename const_iterator::end_tag{}); }

      private:
        RangeMap m_ranges; // coordinate range for each order

    }; // end FromRanges class


    /**
     * This class implements CoordinateSets for any set of coordinates.
     * It is more powerful than the class FromRange, but
     * slightly slower (by 15--20% according to empirical tests).
     */
    class Subsets {

      public:
      /**
       * Constructs a set of all subsets of `coords` with minimum and maximum
       * cardinality specified by `minOrder` and `maxOrder`.
       * For example, to select  all 1, 2, and 3-tuples over coordinates 2, 4, 6,
       * one may use the declaration `Subsets tousens(ens, 1, 3)`,
       * where  set `ens` is `{2,4,6}`; this gives the sets
       * `tousens = {{2}, {4}, {6}, {2, 4}, {2, 6}, {4, 6}, {2, 4, 6}}`.
       */
        Subsets (const Coordinates & coords,
            Coordinates::size_type minOrder,
            Coordinates::size_type maxOrder);

        /**
         * Returns the coordinates, as passed to the constructor \ref
         * Subsets(const Coordinates&,Coordinates::size_type,Coordinates::size_type)
         * */
        const Coordinates& coords() const { return  m_coords; }
        /**
         * Returns `minOrder`, as passed to the constructor \ref
         * Subsets(const Coordinates&,Coordinates::size_type,Coordinates::size_type)
         * */
        Coordinates::size_type minOrder() const { return m_minOrder; }
        /**
         * Returns `maxOrder`, as passed to the constructor \ref
         * Subsets(const Coordinates&,Coordinates::size_type,Coordinates::size_type)
         * */
        Coordinates::size_type maxOrder() const { return m_maxOrder; }

      public:
        /**
         * An iterator class used internaly by the `Subsets` class. Given an
         * object of this class, it is possible to cycle through the element it 
         * contains with the increment (`++`) or decrement (`--`) operators.
         * */
        class const_iterator : public boost::iterators::iterator_facade<const_iterator,
        const Coordinates, // value
        boost::iterators::forward_traversal_tag> // traversal
      {
        public:
          struct end_tag {};

          /**
           * Constructor for an iterator at the begining of the list of sets that
           * `seq` contains.
           * */
          explicit const_iterator(const Subsets& seq):
            m_seq(&seq), m_atEnd(false)
        { resetToOrder (m_seq->minOrder()); }

          /**
           * Constructor for an iterator at the end of the list of sets that
           * `seq` contains.
           * */
          const_iterator(const Subsets& seq, end_tag):
            m_seq(&seq), m_atEnd(true) // end
        {}

          /**
           * Empty constructor.
           * */
          const_iterator():
            m_seq(nullptr), m_atEnd(false)
        { }

        private:
          friend class boost::iterators::iterator_core_access;

          /**
           * The iterator advances of one value. It generates a new subset of
           * coordinates and stores it in m_value that can be accessed with 
           * deference().
           * */
          void increment();

          /**
           * Compares this instance with `other`, returning true if they are 
           * associated with the same FromRange object and if 
           * */
          bool equal(const const_iterator& other) const
          { return m_seq == other.m_seq && (other.m_atEnd ? m_atEnd : m_value == other.m_value); }

          /**
           * Returns the coordinates set at which the iterator is.
           * */
          const Coordinates& dereference() const
          { return m_value; }

          /**
           * Resets the iterator at the beginning of order.
           * */
          void resetToOrder (Coordinates::size_type order);

        private:
          const Subsets* m_seq;
          bool m_atEnd;
          Coordinates m_value;

      }; // end Subsets::const_iterator class

        /**
         * Returns an iterator pointing to the first element of `*this`.
         */
        const_iterator begin() const
        { return const_iterator(*this); }

        /**
         * Returns an iterator pointing past the last element of `*this`.
         */
        const_iterator end() const
        { return const_iterator(*this, typename const_iterator::end_tag{}); }

      private:
        const Coordinates m_coords;
        Coordinates::size_type m_minOrder;
        Coordinates::size_type m_maxOrder;

    }; // end Subsets class



    /**
     * This template class wraps any implementation of CoordinateSets and
     * adds a specific coordinate to each coordinate sets.
     *
     * \tparam BASE Type of coordinate sets that serves as a base.
     */
    template <typename BASE>
      class AddCoordinate {

        /**
         * Constructs a set based on <tt>baseSets</tt> by adding coordinate
         * \c coord to each of the sets in <tt>baseSets</tt>.
         * Warning: do not use a temporary object as <tt>baseSets</tt>.
         */
        public:
          typedef BASE BMat;

          /**
           * Constructs a sequence of coordinate sets by adding the coordinate \c coord
           * to each element in the base sequence \c base.
           */
          AddCoordinate (const BMat& base, Coordinates::value_type coord):
            m_base(base), m_coord(coord)
        { }

#ifdef BOOST_HAS_RVALUE_REFS
          AddCoordinate (BMat&& base, Coordinates::value_type coord):
            m_base(std::move(base)), m_coord(coord)
        { }
#endif

          /**
           * Returns the base sequence.
           */
          const BMat& base() const
          { return  m_base; }

          /**
           * Returns the added coordinate.
           */
          const Coordinates::value_type& coord() const
          { return m_coord; }

        public:
          class const_iterator : public boost::iterators::iterator_adaptor<const_iterator,
          typename BMat::const_iterator>
        {
          public:
            struct end_tag {};

            explicit const_iterator(const AddCoordinate& seq):
              const_iterator::iterator_adaptor_(seq.base().begin()),
              m_seq(&seq)
          { updateValue(); }

            const_iterator(const AddCoordinate& seq, end_tag):
              const_iterator::iterator_adaptor_(seq.base().end()),
              m_seq(&seq)
          {}

            const_iterator():
              const_iterator::iterator_adaptor_(),
              m_seq(nullptr)
          { }

          private:
            friend class boost::iterators::iterator_core_access;

            void increment()
            { ++this->base_reference(); updateValue(); }

            const Coordinates& dereference() const
            { return m_value; }

            void updateValue()
            { m_value = *this->base_reference(); m_value.insert(m_seq->coord()); }

          private:
            const AddCoordinate* m_seq;
            Coordinates m_value;
        }; // end iterator_adaptator class

          /**
           * Returns an iterator pointing to the first element in the seq.
           */
          const_iterator begin() const
          { return const_iterator(*this); }

          /**
           * Returns an iterator pointing past the last element in the seq.
           */
          const_iterator end() const
          { return const_iterator(*this, typename const_iterator::end_tag{}); }

        private:
          BMat m_base;
          Coordinates::value_type m_coord;

      }; // end AddCoordinate class


  } // end CoordinateSets namespace
} // end LatticeTester namspace

#endif // COORDINATE_SETS_H
