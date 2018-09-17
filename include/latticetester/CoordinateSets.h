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

#include <iterator>
#include <map>

#include "latticetester/Coordinates.h"

namespace LatticeTester {

  /**
   * A namespace containing different implementation of sets of coordinates.
   *
   * These can be used to specify a set of projections of lattices or point sets.
   * The classes in this namespace are iterable but are not containers, so they
   * require very little storage. They virtually contain objects of type
   * LatticeTester::Coordinates.
   *
   * \todo present/make a template class for reimplementations so that the
   * classes in this namespace all work the same (they already do) and can be 
   * reimplemented easily.
   * \code
   * class CoordinateSets {
   *   class const_iterator: public std::iterator<std::forward_iterator_tag,
   *     const Coordinates> {
   *     public:
   *       struct end_tag {};
   *
   *       // Constructor that sets the iterator at the beginning
   *       explicit const_iterator(const CoordinateSets& seq):
   *         m_seq(&seq), m_atEnd(false){
   *           resetToOrder (m_seq->ranges().begin());
   *       }
   *
   *       // Constructor that sets the iterator at the end
   *       const_iterator(const FromRanges& seq, end_tag):
   *         m_seq(&seq), m_atEnd(true) {}
   *
   *       // Copy constructor.
   *       const_iterator(const const_iterator& other) {
   *         this->m_seq = other.m_seq;
   *         this->m_atEnd = other.m_atEnd;
   *         this->m_value = other.m_value;
   *       }
   *
   *       // Default constructor that does nothing
   *       const_iterator(): m_seq(nullptr), m_atEnd(false) {}
   *
   *       // Assignment operator. Normally this should do the same as the copy
   *       // constructor.
   *       const_iterator& operator=(const const_iterator& other) {
   *         this->m_seq = other.m_seq;
   *         this->m_atEnd = other.m_atEnd;
   *         this->m_value = other.m_value;
   *         return *this;
   *       }
   *
   *       // Comparison operator.
   *       bool operator==(const const_iterator& other) {
   *         return m_seq == other.m_seq 
   *           && (other.m_atEnd ? m_atEnd : m_value == other.m_value);
   *       }
   * 
   *       // Second comparison operator.
   *       bool operator!=(const const_iterator& other) {
   *         return !(*this == other);
   *       }
   *
   *       // dereference operator. This should return the Coordinates this
   *       // iterator is at.
   *       const Coordinates& operator*() const {
   *         return m_value;
   *       }
   *
   *       // This is definitely the operator that should be reimplemented.
   *       // Prefix increment operator
   *       const_iterator& operator++();
   *
   *       // Same as before
   *       // Postfix iteration operator
   *       const_iterator operator++(int);
   *
   *     private:
   *       // Maybe this is not that usefull
   *       void resetToOrder (const RangeMap::const_iterator& it);
   *
   *       // This is the object this iterator iterates over.
   *       const CoordinateSets* m_seq;
   *
   *       // Maybe this is not usefull
   *       bool m_atEnd;
   *
   *       // The set of values the iterator is at
   *       Coordinates m_value;
   *   };
   *
   *   // Gives an iterator that is at the beginning of the sequence of the
   *   // object
   *   const_iterator begin() const {
   *     return const_iterator(*this);
   *   }
   *
   *   // Gives an iterator that is at the end of the sequence of the object
   *   const_iterator end() const {
   *     return const_iterator(*this, typename const_iterator::end_tag{});
   *   }
   * }
   * \endcode
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
     * */
    class FromRanges {

      private:
        // Represents a range of integers. The first coordinate is the lower
        // bound and the second the upper bound.
        typedef std::pair<Coordinates::value_type, Coordinates::value_type> Range;

        // Maps Range objects to a size of subsets.
        typedef std::map<Coordinates::size_type, Range> RangeMap;

        // For every size of subsets, this contains the the range of integer
        // from which the subset can be drawn.
        RangeMap m_ranges;

        // Having a private getter is not that usefull.
        const RangeMap& ranges() const {return m_ranges;}

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
        FromRanges (Coordinates::size_type minOrder, Coordinates::size_type maxOrder,
            Coordinates::value_type minCoord, Coordinates::value_type maxCoord);

        /**
         * Constructs an empty set of coordinate sets. When using this
         * constructor, it is possible to handpick coordinate sets with the
         * `includeOrder()` method.
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

        /**
         * An iterator class used internaly by the `FromRange` class. Given an
         * object of this class, it is possible to cycle through the element it 
         * contains with the increment (`++`) operator.
         * */
        class const_iterator : public std::iterator<std::forward_iterator_tag,
        const Coordinates>
      {
        public:
          struct end_tag {};

          /**
           * Constructor for an iterator at the begining of the list of sets that
           * `seq` contains.
           * */
          explicit const_iterator(const FromRanges& seq):
            m_seq(&seq), m_atEnd(false){
              resetToOrder (m_seq->ranges().begin());
            }

          /**
           * Constructor for an iterator at the end of the list of sets that
           * `seq` contains.
           * */
          const_iterator(const FromRanges& seq, end_tag):
            m_seq(&seq), m_atEnd(true) {}

          /**
           * Copy constructor.
           * */
          const_iterator(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_atEnd = other.m_atEnd;
            this->m_value = other.m_value;
          }

          /**
           * Empty FromRange::const_iterator constructor
           * */
          const_iterator(): m_seq(nullptr), m_atEnd(false) {}

          /**
           * Assignment operator. This copies the left hand side object in the
           * right hand side one.
           * */
          const_iterator& operator=(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_atEnd = other.m_atEnd;
            this->m_value = other.m_value;
            return *this;
          }

          /**
           * Compares this instance with `other`, returning true if they are 
           * associated with the same FromRange object and if they are at the
           * same point in their enumeration cycle.
           * */
          bool operator==(const const_iterator& other) {
            return m_seq == other.m_seq 
              && (other.m_atEnd ? m_atEnd : m_value == other.m_value);
          }

          /**
           * Compares this instance with `other`, returning false if they are 
           * associated with the same FromRange object and if they are at the
           * same point in their enumeration cycle.
           * */
          bool operator!=(const const_iterator& other) {
            return !(*this == other);
          }

          /**
           * Dereference operator, when dereferencing this object, you get the
           * set of coordinates this iterator points to.
           * */
          const Coordinates& operator*() const {
            return m_value;
          }


          /**
           * Prefix incrementation operator. Increments this iterator and
           * returns the iterator in its new state.
           * */
          const_iterator& operator++();

          /**
           * Postfix incrementation operator. Increments this iterator and
           * returns a copy of the old object.
           * */
          const_iterator operator++(int);

        private:

          /**
           * Resets the iterator at the beginning of the order the iterator `it`
           * is at. `it` should be one of the states of
           *  `m_seq->ranges()::const_iterator` because the values this function
           *  will take to build the first subset will be the one passed by the
           *  iterator.
           * */
          void resetToOrder (const RangeMap::const_iterator& it);

          /**
           * The FromRanges object through which this iterator cycles.
           * */
          const FromRanges* m_seq;

          /**
           * Is true if this iterator is at the end of the sequence and false
           * otherwise.
           * */
          bool m_atEnd;

          /**
           * The coordinates of the current set in the sequence of coordinates
           * this object contains.
           * */
          Coordinates m_value;

      }; // end FromRanges::const_iterator class

        /**
         * Returns a const_iterator pointing to the first element in the sequence
         * of coordinates sets that the object contains. It can then be used to
         * cycle through all the sets with `++`.
         */
        const_iterator begin() const {
          return const_iterator(*this);
        }

        /**
         * Returns a const_iterator pointing past the last element in the seq. 
         * This iterator cannot be used to cycle through the sets, but it can be
         * used as a point of comparison to know when you reached the end of an
         * iterator obtained with `begin()`.
         */
        const_iterator end() const {
          return const_iterator(*this, typename const_iterator::end_tag{});
        }

    }; // End class FromRanges

    /**
     * This class implements a CoordinateSets that will build all the subsets
     * of a Coordinates object that are of a cardinality in a certain range.
     * This is a less flexible class than `FromRanges` but, in the
     * example `Subsets.cc` it has been slightly faster in our tests.
     *
     * \remark it is somewhat hard to build an object of this class, maybe this
     * could be done in a better way.
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
         * Returns the coordinates, as passed to the constructor.
         * */
        const Coordinates& coords() const { return  m_coords; }

        /**
         * Returns `minOrder`, as passed to the constructor.
         * */
        Coordinates::size_type minOrder() const { return m_minOrder; }

        /**
         * Returns `maxOrder`, as passed to the constructor.
         * */
        Coordinates::size_type maxOrder() const { return m_maxOrder; }

        /**
         * An iterator class used internaly by the `Subsets` class. Given an
         * object of this class, it is possible to cycle through the element it 
         * contains with the increment (`++`) operator.
         * */
        class const_iterator : public std::iterator<std::forward_iterator_tag,
        const Coordinates>
      {
        public:
          struct end_tag {};

          /**
           * Constructor for an iterator at the begining of the list of sets that
           * `seq` contains.
           * */
          explicit const_iterator(const Subsets& seq):
            m_seq(&seq), m_atEnd(false){
              resetToOrder (m_seq->minOrder());
            }

          /**
           * Constructor for an iterator at the end of the list of sets that
           * `seq` contains.
           * */
          const_iterator(const Subsets& seq, end_tag):
            m_seq(&seq), m_atEnd(true) {}

          /**
           * Copy constructor.
           * */
          const_iterator(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_atEnd = other.m_atEnd;
            this->m_value = other.m_value;
          }

          /**
           * Default constructor.
           * */
          const_iterator():
            m_seq(nullptr), m_atEnd(false) {}

          /**
           * Assignment operator that does the same as the copy constructor.
           * */
          const_iterator& operator=(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_atEnd = other.m_atEnd;
            this->m_value = other.m_value;
            return *this;
          }

          /**
           * Compares this instance with `other`, returning true if they are 
           * associated with the same Subsets object and if they are at the
           * same point in their enumeration cycle.
           * */
          bool operator==(const const_iterator& other) {
            return m_seq == other.m_seq &&
              (other.m_atEnd ? m_atEnd : m_value == other.m_value);
          }

          /**
           * Compares this instance with `other`, returning false if they are 
           * associated with the same Subsets object and if they are at the
           * same point in their enumeration cycle.
           * */
          bool operator!=(const const_iterator& other) {
            return !(*this == other);
          }

          /**
           * Dereference operator, when dereferencing this object, you get the
           * set of coordinates this iterator points to.
           * */
          const Coordinates& operator*() const {
            return m_value;
          }

          /**
           * Prefix incrementation operator. Increments this iterator and
           * returns the iterator in its new state.
           * */
          const_iterator& operator++();

          /**
           * Postfix incrementation operator. Increments this iterator and
           * returns a copy of the old object.
           * */
          const_iterator operator++(int);

        private:

          /**
           * Resets the iterator at the beginning of the enumeration of subsets
           * of cardinality order.
           * */
          void resetToOrder (Coordinates::size_type order);

          /**
           * The Subsets object through which this iterator cycles.
           * */
          const Subsets* m_seq;

          /**
           * Is true if this iterator is at the end of the sequence and false
           * otherwise.
           * */
          bool m_atEnd;

          /**
           * The coordinates of the current set in the sequence of coordinates
           * this object contains.
           * */
          Coordinates m_value;

      }; // end Subsets::const_iterator class

        /**
         * Returns an iterator pointing to the first element of `*this`.
         */
        const_iterator begin() const {
          return const_iterator(*this);
        }

        /**
         * Returns an iterator pointing past the last element of `*this`.
         */
        const_iterator end() const {
          return const_iterator(*this, typename const_iterator::end_tag{});
        }

      private:

        /**
         * The coordinates from which this object builds subsets.
         * */
        const Coordinates m_coords;

        /**
         * The minimal cardinality of the subsets this object will build.
         * */
        Coordinates::size_type m_minOrder;

        /**
         * The maximal cardinality of the subsets this object will build.
         * */
        Coordinates::size_type m_maxOrder;

    }; // end Subsets class

    //==========================================================================

    /**
     * This template class wraps any implementation of a CoordinateSets and
     * adds a specific coordinate to each coordinate sets. When iterating
     * through the subsets of this class, they won't be order by cardinality
     * because the sets are generated on the fly and the coordinate is added if
     * it is not already contained.
     *
     * \tparam BASE Type of coordinate sets that serves as a base. This should
     * be one of the classes defined in the `CoordinateSets` namespace.
     *
     * \todo Should this class implement the CoordinateSets template? Because it
     * would need a few more attributes. This class also needs to be benchmarked
     * to see how much slower it is (because it will be).
     */
    template <typename BASE> class AddCoordinate {

      public:
        /**
         * Constructs a sequence of coordinate sets by adding the coordinate
         * \c coord to each element in the base sequence \c base.
         */
        AddCoordinate (const BASE& base, Coordinates::value_type coord):
          m_base(base), m_coord(coord) {}

        /**
         * Returns the base object used to produce the subsets.
         */
        const BASE& base() const {
          return  m_base;
        }

        /**
         * Returns the coordinate this object adds.
         */
        const Coordinates::value_type& coord() const {
          return m_coord;
        }

        /**
         * An iterator class used internaly by the `AddCoordinate` class.
         * Given an object of this class, it is possible to cycle through the
         * element it contains with the increment (`++`) operator.
         *
         * Since this iterator is built upon another object, it internally keeps
         * an iterator that cycles through the sets generated by the
         * `CoordinateSets` implementation the associated to the `AddCoordinate`
         * object used for construction.
         * */
        class const_iterator : public std::iterator<std::forward_iterator_tag,
        const Coordinates>
      {
        public:
          struct end_tag {};

          /**
           * Constructor for an iterator at the beginning of the list of sets
           * that `seq` contains.
           * */
          explicit const_iterator(const AddCoordinate& seq):
            m_seq(&seq) {
              this->underlying = m_seq->base().begin();
              updateValue();
            }

          /**
           * Constructor for an iterator at the end of the list of sets that
           * `seq` contains.
           * */
          const_iterator(const AddCoordinate& seq, end_tag):
            m_seq(&seq) {
              this->underlying = m_seq->base().end();
              updateValue();
            }

          /**
           * Copy constructor.
           * */
          const_iterator(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_value = other.m_value;
          }

          /**
           * Default constructor. Does nothing.
           * */
          const_iterator(): m_seq(nullptr) {}

          /**
           * Assignment operator. This copies the left hand side object in the
           * right hand side object.
           * */
          const_iterator& operator=(const const_iterator& other) {
            this->m_seq = other.m_seq;
            this->m_value = other.m_value;
            return *this;
          }

          /**
           * Compares this instance with `other`, returning true if they are 
           * associated with the same AddCoordinate object and if they are at
           * the same point in their enumeration cycle.
           * */
          bool operator==(const const_iterator& other) {
            return m_seq == other.m_seq && this->underlying == other.underlying;
          }

          /**
           * Compares this instance with `other`, returning false if they are 
           * associated with the same AddCoordinate object and if they are at
           * the same point in their enumeration cycle.
           * */
          bool operator!=(const const_iterator& other) {
            return !(*this == other);
          }

          /**
           * Prefix increment operator. This increments the underlying subset
           * and updates `m_value`.
           * */
          const_iterator& operator++() {
            ++(this->underlying);
            updateValue();
            return *this;
          }

          /**
           * Postfix increment operator. This increments the underlying subset
           * and updates `m_value`.
           * */
          const_iterator operator++(int) {
            const_iterator old(*this);
            ++(this->underlying);
            updateValue();
            return old;
          }

          /**
           * Dereference operator that returns the `Coordinates` this iterator
           * is at.
           * */
          const Coordinates& operator*() const {
            return m_value;
          }

        private:
          /**
           * Sets `m_value` to the same set as the underlying iterator and adds
           * the coordinate to it if it is not already contained.
           * */
          void updateValue() {
            m_value = *(this->underlying);
            m_value.insert(m_seq->coord());
          }

          /**
           * Pointer to the `AddCoordinate` object this iterator uses to build
           * the sets.
           * */
          const AddCoordinate* m_seq;

          /**
           * The iterator of the underlying `CoordinateSets`.
           * */
          typename BASE::const_iterator underlying;

          /**
           * The current `Coordinates` in the enumeration.
           * */
          Coordinates m_value;

      }; // end iterator_adaptator class

        /**
         * Returns an iterator pointing to the first element in the seq.
         */
        const_iterator begin() const {
          return const_iterator(*this);
        }

        /**
         * Returns an iterator pointing past the last element in the seq.
         */
        const_iterator end() const {
          return const_iterator(*this, typename const_iterator::end_tag{});
        }

      private:
        /**
         * The underlying `CoordinateSets` that will be augmented.
         * */
        BASE m_base;

        /**
         * The coordinate that will be added.
         * */
        Coordinates::value_type m_coord;

    }; // end AddCoordinate class

  } // End namespace CoordinateSets

} // End namespace LatticeTester

#endif
