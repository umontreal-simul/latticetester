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

#ifndef LATTICETESTER_LACUNARY_H
#define LATTICETESTER_LACUNARY_H
#include "latticetester/Util.h"
#include "latticetester/ntlwrap.h"

#include <string>
#include <cstdint>

namespace LatticeTester {

  /**
   * This class represents a set of indices with lacunary values. This class
   * stores the values of those indices as a vector, and they can be accessed
   * via the `[]` operator. There also is a method to construct a set with
   * values spaced in a certain way. This class is present here to be used with
   * the subclasses of IntLattice.
   *
   * \todo Test the implementation of this class. **Marc-Antoine**: I feel like
   * the current implementation of this class is obsolete : it contains too few
   * methods and functions to be of more use than the vector classes.
   */
  template<typename BasInt>
    class Lacunary {

      private:
        typedef NTL::vector<BasInt> BasIntVec;

      public:

        /**
         * Constructor for a set of \f$t\f$ indices given in the vector `C`.
         * \todo This constructor might not work as intended because BasIntVec
         * is an NTL::vector type and this function dynamically allocates an
         * array to it. There is no real reason to do that since there is a
         * setSize method available.
         */
        Lacunary (const BasIntVec & C, int t)
        {
          m_dim = t; 
          CreateVect (m_lac, t);
          CopyVect (m_lac, C, t);
        }

        /**
         * Constructor for an empty set of \f$t\f$ lacunary indices. This
         * constructor set the value of the indices to `0`. To set indices, one
         * should use either one `[]` or `getLac`, or the `calcIndicesStreams()`
         * method.
         */
        explicit Lacunary (int t = 0) 
        {
          m_dim = t; 
          LatticeTester::CreateVect (m_lac, t); 
        }

        /**
         * Destructor.
         */
        ~Lacunary () { LatticeTester::DeleteVect (m_lac); }

        /**
         * Returns a reference to the value of <tt>m_lac[i]</tt>. (`m_lac` is
         * the underlying vector storing the set of indices.)
         */
        BasInt & operator[] (int i) {return m_lac[i];}

        /**
         * Calling the `object.getLac(i)` gives the same result than calling
         * `object[i]`.
         */
        BasInt & getLac (int i) {return m_lac[i];}

        /**
         * Returns the size of the set, that is the number of elements in the
         * vector of indices.
         */
        int getSize () const { return m_dim; }

        /**
         * Fills the values of this object with `maxDim` indices starting from
         * `0` by groups of \f$s\f$, spaced apart by \f$2^w\f$. If \f$w=0\f$,
         * this is the case of non-lacunary indices. Returns `true` in the
         * lacunary case, and `false` otherwise.
         */
        bool calcIndicesStreams (int s, int w, int maxDim);

        /**
         * Returns a string that describes this object. This string will contain
         * the dimension and all the indices stored.
         */
        std::string toString () const;

      private:

        /**
         * The dimension (or size) of <tt>m_lac</tt>.
         */
        int m_dim;

        /**
         * The set of lacunary indices is <tt>m_lac[j]</tt> for \f$j = 0, 1,
         * \ldots,\f$ tt>m_dim</tt>.
         */
        BasIntVec m_lac;

    }; // End class Lacunary

  /*=========================================================================*/

  template<typename BasInt>
    std::string Lacunary<BasInt>::toString () const
    {
      std::ostringstream out;
      out << "dim = " << m_dim;
      out << "\nLac = {\n   ";
      for (int i = 0; i < m_dim; i++)
        out << m_lac[i] << "\n   ";
      out << "}\n";
      return out.str ();
    }


  /*=========================================================================*/

  template<typename BasInt>
    bool Lacunary<BasInt>::calcIndicesStreams (int s, int w, int maxDim)
    {

      if (m_dim < maxDim) {
        DeleteVect (m_lac);
        CreateVect (m_lac, maxDim);
        m_dim = maxDim;
      }

      BasInt t1;
      NTL::power2 (t1, (std::int64_t) w);
      BasInt t;
      t = 0;
      int i = 0;
      while (true) {
        for (int j = 0; j < s; j++) {
          if (i < maxDim) {
            m_lac[i] = t + j;
            i++;
          } else {
            if (w == 0) return false;
            return true;
          }
        }
        t += t1;
      }
    }

  extern template class Lacunary<std::int64_t>;
  extern template class Lacunary<NTL::ZZ>;

} // End namespace LatticeTester

#endif
