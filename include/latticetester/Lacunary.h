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

#ifndef LATTICETESTER__LACUNARY_H
#define LATTICETESTER__LACUNARY_H
#include "latticetester/Util.h"
#include "latticetester/ntlwrap.h"

#include <string>
#include <cstdint>

namespace LatticeTester {

  /**
   * This class implements sets of lacunary indices.
   */
  template<typename BasInt>
    class Lacunary {
      private:
        typedef NTL::vector<BasInt> BasIntVec;
      public:

        /**
         * Constructor for a set of \f$t\f$ lacunary indices given in \f$C[j]\f$, for
         * \f$j = 1, 2, …, t\f$.
         */
        Lacunary (const BasIntVec & C, int t)
        {
          m_dim = t; 
          LatticeTester::CreateVect (m_lac, t);
          LatticeTester::CopyVect (m_lac, C, t);
        }

        /**
         * Constructor for a set of \f$t\f$ lacunary indices. The lacunary
         * indices can be read later by `ReadLac`.
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
         * Returns the lacunary index <tt>m_lac[j]</tt>.
         */
        BasInt & operator[] (int i) {return m_lac[i];}

        /**
         * Same as above.
         */
        BasInt & getLac (int i) {return m_lac[i];}

        /**
         * Returns the size of the lacunary set (the number of lacunary
         * indices).
         */
        int getSize () const { return m_dim; }

        /**
         * Computes lacunary indices by groups of \f$s\f$, spaced apart by
         * \f$2^w\f$. If \f$w=0\f$, this is the case of non-lacunary indices.
         * If `MinDim` is smaller than `Order`, it is reset to `Order` + 1.
         * Returns `true` in the lacunary case, and `false` otherwise.
         */
        bool calcIndicesStreams (int s, int w, int & minDim, int maxDim, 
            int order);

        /**
         * Returns this object as a string.
         */
        std::string toString () const;
      private:

        /**
         * The dimension (or size) of <tt>m_lac</tt>.
         */
        int m_dim;

        /**
         * The set of lacunary indices is <tt>m_lac[j]</tt> for \f$j = 1, 2,
         * …,\f$ <tt>m_dim</tt>.
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
    bool Lacunary<BasInt>::calcIndicesStreams (int s, int w, 
        int & minDim, int maxDim, int order)
    {
      if (w == 0) {
        if (minDim <= order)
          minDim = order + 1;
        return false;
      }

      if (m_dim < maxDim) {
        LatticeTester::DeleteVect (m_lac);
        LatticeTester::CreateVect (m_lac, maxDim);
        m_dim = maxDim;
      }

      BasInt t1;
      LatticeTester::power2 (t1, (std::int64_t) w);
      BasInt t;
      t = 0;
      int i = 1;
      while (true) {
        for (int j = 0; j < s; j++) {
          if (i <= maxDim) {
            m_lac[i] = t + j;
            i++;
          } else
            return true;
        }
        t += t1;
      }
    }
} // End namespace LatticeTester
#endif
