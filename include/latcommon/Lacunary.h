// This file is part of LatCommon.
//
// LatCommon
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
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

#ifndef LACUNARY_H
#define LACUNARY_H
#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include <string>


namespace LatCommon {

/**
 * This class implements sets of lacunary indices.
 *
 */
class Lacunary {
public:

   /**
    * Constructor for a set of \f$t\f$ lacunary indices given in \f$C[j]\f$, for
    * \f$j = 1, 2, …, t\f$.
    */
   Lacunary (const BVect & C, int t) { m_dim = t; CreateVect (m_lac, t);
        CopyVect (C, m_lac, t); }

   /**
    * Constructor for a set of \f$t\f$ lacunary indices. The lacunary
    * indices can be read later by `ReadLac`.
    */
   explicit Lacunary (int t = 0) { m_dim = t; CreateVect (m_lac, t); }

   /**
    * Destructor.
    */
   ~Lacunary () { DeleteVect (m_lac); }

   /**
    * Returns the lacunary index <tt>m_lac[j]</tt>.
    */
   BScal & operator[] (int i) {return m_lac[i];}

   /**
    * Same as above.
    */
   BScal & getLac (int i) {return m_lac[i];}

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
   bool calcIndicesStreams (int s, int w, int & minDim, int maxDim, int order);

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
   BVect m_lac;
};
}
#endif
