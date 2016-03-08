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

#ifndef LATCOMMON__RANK1LATTICE_H
#define LATCOMMON__RANK1LATTICE_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latcommon/IntLattice.h"


namespace LatCommon {

/**
 * This class implements a general rank 1 lattice basis. For the values
 * \f$a_1, a_2, …, a_d\f$ given, the \f$d\f$-dimensional lattice basis is
 * formed as:
 * \f[
 * \mathbf{b_1} = (a_1, a_2, …, a_d),\quad\mathbf{b_2} = (0, n, 0, …, 0),\quad…, \quad\mathbf{b_d} = (0, …, 0, n)
 * \f]
 * Without loss of generality, one may choose \f$a_1 = 1\f$.
 *
 */
class Rank1Lattice: public IntLattice {
public:

   /**
    * Constructor. \f$d\f$ represents the number of multipliers in the array
    * `a`.
    */
   Rank1Lattice (const MScal & n, const MVect & a, int d,
                    NormType norm = L2NORM);

   /**
    * Copy constructor.
    */
   Rank1Lattice (const Rank1Lattice & Lat);

   /**
    * Assigns `Lat` to this object.
    */
   Rank1Lattice & operator= (const Rank1Lattice & Lat);

   /**
    * Destructor.
    */
   ~Rank1Lattice();

   /**
    * Builds the basis in dimension \f$d\f$.
    */
   void buildBasis (int d);

protected:

   /**
    * Initializes the rank 1 lattice.
    */
   void init();

   /**
    * The multipliers of the rank 1 lattice rule.
    */
   MVect m_a;
};

}
#endif
