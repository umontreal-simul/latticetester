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

/* NormaLaminated.h for ISO C++ */
#ifndef LATCOMON__NORMALAMINATED_H
#define LATCOMON__NORMALAMINATED_H
#include "latcommon/Normalizer.h"
#include <stdexcept>


namespace LatCommon {

/**
 * This class implements theoretical bounds on the length of the shortest
 * nonzero vector in a lattice, based on the densest sphere packing in
 * *laminated* lattices. The length of vectors is computed using the
 * \f${\mathcal{L}}_2\f$ norm. The bounding lengths, for a lattice of rank
 * \f$k\f$ containing \f$m\f$ points per unit volume in dimension \f$t\f$,
 * are given by \f$\ell_t^* = \gamma_t m^{k/t}\f$ for \f$t \ge k\f$, where
 * the \f$\gamma_t\f$ are the lattice constants for the best *laminated*
 * lattices \cite mCON99a&thinsp;.
 *
 */
class NormaLaminated : public Normalizer {
public:

   /**
    * Constructor for the bounds obtained for laminated lattices. The lattices
    * are those of rank \f$k\f$, with \f$m\f$ points per unit volume, in all
    * dimensions \f$\le t\f$. The bias factor `beta` \f$= \beta\f$ gives more
    * weight to some of the dimensions. Restriction: \f$t \le48\f$.
    */
   NormaLaminated (const MScal & m, int k, int t, double beta = 1);

   /**
    * Returns the value of the lattice constant \f$\gamma_j\f$ in
    * dimension \f$j\f$.
    */
   double getGamma (int j) const throw (std::out_of_range);
private:

   /**
    * Lattice constants \f$\gamma_j\f$ for the laminated lattices in each
    * dimension \f$j\f$.
    */
   static const double m_gamma[1 + Normalizer::MAX_DIM];
};

}
#endif
