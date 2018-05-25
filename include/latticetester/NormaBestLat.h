// This file is part of LatticeTester.
//
// LatticeTester
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

/* NormaBestLat.h for ISO C++ */
#ifndef LATTICETESTER__NORMABESTLAT_H
#define LATTICETESTER__NORMABESTLAT_H
#include "latticetester/Normalizer.h"
#include <stdexcept>


namespace LatticeTester {

/**
 * This class implements the *best* theoretical bounds on the length of the
 * shortest nonzero vector in a lattice, based on the densest sphere packing
 * in lattices. The length of a vector is computed using the \f${\mathcal{L}}_2\f$
 * norm. The bounding lengths for a lattice containing \f$n\f$ points per 
 * unit volume in dimension \f$t\f$, are given by 
 * \f$\ell_t^* = \gamma_t^{1/2} n^{-1/t}\f$, where the \f$\gamma_t\f$ are 
 * the lattice constants for the *best* known lattices \cite mCON99a&thinsp;.
 * Note this class stores the log value of the density to handle larger values.
 */
class NormaBestLat : public Normalizer {
public:

   /**
    * Constructor for the best bounds obtained for lattices. The lattices have
    * \f$Density\f$ points per unit volume, for all dimensions \f$\le t\f$. The bias 
    * factor `beta` \f$= \beta\f$ gives more weight to some of the dimensions. 
    * Restriction: \f$t \le48\f$.
    */
   NormaBestLat (RScal & logDensity, int t, double beta = 1);

   /**
    * Returns the value of the lattice constant \f$\gamma_j\f$ in
    * dimension \f$j\f$.
    */
   double getGamma (int j) const throw (std::out_of_range);
private:

   /**
    * Lattice constants \f$\gamma_j\f$ for the most general lattices in each
    * dimension \f$j\f$.
    */
   static const double m_gamma[1 + Normalizer::MAX_DIM];
};

}
#endif
