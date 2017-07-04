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

/* NormaRogers.h for ISO C++ */
#ifndef LATTICETESTER__NORMAROGERS_H
#define LATTICETESTER__NORMAROGERS_H
#include "latticetester/Normalizer.h"
#include <stdexcept>


namespace LatticeTester {

/**
 * This class implements the *Rogers* bounds on the density of sphere
 * packing. The length of vectors is computed using the
 * \f${\mathcal{L}}_2\f$ norm. The bounding lengths, for a lattice containing 
 * \f$n\f$ points per unit volume in dimension \f$t\f$, are given by 
 * \f$\ell_t^* = \gamma_t^{1/2} n^{-1/t}\f$, where the \f$\gamma_t\f$ are 
 * the *Rogers* lattice constants.
 */
class NormaRogers : public Normalizer {
public:

   /**
    * Constructor for the Rogers bounds. The lattices have \f$n\f$ points per 
    * unit volume, in all dimensions \f$\le t\f$. The bias factor `beta` 
    * \f$= \beta\f$ gives more weight to some of the dimensions.
    * Note this class stores the log value of the density to handle larger values.
    * There is no restriction on the dimension \f$t\f$ which can be larger than 48.
    */
   NormaRogers (const RScal & logDensity, int t, double beta = 1);

   /**
    * Destructor.
    */
   ~NormaRogers();

   /**
    * Returns the value of the Rogers lattice constant \f$\gamma_j\f$ in
    * dimension \f$j\f$.
    */
   double getGamma (int j) const throw (std::out_of_range);
private:

   /**
    * The lattice constants \f$\gamma_j\f$ are the Rogers bounds in each
    * dimension \f$j\f$.
    */
   double *m_gamma;

   /**
    * Precomputed lattice constants \f$\gamma_j\f$ for the Rogers bounds
    * in each dimension \f$j \le48\f$.
    */
   static const double m_gamma0[1 + Normalizer::MAX_DIM];

   /**
    * Computes the Rogers bound in dimension \f$d\f$.
    */
   double calcGamma (int d);
};

}
#endif
