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

/* Normalizer.h for ISO C++ */
#ifndef LATCOMMON__NORMALIZER_H
#define LATCOMMON__NORMALIZER_H 
#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include "latcommon/Const.h"
#include <string>


namespace LatCommon {

/**
 * Classes which inherit from this base class are used in implementing bounds
 * on the length of the shortest nonzero vector in a lattice
 * \cite mCON99a&thinsp;. These bounds are used to normalize the length of
 * the shortest vectors. Tight lower bounds are available for all dimensions
 * for many important cases. In most cases, the \f${\mathcal{L}}_2\f$ norm is
 * used to compute the length of vectors.
 *
 * For some figures of merit, no useful bounds are known to normalize the
 * length of the shortest vector. In these cases, this base class will be
 * used as normalizer since it simply sets all normalization constants to 1.
 * This is necessary because the tests compare the normalized values of the
 * merit when searching for good lattices.
 *
 */
class Normalizer {
public:

   static const int MAX_DIM = 48;

   /**
    * Constructor for the bounds. Deals with lattices of rank \f$k\f$, having
    * \f$m\f$ points per unit volume, in all dimensions \f$\le t\f$. `Name` is
    * the name of the Normalizer. The bias factor `beta` \f$= \beta\f$ gives
    * more weight to some of the dimensions: taking \f$\beta< 1\f$ inflates the
    * figure of merit by \f$(1/\beta)^t\f$, thus weakening the requirements for
    * large \f$t\f$ in a worst-case figure of merit. One normally uses
    * \f$\beta= 1\f$.
    * \remark **Richard:** Je crois que ce facteur `beta` devrait
    * disparaître car des poids beaucoup plus généraux sont maintenant
    * implantés dans les classes `*Weights`.
    */
   Normalizer (const MScal & m, int k, int t, std::string Name,
                  NormType norm = L2NORM, double beta = 1);

   /**
    * Destructor.
    */
   virtual ~Normalizer () { delete [] m_cst; }

   /**
    * Initializes the bounds on the length of the shortest vector. The
    * lattices have \f$m\f$ points per unit volume, are of rank \f$k\f$,
    * and the bias factor is `beta` for all dimensions \f$j
    * \le\f$ `maxDim`.
    */
   void init (const MScal & m, int k, double beta);

   /**
    * Returns this object as a string.
    */
   std::string ToString () const;

   /**
    * Returns the norm associated with this object.
    */
   NormType getNorm () const { return m_norm; }

   /**
    * Sets the norm associated with this object to `norm`.
    */
   void setNorm (NormType norm) { m_norm = norm; }

   /**
    * Returns the maximal dimension for this object.
    */
   int getDim () const { return m_maxDim; }

   /**
    * Returns the bound on the length of the shortest nonzero vector in
    * dimension \f$j\f$.
    */
   double & getCst (int j);

   /**
    * Returns the value of the lattice constant \f$\gamma_j\f$ in
    * dimension \f$j\f$. For this base class, always returns 1.
    */
   virtual double getGamma (int j) const;
protected:

   /**
    * Name of the normalizer.
    */
   std::string m_name;

   /**
    * Norm associated with this object.
    */
   NormType m_norm;

   /**
    * Number of points of the lattice per unit volume.
    */
   MScal m_m;

   /**
    * Rank of the lattice.
    */
   int m_rank;

   /**
    * Only elements 1 to <tt>m_maxDim</tt> (inclusive) of arrays are
    * defined.
    */
   int m_maxDim;

   /**
    * Beta factor.
    */
   double m_beta;

   /**
    * Contains the bounds on the length of the shortest nonzero vector in
    * the lattice in each dimension.
    */
   double *m_cst;
private:

   /**
    * Use of the copy-constructor is forbidden.
    */
   Normalizer (const Normalizer &);

   /**
    * Use of assigment is forbidden.
    */
   Normalizer & operator= (const Normalizer &);
};

}
#endif
