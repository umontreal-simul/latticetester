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

#ifndef LATTICETESTER__BASIS_H
#define LATTICETESTER__BASIS_H

#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include <string>


namespace LatticeTester {

/**
 * This class represents a basis for a lattice. To compute the length of
 * vectors, one may use either the \f$L_1\f$, the \f$L_2\f$ or the
 * \f$L_{\infty}\f$ norms. If one uses the \f$L_2\f$ norm and if `BScal` and
 * `NScal` are of type `double`, then the norm could overflow if the
 * components of the basis vectors are larger than \f$\approx2^{500}\f$. In
 * that case, `NScal` may be chosen as `RR`. If one uses the \f$L_1\f$ or the
 * \f$L_{\infty}\f$ norms, `NScal` may be chosen the same as `BScal`.
 *
 */

/* Erwan
Cette classe représente une Basis pour un réseau. Elle hérite de BMat (aka
NTL::matrix<long> ou boost::numeric::ublas::matrix<long>) définit dans Types.h.
Elle hérite donc d'une structure de matrice : une Basis A contient des éléments
A(i,j) (j-ème composante du i-ème vecteur).

Ces éléments sont modifiés grâces aux méthodes fournies dans IntLattice.
*/

class Basis : public BMat {
public:

   /**
    * Constructor. Builds a basis of actual dimension \f$d\f$, maximum dimension
    * `maxDim` and with norm `norm`.
    */
   Basis (int dim, int maxDim, NormType norm = SUPNORM);

   /** Erwan
    * Constructor. Builds a basis of actual dimension \f$d\f$,
    * and with norm `norm`.
    */
   Basis (int dim, NormType norm = SUPNORM);

   /**
    * Copy constructor.
    */
   Basis (const Basis &);

   /**
    * Destructor.
    */
   ~Basis ();

   /**
    * Cleans and releases all memory used by the basis.
    */
   void kill ();

   /**
    * Assignment operator.
    */
   Basis & operator= (const Basis &);

   /**
    * Swaps this basis with the basis `b`.
    */
   void swap (Basis & b);

   /**
    * Exchanges vectors \f$i\f$ and \f$j\f$ in the basis.
    */
   void permute (int i, int j);

   /**
    * Sets the actual dimension of the basis to `d`.
    */
   void setDim (int d);

   /**
    * Sets the norm used by the basis.
    */
   void setNorm (NormType norm);

   /**
    * Sets the \f$i\f$-th vector’s norm to `value`. The negative flag for
    * this vector is set to false and no data integrity is made.
    */
   void setVecNorm (NScal & value, int i);

   /**
    * Recalculates the norm of each vector in the basis.
    */
   void updateVecNorm ();

   /**
    * Same as above, except that the recalculation begins at dimension `d`
    *
    */
   void updateVecNorm (int d);

   /**
    * Updates the norm of vector at dimension `d` using the `L2NORM`.
    */
   void updateScalL2Norm (int d);

   /**
    * Updates the norm of all basis vectors from dimensions `d1` to `d2`
    * (inclusive) using the `L2NORM`.
    */
   void updateScalL2Norm (int d1, int d2);

   /**
    * Sets the dirty flag of the norms to `flag`.
    */
   void setNegativeNorm (bool flag);

   /**
    * Writes the basis in a string and returns it.
    */
   std::string toString () const;

   /**
    * Writes the \f$i\f$-th basis vector in a string and returns it.
    */
   std::string toString (int i) const;

   /**
    * Writes the basis on standard output.
    */
   void write () const;

   /**
    * Writes \f$i\f$-th basis vector on standard output.
    */
   void write (int i) const;

   /**
    * Returns the actual dimension of the basis.
    */
   int getDim () const { return m_dim; }

   /**
    * Returns the maximal dimension of the basis.
    */
   int getMaxDim () const { return m_maxDim; }

   /**
    * Returns the norm used by the basis.
    */
   NormType getNorm () const { return m_norm; }

   /**
    * Returns `true` if the \f$i\f$-th vector’s norm is not updated and
    * erroneous. Returns `false` otherwise.
    */
   bool isNegativeNorm (int i) const { return m_negFlag[i]; }

   /**
    * Sets the negative flag for the \f$j\f$-th vector to `flag`.
    */
   void setNegativeNorm (bool flag, int j) { m_negFlag[j] = flag; }

   /**
    * Returns the \f$i\f$-th vector’s norm.
    */
   NScal getVecNorm (int i) const { return m_vecNorm[i]; }

protected:

   /**
    * Actual dimension of the basis.
    */
   int m_dim;

   /**
    * Maximum dimension of the basis.
    */
   int m_maxDim;

   /**
    * Norm used to calculate the norm of the vectors.
    */
   NormType m_norm;

   /**
    * The norm of each vector in the basis.
    */
   NVect m_vecNorm;

   /**
    * Indicates whether a vector norm must be recalculated or not.
    */
   bool* m_negFlag;
};

}
#endif
