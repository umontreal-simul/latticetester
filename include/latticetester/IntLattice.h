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

#ifndef LATTICETESTER__INTLATTICE_H
#define LATTICETESTER__INTLATTICE_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latticetester/Basis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/Lacunary.h"
#include <string>


namespace LatticeTester {

/**
 * This class offers tools to manipulate lattice bases. Each lattice is
 * represented by a basis \f$V\f$ and its dual \f$W\f$. It is sometimes
 * possible, as in the case of lattices associated with LCGs (linear
 * congruential generators) or MRGs (multiple recursive generators), to
 * multiply a lattice (and its dual) by a constant factor \f$m\f$ in such a
 * way that they are included in \f$\mathbb Z^t\f$, allowing exact
 * representation of basis vector coordinates. The duality relation will then
 * read \f$V_i\cdot W_j = m\delta_{ij}\f$ for some integer constant \f$m\f$.
 *
 */

/*
* Cette classe représente une lattice. Une lattice est caractérisée par
* une base et un modulo (m).



*/
class IntLattice {
public:

   /** Erwan
    * Constructor. The modulus is initialized to \f$modulus\f$, and the primal
    * and dual vectors are initialize to \f$base_primal\f$.
    */
   IntLattice (const Base base_primal, int modulus);

   /** Erwan
    * Constructor. The modulus is initialized to \f$modulus\f$, the primal
    * vectors are initialize to \f$base_primal\f$, the dual vectors are
    * initialize to \f$base_dual\f$.
    */
   IntLattice (const Base base_primal, const Base base_dual, int modulus);

   /**
    * Constructor. The modulus is initialized to \f$m\f$, the maximal dimension
    * of the bases is set to `maxDim`, and the order is set to \f$k\f$.
    */
   IntLattice (const MScal & m, int k, int maxDim, NormType norm = L2NORM);

   /**
    * Copy constructor. The maximal dimension of the created basis is set
    * equal to <tt>Lat</tt>’s current dimension.
    */
   IntLattice (const IntLattice & Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set equal to <tt>Lat</tt>’s current dimension.
    */
   IntLattice & operator= (const IntLattice & Lat);

   /**
    * Destructor.
    */
   virtual ~IntLattice ();

   /**
    * Same as assignment operator above.
    */
   void copy (const IntLattice & lattice);

   /**
    * Returns actual dimension `Dim`.
    */
   int getDim () const { return m_v.getDim(); }

   /**
    * Sets actual dimension to \f$d\f$.
    */
   void setDim (int d);

   /**
    * Returns the maximal dimension of the lattice.
    */
   int getMaxDim () const { return m_v.getMaxDim(); }

   /**
    * Increases dimension `Dim` by 1.
    */
   virtual void incDim();

   /**
    * Builds the basis for the lattice in dimension `d` (the basis is m_v)
    */
   virtual void buildBasis (int d);

   /**
    * Returns the norm used by the lattice.
    */
   NormType getNorm () const { return m_v.getNorm(); }

   /**
    * Returns the value of the modulus \f$m\f$ of the recurrence (the
    * number of points of the lattice).
    */
   MScal getM () const { return m_m; }

   /**
    * Returns the square number of points in the lattice.
    */
   NScal getM2 () const { return m_m2; }

   /**
    * Returns the order.
    */
   int getOrder () const { return m_order; }

   /**
    * Creates and returns the normalizer corresponding to criterion
    * `norma`. In the case of the \f$P_{\alpha}\f$ test, the argument
    * `alpha` = \f$\alpha\f$. In all other cases, it is unused.
    */
   Normalizer * getNormalizer (NormaType norma, int alpha);

   /**
    * Returns the vector of multipliers (or coefficients) \f$A\f$ as a
    * string.
    */
   virtual std::string toStringCoef() const;

#ifdef WITH_NTL
   /**
    * Must be implemented in subclasses.
    */
   virtual const MVect & getCoef() const {
      MyExit(1, "IntLattice.getCoef() is empty");  return m_dummyCoef; }
#endif

   /**
    * Returns the primal basis `V`.
    */
   Basis & getPrimalBasis () { return m_v; }

   /**
    * Returns the dual basis `W`.
    */
   Basis & getDualBasis () { return m_w; }
   double getLgVolDual2 (int i) const { return m_lgVolDual2[i]; }

   /**
    * This function is called to fix the normalization constants to get
    * the normalized merit from the shortest distance in the lattice. If
    * `dualF` is `true`, the normalization constant is reset for the dual
    * lattice, otherwise it is reset for the primal lattice.
    */
   void fixLatticeNormalization (bool dualF);
   bool getXX (int i) const { return m_xx[i]; } //???
   void setXX (bool val, int i) { m_xx[i] = val; } //??

   /**
    * Sorts the basis vectors with indices from \f$d+1\f$ to the dimension
    * of the basis by increasing length. The dual vectors are permuted
    * accordingly. Assumes that the lengths of the corresponding basis
    * vectors are up to date.
    */
   void sort (int d);

   /**
    * Exchanges vectors \f$i\f$ and \f$j\f$ in the basis and in the dual
    * basis.
    */
   void permute (int i, int j);

   /**
    * Writes this basis in file named `filename`.
    */
   void write (const char *filename) const;

   /**
    * For debugging purposes.
    */
   void trace (char* msg);

   /**
    * Writes this basis on standard output. If `flag` \f$ > 0\f$, the norm
    * of the basis vectors will be recomputed; otherwise not.
    */
   void write (int flag);

   /**
    * Exchanges basis \f$V\f$ and its dual \f$W\f$.
    */
   void dualize ();

   /**
    * Checks that the bases satisfy the duality relation
    * \f$V[i]\cdot W[j] = m \delta_{ij}\f$. If so, returns `true`,
    * otherwise returns `false`.
    */
   bool checkDuality () const;

   /**
    * Checks that <tt>Lat</tt>’s basis and this basis are equivalent. If
    * so, returns `true`, otherwise returns `false`.
    */
   bool baseEquivalence (const IntLattice & Lat) const;

   /**
    * Builds the basis (and dual basis) of the projection `proj` for this
    * lattice. The result is placed in the `lattice` lattice. The basis is
    * triangularized to form a proper basis.
    */
   void buildProjection (IntLattice* lattice, const Coordinates & proj);

   /**
    * Does nothing is this base class.
    */
   virtual MScal getRho() const {
       MyExit(1, "IntLattice.getRho is empty");  return m_t1; }

   /**
    * Does nothing is this base class.
    */
   virtual MScal getLossRho() const {
       MyExit(1, "IntLattice.getLossRho is empty");   return m_t1; }

   /**
    * Does nothing is this base class.
    */
   virtual void setLac (const Lacunary &) {
          MyExit(1, "IntLattice.setLac is empty"); }

protected:

   void init ();

   /**
    * Cleans and releases all the memory allocated for this lattice.
    */
   void kill ();

   /**
    * The order of the basis.
    */
   int m_order;

   /**
    * Number of points per unit volume.
    */
   MScal m_m;

   /**
    * Square of the number of points per unit volume.
    */
   NScal m_m2;

   /**
    * The logarithm \f$\log_2 (m^2)\f$.
    */
   double m_lgm2;

   /**
    * Primal basis of the lattice.
    */
   Basis m_v;

   /**
    * Dual basis of the lattice.
    */
   Basis m_w;
   double *m_lgVolDual2;

   /**
    * Computes the logarithm of the normalization factor
    * (<tt>m_lgVolDual2</tt>) in all dimensions \f$\le\f$ `MaxDim` for
    * the lattice. `lgm2` is the logarithm in base 2 of \f$m^2\f$.
    */
   void calcLgVolDual2 (double lgm2);

   bool *m_xx;
   MScal m_t1, m_t2, m_t3;
   BMat m_vSI;

   /**
    * Work variables.
    */
   Basis m_vTemp;

   /**
    * Default return value for getCoef() (empty vector).
    */
   MVect m_dummyCoef;
};

}
#endif
