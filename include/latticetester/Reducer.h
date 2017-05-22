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

#ifndef LATTICETESTER__REDUCER_H
#define LATTICETESTER__REDUCER_H

#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latticetester/Util.h"
//#include "latticetester/Basis.h"
#include "latticetester/IntLatticeBasis.h"

#include <fstream>
#include <sstream>
#include <vector>


namespace LatticeTester {

/**
 * For a given lattice, this class implements methods to reduce its basis in
 * the sense of Minkowski and to find the shortest non-zero vector of the
 * lattice using pre-reductions and a branch-and-bound (BB) algorithm
 * \cite rLEC97c&thinsp;. It also implements the method of Lenstra, Lenstra
 * and Lovasz (LLL) \cite mLEN82a&thinsp; as well as the method of Dieter
 * \cite rDIE75a&thinsp; to reduce a lattice basis. The reduction algorithms
 * in this class use both the primal and the dual lattices, so both lattices
 * must be defined.
 *
 */
class Reducer {
public:

   /**
    * Whenever the number of nodes in the branch-and-bound tree exceeds
    * <tt>SHORT_DIET</tt> in the method `ShortestVector`, `PreRedDieterSV` is
    * automatically set to `true` for the next call; otherwise it is set to
    * `false`.
    */
   static const long SHORT_DIET = 1000;

   /**
    * Whenever the number of nodes in the branch-and-bound tree exceeds
    * <tt>SHORT_LLL</tt> in the method `ShortestVector`, `PreRedLLLSV` is
    * automatically set to `true` for the next call; otherwise it is set
    * to `false`.
    */
   static const long SHORT_LLL = 1000;

   /**
    * Whenever the number of nodes in the branch-and-bound tree exceeds
    * <tt>MINK_LLL</tt> in the method <tt>reductMinkowski</tt>,
    * `PreRedLLLRM` is automatically set to `true` for the next call;
    * otherwise it is set to `false`.
    */
   static const long MINK_LLL = 500000;

   /**
    * Maximum number of transformations in the method `PreRedDieter`.
    * After <tt>MAX_PRE_RED</tt> successful transformations have been
    * performed, the prereduction is stopped.
    */
   static const long MAX_PRE_RED = 1000000;

   /**
    * The maximum number of nodes in the branch-and-bound tree when
    * calling `ShortestVector` or `reductMinkowski`. When this number is
    * exceeded, the method aborts and returns `false`.
    */
   static long maxNodesBB;

   /**
    * \name Pre-reduction flags
    * @{
    *
    * These boolean variables indicate which type of pre-reduction is to be
    * performed for `ShortestVector` (SV) and for `reductMinkowski` (RM).
    * `Dieter` means the pairwise pre-reduction as in the method `PreRedDieter`.
    * `LLL` means the LLL reduction of Lenstra, Lenstra, and Lovász. The
    * variable `PreRedDieterSV` is originally set to `true` and the two others
    * are originally set to `false`. These variables are reset automatically
    * depending on the thresholds `MinkLLL, ShortDiet, ShortLLL` as explained
    * above.
    */
   static bool PreRedDieterSV;
   static bool PreRedLLLSV;
   static bool PreRedLLLRM;

   /**
    * @}
    */

   /**
    * Constructor. Initializes the reducer to work on lattice `lat`.
    */
   Reducer (IntLatticeBasis & lat);

   /**
    * Copy constructor.
    */
   Reducer (const Reducer & red);

   /**
    * Destructor.
    */
   ~Reducer ();

   /**
    * Assignment operator.
    */
   Reducer & operator= (const Reducer & red);

   /**
    * Copies the reducer `red` into this object.
    * \remark **Richard:** Encore utile?
    */
   void copy (const Reducer & red);

   /**
    * Method used in `reductMinkowski` to perform a transformation of
    * stage 3 described in \cite rAFF85a&thinsp;. Also used in
    * `ShortestVector`. Assumes that \f$\sum_{i=1}^t z_i V_i\f$ is a
    * short vector that will enter the basis. Tries to reduce some vectors
    * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
    * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
    * \f$V_j\f$ when this happens. Returns in \f$k\f$ the last index
    * \f$j\f$ such that \f$|z_j|=1\f$.
    */
   void transformStage3 (std::vector<long> & z, int & k);

   /**
    * Finds a Choleski decomposition of the basis. Returns in `C0` the
    * elements of the upper triangular matrix of the Choleski
    * decomposition that are above the diagonal. Returns in `DC2` the
    * squared elements of the diagonal.
    */
   bool calculCholeski (RVect & DC2, RMat & C0);

   /**
    * Tries to find shorter vectors in `reductMinkowski`.
    */
   bool tryZ  (int j, int i, int Stage, bool & smaller, const BMat & WTemp);

   /**
    * Tries to find a shorter vector in `shortestVector`.
    */
   bool tryZ0 (int j, bool & smaller);

   /**
    * Computes the shortest non-zero vector of this lattice with respect
    * to norm `norm` using branch-and-bound and the algorithm described in
    * \cite rLEC97c&thinsp;. The `Norm` member of this object will be
    * changed to `norm`. If `MaxNodesBB` is exceeded during one of the
    * branch-and-bounds, the method aborts and returns `false`. Otherwise,
    * it returns `true`. Uses the pre-reduction algorithms of Dieter and
    * of Lenstra-Lenstra-Lovasz.
    */
   bool shortestVector (NormType norm);

   /**
    * Similar to `ShortestVector` but uses the algorithm of Dieter
    * \cite rDIE75a, \cite rKNU98a&thinsp;.
    */
   bool shortestVectorDieter (NormType norm);

   /**
    * Tries to shorten the vectors of the primal basis using
    * branch-and-bound, in `reductMinkowski`.
    */
   bool redBB (int i, int d, int Stage, bool & smaller);

   /**
    * Tries to shorten the smallest vector of the primal basis using
    * branch-and-bound, in `ShortestVector`.
    */
   bool redBB0 (NormType norm);

   /**
    * Performs the reductions of the preceding two methods using
    * cyclically all values of \f$i\f$ (only for \f$i > d\f$ in the latter
    * case) and stops after either `MaxPreRed` successful transformations
    * have been achieved or no further reduction is possible. Always use
    * the Euclidean norm.
    */
   void preRedDieter (int d);
   void preRedDieterPrimalOnly (int d);
   void preRedDieterPrimalOnlyRandomized (int d, int seed);


   /**
    * Finds the shortest non-zero vector using norm `norm`. Returns `true`
    * upon success. Uses the algorithm of Dieter \cite rDIE75a&thinsp;
    * given in Knuth \cite rKNU98a&thinsp;.
    */
   bool redDieter (NormType norm);

   /**
    * Performs a LLL (Lenstra-Lenstra-Lovasz) basis reduction up to
    * dimension `dim` with coefficient `fact`, which must be smaller than
    * 1. If `fact` is closer to 1, the basis will be (typically) "more
    * reduced", but that will require more work. The reduction algorithm
    * is applied until `maxcpt` successful transformations have been done.
    * Always uses the Euclidean norm.
    */
   void redLLL (double fact, long maxcpt, int dim);

   /**
    * With NTL
    */
   void redLLLNTLProxy(double fact);
   void redLLLNTLExact(ZZ & det, long a, long b);
   void redBKZ(double fact, long Blocksize);

   /**
    * Reduces the current basis to a Minkowski reduced basis with respect
    * to the Euclidean norm, assuming that the first \f$d\f$ vectors are
    * already reduced and sorted. If `MaxNodesBB` is exceeded during one
    * of the branch-and-bound step, the method aborts and returns `false`.
    * Otherwise it returns `true`, the basis is reduced and sorted by
    * vector lengths (the shortest vector is `V[1]` and the longest is
    * <tt>V[Dim]</tt>).
    */
   bool reductMinkowski (int d);

   /**
    * Performs pairwise reductions. This method tries to reduce each basis
    * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
    * adding to it a multiple of the \f$i\f$-th vector. Always uses the
    * Euclidean norm.
    */
   void pairwiseRedPrimal (int i, int d);

   /**
    * Performs pairwise reductions, trying to reduce every other vector of
    * the *dual* basis by adding multiples of the \f$i\f$-th vector. That
    * may change the \f$i\f$-th vector in the primal basis. Each such dual
    * reduction is actually performed only if that does not increase the
    * length of vector \f$i\f$ in the primal basis. Always uses the
    * Euclidean norm.
    */
   void pairwiseRedDual (int i);

   /**
    * Returns the length of the shortest basis vector in the lattice.
    */
   NScal getMinLength () {
      if (m_lat->getNorm() == L2NORM)
         return sqrt(m_lMin2);
      else return m_lMin; }

   /**
    * Sets the lower bound on the square length of the shortest vector in
    * dimension \f$i\f$ to \f$V[i]\f$, for \f$i\f$ going from `dim1` to
    * `dim2`.
    */
   void setBoundL2 (const NVect & V, int dim1, int dim2);

   /**
    * Debug function that print the primal and dual bases.
    */
   void trace (char *mess);
private:

   /**
    * Lattice on which the reduction will be performed.
    */
   IntLatticeBasis* m_lat;

   /**
    * Permutes the \f$i^{th}\f$ and the \f$j^{th}\f$ line, and the
    * \f$i^{th}\f$ and the \f$j^{th}\f$ column of the scalar product’s
    * matrix.
    */
   void permuteGramVD (int i, int j, int n);

   /**
    * Recalculates the first \f$n\f$ entries of the \f$j^{th}\f$ column of
    * the Choleski’s matrix of order 2.
    */
   void calculCholeski2LLL (int n, int j);

   /**
    * Recalculates the entry (\f$i\f$, \f$j\f$) of the Choleski’s matrix
    * of order 2.
    */
   void calculCholeski2Ele (int i, int j);
   void miseAJourGramVD (int j);
   void calculGramVD ();

   /**
    * Reduce the Choleski matrix with adding a multiple of the i-th vector
    * to the j-th vector. It updates the Gram Schmidt matrix
    */
   void reductionFaible (int i, int j);

   /**
    * Lower bound on the length squared of the shortest vector in each
    * dimension. If any vector of the lattice is shorter than this bound,
    * we stop the reduction immediately and reject this lattice since its
    * shortest vector will be even smaller.
    */
   NVect m_BoundL2;

   /**
    * \name Work variables.
    * @{
    */
   BScal m_bs;
   BVect m_bv;            // Saves shorter vector in primal basis
   BVect m_bw;            // Saves shorter vector in dual basis

   NScal m_lMin, m_lMin2;
   NScal m_ns;
   NVect m_nv;

   //RScal m_rs;
   RVect m_zLR, m_n2, m_dc2;
   RMat m_c0, m_c2, m_cho2, m_gramVD;
   int *m_IC;             // Indices in Cholesky

   std::vector<long> m_zLI;
   std::vector<long> m_zShort;
   long m_countNodes;     // Counts number of nodes in the branch-and-bound tree
   long m_countDieter;    // Counts number of attempts since last successful
                          // Dieter transformation
   long m_cpt;            // Counts number of successes in pre-reduction
                          // transformations
   bool m_foundZero;      // = true -> the zero vector has been found
// long m_BoundL2Count;   // Number of cases for which the reduction stops
                          // because any vector is shorter than L2 Bound.
   /**
    * @}
    */
};

}     // namespace LatticeTester
#endif // REDUCER_H
