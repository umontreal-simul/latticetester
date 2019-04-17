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

#ifndef LATTICETESTER_REDUCER_H
#define LATTICETESTER_REDUCER_H

#include "NTL/LLL.h"

#include "latticetester/Const.h"
#include "latticetester/Util.h"
#include "latticetester/IntLatticeBasis.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>

namespace LatticeTester {

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class Reducer;

  /// \cond specReducerDec
  /*
   * This struct specializes some of the functions in a `Reducer`. This is a
   * workaround needed, implementation wise, since you can't specialize member
   * functions of a class without specializing the whole class. 
   *
   * What this does is that this structure contains specific functions of the 
   * original class Reducer that will act differently depending on the types of
   * the arguments. Instead of specializing all the methods of Reducer for all 
   * our use cases (since some of them depend on the type) we can specialize
   * only this structure.
   * */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      struct specReducer {
        void redLLLNTLExact(Reducer<Int, BasInt, Dbl, RedDbl>& red, double fact);
        void redBKZ(Reducer<Int, BasInt, Dbl, RedDbl>& red, double fact,
            std::int64_t blocksize, PrecisionType precision, int dim);
        void redLLLNTL(Reducer<Int, BasInt, Dbl, RedDbl>& red, double fact,
            PrecisionType precision, int dim);
      };
  /// \endcond

  /**
   * This class implements (or wraps from NTL) all the functions that are
   * needed to reduce a basis.
   * For a given lattice basis, stored in an IntLatticeBasis, this class can
   * reduce it in the sense of Minkowski, as well as find the shortest vector
   * in the lattice with the Branch-and-Bound algorithm. It is also possible to
   * use weaker, but much faster, reductions such as pairwise-reduction 
   * \cite rDIE75a, LLL reduction \cite mLEN82a and BKZ reduction \cite mSCH91a.
   * The theoretical details of these algorithm are presented in \ref a_intro.
   *
   * Unless using Dieter or Minkowski reduction, the methods of this class do not
   * change the dual of the lattice that is being reduced.
   *
   * The LLL reduction has been implemented both by us and in NTL, but since we
   * did not test our algorithm for efficiency, it is recommended to stick with
   * the NTL implementation wrapped by the methods redLLLNTL and redLLLNTLExact.
   * Note that redBKZ is also a wrapper for NTL algorithm for BKZ reduction.
   *
   * To use this class, it suffices to create an instance of it by passing a
   * `IntLatticeBasis` to the constructor. You can then simply call the methods
   * and the `IntLatticeBasis` passed to the Reducer will be modified since it
   * is stored as a pointer. You should note that the methods to compute the 
   * shortest vector apply no pre-reduction. In a real context, you should
   * always reduce the basis separately with LLL or BKZ reduction before
   * searching for the shortest vector since it drastically reduces the size of
   * the branch-and-bound search.
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class Reducer {
        private:
          // Local typedefs for matrix and vector types needed in the class.
          typedef NTL::vector<BasInt> BasIntVec;
          typedef NTL::matrix<BasInt> BasIntMat;
          typedef NTL::vector<Dbl> DblVec;
          typedef NTL::vector<RedDbl> RedDblVec;
          typedef NTL::matrix<RedDbl> RedDblMat;
        public:

          /**
           * The maximum number of nodes in the branch-and-bound tree when
           * calling `ShortestVector` or `reductMinkowski`. When this number is
           * exceeded, the method aborts and returns `false`.
           */
          static std::int64_t maxNodesBB;

          /**
           * \name Pre-reduction flags
           * @{
           * **Marc-Antoine**: it is very possible that we can remove those variables.
           * We just have to make sure that reductMinkowski does not need them 
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
          static bool PreRedDieterSV; // Not used
          static bool PreRedLLLSV; // Not used
          static bool PreRedLLLRM;
          static bool PreRedBKZ; // Not used

          /**
           * @}
           */

          /**
           * Constructor that initializes the reducer to work on lattice basis `lat`.
           */
          Reducer (IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat);

          /**
           * Copy constructor.
           */
          Reducer (const Reducer<Int, BasInt, Dbl, RedDbl> & red);

          /**
           * Destructor.
           */
          ~Reducer ();

          /**
           * Assignment operator that allows copying `red` into the object by 
           * usig copy.
           */
          Reducer<Int, BasInt, Dbl, RedDbl> & operator= (
              const Reducer<Int, BasInt, Dbl, RedDbl> & red);

          /**
           * Copies the reducer `red` into this object.
           */
          void copy (const Reducer<Int, BasInt, Dbl, RedDbl> & red);

          /**
           * Computes the shortest non-zero vector of the IntLatticeBasis stored
           * in this object with respect to norm `norm` using branch-and-bound 
           * and the algorithm described in \cite rLEC97c. The `Norm` member of
           * this object will be changed to `norm`. If `MaxNodesBB` is exceeded
           * during one of the branch-and-bounds, the method aborts and returns
           * `false`. Otherwise, it returns `true`. If the reduction was
           * successful, the new reduced basis can be accessed as desired via
           * getIntLatticeBasis().
           *
           * It is strongly recommended to use some kind of pre-reduction before
           * using this method. We suggest redLLLNTL that should be god for most
           * use cases, but redBKZ can sometimes git better results. None of the
           * Dieter pre-reduction should be called before calling this, unless
           * the user's intention is to benchmark the different methods in his
           * use case.
           */
          bool shortestVector (NormType norm);

          /**
           * This method performs pairwise reduction sequentially on all vectors
           * of the basis whose indice is greater of equal to `d`.
           *
           * `xx[]` is used internally when doing Minkowski reduction. This is
           * somewhat ugly but this fixes a leak.
           */
          void redDieter (int d, bool xx[] = NULL);


          /**
           * This method does the same thing as redDieter(int) but the choice of
           * vectors on which to perform pairwise reduction is randomized.
           */
          void redDieterRandomized (int d, int seed);

          /**
           * Performs a LLL (Lenstra-Lenstra-Lovasz) basis reduction for the 
           * first `dim` vector of the basis with coefficient `fact`.
           * `fact` is a floating point number between 1/2 and 1 that
           * specifies the algorithm's tolerance to floating point arithmetic
           * error. If `fact` is closer to 1,
           * the basis will be (typically) "more reduced" as the
           * algorithm will enforce a tighter condition on the basis,
           * but that will require more work. The reduction algorithm is
           * applied until `maxcpt` successful transformations have been done,
           * or until the basis is correctly reduced.
           * 
           * This algorithm is an implementation of the LLL reduction
           * algorithm made in our lab. It is considerably slower than what is
           * available through NTL for NTL types, but performs well on standard
           * C++ types. This implementation
           * always uses the Euclidean norm.
           */
          void redLLL (double fact = 0.999999, std::int64_t maxcpt = 1000000000, int dim = 0);

          /**
           * This is the NTL implementation of the floating point version of the
           * LLL reduction algorithm presented in \cite mSCH91a.
           * 
           * `1/2 < fact < 1` is a number that specifies the strenght of the 
           * condition on the basis, closer to 1 meaning a stronger condition.
           * `precision` specifies the size of the floating point numbers the 
           * algorithm will use. Const.h provides a list of the possible values,
           * but their description is done in the module LLL of NTL. `dim` 
           * is the number of vectors on which to apply the reduction (zero being
           * the entire basis).
           *
           * There should not really be a reason to change the default
           * parameters suggested in this function. It is nonetheless possible
           * to be in a situation where they do not work. In the occurence of 
           * poor quality results, it might be necessary to choose a better
           * `precision`. Also, in a few rare cases, it is possible that fact 
           * be "too close to 1" and that the algorithm takes too much time.
           */
          void redLLLNTL(double fact = 0.999999, 
              PrecisionType precision = QUADRUPLE, int dim = 0);

          /**
           * This implements an algorithm to perform the original LLL reduction.
           */
          void redLLLNTLExact(double fact);

          /**
           * This is the NTL implementation of the floating point version of the
           * BKZ reduction algorithm presented in \cite mSCH91a.
           * 
           * `1/2 < fact < 1` is a number that specifies the strenght of the 
           * condition on the basis, closer to 1 meaning a stronger condition.
           * `precision` specifies the size of the floating point numbers the 
           * algorithm will use. Const.h provides a list of the possible values,
           * but their description is done in the module LLL of NTL. `dim` 
           * is the number of vectors on which to apply the reduction with zero being
           * the entire basis. `blocksize` is the size of the blocks of the
           * reduction in the BKZ reduction condition. Roughly, larger blocks
           * mean a stronger condition. A `blocksize` of 2 is equivalent to LLL
           * reduction.
           *
           * There should not really be a reason to change the default
           * parameters suggested in this function. It is nonetheless possible
           * to be in a situation where they do not work. In the occurence of 
           * poor quality results, it might be necessary to choose a better
           * `precision`. Also, in a few rare cases, it is possible that fact 
           * be "too close to 1" and that the algorithm takes too much time.
           * A user whishing to get a faster execution at the cost of condition
           * strenght could also choose a smaller `blocksize`.
           */
          void redBKZ(double fact = 0.999999, std::int64_t blocksize = 10,
              PrecisionType precision = QUADRUPLE, int dim = 0);

          /**
           * Reduces the current basis to a Minkowski reduced basis with respect
           * to the Euclidean norm, assuming that the first \f$d\f$ vectors are
           * already reduced and sorted. If `MaxNodesBB` is exceeded during one
           * of the branch-and-bound step, the method aborts and returns `false`.
           * Otherwise it returns `true`, the basis is reduced and sorted by
           * vector lengths (the shortest vector is `V[0]` and the std::int64_test is
           * <tt>V[Dim-1]</tt>).
           */
          bool reductMinkowski (int d);

          /**
           * Returns the length of the shortest basis vector in the lattice.
           */
          Dbl getMinLength () {
            if (m_lat->getNorm() == L2NORM)
              return sqrt(m_lMin2);
            else return m_lMin; }

          /**
           * Returns the length of the longest basis vector in the lattice.
           */
          Dbl getMaxLength () {
            if (m_lat->getNorm() == L2NORM)
              return sqrt(m_lat->getVecNorm(m_lat->getDim() - 1));
            else return m_lat->getVecNorm(m_lat->getDim() - 1); }

          /**
           * Sets a bound on the square length of the shortest vector in the 
           * lattice in dimensions from `dim1+1` to `dim2`. `V[i]` has to 
           * contain a lower bound on the square of the lenght of the shortest
           * vector in the lattice in dimension `i+1`. Such a bound, if it is
           * set, will be used during during the Branch-and-Bound step when
           * searching the sortest vector of a lattice. If the Branch-and-Bound
           * proves that the shortest vector of the lattice is smaller than the
           * bound, the algorithm will stop. This is usefull in the case where
           * the user is searching for a good lattice with the spectral test
           * since this can cut off some work.
           */
          void setBoundL2 (const DblVec & V, int dim1, int dim2);

          /**
           * Returns the basis of this object is working on.
           * */
          IntLatticeBasis<Int, BasInt, Dbl, RedDbl>* getIntLatticeBasis ()
            {
              return m_lat;
            }

        private:

          /**
           * Lattice on which the reduction will be performed.
           */
          IntLatticeBasis<Int, BasInt, Dbl, RedDbl>* m_lat;

          /** 
           * Contains specialized implementations of member methods depending on
           * the types given to the `Reducer` template. This is private for 
           * because it is possible to call specialized method at an higher 
           * level.
           * */
          struct specReducer<Int, BasInt, Dbl, RedDbl> spec;

          /**
           * Method used in `reductMinkowski` to perform a transformation of
           * stage 3 described in \cite rAFF85a. It is also used in
           * `ShortestVector`. Assumes that \f$\sum_{i=1}^t z_i V_i\f$ is a
           * short vector that will enter the basis. Tries to reduce some vectors
           * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
           * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
           * \f$V_j\f$ when this happens. Returns in \f$k\f$ the last index
           * \f$j\f$ such that \f$|z_j|=1\f$.
           */
          void transformStage3 (std::vector<std::int64_t> & z, int & k);

          /**
           * Performs pairwise reductions. This method tries to reduce each basis
           * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
           * adding to it a multiple of the \f$i\f$-th vector. Always uses the
           * Euclidean norm.
           */
          void pairwiseRedPrimal (int i, int d, bool xx[] = NULL);

          /**
           * Performs pairwise reductions, trying to reduce every other vector of
           * the *dual* basis by adding multiples of the \f$i\f$-th vector. That
           * may change the \f$i\f$-th vector in the primal basis. Each such dual
           * reduction is actually performed only if that does not increase the
           * length of vector \f$i\f$ in the primal basis. Always uses the
           * Euclidean norm.
           */
          void pairwiseRedDual (int i, bool xx[] = NULL);

          /**
           * Tries to shorten the vectors of the primal basis using
           * branch-and-bound, in `reductMinkowski`.
           */
          bool redBB (int i, int d, int Stage, bool & smaller, bool xx[] = NULL);

          /**
           * Tries to shorten the smallest vector of the primal basis using
           * branch-and-bound, in `ShortestVector`.
           */
          bool redBB0 (NormType norm);

          /**
           * Tries to find shorter vectors in `reductMinkowski`.
           */
          bool tryZ  (int j, int i, int Stage, bool & smaller, 
              const BasIntMat & WTemp);

          /**
           * Tries to find a shorter vector in `shortestVector`.
           */
          bool tryZ0 (int j, bool & smaller);

          /**
           * Finds a Choleski decomposition of the basis. Returns in `C0` the
           * elements of the upper triangular matrix of the Choleski
           * decomposition that are above the diagonal. Returns in `DC2` the
           * squared elements of the diagonal.
           */
          bool calculCholeski (RedDblVec & DC2, RedDblMat & C0);

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
           * shortest vector will be even smaller. Used only in Seek
           */
          DblVec m_BoundL2;

          /**
           * Whenever the number of nodes in the branch-and-bound tree exceeds
           * <tt>SHORT_DIET</tt> in the method `ShortestVector`, `PreRedDieterSV` is
           * automatically set to `true` for the next call; otherwise it is set to
           * `false`.
           */
          static const std::int64_t SHORT_DIET = 1000;

          /**
           * Whenever the number of nodes in the branch-and-bound tree exceeds
           * <tt>SHORT_LLL</tt> in the method `ShortestVector`, `PreRedLLLSV` is
           * automatically set to `true` for the next call; otherwise it is set
           * to `false`.
           */
          static const std::int64_t SHORT_LLL = 1000;

          /**
           * Whenever the number of nodes in the branch-and-bound tree exceeds
           * <tt>MINK_LLL</tt> in the method <tt>reductMinkowski</tt>,
           * `PreRedLLLRM` is automatically set to `true` for the next call;
           * otherwise it is set to `false`.
           */
          static const std::int64_t MINK_LLL = 500000;

          /**
           * Maximum number of transformations in the method `PreRedDieter`.
           * After <tt>MAX_PRE_RED</tt> successful transformations have been
           * performed, the prereduction is stopped.
           */
          static const std::int64_t MAX_PRE_RED = 1000000;

          /*
           * Work variables.
           */
          BasInt m_bs;
          BasIntVec m_bv;            // Saves shorter vector in primal basis
          BasIntVec m_bw;            // Saves shorter vector in dual basis

          Dbl m_lMin;          // The norm of the shorter vector in the primal basis
          // according to the norm considered
          Dbl m_lMin2;         // The norm of the shorter vector in the primal basis
          // according to the norm L2
          Dbl m_ns;
          DblVec m_nv;

          //RedDbl m_rs;
          RedDblVec m_zLR, m_n2, m_dc2;
          RedDblMat m_c0, m_c2, m_cho2, m_gramVD;
          int *m_IC;             // Indices in Cholesky

          std::vector<std::int64_t> m_zLI;
          std::vector<std::int64_t> m_zShort;
          std::int64_t m_countNodes;     // Counts number of nodes in the branch-and-bound tree
          std::int64_t m_countDieter;    // Counts number of attempts since last successful
          // Dieter transformation
          std::int64_t m_cpt;            // Counts number of successes in pre-reduction
          // transformations
          bool m_foundZero;      // = true -> the zero vector has been found
          // std::int64_t m_BoundL2Count;   // Number of cases for which the reduction stops
          // because any vector is shorter than L2 Bound.

          /**
           * Debug function that print the primal and dual bases.
           */
          void trace (char *mess);

      }; // End class Reducer

  //============================================================================

  /// \cond specReducerSpec
  // Structure specialization

  // Specialization for the case LLXX
  template<typename Dbl, typename RedDbl>
      struct specReducer<std::int64_t, std::int64_t, Dbl, RedDbl> {

      void redLLLNTLExact(Reducer<std::int64_t, std::int64_t, Dbl, RedDbl>& red,
          double fact)
      {
        std::cout << "**** WARNING redLLLNTLExact cannot be used with integers (std::int64_t)\n";
        std::cout << "**** in NTL nor without NTL. LLL reduction is used with our algorithm \n";
        std::cout << "**** redLLL which does not do the same thing.\n";
        red.redLLL(fact, 1000000, red.getIntLatticeBasis()->getDim ());
      }

      void redBKZ(Reducer<std::int64_t, std::int64_t, Dbl, RedDbl>& red,
          double fact, std::int64_t blocksize, PrecisionType precision, int dim)
      {
          IntLatticeBasis<std::int64_t, std::int64_t, Dbl, RedDbl> *lattmp = 0;
          if (dim > 0) {
            lattmp = new IntLatticeBasis<std::int64_t, std::int64_t, Dbl, RedDbl>(
                       dim, red.getIntLatticeBasis()->getNorm());
            lattmp->copyBasis(*red.getIntLatticeBasis(), dim);
          }
          else
            lattmp = red.getIntLatticeBasis();
          std::cout << "\n**** WARNING redBKZ cannot be use with std::int64_t type for integers\n";
          std::cout << "**** and requires ZZ type. Instead, LLL reduction is\n";
          std::cout << "**** performed with the algorithm redLLL which does not do the same thing.\n";
          std::cout << std::endl;
          red.redLLL(fact, 1000000, red.getIntLatticeBasis()->getDim ());
          if (dim > 0) delete lattmp; 
      }

      /* This is a test to see if making the program promote to NTL types and
       * performs NTL LLL reduction is faster than our own LLL.
       * If it is, we should reimplement our LLL to match what is done in NTL.
       * */
      void redLLLNTL(Reducer<std::int64_t, std::int64_t, Dbl, RedDbl>& red,
          double fact, PrecisionType precision, int dim)
       {
         IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> *lattmp = 0;
         NTL::matrix<std::int64_t> basis = red.getIntLatticeBasis()->getBasis();
         NTL::mat_ZZ U;
         U.SetDims(basis.NumRows(), basis.NumCols());
         for (int i = 0; i < basis.NumRows(); i++) {
           for (int j = 0; j < basis.NumCols(); j++) {
             U[i][j] = basis[i][j];
           }
         }
         lattmp = new IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>(U, dim,
             red.getIntLatticeBasis()->getNorm());
         U.kill();

         switch(precision){
           case DOUBLE:
             LLL_FP(lattmp->getBasis(), U, fact, 0, 0);
             break;
           case QUADRUPLE:
             LLL_QP(lattmp->getBasis(), U, fact, 0, 0);
             break;
           case EXPONENT:
             LLL_XD(lattmp->getBasis(), U, fact, 0, 0);
             break;
           case ARBITRARY:
             LLL_RR(lattmp->getBasis(), U, fact, 0, 0);
             break;
           case EXACT:
             break;
           default:
             MyExit(1, "LLL PrecisionType:   NO SUCH CASE");
         }

         for (int i = 0; i < basis.NumRows(); i++) {
           for (int j = 0; j < basis.NumCols(); j++) {
             red.getIntLatticeBasis()->getBasis()[i][j] = NTL::trunc_long(lattmp->getBasis()[i][j], 63);
             red.getIntLatticeBasis()->getBasis()[i][j] *= NTL::sign(lattmp->getBasis()[i][j]);
           }
         }
         red.getIntLatticeBasis()->updateVecNorm();

         delete lattmp;
       }
      /*
       * void redLLLNTL(Reducer<std::int64_t, std::int64_t, Dbl, RedDbl>& red,
       *     double fact, PrecisionType precision, int dim)
       * {
       *   IntLatticeBasis<std::int64_t, std::int64_t, Dbl, RedDbl> *lattmp = 0;
       *   if(dim > 0){
       *     lattmp = new IntLatticeBasis<std::int64_t, std::int64_t, Dbl, RedDbl>(
       *                dim, red.getIntLatticeBasis()->getNorm());
       *     lattmp->copyBasis(*red.getIntLatticeBasis(), dim);
       *   }
       *   else
       *     lattmp = red.getIntLatticeBasis();
       *   std::cout << "\n**** WARNING redLLLNTL cannot be use with std::int64_t type for integers\n";
       *   std::cout << "**** (LLDD) and requires ZZ type. Instead, LLL reduction is performed\n";
       *   std::cout << "**** with our algorithm which can be much slower.\n";
       *   std::cout << std::endl;
       *   red.redLLL(fact, 1000000, red.getIntLatticeBasis()->getDim ());
       *   if (dim>0) delete lattmp;
       * }
       * */
    };

  // Specialization for the case ZZXX
  template<typename Dbl, typename RedDbl>
      struct specReducer<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> {

      void redLLLNTLExact(Reducer<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> & red,
          double fact)
      {
        /*if (red.getIntLatticeBasis()->withDual()) {
          NTL::mat_ZZ U;
          U.SetDims(red.getIntLatticeBasis()->getBasis().NumRows(), 
              red.getIntLatticeBasis()->getBasis().NumCols());
          NTL::ZZ det(0);
          NTL::LLL(det, red.getIntLatticeBasis()->getBasis(), U, 99999, 100000);
          red.getIntLatticeBasis()->getDualBasis() = 
            transpose(inv(U)) * red.getIntLatticeBasis()->getDualBasis();
        } else{*/
          NTL::ZZ det(0);
          NTL::LLL(det, red.getIntLatticeBasis()->getBasis(), 99999, 100000);
        //}
      }


      void redBKZ(Reducer<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>& red, double fact,
          std::int64_t blocksize, PrecisionType precision, int dim)
      {
        IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> *lattmp = 0;
        if(dim >0){
          lattmp = new IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>(
                     dim, red.getIntLatticeBasis()->getNorm());
          lattmp->copyBasis(*red.getIntLatticeBasis(), dim);
        }
        else
          lattmp = red.getIntLatticeBasis();

        //bool withDual = red.getIntLatticeBasis()->withDual();
        NTL::mat_ZZ U;
        U.SetDims(lattmp->getBasis().size1(), lattmp->getBasis().size2());

        switch(precision){
          case DOUBLE:
            NTL::BKZ_FP(lattmp->getBasis(), U, fact, blocksize);
            break;
          case QUADRUPLE:
            NTL::BKZ_QP(lattmp->getBasis(), U, fact, blocksize);
            break;
          case EXPONENT:
            NTL::BKZ_XD(lattmp->getBasis(), U, fact, blocksize);
            break;
          case ARBITRARY:
            NTL::BKZ_RR(lattmp->getBasis(), U, fact, blocksize);
            break;
          default:
            MyExit(1, "BKZ PrecisionType:   NO SUCH CASE");
        }

        //if (withDual)
          //lattmp->getDualBasis() = transpose(inv(U)) * lattmp->getDualBasis();

        red.getIntLatticeBasis()->copyBasis(*lattmp, dim);

        if (dim>0) delete lattmp;
      }

      void redLLLNTL(Reducer<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>& red, double fact,
          PrecisionType precision, int dim)
      {
        IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> *lattmp = 0;
        if(dim >0){
          lattmp = new IntLatticeBasis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>(
                     dim, red.getIntLatticeBasis()->getNorm());
          lattmp->copyBasis(*red.getIntLatticeBasis(), dim);
        }
        else
          lattmp = red.getIntLatticeBasis();
        if(precision == EXACT)
          red.redLLLNTLExact(fact);
        else{
          //bool withDual = red.getIntLatticeBasis()->withDual();
          NTL::mat_ZZ U;
          U.SetDims(lattmp->getBasis().NumRows(), lattmp->getBasis().NumCols());

          switch(precision){
            case DOUBLE:
              LLL_FP(lattmp->getBasis(), U, fact, 0, 0);
              break;
            case QUADRUPLE:
              LLL_QP(lattmp->getBasis(), U, fact, 0, 0);
              break;
            case EXPONENT:
              LLL_XD(lattmp->getBasis(), U, fact, 0, 0);
              break;
            case ARBITRARY:
              LLL_RR(lattmp->getBasis(), U, fact, 0, 0);
              break;
            case EXACT:
              break;
            default:
              MyExit(1, "LLL PrecisionType:   NO SUCH CASE");
          }

          //if (withDual)
          //  lattmp->getDualBasis() = transpose(inv(U)) * lattmp->getDualBasis();

          red.getIntLatticeBasis()->copyBasis(*lattmp, dim);

        }
        if (dim>0) delete lattmp;

      }
    };
  /// \endcond


  //===========================================================================

  // Initialization of non-const static members
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::PreRedDieterSV = false;
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::PreRedLLLSV = false;
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::PreRedLLLRM = false;
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::PreRedBKZ = true;
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      std::int64_t Reducer<Int, BasInt, Dbl, RedDbl>::maxNodesBB = 1000000000;


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Reducer<Int, BasInt, Dbl, RedDbl>::Reducer (
          IntLatticeBasis<Int, BasInt, Dbl, RedDbl> & lat)
    {
      // Squared length of vectors must not overflow max(double)
      m_lat = &lat;
      const int dim1 = m_lat->getDim ();
      int dim2 = dim1;
      if (dim2 <= 2)
        dim2++;

      m_c0.resize (dim1, dim1);
      m_c2.resize (dim1, dim1);
      m_cho2.resize (dim2, dim2);
      m_gramVD.resize (dim2, dim2);
      // Indices in Cholesky go as high as m_lat->getMaxDim() + 2
      m_IC = new int[2 + dim1];

      m_nv.resize (dim1);
      m_bv.resize (dim1);
      m_bw.resize (dim1);
      m_n2.resize (dim1);
      m_zLR.resize (dim1);
      m_zLI.resize (dim1);
      m_zShort.resize (dim1);
      m_dc2.resize (dim1);
      m_BoundL2.resize (dim1);

      m_lMin = std::numeric_limits <double>::max ();
      m_lMin2 = m_lMin;
      for (int i = 0; i < dim1; i++) {
        m_zLI[i] = -1;
        m_zShort[i] = -1;
        m_BoundL2[i] = -1;
        m_IC[i] = -1;
      }
      m_countNodes = 0;
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Reducer<Int, BasInt, Dbl, RedDbl>::Reducer (const Reducer<Int, BasInt, 
        Dbl, RedDbl> & red)
    {
      copy (red);
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Reducer<Int, BasInt, Dbl, RedDbl> &
      Reducer<Int, BasInt, Dbl, RedDbl>::operator= (const Reducer<Int, BasInt, 
        Dbl, RedDbl> & red)
    {
      if (this != &red)
        copy (red);
      return *this;
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::copy (const Reducer<Int, BasInt, 
        Dbl, RedDbl> & red)
    {
      m_lat = red.m_lat;
      m_c0 = red.m_c0;
      m_c2 = red.m_c2;
      m_dc2 = red.m_dc2;
      m_nv = red.m_nv;
      m_bv = red.m_bv;
      m_bw = red.m_bw;
      m_n2 = red.m_n2;
      m_zLR = red.m_zLR;
      m_zLI = red.m_zLI;
      m_zShort = red.m_zShort;
      m_cho2 = red.m_cho2;
      m_gramVD = red.m_gramVD;
      m_lMin = red.m_lMin;
      m_lMin2 = red.m_lMin2;
      m_BoundL2 = red.m_BoundL2;
      if (m_IC != 0) delete[] m_IC;
      m_IC = new int[3 + m_lat->getDim ()];
      for (int i = 0; i < 2 + m_lat->getDim (); i++)
        m_IC[i] = red.m_IC[i];
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Reducer<Int, BasInt, Dbl, RedDbl>::~Reducer ()
    {
      m_c0.clear();
      m_c2.clear();
      m_cho2.clear();
      m_gramVD.clear();
      m_nv.clear();
      m_bv.clear();
      m_bw.clear();
      m_n2.clear();
      m_zLR.clear();
      m_zLI.clear();
      m_zShort.clear();
      m_dc2.clear();
      m_BoundL2.clear();
      delete[] m_IC;
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::setBoundL2 (const DblVec & V,
          int dim1, int dim2)
    {
      m_BoundL2.resize(dim2);
      for (int i = dim1; i < dim2; i++)
        m_BoundL2[i] = V[i];
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::permuteGramVD (int i, int j, int n)
      /*
       * Permute la i-ieme et la j-ieme ligne, et la i-ieme et j-ieme colonne
       * de la matrice des produits scalaires.   Pour redLLL.
       */
    {
      int k;
      for (k = 0; k < n; k++)
      {
        std::swap(m_gramVD(i,k), m_gramVD(j,k));
      }
      for (k = 0; k < n; k++)
      {
        std::swap(m_gramVD(k,i), m_gramVD(k,j));
      }
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::calculCholeski2LLL (int n, int j)
      /*
       * Recalcule les n premieres entrees de la colonne j de la matrice de
       * Choleski d'ordre 2.   Pour redLLL.
       */
    {
      m_cho2(0,j) = m_gramVD(0,j);
      for (int i = 1; i <= n; i++) {
        m_cho2(i,j) = m_gramVD(i,j);
        for (int k = 0; k < i; k++) {
          m_cho2(i,j) -= (m_cho2(k,j) / m_cho2(k,k)) * m_cho2(k,i);
        }
      }
    }


  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      inline void Reducer<Int, BasInt, Dbl, RedDbl>::calculCholeski2Ele (int i, int j)
      // Recalcule l'entree (i,j) de la matrice de Choleski d'ordre 2
    {
      m_cho2(i,j) = m_gramVD(i,j);
      for (int k = 0; k < i; k++) {
        m_cho2(i,j) -= m_cho2(k,i) * (m_cho2(k,j) / m_cho2(k,k));
      }
    }

  //=========================================================================

  void negativeCholeski ();

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      inline void Reducer<Int, BasInt, Dbl, RedDbl>::calculGramVD ()
      /*
       * Retourne dans m_gramVD la matrice des produits scalaires m_lat->V[i]*m_lat->V[j].
       * Rem.: m_lat->V.vecNorm ne contient que les m_lat->V[i]*m_lat->V[i].
       * Reviens à faire m_lat->V * transpose(m_lat->V).
       * Utilise pour redLLL.
       * */
    {
      const int dim = m_lat->getDim ();

      for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
          NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), i);
          NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), j);
          ProdScal<Int> (row1, row2, dim, m_gramVD(i,j));
          m_gramVD(j,i) = m_gramVD(i,j);
        }
      }
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      inline void Reducer<Int, BasInt, Dbl, RedDbl>::miseAJourGramVD (int j)
      /*
       * Recalcule la j-ieme ligne et la j-ieme colonne de la matrice des
       * produits scalaires.  Pour redLLL.
       */
    {
      const int dim = m_lat->getDim ();
      for (int i = 0; i < dim; i++) {
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), i);
        NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), j);
        ProdScal<Int> (row1, row2, dim, m_gramVD(i,j));
        m_gramVD(j,i) = m_gramVD(i,j);
      }
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::calculCholeski (RedDblVec & DC2,
          RedDblMat & C0)
      /*
       * Returns in C0 the elements of the upper triangular matrix of the
       * Choleski decomposition that are above the diagonal. Returns in DC2 the
       * squared elements of the diagonal. These elements are rescaled by EchVV
       * when SISquares= true.
       */
    {
      const int dim = m_lat->getDim ();

      int k, j, i;
      RedDbl m2;
      // C2(i,j) = C0(i,j) * C2(i,i) if i != j.
      // C2(i,i) = DC2[i].
      NTL::conv (m2, m_lat->getModulo ());
      m2 = m2*m2;
      int d = dim;
      //if(m_lat->withDual())
      //  d = dim / 2; // If we use the Dual, we compute Choleski
      // with the Dual

      // Compute the d first lines of C0 with the primal Basis.
      for (i = 0; i < d; i++) {
        m_lat->updateScalL2Norm (i);
        for (j = i; j < dim; j++) {
          if (j == i)
            NTL::conv (m_c2(i,i), m_lat->getVecNorm (i));
          else {
            NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), i);
            NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), j);
            ProdScal<Int> (row1, row2, dim, m_c2(i,j));
          }
          for (k = 0; k < i; k++)
            m_c2(i,j) -= C0(k,i) * m_c2(k,j);
          if (i == j) {
            DC2[i] = m_c2(i,i);
            if (DC2[i] < 0.0) {
              negativeCholeski ();
              return false;
            }
          } else
            C0(i,j) = m_c2(i,j) / DC2[i];
        }
      }

      // Compute the d last lines of C0 with the dual Basis.
      /* This operation with the dual is needed in case of high dimension
       * and large number (30 bits). The choleski decomposition
       * convert numbers in double which is not sufficient in that case.
       * You need to use RR of the NTL library for this calculation.
       */
      //if(m_lat->withDual()){
      //  for (i = dim-1; i >= d; i--)
      //  {
      //    m_lat->updateDualScalL2Norm (i);
      //    for (j = i; j >= 0; j--) {
      //      if (j == i)
      //        NTL::conv (m_c2(i,i), m_lat->getDualVecNorm (i));
      //      else {
      //        NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), i);
      //        NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), j);
      //        ProdScal<Int> (row1, row2, dim, m_c2(i,j));
      //      }
      //      for (k = i + 1; k < dim; k++){
      //        m_c2(i,j) -= C0(k,i) * m_c2(k,j);
      //      }
      //      if (i != j)
      //        C0(i,j) = m_c2(i,j) / m_c2(i,i);
      //    }
      //    DC2[i] = m2 / m_c2(i,i);
      //    if (DC2[i] < 0.0) {
      //      negativeCholeski();
      //      return false;
      //    }
      //    for (j = i + 1; j < dim; j++) {
      //      C0(i,j) = -C0(j,i);
      //      for (k = i + 1; k < j; k++) {
      //        C0(i,j) -= C0(k,i) * C0(k,j);
      //      }
      //    }
      //  }
      //}
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::trace (char *mess)
    {
      std::cout << std::endl << "================================= " << mess << std::endl;
      //std::cout << "dim = " << m_lat->getDim () << std::endl;
      m_lat->setNegativeNorm();
      //m_lat->setDualNegativeNorm();
      m_lat->updateVecNorm ();
      //m_lat->updateDualVecNorm();
      m_lat->sort(0);
      m_lat->write();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::pairwiseRedPrimal (int i, int d, bool xx[])
    {
      const int dim = m_lat->getDim ();
      ++m_countDieter;
      m_lat->updateScalL2Norm (i);
      bool modifFlag;

      for (int j = d; j < dim; j++) {
        if (i == j)
          continue;
        modifFlag = false;
        {
          NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), i);
          NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), j);
          ProdScal<Int> (row1, row2, dim, m_ns);
        }
        DivideRound <Dbl> (m_ns, m_lat->getVecNorm (i),
            m_ns);
        if (m_ns == 0)
          continue;
        NTL::conv (m_bs, m_ns);
        if (m_ns < 1000 && m_ns > -1000) {
          m_lat->updateScalL2Norm (j);
          {
            NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
            NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
            ModifVect (row1, row2, -m_bs, dim);
          }

          // Verify that m_lat->getBasis()[j] is really shorter
          {
            NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
            ProdScal<Int> (row1, row1, dim, m_ns);
          }
          if (m_ns >= m_lat->getVecNorm (j)) {
            NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
            NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
            ModifVect (row1, row2, m_bs, dim);
          } else {
            modifFlag = true;
            m_lat->setVecNorm (m_ns, j);
          }
        } else {
          NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
          NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
          ModifVect (row1, row2, -m_bs, dim);
          m_lat->setNegativeNorm (j);
          modifFlag = true;
        }

        if (modifFlag) {
          m_countDieter = 0;
          ++m_cpt;
          if(m_lat->withDual()){
            NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), i);
            NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), j);
            ModifVect (row1, row2, m_bs, dim);
            m_lat->setDualNegativeNorm(i);

          }
          if(xx) {
            xx[i] = false;
            xx[j] = false;
          }
        }
      }
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::pairwiseRedDual (int i, bool xx[])
    {
      int j;
      const int dim = m_lat->getDim ();

      ++m_countDieter;
      m_lat->updateDualScalL2Norm (i);
      NTL::matrix_row<BasIntMat> row9(m_lat->getBasis(), i);
      m_bv = row9;
      for (j = 0; j < dim; j++) {
        if (i != j) {
          NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), i);
          NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), j);
          ProdScal<Int> (row1, row2, dim, m_ns);
          // ProdScal<Int> (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
          //           dim, m_ns);
          DivideRound <Dbl> (m_ns, m_lat->getDualVecNorm (i),
              m_nv[j]);
          if (m_nv[j] != 0) {
            NTL::conv (m_bs, m_nv[j]);
            NTL::matrix_row<BasIntMat> row7(m_lat->getBasis(), j);
            ModifVect (m_bv, row7, m_bs, dim);
          }
        }
      }

      m_lat->updateScalL2Norm (i);
      ProdScal<Int> (m_bv, m_bv, dim, m_ns);
      if (m_ns < m_lat->getVecNorm (i)) {
        ++m_cpt;
        m_countDieter = 0;
        NTL::matrix_row<BasIntMat> row6(m_lat->getBasis(), i);
        for (j = 0; j < dim; j++)
          row6(j) = m_bv[j];
        m_lat->setNegativeNorm (i);
        if (xx) xx[i] = false;
        m_lat->setVecNorm (m_ns, i);
        for (j = 0; j < dim; j++) {
          if (i != j && m_nv[j] != 0) {
            NTL::conv (m_bs, -m_nv[j]);
            NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), j);
            NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), i);
            ModifVect (row1, row2, m_bs, dim);
            //  ModifVect (m_lat->getDualBasis ()[j], m_lat->getDualBasis ()[i],
            //            m_bs, dim);
            m_lat->setDualNegativeNorm (j);
            if (xx) xx[j] = false;
          }
        }
      }
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redDieter(int d, bool xx[])
    {
      std::int64_t BoundCount;
      const int dim = m_lat->getDim ();
      bool withDual = m_lat->withDual();

      m_lat->updateScalL2Norm (d, dim);
      m_lat->sort (d);
      int i = dim-1;
      m_cpt = 0;
      m_countDieter = 0;
      BoundCount = 2 * dim - d;
      do {
        pairwiseRedPrimal (i, d, xx);
        if (i > d && withDual)
          pairwiseRedDual (i, xx);
        if (i < 1)
          i = dim-1;
        else
          --i;
      } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redDieterRandomized (int d,
          int seed)
    {
      std::int64_t BoundCount;
      const int dim = m_lat->getDim ();
      bool withDual = m_lat->withDual();

      m_lat->updateScalL2Norm (d, dim);
      //m_lat->getDualBasis ().updateScalL2Norm (d, dim);
      m_lat->sort (d);
      int i = dim-1;
      m_cpt = 0;
      m_countDieter = 0;
      BoundCount = 2 * dim - d;
      srand(seed);
      do {
        pairwiseRedPrimal (rand() % dim, d);
        if (i > d && withDual)
          pairwiseRedDual (i);
        if (i < 1)
          i = dim-1;
        else
          --i;
      } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
    }

  //=========================================================================

  /**
   * Reduce the Choleski matrix with adding a multiple of the i-th vector
   * to the j-th vector. It updates the Gram Schmidt matrix
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::reductionFaible (int i, int j)
      /*
       * Reduit la matrice de Choleski (d'ordre 2) en ajoutant un multiple du
       * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
       * consequence et remet a jour la matrice des produits scalaires.
       * Utilise par redLLL.
       */
    {
      RedDbl cte;
      std::int64_t cteLI;
      //bool withDual = m_lat->withDual();
      cte = m_cho2(i,j) / m_cho2(i,i);

      const int dim = m_lat->getDim ();

      if (abs (cte) < std::numeric_limits<double>::max()) {
        // On peut representer cte en LONGINT.
        if (abs (cte) > 0.5) {
          NTL::conv (cteLI, Round (cte));
          NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
          NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
          ModifVect (row1, row2, -cteLI, dim);
          //if(withDual){
          //  NTL::matrix_row<BasIntMat> row3(m_lat->getDualBasis(), i);
          //  NTL::matrix_row<BasIntMat> row4(m_lat->getDualBasis(), j);
          //  ModifVect (row3, row4, cteLI, dim);
          //}

        } else
          return;

      } else {
        // On represente cte en double.
        if (abs (cte) < std::numeric_limits<long double>::max())
          cte = Round (cte);
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
        NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
        ModifVect (row1, row2, -cte, dim);
        //if(withDual){
        //  NTL::matrix_row<BasIntMat> row3(m_lat->getDualBasis(), i);
        //  NTL::matrix_row<BasIntMat> row4(m_lat->getDualBasis(), j);
        //  ModifVect (row3, row4, cte, dim);
        //}

      }
      m_lat->setNegativeNorm (j);
      m_lat->updateVecNorm(j);
      //if(withDual){
      //  m_lat->setDualNegativeNorm (i);
      //  m_lat->updateDualVecNorm(i);
      //}


      miseAJourGramVD (j);
      calculCholeski2LLL (i, j);
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redLLLNTLExact(double fact) {
      spec.redLLLNTLExact(*this, fact);
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redBKZ(
        double fact, std::int64_t blocksize, PrecisionType precision, int dim)
    {
      spec.redBKZ(*this, fact, blocksize, precision, dim);
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redLLLNTL(
          double fact, PrecisionType precision, int dim) {
      spec.redLLLNTL(*this, fact, precision, dim);
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::redLLL (
          double fact, std::int64_t maxcpt, int Max)
      /*
       * Effectue la pre-reduction de B au sens de Lenstra-Lenstra-Lovasz. N'utilise
       * pas les vecteurs m_lat->getBasis().vecNorm, Wm_lat->getDualBasis().
       */
    {
      //bool withDual = m_lat->withDual();
      const int REDBAS_e = 40;
      int i, j, k, h;
      RedDbl Cho0ij;
      RedDbl limite;
      std::int64_t cpt;

      const int dim = m_lat->getDim ();
      if (Max == 0) Max = dim;

      cpt = 0;
      calculGramVD ();
      limite = 1.0;
      for (k = 1; k <= REDBAS_e; k++)
        limite *= 2.0;
      limite *= dim;
      m_cho2(0,0) = m_gramVD(0,0);
      m_cho2(0,1) = m_gramVD(0,1);
      m_IC[0] = 1;
      m_cho2(1,1) = m_gramVD(1,1) - m_cho2(0,1) * (m_cho2(0,1) / m_cho2(0,0));
      m_IC[1] = 1;
      for (i = 2; i < dim; i++)
        m_IC[i] = -1;
      h = 0;

      while (h < Max-1 && cpt < maxcpt) {
        if (m_gramVD(h + 1,h + 1) > limite) {
          for (i = h; i >= 0; i--)
            reductionFaible (i, h + 1);
        } else
          reductionFaible (h, h + 1);

        calculCholeski2Ele (h + 1, h + 1);
        if (m_IC[h + 1] == -1)
          m_IC[h + 1] = h + 1;

        if (m_cho2(h + 1,h + 1)/m_cho2(h,h) + (m_cho2(h,h + 1)/m_cho2(h,h))
            * (m_cho2(h,h + 1) / m_cho2(h,h)) < fact) {
          ++cpt;

          m_lat->permute (h, h + 1);
          permuteGramVD (h, h + 1, dim);
          m_cho2(h,h) = m_gramVD(h,h);
          for (i = 0; i < h; i++) {
            std::swap(m_cho2(i,h), m_cho2(i,h + 1));
            m_cho2(h,h) -= m_cho2(i,h) * (m_cho2(i,h) / m_cho2(i,i));
          }
          if (h == 0) {
            Cho0ij = m_cho2(0,1) / m_cho2(0,0);
            if (abs (Cho0ij) > 0.5) {
              m_IC[0] = 1;
              m_IC[1] = -1;
              h = 0;
            } else {
              m_cho2(1,1) = m_gramVD(1,1) -
                m_cho2(0,1) * m_cho2(0,1) / m_cho2(0,0);
              calculCholeski2LLL (2, 2);
              m_IC[0] = 2;
              m_IC[1] = 2;
              m_IC[2] = 2;
              h = 1;
            }
          } else {
            m_IC[h] = h + 1;
            m_IC[h + 1] = -1;
            --h;
          }

        } else {
          for (i = 0; i <= h + 2; i++) {
            if (h + 2 > m_IC[i]) {
              if (h + 2 < dim)
                calculCholeski2Ele (i, h + 2);
              m_IC[i] = h + 2;
            }
          }
          ++h;
        }
      }

      if(cpt == maxcpt){
        std::cout << "***** in redLLL cpt > maxcpt = " << maxcpt << std::endl;
      }

      for (j = 2; j < Max; j++) {
        for (i = j - 2; i >= 0; i--)
          reductionFaible (i, j);
      }
      m_lat->setNegativeNorm ();
      //if(withDual){
      //  m_lat->setDualNegativeNorm ();
      //}
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void Reducer<Int, BasInt, Dbl, RedDbl>::transformStage3 (
          std::vector<std::int64_t> & z, int & k)
    {
      int j, i;
      std::int64_t q;
      const int dim = m_lat->getDim ();
      //bool withDual = m_lat->withDual();

      j = dim-1;
      while (z[j] == 0)
        --j;
      while (abs (z[j]) > 1) {
        i = j - 1;
        while (z[i] == 0)
          --i;
        // On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
        while (z[j]) {
          // Troncature du quotient vers 0
          q = z[i] / z[j];
          if (q) {
            // On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
            z[i] -= q * z[j];
            NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), i);
            NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), j);
            //    ModifVect (m_lat->getBasis ()[j], m_lat->getBasis ()[i],
            //            q, dim);
            ModifVect (row1, row2, q, dim);
            m_lat->setNegativeNorm (j);

            //if(withDual){
            //  NTL::matrix_row<BasIntMat> row3(m_lat->getDualBasis(), i);
            //  NTL::matrix_row<BasIntMat> row4(m_lat->getDualBasis(), j);
            //  //    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
            //  //             -q, dim);
            //  ModifVect (row3, row4, -q, dim);
            //  m_lat->setDualNegativeNorm (i);
            //}
          }
          // Permutation.
          std::swap <std::int64_t>(z[i], z[j]);
          m_lat->permute (i, j);
        }
        j = i;
      }
      k = j;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::tryZ (
          int j, int i, int Stage, bool & smaller, const BasIntMat & WTemp)
      // Si m_countNodes > MaxNodes retourne false, sinon retourne true.
    {
      std::int64_t max0, min0;
      RedDbl x, dc;
      RedDbl center;
      std::int64_t zhigh, zlow, h;
      bool high;
      int k;
      RedDbl S1, S2, S3, S4, mR;

      const int dim = m_lat->getDim ();
      NTL::conv (mR, m_lat->getModulo());


      ++m_countNodes;
      if (m_countNodes > maxNodesBB) {
        std::cerr << "-------- m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
        return false;
      }


      // Calcul d'un intervalle contenant les valeurs admissibles de zj.
      center = 0.0;
      if (j < dim-1) {
        // Calcul du centre de l'intervalle.
        for (k = j + 1; k < dim; k++)
          center = center - m_c0(j,k) * m_zLR[k];

        // Distance du centre aux extremites de l'intervalle.
        dc = sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j]);

        /* Calcul de deux entiers ayant la propriete qu'un entier */
        /* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
        /* n'appartient pas a l'intervalle.  */
        x = center - dc;
        NTL::conv (min0, trunc (x));
        if (x > 0.0)
          ++min0;

        x = center + dc;
        NTL::conv (max0, trunc (x));
        if (x < 0.0)
          --max0;

        // En vue du choix initial de zj. On determine zlow et zhigh.
        if (min0 > max0){
          return true;
        }
        if (min0 == max0) {
          zlow = min0;
          zhigh = max0 + 1;
          high = false;
        } else if (center >= 0.0) {
          NTL::conv (zlow, trunc (center));
          zhigh = zlow + 1;
          NTL::conv (h, trunc (2.0 * center));
          high = h & 1;
        } else {
          NTL::conv (zhigh, trunc (center));
          zlow = zhigh - 1;
          NTL::conv (h, -trunc (2.0 * center));
          high = (h & 1) == 0;
        }

      } else {  // j = dim-1

        zlow = 0;
        high = true;
        if (Stage == 2) {
          min0 = 1;
          max0 = 1;
          zhigh = 1;
        } else {
          min0 = 2;
          zhigh = 2;
          NTL::conv (max0, trunc (sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j])));
        }
      }


      Dbl temp;
      /* On essaie maintenant chacun des z[j] dans l'intervalle, en */
      /* commencant par le centre puis en alternant d'un cote a l'autre. */
      while (zlow >= min0 || zhigh <= max0) {

        if (high) {
          m_zLI[j] = zhigh;
        } else {
          m_zLI[j] = zlow;
        }
        m_zLR[j] = m_zLI[j];

        // Calcul de m_n2[j-1].
        x = m_zLR[j] - center;

        if (j == 0) {
          RedDbl tmps_n2 = m_n2[0] + x * x * m_dc2[0];
          if (tmps_n2 < m_lMin2) {
            // On verifie si on a vraiment trouve un vecteur plus court
            NTL::matrix_row<const BasIntMat> row1(m_lat->getBasis(), dim-1);
            m_bv = row1;
            for (k = 0; k < dim-1; k++) {
              if (m_zLI[k] != 0) {
                NTL::matrix_row<const BasIntMat> row1(m_lat->getBasis(), k);
                ModifVect (m_bv, row1, m_zLI[k], dim);
              }
            }
            if (Stage == 3) {
              NTL::matrix_row<const BasIntMat> row1(m_lat->getBasis(), dim-1);
              ModifVect (m_bv, row1, m_zLR[dim-1] - 1.0, dim);
            }

            ProdScal<Int> (m_bv, m_bv, dim, S1);
            NTL::conv (S4, m_lat->getVecNorm (dim-1));
            if (S1 < S4) {
              if (Stage == 2) {
                smaller = true;
                if (!PreRedLLLRM)
                  m_zShort = m_zLI;
                else {
                  for (k = 1; k < dim; k++) {
                    NTL::matrix_row<const BasIntMat> row1(WTemp, k);
                    ProdScal<Int> (m_bv, row1, dim, S2);
                    Quotient (S2, mR, S3);
                    NTL::conv (m_zShort[k], S3);
                  }
                  m_zShort[dim-1] = 1;
                }
              } else if (Stage == 3 && !PreRedLLLRM) {
                if (GCD2vect (m_zLI, i, dim) == 1) {
                  m_zShort = m_zLI;
                  smaller = true;
                }
              } else {
                for (k = 0; k < dim; k++) {
                  NTL::matrix_row<const BasIntMat> row1(WTemp, k);
                  ProdScal<Int> (m_bv, row1, dim, S2);
                  Quotient (S2, mR, S3);
                  NTL::conv (m_zShort[k], S3);
                }
                if (GCD2vect (m_zShort, i, dim) == 1) {
                  smaller = true;
                }
              }
              if (smaller) {
                NTL::conv (temp, S1);
                m_lat->setVecNorm (temp, dim-1);
                return true;
              }
            }
          }
        } else { // j > 0
          m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
          if (m_lMin2 >= m_n2[j - 1]) {
            if (!tryZ (j - 1, i, Stage, smaller, WTemp))
              return false;
            // Des qu'on a trouve quelque chose, on sort de la recursion
            // et on retourne dans reductMinkowski.
            if (smaller)
              return true;
          }
        }
        if (high) {
          ++zhigh;
          if (zlow >= min0)
            high = false;
        } else {
          --zlow;
          if (zhigh <= max0)
            high = true;
        }
      }
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::redBB (
          int i, int d, int Stage, bool & smaller, bool xx[])
      /*
       * Tries to shorten m_lat->getBasis()[i] using branch-and-bound.
       * Used in Minkowski Reduction.
       * Stage is 2 or 3.
       * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
       * Stops and returns false if not finished after examining MaxNodes
       * nodes in the branch-and-bound tree.  When succeeds, returns true.
       * Assumes that the norm is Euclidean.
       */
    {
      bool withDual = m_lat->withDual();
      const int dim = m_lat->getDim ();
      BasIntMat VTemp (dim, dim), WTemp (dim, dim);
      bool XXTemp[dim];
      Dbl tmp;
      // trace( "AVANT redBB");
      smaller = false;

      // Approximation du carre de la std::int64_tueur de Vi.
      if (m_lat->getVecNorm(i) < 0) {
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), i);
        ProdScal<Int> (row1, row1, dim, tmp);
        //  ProdScal<Int> (m_lat->getBasis()[i], m_lat->getBasis()[i],
        //            dim, tmp);
        m_lat->setVecNorm (tmp, i);
      }
      NTL::conv (m_lMin2, m_lat->getVecNorm (i));

      if (Stage == 3 && withDual)
      {
        if (m_lat->getDualVecNorm (i) < 0) {
          NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), i);
          ProdScal<Int> (row1, row1, dim, tmp);
          //   ProdScal<Int> (m_lat->getDualBasis()[i], m_lat->getDualBasis()[i],
          //            dim, tmp);
          m_lat->setDualVecNorm (tmp, i);
        }
        Dbl m2;
        NTL::conv(m2, m_lat->getModulo());
        m2 = m2*m2;
        if (m_lMin2 * m_lat->getDualVecNorm (i) < 4*m2)
          return true; // if the angle between the basis vector i and the dual
        // basis vector i is between -Pi/3 and Pi/3
      }
      if(withDual){
        m_lat->updateDualVecNorm ();
      }
      m_lat->updateVecNorm ();
      m_lat->permute (i, dim-1);

      int k, h;

      if (PreRedLLLRM)
      {
        // On memorise la base courante.
        VTemp = m_lat->getBasis ();
        if(withDual){
          WTemp = m_lat->getDualBasis ();
        }
        for (h = 0; h < dim; h++)
          XXTemp[h] = true;
        redLLL (1.0, 1000000, dim - 1);
        m_lat->updateVecNorm ();
        if(withDual){
          m_lat->updateDualVecNorm ();
        }
      }
      if (!calculCholeski (m_dc2, m_c0))
        return false;
      m_countNodes = 0;
      m_n2[dim-1] = 0.0;
      if (!tryZ (dim-1, i, Stage, smaller, WTemp))
        return false;

      if (PreRedLLLRM)
      {
        /* On remet l'anciennne base, celle d'avant LLL, avant de considerer
           la prochaine m_lat->dimension.  */
        m_lat->getBasis () = VTemp;
        m_lat->updateVecNorm ();
        if(withDual){
          m_lat->getDualBasis () = WTemp;
          m_lat->updateDualVecNorm ();
        }
        for (h = 0; h < dim; h++)
          xx[h] = XXTemp[h];
      }
      if (smaller)
      {
        /* On a trouve un plus court vecteur.  On ameliore
           m_lat->getBasis()[k].  */
        if (Stage == 2)
          k = dim-1;
        else
          transformStage3 (m_zShort, k);
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), k);
        for (h = 0; h < dim; h++)
          row1(h) = m_bv[h];
        //  m_lat->getBasis ()[k] = m_bv;
        m_lat->setNegativeNorm (k);
        if (m_zShort[k] < 0 && withDual) {
          NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), k);
          ChangeSign (row2, dim);
        }
        /* Mise a jour des vecteurs de la base duale selon le nouveau
           m_lat->getBasis()[k] */
        for (h = 0; h < dim; h++) {
          if ((m_zShort[h] != 0) && (h != k)) {
            if(withDual){
              NTL::matrix_row<BasIntMat> row1(m_lat->getDualBasis(), h);
              NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), k);
              ModifVect (row1, row2, -m_zShort[h], dim);
              m_lat->setDualNegativeNorm (h);
            }
            if (Stage == 2) {
              if (h >= d)
                xx[h] = false;
            }
          }
        }
      } else if (Stage == 2)
        xx[dim-1] = true;

      m_lat->permute (i, dim-1);
      // trace( "APRES redBB");
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::tryZ0 (int j, bool & smaller)
      /*
       * Procedure recursive implantant le BB pour shortestVector (redBB0).
       * Si m_countNodes > MaxNodes, retourne false.
       * Sinon, retourne true: on a trouve le plus court vecteur.
       * Le carre de la std::int64_tueur du plus court vecteur a date est dans m_lMin2.
       * Si j=1 et on trouve un plus court vecteur, on retrouve les zj
       * correspondants dans m_zShort.
       */
    {
      /* Note: Pour une implantation non recursive, ces variables devraient
         etre des tableaux indices par j. */
      RedDbl dc, x, center;
      std::int64_t min0, max0;            // Bornes de l'intervalle pour les z_j.
      std::int64_t zlow, zhigh;           // Valeur courante examinee a gauche et a
      // droite du centre.
      bool high;                  // Indique si on est a d. ou g. du centre.
      int k;
      std::int64_t temp;

      const int dim = m_lat->getDim ();

      ++m_countNodes;
      if (m_countNodes > maxNodesBB) {
        std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
        return false;
      }
      /* Calcul d'un intervalle contenant les valeurs admissibles de zj. */
      /* Calcul du centre de l'intervalle.  */
      center = 0.0;
      for (k = j + 1; k < dim; ++k)
        center -= m_c0(j,k) * m_zLR[k];

      // Distance du centre aux extremites de l'intervalle.

      dc = sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j]);

      /* Calcul de deux entiers ayant la propriete qu'un entier */
      /* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
      /* n'appartient pas a l'intervalle.  */
      /* Si NOT m_foundZero on pose min egal a zero par symetrie.  */
      if (!m_foundZero)
        min0 = 0;
      else {
        x = center - dc;
        NTL::conv (min0, trunc (x));
        if (x > 0.0)
          ++min0;
      }
      x = center + dc;
      NTL::conv (max0, trunc (x));
      if (x < 0.0)
        --max0;

      // En vue du choix initial de zj. On determine zlow et zhigh.
      if (min0 > max0)
        return true;
      if (min0 == max0) {
        zlow = min0;
        zhigh = max0 + 1;
        high = false;
      } else if (center >= 0.0) {
        NTL::conv (zlow, trunc (center));
        zhigh = zlow + 1;
        NTL::conv (temp, trunc (2.0 * center));
        high = (temp & 1);
      } else {
        NTL::conv (zhigh, trunc (center));
        zlow = zhigh - 1;
        NTL::conv (temp, -trunc (2.0 * center));
        high = (temp & 1) == 0;
      }

      // On essaie maintenant chacun des z[j] dans l'intervalle, en
      // commencant par le centre puis en alternant d'un cote a l'autre.
      while (zlow >= min0 || zhigh <= max0) {
        if (high)
          m_zLI[j] = zhigh;
        else
          m_zLI[j] = zlow;
        m_zLR[j] = m_zLI[j];

        // Calcul de m_n2[j-1], qui est normalise dans le cas SISquares.
        x = m_zLR[j] - center;


        if (j == 0) {
          // Tous les zj sont choisis: on a un vecteur a tester.
          if (m_lMin2 > m_n2[0] + x * x * m_dc2[0]) {
            /* On verifie si on a vraiment trouve un vecteur */
            /* non-nul et plus court.  */
            if (!m_foundZero) {
              // Le premier vecteur trouve sera zero.
              m_foundZero = true;
            } else {
              SetZero (m_bv, dim);
              for (k = 0; k < dim; k++) {
                if (m_zLI[k] != 0) {
                  NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), k);
                  ModifVect (m_bv, row1, m_zLI[k], dim);
                }
              }
              // Le nouveau vecteur trouve est maintenant dans m_bv.
              if (m_lat->getNorm () == L2NORM) {
                ProdScal<Int> (m_bv, m_bv, dim, x);
              } else {
                CalcNorm <BasIntVec, RedDbl> (m_bv, dim, x, m_lat->getNorm ());
                x = x * x;
              }
              if (x < m_lMin2) {
                /* La condition suivante ralentit le programme; je l'ai mise donc en
                   commentaire. Il est très rare qu'elle soit effective et permette de sortir
                   prématurément, et seulement pour de grandes dimensions dans le cas L2NORM.
                   Mais on doit la tester à chaque passage ici, et la std::int64_tueur des candidats
                   plus courts diminue lentement.
                   Il se pourrait que ce soit parfois plus rapide dans le cas L1NORM,
                   dépendant de la dimension. Mais L2NORM est primordial. */
                /*
                   if (x <= m_BoundL2[dim]) {
                   return false;
                   }
                   */

                // Il est plus court!
                smaller = true;
                NTL::conv (m_lMin2, x);
                m_zShort = m_zLI;
                m_bw = m_bv;
              }
            }
          }
        } else if (m_lMin2 > m_n2[j] + x * x * m_dc2[j]) {
          // Encore de l'espoir: recursion.
          m_n2[j-1] = m_n2[j] + x * x * m_dc2[j];
          if (!tryZ0 (j-1, smaller))
            return false;
        } else
          m_n2[j-1] = m_n2[j] + x * x * m_dc2[j];

        if (high) {
          ++zhigh;
          if (zlow >= min0)
            high = false;
        } else {
          --zlow;
          if (zhigh <= max0)
            high = true;
        }
      }
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::redBB0 (NormType norm)
      /*
       * Finds shortest non-zero vector, using branch-and-bound. Stops and returns
       * false if not finished after examining MaxNodes nodes in the
       * branch-and-bound tree.  When succeeds, returns true, and the shortest
       * vector square length will be in m_lMin2.
       */
    {

      //bool withDual = m_lat->withDual();
      bool smaller = false;
      int k, h;
      const int dim = m_lat->getDim ();

      RedDbl x;
      /* We sort here to get same results as in xds98 version. Otherwise,
         Cholesky will fail more rapidly due to floating-point errors. We do
         not sort after redLLL because this greatly slows down the program in
         the case of the L2 Norm. */

      m_lat->updateScalL2Norm (0, dim);
      //if(withDual){
      //  m_lat->updateDualScalL2Norm (0, dim);
      //}
      if (m_countNodes < SHORT_LLL)
        m_lat->sort (0);

      /* Approximation de la norme du plus court vecteur. */
      if (norm == L2NORM) {
        NTL::conv (m_lMin2, m_lat->getVecNorm (0));
      } else {
        // Looking for the min of of the vecteur according to the considered norm
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), 0);
        CalcNorm <BasIntVec, Dbl> (row1, dim, m_lMin, norm);
        for (k = 1; k < dim; k++) {
          NTL::matrix_row<BasIntMat> row2(m_lat->getBasis(), k);
          CalcNorm <BasIntVec, Dbl> (row2, dim, x, norm);
          if (x < m_lMin)
            m_lMin = x;
        }
        m_lMin2 = m_lMin * m_lMin;
      }

      if (m_lMin2 <= m_BoundL2[dim-1]) {

        /* If there is in this lattice a vector that are shorter than the shortest
           of the current lattice, we don't need to find the shortest vector in this
           lattice. Use only in seek function

           S'il existe dans ce réseau un vecteur de std::int64_tueur plus courte que celle
           du meilleur réseau trouvé à date, il n'est pas nécessaire de trouver le
           plus court vecteur de ce réseau: on peut l'éliminer immédiatement. */

        return false;
      }

      if (!calculCholeski (m_dc2, m_c0))
        return false;

      /* On effectue le branch and bound.  */
      /* m_n2[j] sera la somme des termes |z*k|^2 ||v*k||^2 pour k > j.  */
      m_n2[dim-1] = 0.0;
      m_countNodes = 0;
      smaller = false;
      m_foundZero = false;
      if (!tryZ0 (dim-1, smaller))
        return false;

      if (smaller) {
        // We have found a shorter vector. Its square length is in m_lMin2.
        transformStage3 (m_zShort, k);
        NTL::matrix_row<BasIntMat> row1(m_lat->getBasis(), k);
        for (h = 0; h < dim; h++)
          row1(h) = m_bw[h];
        m_lat->setNegativeNorm (k);

        /* The new shorter vector is now in m_lat->getBasis()[k].  */
        /*  We update the vectors of the dual basis. */
        //if (withDual){
        //  if (m_zShort[k] < 0) {
        //    NTL::matrix_row<BasIntMat> row2(m_lat->getDualBasis(), k);
        //    ChangeSign (row2, dim);
        //  }
        //  for (h = 0; h < dim; h++) {
        //    if (m_zShort[h] && h != k) {
        //      NTL::matrix_row<BasIntMat> row3(m_lat->getDualBasis(), h);
        //      NTL::matrix_row<BasIntMat> row4(m_lat->getDualBasis(), k);
        //      ModifVect (row3, row4, -m_zShort[h], dim);
        //      m_lat->setDualNegativeNorm (h);
        //    }
        //  }
        //}

        /* The new candidate for a shortest vector will be in
           m_lat->getBasis()(0). */
        /* In the case of L1NORM or others, check that it is really smaller.  */
        if (norm == L2NORM)
          m_lat->permute (k, 0);
        else {
          NTL::matrix_row<BasIntMat> row5(m_lat->getBasis(), k);
          CalcNorm (row5, dim, x, norm);
          if (x < m_lMin) {
            m_lMin = x;
            m_lMin2 = m_lMin * m_lMin;
            m_lat->permute (k, 0);
          }
        }
      }
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::reductMinkowski (int d)
    {
      bool withDual = m_lat->withDual();
      const int dim = m_lat->getDim ();
      int i;
      std::int64_t totalNodes = 0;
      bool found;
      bool smaller;               // A smaller vector has been found
      bool xx[dim];               // This has a use, I just don't know which

      do {
        // The first d vectors should not be modified.
        for (i = 0; i < d; i++)
          xx[i] = true;
        for (i = d; i < dim; i++)
          xx[i] = false;


        found = false; 

        do {
          redDieter (d, xx);

          m_lat->setNegativeNorm(d);
          m_lat->updateVecNorm (d);
          if(withDual){
            m_lat->setDualNegativeNorm(d);
            m_lat->updateDualVecNorm (d);
          }
          m_lat->sort (d);
          found = false;

          for (i = 0; i < dim; i++) {
            if (!xx[i]) {
              // On essaie de reduire le i-eme vecteur.
              if (!redBB (i, d, 2, smaller, xx))
                return false;
              totalNodes += m_countNodes;
              if (smaller){
                found = true;
              }
            }
          }

        } while (found);
        // Stage 3
        if (dim > 7) {
          for (i = d; i < dim; i++) {
            if (!redBB (i, d, 3, smaller, xx))
              return false;
            totalNodes += m_countNodes;
            if (smaller)
              found = true;
          }
        }
      } while (found);

      /*
         if (totalNodes > MINK_LLL) {
         PreRedLLLRM = true;
         }
         */

      m_lat->setNegativeNorm();
      m_lat->updateScalL2Norm (0, dim);
      if(withDual){
        m_lat->setDualNegativeNorm();
        m_lat->updateDualScalL2Norm (0, dim);
      }
      m_lMin2 = m_lat->getVecNorm(0);
      return true;
    }

  //=========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool Reducer<Int, BasInt, Dbl, RedDbl>::shortestVector (NormType norm)
      // Square length of shortest vector can be recovered in m_lMin2
    {

      if (norm != L2NORM) {
        m_lat->setNegativeNorm ();
        //if(m_lat->withDual())
        //  m_lat->setDualNegativeNorm ();
      }

      /* Find the shortest vector for the selected norm.  */
      /* The L2 norm is used for the Choleski decomposition and BB bounds. */
      bool ok;
      if (norm == L1NORM || norm == L2NORM || norm == ZAREMBANORM) {
        ok = redBB0 (norm);
      } else {
        ok = false;
        std::cerr << "RedLattice::shortestVector:   wrong norm";
        exit (3);
      }

      m_lat->updateVecNorm();
      m_lat->sort(0);
      //if(m_lat->withDual())
      //  m_lat->updateDualVecNorm();

      return ok;
    }

  //============================================================================

  extern template class Reducer<std::int64_t, std::int64_t, double, double>;
  extern template class Reducer<NTL::ZZ, NTL::ZZ, double, double>;
  extern template class Reducer<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>;

}     // namespace LatticeTester

#endif // REDUCER_H
