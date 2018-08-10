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

#ifndef LATTICEANALYSIS_H
#define LATTICEANALYSIS_H

#include <string>
#include <list>
#include <dirent.h>
#include <fnmatch.h>
#include <typeinfo>
#include <cstdint>

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/Reducer.h"
#include "latticetester/Writer.h"
#include "latticetester/WriterRes.h"
#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/ParamReader.h"

namespace LatticeTester {

  // Declaration is needed for specializer structure definition
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class LatticeAnalysis;

  /**
   * \cond
   * This structure contains methods that also are in LatticeAnalysis.
   * It can then be specialized to specialize these methods instead of the whole
   * class.
   * This structure is invisible in the documentation to keep it light (it is
   * not usefull to know this implementation detail when using the API).
   * */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      struct specLatticeAnalysis {
        Normalizer<RedDbl>* initNormalizer(LatticeAnalysis<Int, BasInt, Dbl, RedDbl>& latanal,
            NormaType norma, int alpha);
      };
  /// \endcond

  /**
   * Objects of this class can perform various tests on lattices. There are two
   * intended usages of this class.
   * - The first one is the one used in the 
   * **lattest** executable, it is to simply pass a directory name (resp. a file
   * name) to the program to perform the tests specified in the `.dat` files in
   * the directory (resp. the file itself).
   * - The other possibility is to instanciate all the needed fields of the
   *   class independently and to call the doTest() method. This approach is
   *   more flexible and allows easy implementation of small scripts.
   *
   * The tests consist on the computation of one of the figures of merit
   * enumerated in CriterionType. If the user intend to simply reduce a lattice
   * basis, he should look into the Reducer module since it offers much more
   * flexibility. This class curently only implements the spectral and the Beyer
   * test. \todo implement \f$\mathcal{P}_\alpha\f$ and implement the other one
   * liste in the enum.
   *
   * To use this class it is imperative to instanciate a Reducer object with the
   * basis of the lattice to analyse. After that, you have to pass, either by
   * the constructor or by the various methods available, all the parameters 
   * necessary to do a test (these are listed under the empty constructor).
   * Finally, you simply have to call the doTest() method and the `m_merit`
   * field will contain the figure of merit that was asked for.
   */

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class LatticeAnalysis {
        private:
          typedef NTL::vector<Int> IntVec;
          typedef NTL::matrix<Int> IntMat;
          typedef NTL::vector<BasInt> BasIntVec;
          typedef NTL::matrix<BasInt> BasIntMat;
          typedef NTL::vector<Dbl> DblVec;
          typedef NTL::vector<RedDbl> RedDblVec;
          typedef NTL::matrix<RedDbl> RedDblMat;
        public:

          /**
           * Base constructor for the case where the data will come from files.
           * In that case, the fields of the class will be instanciated when
           * read. This constructor can also be used to create a default object
           * for the class that will be populated one field at a time by the user.
           *
           * If the user choose to fill by hand the fields that are needed, he
           * needs to fill the following fields:
           * - `m_reducer` with a valid reducer containing a basis.
           * - `m_criterion` with the criterion for the figure of merit he wants
           * - `m_preRed` with one of the pre-reduction strategies.
           * - `m_norm` with one of the norms for which the selected test has
           *   been implemented.
           * - `m_normalizer` and `m_normalizerType` to describe the normalizer
           *   to use. To do that, one should use the
           *   initNormalizer() method that does dynamic memory management.
           * */
          LatticeAnalysis ();

          /**
           * Constructor containing all the fields needed to perform a test.
           * - `reducer` is the reducer that will be used to perform the tests.
           *   It will be stored as a pointer, so the object can be retrieved
           *   even after the tests have been done. Also note that the reducer
           *   contains the basis and that all interactions with the lattice
           *   basis will be done through this object.
           * - `criterion` is the criterion for which the tests will be applied.
           *   For the various options, see the CriterionType enum.
           * - `preRed` is the chose of pre-reduction algorithm. See
           *   PreReductionType for a list of option.
           * - `norm` is the norm that will be used to compute the different
           *   vector lengths in the algorithms. For a list of options, look at
           *   NormType.
           * - `MaxNodesBB` is the maximum number of nodes that the program will
           *   use when doing the branch-and-bound algorithm to find the shortest
           *   vector.
           * - `normaType` is the normalization type to be used. See NormaType
           *   for a list of options.
           * - `alpha` is the \f$\alpha\f$ parameter for the \f$\mathcal{P}_\alpha\f$
           *   test. As of now, this test is not implemented.
           * */
          LatticeAnalysis (Reducer<Int, BasInt, Dbl, RedDbl> & reducer,
              CriterionType criterion,
              PreReductionType preRed = BKZ,
              NormType norm = L2NORM,
              std::int64_t maxNodesBB = 10000000,
              NormaType normaType = NORMA_GENERIC, 
              int alpha = 0);

          /**
           * Destructor. This will deallocate memory to `m_normalizer` if it is
           * not null.
           */
          ~LatticeAnalysis ();

          /**
           * This will perform the tests specified in `m_criterion` on the basis
           * contained in `m_reducer`. Before searching for the shortest vector
           * or performing Minkowski reduction, it will use the reduction in
           * `m_preRed`. Depending on the criterion, this will calculate a
           * figure of merit normalized with `m_normalizer` or simply a figure
           * of merit that is not standardized and store it in `m_merit`.
           *
           * This method will return `true` upon completion if the test has been
           * successful, and `false` if there has been some kind of problem in
           * the execution.
           * 
           * To access the results of the test, it is mandatory to call
           * getMerit() const before calling this method again.
           *
           * The parameters of this method are used by the pre-reduction algorithms
           * if BKZ or LLL pre-reduction has been choosen. Most of the time,
           * there will be no need to touch them. If you want to modify them,
           * please look at the redBKZ and redLLLNTL methods in Reducer.
           */
          bool doTest (double fact = 0.999999, PrecisionType precision = QUADRUPLE,
              int blocksize = 10);

          /**
           * This prints the results of the tests in a standard format on
           * standard output. This will print the criterion and the
           * pre-reduction used, the length of the shortest vector, the figure
           * of merit obtained and the normalizer used.
           * */
          void printTestResults ();

          /**
           * Reads the parameters of the test in input text file `datafile`; then 
           * do the test. The format of the data file must be as presented in 
           * \ref detailed_usage.
           *
           * The data file must always have the extension ".dat", 
           * but must be given as argument here *without extension*. For example, 
           * if the data file is named `myLattice.dat`, then the method must be 
           * called as `doTest("myLattice")`.
           *
           * This method returns 0 if the test completed successfully and
           * returns a negative integer if there was an error.
           */
          int doTestFromInputFile (const char *datafile);

          /**
           * Applies the method `doTestFromInputFile` to all the files with
           * extension `".dat"` in directory named `dirname`. Returns 0 if all
           * the tests completed successfully; returns a non-zero integer if
           * there was an error. Even if there is an error with one of the files,
           * this method will try to do the tests on all the files.
           */
          int doTestFromDirectory (const char *dirname);


          /**
           * Initialize `m_normalizer` to a pointer on a normalizer object of
           * type `norma`. This will set `m_normalizerType` to `norma`,
           * deallocate memory of `m_normalizer` (if it is not null) and make
           * `m_normalizer` point to a new dynamically allocated normalizer
           * of type `norma`.
           */
          void initNormalizer (NormaType norma, int alpha = 0);

          /**
           * This will set `m_reducer` to point to `red`.
           */
          void setReducer (Reducer<Int, BasInt, Dbl, RedDbl> & red)
          {
            m_reducer = &red;
          }

          /**
           * Sets `m_criterion` to `criterion`. This changes the critertion
           * for the test.
           * */
          void setCriterion (CriterionType criterion) 
          {
            m_criterion = criterion;
          }

          /**
           * Sets `m_preRed` to `preRed`. This changes the pre-reduction for the
           * test.
           * */
          void setPreReduction (PreReductionType preRed) { m_preRed = preRed; }

          /**
           * Sets `m_norm` to `norm`. This changes the norm used in the
           * reduction for the test.
           * */
          void setNorm (NormType norm) { m_norm = norm; }


          /**
           * Sets `m_maxNodesBB` to `maxNodesBB`. This changes the maximum
           * number of nodes there can be in the branch and bound tree when doing
           * the shortest vector problem.
           * */
          void setMaxNodesBB (std::int64_t maxNodesBB)
          {
            Reducer<Int, BasInt, Dbl, RedDbl>::maxNodesBB = m_maxNodesBB 
              = maxNodesBB;
          }

          /**
           * Returns a pointer to the reducer used for the analysis.
           * */
          Reducer<Int, BasInt, Dbl, RedDbl>* getReducer() 
            {
              return m_reducer;
            }

          /**
           * Returns the figure of merit calculated by doTest() as a `double`
           */
          double getMerit () const { return m_merit; }

        protected:
          /* Making these fields protected instead of private has two effects
           * that are both intended.
           *
           * First, it is possible to expand this class or overload some of its
           * methods by implementing a subclass and have access to those fields
           * directly while in this class.
           *
           * Second, it makes those variables visible in the documentation, but
           * prohibits a direct interraction with the variables. Since it is
           * possible that a user will interact with those fields (via the
           * setters and getters) it is important that they are documented.
           * */

          /**
           * Sets `m_normalizerType` to `normalizerType`. There should be no
           * need to call this method directly because initNormalizer() calls it.
           * */
          void setNormalizerType (NormaType normalizerType)
          {
            m_normalizerType = normalizerType;
          }

          /**
           * Deallocates the memory allocated to `m_normalizer` and makes it
           * point to `normalizer`. Since this method calls a delete on
           * `m_normalizer`, it is advisable to never call it directly and to
           * let initNormalizer() manage the memory.
           * */
          void setNormalizer (Normalizer<RedDbl>* normalizer) 
          {
            // This feels like a bad idea, but this is definitely the only way
            // to really make the object pointed to by m_normalizer is freed.
            // The implementation of nomalizers suggests that they should always
            // be dynamically allocated.
            if (m_normalizer) delete m_normalizer;
            m_normalizer = normalizer;
          }

          /**
           * Returns a pointer to the normalizer used in this class. Using this
           * is unsafe because methods of this class might deallocate the memory
           * emplacement this pointer points to and this will point to garbage.
           * */
          Normalizer<RedDbl>* getNormalizer() { return m_normalizer;}

          /**
           * Pointer to the reducer class used to perform pre-reductions and 
           * Branch-and-Bound. The object this points to also contains the basis
           * on which the tests are applied.
           */
          Reducer<Int, BasInt, Dbl, RedDbl>* m_reducer;

          /**
           * Type of pre-reduction that will be used before the Branch-and-Bound.
           */
          PreReductionType m_preRed;

          /**
           * Type of test applied on the basis.
           */
          CriterionType m_criterion;

          /**
           * Pointer to the Normalizer class used to normalize the results. To
           * interract with this data field, you should only use the
           * initNormalizer() method.
           */
          Normalizer<RedDbl>* m_normalizer;

          /**
           * Type of normalization chosen for the test.
           */
          NormaType m_normalizerType;

          /**
           * Norm used to evaluate the vector lengths in the reduction and
           * pre-reduction.
           */
          NormType m_norm;

          /**
           * Contains the results of the test. You can get this after doTest() has
           * been called use the getMerit() const method.
           */
          double m_merit;

          /**
           * Contains the maximum number of nodes visited in the Branch-and-Bound 
           * procedure.
           */
          std::int64_t m_maxNodesBB;

        private:

          /**
           * This creates and returns a Writer object
           * Returns a `Writer` created from the input file `infile` and the 
           * given `OutputType`.
           */
          Writer<Int>* createWriter (const char *infile,
              OutputType ot);

          /**
           * This contains methods that need to be specialized
           * */
          struct specLatticeAnalysis<Int, BasInt, Dbl, RedDbl> spec;

      }; // End class LatticeAnalysis

  //===========================================================================

  // specLatticeAnalysis specialization

  ///\cond
  // LLXX case specialization
  template<typename Dbl, typename RedDbl>
      struct specLatticeAnalysis<std::int64_t, std::int64_t, Dbl, RedDbl> {
      Normalizer<RedDbl>* initNormalizer(
          LatticeAnalysis<std::int64_t, std::int64_t, Dbl, RedDbl>& latanal,
          NormaType norma, int alpha) {
        RedDbl logDensity;
        // We have to cast to NTL::matrix<ZZ> because it does not compute the 
        // determinant of a general matrix by default.
        // This may fix some bugs at places where the old implementation did not
        // because boost did not compute the good determinant for big matrixes
        // somehow.
        NTL::mat_ZZ temp;
        int dim = latanal.getReducer()->getIntLatticeBasis()->getDim();
        temp.SetDims(dim, dim);
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++){
            temp[i][j] = latanal.getReducer()->getIntLatticeBasis()->
              getBasis()(i,j);
          }
        }
        logDensity = -log(abs(NTL::determinant(temp)));
        switch (norma) {
          case BESTLAT:
            return new NormaBestLat<RedDbl> (logDensity, dim);
            break;
          case LAMINATED:
            return new NormaLaminated<RedDbl> (logDensity, dim);
            break;
          case ROGERS:
            return new NormaRogers<RedDbl> (logDensity, dim);
            break;
          case MINKL1:
            return new NormaMinkL1<RedDbl> (logDensity, dim);
            break;
          case MINKOWSKI:
            return new NormaMinkowski<RedDbl> (logDensity, dim);
            break;
          case NORMA_GENERIC:
            return new Normalizer<RedDbl> (logDensity, dim, "Norma_generic");
            break;
          case PALPHA_N:
            return new NormaPalpha<std::int64_t, RedDbl> (
                latanal.getReducer()->getIntLatticeBasis()->getModulo(), alpha, dim);
            break;
          case L1:
            break;
          case L2:
            break;
          default:
            std::cout << "LatticeAnalysis::initNormalizer:   no such case";
            exit (2);
        }
        return NULL;
      } 
    };

  // ZZXX case specialization
  template<typename Dbl, typename RedDbl>
      struct specLatticeAnalysis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> {
      Normalizer<RedDbl>* initNormalizer(LatticeAnalysis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>&
          latanal, NormaType norma, int alpha) {
        RedDbl logDensity;
        // We suppose we already have NTL::ZZ as integer type
        // Could have some kind of error detection here
        logDensity = -log(abs(NTL::determinant(
                latanal.getReducer()->getIntLatticeBasis()->getBasis())));
        int dim = latanal.getReducer()->getIntLatticeBasis()->getDim();
        switch (norma) {
          case BESTLAT:
            return new NormaBestLat<RedDbl> (logDensity, dim);
            break;
          case LAMINATED:
            return new NormaLaminated<RedDbl> (logDensity, dim);
            break;
          case ROGERS:
            return new NormaRogers<RedDbl> (logDensity, dim);
            break;
          case MINKL1:
            return new NormaMinkL1<RedDbl> (logDensity, dim);
            break;
          case MINKOWSKI:
            return new NormaMinkowski<RedDbl> (logDensity, dim);
            break;
          case NORMA_GENERIC:
            return new Normalizer<RedDbl> (logDensity, dim, "Norma_generic");
            break;
          case PALPHA_N:
            return new NormaPalpha<NTL::ZZ, RedDbl> (
                  latanal.getReducer()->getIntLatticeBasis()->getModulo(),
                  alpha, dim);
            break;
          case L1:
            break;
          case L2:
            break;
          default:
            std::cout << "LatticeAnalysis::initNormalizer:   no such case";
            exit (2);
        }
        return NULL;
      }
    };
  /// \endcond


  //===========================================================================

  // Utility functions for this class
  namespace {

    int getDir (std::string dir, std::vector <std::string> & files)
    {
      DIR *dp;
      struct dirent *dirp;
      if ((dp = opendir (dir.c_str())) == NULL) {
        std::cerr << "Directory: " << dir << std::endl;
        perror ("Couldn't open the directory");
        return errno;
      }

      // Does directory name ends with /
      size_t j = dir.rfind('/');
      std::string SEP("");
      // if not, add one /
      if (dir.size() != (1 + j))
        SEP += "/";

      while ((dirp = readdir (dp)) != NULL) {
        if (0 == fnmatch("*.dat", dirp->d_name, 0))
          // keeps full name including directory name
          files.push_back (std::string (dir + SEP + dirp->d_name));
      }
      closedir (dp);
      return 0;
    }

    //==========================================================================

    void eraseExtension (std::vector <std::string> & files)
    {
      for (unsigned int i = 0; i < files.size (); i++) {
        size_t j = files[i].rfind(".dat");
        if (j != std::string::npos)
          files[i].erase(j);
      }
    }

    //==========================================================================

    void printFileNames (std::vector <std::string> & files)
    {
      std::cout << "----------------------------------------------" << std::endl;
      for (unsigned int i = 0; i < files.size (); i++) {
        std::cout << files[i] << std::endl;
      }
    }
  }

  //===========================================================================
  // Class implementation

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::LatticeAnalysis ()
    {
      m_normalizer = 0;
      m_reducer = 0;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::LatticeAnalysis (
          Reducer<Int, BasInt, Dbl, RedDbl> & reducer,
          CriterionType criterion,
          PreReductionType preRed,
          NormType norm,
          std::int64_t maxNodesBB,
          NormaType normaType,
          int alpha)
    {
      m_reducer = &reducer;
      m_criterion = criterion;
      m_normalizerType = normaType;
      initNormalizer(normaType, alpha);
      m_preRed = preRed;
      m_norm = norm;
      Reducer<Int, BasInt, Dbl, RedDbl>::maxNodesBB = m_maxNodesBB = maxNodesBB;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::~LatticeAnalysis ()
    {
      if (m_normalizer != 0)
        delete m_normalizer;
    }

  //===========================================================================

  // This is the old implementation (templatized)
  /* 
   *    template<typename Int, typename IntVec, typename IntMat, typename BasInt,
   *    typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec,
   *    typename RedDbl, typename RedDblVec, typename RedDblMat>
   *    void LatticeAnalysis<Int, IntVec, IntMat, BasInt, BasIntVec, BasIntMat,
   *    Dbl, DblVec, RedDbl, RedDblVec, RedDblMat>::initNormalizer (
   *    NormaType norma, int alpha)
   *    {
   *    RedDbl logDensity;
   * #if NTL_TYPES_CODE > 1
   * RedDbl logDensity;
   * logDensity = - log( abs( NTL::determinant(
   * m_reducer->getIntLatticeBasis()->getBasis()) ) );
   * #else
   * // As NTL library does not support matrix with std::int64_t
   * // we compute the determinant with the boost library
   * boost::numeric::ublas::matrix<std::int64_t>  mat_tmps;
   * mat_tmps.resize(m_dim, m_dim);
   * for(unsigned int i = 0; i < m_dim; i++){
   * for(unsigned int j = 0; j < m_dim; j++){
   * mat_tmps(i,j) = m_reducer->getIntLatticeBasis()->getBasis()(i,j);
   * }
   * }
   * RedDbl logDensity(-log( abs( det_double(mat_tmps) ) ) );
   * #endif
   * switch (norma) {
   * case BESTLAT:
   * m_normalizer = new NormaBestLat<RedDbl> (logDensity, m_dim);
   * break;
   * case LAMINATED:
   * m_normalizer = new NormaLaminated<RedDbl> (logDensity, m_dim);
   * break;
   * case ROGERS:
   * m_normalizer = new NormaRogers<RedDbl> (logDensity, m_dim);
   * break;
   * case MINKL1:
   * m_normalizer = new NormaMinkL1<RedDbl> (logDensity, m_dim);
   * break;
   * case MINKOWSKI:
   * m_normalizer = new NormaMinkowski<RedDbl> (logDensity, m_dim);
   * break;
   * case NORMA_GENERIC:
   * m_normalizer = new Normalizer<RedDbl> (logDensity, m_dim, "Norma_generic");
   * break;
   * case PALPHA_N:
   * m_normalizer = new NormaPalpha<Int, RedDbl> (
   * m_reducer->getIntLatticeBasis()->getModulo(), alpha, m_dim);
   * break;
   * case L1:
   * break;
   * case L2:
   * break;
   * default:
   * std::cout << "LatticeAnalysis::initNormalizer:   no such case";
   * exit (2);
   * }
   * }
   * */

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::initNormalizer (
          NormaType norma, int alpha)
    {
      setNormalizerType(norma);
      setNormalizer(spec.initNormalizer(*this, norma, alpha));
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      bool LatticeAnalysis<Int, BasInt,Dbl, RedDbl>::doTest (
        double fact, PrecisionType precision, int blocksize)
    {
      bool result = false;

      // performing pre-reduction
      switch (m_preRed) {
        case BKZ:
          m_reducer->redBKZ(fact, blocksize, precision);
          break;
        case LenstraLL:
          m_reducer->redLLLNTL(fact, precision);
          break;
        case PreRedDieter:
          m_reducer->redDieter(0);
          break;
        case NOPRERED:
          break;
        default:
          MyExit(1, "LatticeLatticeAnalysis::doTest:   no such case");
          exit(1);
      }

      // performing the Branch-and-Bound procedure to find the shortest 
      // non-zero vector
      switch (m_criterion) {
        case SPECTRAL:
          // performing the Branch-and-Bound procedure to find the shortest 
          // non-zero vector
          result = m_reducer->shortestVector(m_norm);
          // calculating the Figure of Merit
          if(m_normalizerType == L1 || m_normalizerType == L2) {
            m_merit = NTL::conv<double>(m_reducer->getMinLength());
          } else {
            m_merit = NTL::conv<double>(m_reducer->getMinLength())
              / m_normalizer->getPreComputedBound(
                  m_reducer->getIntLatticeBasis()->getDim());
          }
          break;
        case BEYER:
          //performing the Branch-and-Bound procedure to find the 
          //Minkowski-reduced matrix
          result = m_reducer->reductMinkowski(0);
          // calculating the Figure of Merit
          m_merit = NTL::conv<double>(m_reducer->getMinLength())
            /NTL::conv<double>(m_reducer->getMaxLength());
          break;
        case PALPHA:
          MyExit(1, "PALPHA:   to be implemented");
          break;
        case BOUND_JS:
          MyExit(1, "BOUND_JS:   NOT YET");
          break;
        default:
          MyExit(1, "LatticeAnalysis::doTest:   NO SUCH CASE");
          exit(1);
      }

      return result;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::printTestResults ()
    {
      std::cout << "\n----------------------------------------------------------" 
        << std::endl;
      std::cout << "Criterion: " << toStringCriterion(m_criterion) << std::endl;
      std::cout << "Prereduction used: " << toStringPreRed(m_preRed) << std::endl;
      std::cout << "Length of shortest non-zero vector = " 
        << NTL::conv<double>(m_reducer->getMinLength());
      std::cout << " (" << toStringNorm(m_norm) << ")" << std::endl;
      std::cout << "Figure of Merit = " << m_merit;
      std::cout << " (" << toStringNorma(m_normalizerType) << " normalization)" 
        << std::endl;
      std::cout << "----------------------------------------------------------\n" 
        << std::endl;
    }

  //===========================================================================

  /*
   * Reads the test parameters in infile; then do the test.
   * infile is the data file name without extension: if the data file is named
   * "poil.dat", then infile is "poil".
   * Data files must always have the extension "dat".
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      int LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::doTestFromInputFile (
        const char *infile)
    {   
      // parameters reading
      std::string fname (infile);
      fname += ".dat";
      ParamReader<Int, BasInt, RedDbl> paramRdr (fname.c_str ());
      fname.clear ();

      LatticeTesterConfig<Int, BasIntMat> config;
      paramRdr.read (config);
      //config.write();

      // creating the Reducer object from input
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl> basis (config.basis,
          config.dim, config.norm);
      Reducer<Int, BasInt, Dbl, RedDbl> red (basis);
      // update parameters
      setReducer(red);
      setCriterion(config.test);
      setNormalizerType(config.normalizer);
      initNormalizer (config.normalizer);
      setPreReduction(config.prereduction);
      setNorm(config.norm);
      setMaxNodesBB(config.maxNodesBB);

      if (!doTest(config.fact, config.precision, config.blocksize)) {
        MyExit(1, "error in LatticeAnalysis::doTestFromInputFile");
        exit(1);
      }

      // putting the results in the output stream
      Writer<Int>* rw = createWriter (infile, config.outputType);

      rw->writeString(
          "\n----------------------------------------------------------");
      rw->newLine();
      rw->writeString("Criterion: "); 
      rw->writeString(toStringCriterion(m_criterion));
      rw->newLine();
      rw->writeString("Prereduction used: ");
      rw->writeString(toStringPreRed(m_preRed));
      rw->newLine();
      rw->writeString("Length of shortest non-zero vector = ");
      rw->writeDouble(NTL::conv<double>(m_reducer->getMinLength()));
      rw->writeString(" (");
      rw->writeString(toStringNorm(m_norm));
      rw->writeString(")");
      rw->newLine();
      rw->writeString("Figure of Merit = ");
      rw->writeDouble(m_merit);
      rw->writeString(" (");
      rw->writeString(toStringNorma(m_normalizerType));
      rw->writeString(" normalization)");
      rw->newLine();
      rw->writeString(
          "----------------------------------------------------------\n");
      rw->newLine();

      delete rw;
      return 0;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      int LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::doTestFromDirectory (
        const char *dirname)
    {
      std::string dir = std::string (dirname);
      std::vector <std::string> files = std::vector <std::string> ();

      getDir (dir, files);
      printFileNames (files);
      eraseExtension (files);

      int flag = 0;
      for (unsigned int i = 0; i < files.size (); i++)
        flag |= doTestFromInputFile (files[i].c_str());

      return flag;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Writer<Int>* LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::
      createWriter (const char *infile, OutputType ot)
    {
      Writer<Int> *rw = 0;
      std::string fname;

      switch (ot) {
        case RES:
          fname = infile;
          fname += ".res";
          rw = new WriterRes<Int> (fname.c_str ());
          break;

        case TEX:
          fname = infile;
          fname += ".tex";
          //rw = new WriterTex(fname.c_str()); //EB Ne permet pas d'écrire en Tex
          std::cerr << "\n*** outputType:   TEX not implemented" << std::endl;
          return 0;
          break;

        case TERMINAL:
          rw = new WriterRes<Int> (&std::cout);
          break;

        default:
          std::cerr << "\n*** outputType:   no such case" << std::endl;
          return 0;
      }
      return rw;
    }

} // end namespace

#endif
