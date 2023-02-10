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

#ifndef LATTICETESTER_LATTICEANALYSIS_H
#define LATTICETESTER_LATTICEANALYSIS_H

#include <string>
#include <list>
#include <dirent.h>
#include <fnmatch.h>
#include <typeinfo>
#include <cstdint>

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkL2.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/Reducer.h"
#include "latticetester/Writer.h"
#include "latticetester/WriterRes.h"
#include "../examples/Config.h"
//#include "Config.h"
#include "latticetester/ParamReader.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"

namespace LatticeTester {

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
   * test.
   *
   * To use this class it is imperative to instanciate a Reducer object with the
   * basis of the lattice to analyse. After that, you have to pass, either by
   * the constructor or by the various methods available, all the parameters 
   * necessary to do a test (these are listed under the empty constructor).
   * Finally, you simply have to call the doTest() method and the `m_merit`
   * field will contain the figure of merit that was asked for.
   */

  template<typename Int, typename Real, typename RealRed>
     class LatticeAnalysis {
        private:
          typedef NTL::vector<Int> IntVec;
          typedef NTL::matrix<Int> IntMat;
          typedef NTL::vector<Real> RealVec;
          typedef NTL::vector<RealRed> RealRedVec;
          typedef NTL::matrix<RealRed> RealRedMat;
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
          LatticeAnalysis (Config<Int, IntMat>& config);

          /**
           * Destructor. This will deallocate memory to `m_normalizer` if it is
           * not null.
           */
          ~LatticeAnalysis (){}

          /**
           * Launches the computation specified in config. This will not print
           * the results, but they can be retrieved via the printTestResults()
           * method.
           * This will store the `config` passed to this object for printing
           * latter.
           * */
          bool doTest(Config<Int, IntMat>& config);

          /**
           * This prints the result of the last computation that was asked to
           * this object, as well as the configuration used.
           * */
          void printTestResults (const char* infile);

          /**
           * Reads the parameters of the test in input text file `datafile`;
           * then does the specified computation. The format of the data file
           * must be as presented in \ref usage.
           * This will also print the results of the test in the way specified
           * by the configuration file.
           *
           * The data file must always have the extension ".dat", 
           * but must be given as argument here *without extension*. For example, 
           * if the data file is named `myLattice.dat`, then the method must be 
           * called as `doTestFromInputFile("myLattice")`.
           *
           * This method returns 0 if the test completed successfully and
           * returns a negative integer if there was an error.
           */
          int doTestFromInputFile (const char *datafile);

          /**
           * Applies the method `doTestFromInputFile` to all the files with
           * extension `".dat"` in directory named `dirname`. This will print
           * all the results, but will only keep in memory the results of the
           * last computation done.
           *
           * Returns 0 if all
           * the tests completed successfully; returns a non-zero integer if
           * there was an error. Even if there is an error with one of the files,
           * this method will try to do the tests on all the files.
           */
          int doTestFromDirectory (const char *dirname);


          Config<Int, IntMat>* getConfig() {
            return m_config;
          }

        protected:
          /**
           * This will perform a basis construction according to the
           * configuration and the basis of `config`.
           * */
          bool performBasis(Config<Int, IntMat>& config);

          /**
           * This will perform a dual construction according to the
           * configuration and the basis of `config`. `config` is stored in this
           * object so that the results can be printed with printTestResults().
           * */
          bool performDual(Config<Int, IntMat>& config);

          /**
           * This will perform a lattice reduction according to the
           * configuration and the basis of `config`. `config` is stored in this
           * object so that the results can be printed with printTestResults().
           * */
          bool performReduction(Config<Int, IntMat>& config);

          /**
           * This will solve the shortest vector problem according to the
           * configuration and the basis of `config`. `config` is stored in this
           * object so that the results can be printed with printTestResults().
           * */
          bool performShortest(Config<Int, IntMat>& config);

          /**
           * This will compute a figure of merit according to the
           * configuration and the basis of `config`. `config` is stored in this
           * object so that the results can be printed with printTestResults().
           * */
          bool performMerit(Config<Int, IntMat>& config);

          /**
           * This contains the last figure of merit computed. This is not erased
           * when a new computation not requiring a figure of merit is done.
           */
          double m_merit;

          /**
           * This contains the length of the shortest vector that was last
           * computed. This is not erased when a new computation not requiring
           * the shortest vector is done.
           */
          double m_shortest;

        private:

          /**
           * This creates and returns a Writer object
           * Returns a `Writer` created from the input file `infile` and the 
           * given `OutputType`.
           */
          Writer<Int>* createWriter (const char *infile,
              OutputType ot);

          /**
           * The configuration of the last test this object has done.
           * */
          Config<Int, IntMat>* m_config;

      }; // End class LatticeAnalysis

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

  template<typename Int, typename Real, typename RealRed>
      LatticeAnalysis<Int, Real, RealRed>::LatticeAnalysis ()
    {
      m_config = NULL;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
      LatticeAnalysis<Int, Real, RealRed>::LatticeAnalysis (
          Config<Int, IntMat>& config)
    {
      m_config = &config;
    }

  //===========================================================================

  /*
   * This is a big ugly function that will print the results after the program
   * execution.
   * */
  template<typename Int, typename Real, typename RealRed>
      void LatticeAnalysis<Int, Real, RealRed>::printTestResults (const char* infile)
    {
      // putting the results in the output stream
      Writer<Int>* rw = createWriter (infile, m_config->outputType);

      rw->writeString(
          "----------------------------------------------------------\n");
      rw->newLine();
      rw->writeString("Problem: "); 
      rw->writeString(toStringProblem(m_config->prob));
      rw->newLine();
      if (m_config->prob == BASIS) {
        rw->writeString("Construction method: ");
        if (m_config->config.basis.method) {
          rw->writeString("GCD");
        } else {
          rw->writeString("LLL");
        }
        rw->newLine();
        rw->writeString("Resulting matrix:\n");
        rw->writeMMat(m_config->basis);
      } else if (m_config->prob == DUAL) {
        rw->writeString("Rescaling factor (m): ");
        rw->getStream() << m_config->m;
        rw->newLine();
        rw->writeString("Dual matrix:\n");
        rw->writeMMat(m_config->dual_basis);
      } else if (m_config->prob == REDUCTION) {
        rw->writeString("Reduction method: ");
        rw->writeString(toStringPreRed(m_config->config.reduct.method));
        rw->newLine();
        rw->writeString("Resulting matrix:\n");
        rw->writeMMat(m_config->basis);
      } else if (m_config->prob == SHORTEST) {
        rw->writeString("Reduction used?");
        if (m_config->config.shortest.reduction) {
          rw->writeString("(y)\n");
          rw->writeString("Reduction method: ");
          rw->writeString(toStringPreRed(m_config->config.shortest.method));
          rw->newLine();
        } else rw->writeString("(n)\n");
        rw->writeString("Resulting basis:\n");
        rw->writeMMat(m_config->basis);
        rw->newLine();
        rw->writeString("Length of the shortest vector: ");
        rw->writeDouble(m_shortest);
      } else if (m_config->prob == MERIT) {
        rw->writeString("Figure of merit computed: ");
        rw->writeString(toStringCriterion(m_config->config.merit.figure));
        rw->newLine();
        rw->writeString("Reduction used? ");
        if (m_config->config.merit.reduction) {
          rw->writeString("(y)\n");
          rw->writeString("Reduction method: ");
          rw->writeString(toStringPreRed(m_config->config.merit.method));
          rw->newLine();
        } else rw->writeString("(n)\n");
        rw->writeString("Normalization used: ");
        rw->writeString(toStringNorma(m_config->config.merit.norma));
        rw->newLine();
        rw->writeString("Resulting basis:\n");
        rw->writeMMat(m_config->basis);
        rw->newLine();
        rw->writeString("Length of the shortest vector: ");
        rw->writeDouble(m_shortest);
        rw->newLine();
        rw->writeString("Figure of merit computed: ");
        rw->writeDouble(m_merit);
      } 
      rw->newLine();
      rw->newLine();
      rw->writeString(
          "----------------------------------------------------------\n");
      delete rw;
    }

  //===========================================================================

  /*
   * Reads the test parameters in infile; then do the test.
   * infile is the data file name without extension: if the data file is named
   * "poil.dat", then infile is "poil".
   * Data files must always have the extension "dat".
   */
  template<typename Int, typename Real, typename RealRed>
      int LatticeAnalysis<Int, Real, RealRed>::doTestFromInputFile (
        const char *infile)
    {   
      bool result = false;
      // parameters reading
      std::string fname (infile);
      fname += ".dat";
      ParamReader<Int, RealRed> paramRdr (fname.c_str ());
      fname.clear ();

      Config<Int, IntMat> config;
      paramRdr.read (config);
      doTest(config);

      printTestResults(infile);

      return result;
    }

  //============================================================================

  template<typename Int, typename Real, typename RealRed>
    bool LatticeAnalysis<Int, Real, RealRed>::doTest(
        Config<Int, IntMat>& config) {
      m_config = &config;
      bool result = false;
      if (config.prob == BASIS) {
        result = performBasis(config);
      } else if (config.prob == DUAL) {
        result = performDual(config);
      } else if (config.prob == REDUCTION) {
        result = performReduction(config);
      } else if (config.prob == SHORTEST) {
        result = performShortest(config);
      } else if (config.prob == MERIT) {
        result = performMerit(config);
      }
      return result;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
    bool LatticeAnalysis<Int, Real, RealRed>::performBasis(
        Config<Int, IntMat>& config) {
          Int m(1021);
      BasisConstruction<Int> basis;
      if (config.config.basis.method) {
        basis.GCDTriangularBasis(config.basis,m);
      } else {
        basis.LLLConstruction(config.basis);
      }
      return true;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
    bool LatticeAnalysis<Int, Real, RealRed>::performDual(
        Config<Int, IntMat>& config) {
      config.m = 1;
      BasisConstruction<Int> basis;
      basis.mDualTriangular(config.basis, config.dual_basis, config.m);
      return true;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
    bool LatticeAnalysis<Int, Real, RealRed>::performReduction(
        Config<Int, IntMat>& config) {

      IntLatticeBase<Int, Real, RealRed>
        Basis(config.basis, config.NumCols);
      Reducer<Int, Real, RealRed> Red(Basis);
      PreReductionType reduction = config.config.reduct.method;

      if (reduction == BKZ) {
        Red.redBKZ();
      } else if (reduction == LLL) {
        Red.redLLLNTL();
      } else if (reduction == DIETER) {
        Red.redDieter(0);
      }
      config.basis = Basis.getBasis();
      return true;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
    bool LatticeAnalysis<Int, Real, RealRed>::performShortest(
        Config<Int, IntMat>& config) {
      bool result = false;
      IntLatticeBase<Int, Real, RealRed>
        Basis(config.basis, config.NumCols);
      Reducer<Int, Real, RealRed> Red(Basis);

      if (config.config.shortest.reduction) {
        PreReductionType reduction = config.config.shortest.method;
        if (reduction == BKZ) {
          Red.redBKZ();
        } else if (reduction == LLL) {
          Red.redLLLNTL();
        } else if (reduction == DIETER) {
          Red.redDieter(0);
        }
      }
      std::string ch("cholesky");
      result = Red.shortestVector(L2NORM,ch);
      m_shortest = NTL::conv<double>(Red.getMinLength());
      config.basis = Basis.getBasis();

      return result;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
      bool LatticeAnalysis<Int, Real, RealRed>::performMerit (
          Config<Int, IntMat>& config)
    {
      bool result = false;
      IntLatticeBase<Int, Real, RealRed>
        Basis(config.basis, config.NumCols);
      Reducer<Int, Real, RealRed> Red(Basis);

      if (config.config.merit.reduction) {
        PreReductionType reduction = config.config.merit.method;
        if (reduction == BKZ) {
          Red.redBKZ();
        } else if (reduction == LLL) {
          Red.redLLLNTL();
        } else if (reduction == DIETER) {
          Red.redDieter(0);
        }
      }

      // performing the Branch-and-Bound procedure to find the shortest 
      // non-zero vector
      if (config.config.merit.figure == SPECTRAL) {
        // performing the Branch-and-Bound procedure to find the shortest 
        // non-zero vector
         std::string ch("cholesky");
        result = Red.shortestVector(L2NORM,ch);
        // calculating the Figure of Merit
        NormaType norma(config.config.merit.norma);
        //RealRed density = RealRed(-log(abs(NTL::determinant(config.basis))));
        double density = (double)(-log(abs(NTL::determinant(config.basis))));
        Normalizer* normalizer = NULL;
        m_shortest = NTL::conv<double>(Red.getMinLength());
        if (norma == NONE) {
          m_merit = m_shortest;
        } else if (norma == BESTLAT) {
          normalizer = new NormaBestLat(density, config.NumCols);
          m_merit = NTL::conv<double>(Red.getMinLength())
            / normalizer->getBound(config.NumCols);
        } else if (norma == BESTBOUND) {
          normalizer = new NormaBestBound(density, config.NumCols);
          m_merit = NTL::conv<double>(Red.getMinLength())
            / normalizer->getBound(config.NumCols);
        } else if (norma == MINKL2) {
          normalizer = new NormaMinkL2(density, config.NumCols);
          m_merit = NTL::conv<double>(Red.getMinLength())
            / normalizer->getBound(config.NumCols);
        } else if (norma == LAMINATED) {
          normalizer = new NormaLaminated(density, config.NumCols);
          m_merit = NTL::conv<double>(Red.getMinLength())
            / normalizer->getBound(config.NumCols);
        } 
        if (normalizer != NULL) delete normalizer;
      } else if (config.config.merit.figure == BEYER) {
        //performing the Branch-and-Bound procedure to find the 
        //Minkowski-reduced matrix
        result = Red.reductMinkowski(0);
        // calculating the Figure of Merit
        m_merit = NTL::conv<double>(Red.getMinLength())
          /NTL::conv<double>(Red.getMaxLength());
      }
      config.basis = Basis.getBasis();

      return result;
    }

  //===========================================================================

  template<typename Int, typename Real, typename RealRed>
      int LatticeAnalysis<Int, Real, RealRed>::doTestFromDirectory (
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

  template<typename Int, typename Real, typename RealRed>
      Writer<Int>* LatticeAnalysis<Int, Real, RealRed>::
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

        case TERM:
          rw = new WriterRes<Int> (&std::cout);
          break;

        default:
          std::cerr << "\n*** outputType:   no such case" << std::endl;
          return 0;
      }
      return rw;
    }

  extern template class LatticeAnalysis<std::int64_t, double, double>;
  extern template class LatticeAnalysis<NTL::ZZ, double, double>;
  extern template class LatticeAnalysis<NTL::ZZ, NTL::RR, NTL::RR>;

} // end namespace

#endif
