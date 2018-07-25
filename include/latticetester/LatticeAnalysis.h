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
#include "latticetester/LatTestWriter.h"
#include "latticetester/LatTestWriterRes.h"
#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/ParamReader.h"

namespace LatticeTester {

  // Declaration is needed for specializer structure definition
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class LatticeAnalysis;

  /**
   * This structure specializes certain members of LatticeAnalysis. 
   * \todo Render this struct invisible outside of this file with an unnamed 
   * namespace or something like that.
   * */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      struct specLatticeAnalysis {
        void initNormalizer(LatticeAnalysis<Int, BasInt, Dbl, RedDbl>& latanal,
            NormaType norma, int alpha);
      };

  /**
   * This class gathers other classes of LatticeTester to create an object
   * performing tests on lattices. These tests are applied on lattices to
   * assess their structural properties and their qualities with respect to
   * different criteria. Included are well-known tests such as the *spectral*
   * test, the *Beyer* test, the \f$P_{\alpha}\f$ test. The corresponding
   * figures of merit for the lattice are the length of the shortest vector
   * in the lattice computed with different norms, the Beyer quotient, or the
   * \f$P_{\alpha}\f$ criterion. For the standard spectral test, the figure of 
   * merit is based on the length of the shortest non-zero vector in the 
   * lattice, using the \f${\mathcal{L}}_2\f$ norm to compute the length of 
   * vectors, and the inverse of this length gives the maximal distance between 
   * successive hyperplanes covering all the points in the *primal* lattice. If 
   * one computes the length of the shortest non-zero vector in the *dual* 
   * lattice using the \f${\mathcal{L}}_1\f$ norm, one obtains the minimal 
   * number of hyperplanes covering all the points of the *primal* lattice.
   *
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
           * Base constructor. The test will be applied on `lattice`, with the 
           * selected `normalizer`.
           */
          LatticeAnalysis ();

          /**
           * Constructor. The test will be applied on `lattice`, with the 
           * selected `normalizer` and `criterion`.
           */
          LatticeAnalysis (Reducer<Int, BasInt, Dbl, RedDbl> & reducer,
              CriterionType criterion, NormaType normaType, 
              PreReductionType preRed, NormType norm, int alpha = 0, 
              std::int64_t maxNodesBB = 10000000);

          /**
           * Destructor.
           */
          ~LatticeAnalysis ();

          /**
           * Performs the test in dimension `dim`.
           * The method returns `false` if the test was interrupted for any 
           * reason before completion, and it returns `true` upon success. The 
           * result of the test is kept in <tt>m_merit</tt>.
           */
          bool doTest (double fact, PrecisionType precision,
              int blocksize = 20);

          void printTestResults ();

          /**
           * Reads the parameters of the test in input text file `datafile`; then 
           * do the test. The data file must always have the extension `".dat"`, 
           * but must be given as argument here *without extension*. For example, 
           * if the data file is named `myLattice.dat`, then the method must be 
           * called as `doTest("myLattice")`. Returns 0 if the test completed 
           * successfully; returns a negative integer if there was an error.
           */
          int doTestFromInputFile (const char *datafile);

          /**
           * Applies the method `doTest` to all the files with extension `".dat"`
           * in directory named `dirname`. Returns 0 if all the tests completed
           * successfully; returns a non-zero integer if there was an error.
           */
          int doTestFromDirectory (const char *dirname);


          /**
           * Initialize m_normalizer to a pointer on a normalizer object
           * of type norma.
           */
          void initNormalizer (NormaType norma, int alpha = 0);

          /**
           * Gets the results of the applied test.
           */
          double getMerit () const { return m_merit; }

          /**
           * Set functions
           */
          void setReducer (Reducer<Int, BasInt, Dbl, RedDbl> & red)
          {
            m_reducer = &red;
          }

          void setCriterion (CriterionType criterion) 
          {
            m_criterion = criterion;
          }

          void setNormalizerType (NormaType normalizerType)
          {
            m_normalizerType = normalizerType;
          }

          void setNormalizer (Normalizer<RedDbl>* normalizer) 
          {
            m_normalizer = normalizer;
          }

          void setPreReduction (PreReductionType preRed) { m_preRed = preRed; }

          void setNorm (NormType norm) { m_norm = norm; }

          void setDim (int dim) { m_dim = dim; }

          void setMaxNodesBB (std::int64_t maxNodesBB)
          {
            Reducer<Int, BasInt, Dbl, RedDbl>::maxNodesBB = m_maxNodesBB 
              = maxNodesBB;
          }

          /**
           * Get functions.
           * */
          int getDim() { return m_dim;}

          Reducer<Int, BasInt, Dbl, RedDbl>* getReducer() 
            {
              return m_reducer;
            }

          Normalizer<RedDbl>* getNormalizer() { return m_normalizer;}

        private:
          /**
           * Pointer to the reducer class used to perform pre-reductions and 
           * Branch-and-Bound
           */
          Reducer<Int, BasInt, Dbl, RedDbl>* m_reducer;

          /**
           * Type of pre-reduction
           */
          PreReductionType m_preRed;

          /**
           * Type of test applied: SPECTRAL, BEYER, ...
           */
          CriterionType m_criterion;

          /**
           * Pointer to the Normalizer class used to normalize the results
           */
          Normalizer<RedDbl>* m_normalizer;

          /**
           * Type of normalization chosen for the test
           */
          NormaType m_normalizerType;

          /**
           * Norm used
           */
          NormType m_norm;

          /**
           * the dimension of the test
           */
          int m_dim;

          /**
           * Contains the results of the test
           */
          double m_merit;

          /**
           * Contains the maximum number of nodes visited in the Branch-and-Bound 
           * procedure
           */
          std::int64_t m_maxNodesBB;

          /**
           * Returns a `Writer` created from the input file `infile` and the 
           * given `OutputType`.
           */
          LatTestWriter<Int>* createLatTestWriter (const char *infile,
              OutputType ot);

          /**
           * This contains methods that need to be specialized
           * */
          struct specLatticeAnalysis<Int, BasInt, Dbl, RedDbl> spec;

      }; // End class LatticeAnalysis

  //===========================================================================

  // specLatticeAnalysis specialization

  // LLXX case specialization
  template<typename Dbl, typename RedDbl>
      struct specLatticeAnalysis<std::int64_t, std::int64_t, Dbl, RedDbl> {
      void initNormalizer(
          LatticeAnalysis<std::int64_t, std::int64_t, Dbl, RedDbl>& latanal,
          NormaType norma, int alpha) {
        RedDbl logDensity;
        // We have to cast to NTL::matrix<ZZ> because it does not compute the 
        // determinant of a general matrix by default.
        // This may cause bugs at places where the old implementation did not
        // because boost did not compute the good determinant for big matrixes
        // somehow.
        NTL::mat_ZZ temp;
        temp.SetDims(latanal.getDim(), latanal.getDim());
        for (int i = 0; i < latanal.getDim(); i++) {
          for (int j = 0; j < latanal.getDim(); j++){
            temp[i][j] = latanal.getReducer()->getIntLatticeBasis()->
              getBasis()(i,j);
          }
        }
        logDensity = -log(abs(NTL::determinant(temp)));
        switch (norma) {
          case BESTLAT:
            latanal.setNormalizer(new NormaBestLat<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case LAMINATED:
            latanal.setNormalizer(new NormaLaminated<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case ROGERS:
            latanal.setNormalizer(new NormaRogers<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case MINKL1:
            latanal.setNormalizer(new NormaMinkL1<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case MINKOWSKI:
            latanal.setNormalizer(new NormaMinkowski<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case NORMA_GENERIC:
            latanal.setNormalizer(new Normalizer<RedDbl> (logDensity,
                  latanal.getDim(), "Norma_generic"));
            break;
          case PALPHA_N:
            latanal.setNormalizer(new NormaPalpha<std::int64_t, RedDbl> (
                  latanal.getReducer()->getIntLatticeBasis()->getModulo(),
                  alpha, latanal.getDim()));
            break;
          case L1:
            break;
          case L2:
            break;
          default:
            std::cout << "LatticeAnalysis::initNormalizer:   no such case";
            exit (2);
        }
      } 
    };

  // ZZXX case specialization
  template<typename Dbl, typename RedDbl>
      struct specLatticeAnalysis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl> {
      void initNormalizer(LatticeAnalysis<NTL::ZZ, NTL::ZZ, Dbl, RedDbl>&
          latanal, NormaType norma, int alpha) {
        RedDbl logDensity;
        // We suppose we already have NTL::ZZ as integer type
        // Could have some kind of error detection here
        logDensity = -log(abs(NTL::determinant(
                latanal.getReducer()->getIntLatticeBasis()->getBasis())));
        switch (norma) {
          case BESTLAT:
            latanal.setNormalizer(new NormaBestLat<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case LAMINATED:
            latanal.setNormalizer(new NormaLaminated<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case ROGERS:
            latanal.setNormalizer(new NormaRogers<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case MINKL1:
            latanal.setNormalizer(new NormaMinkL1<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case MINKOWSKI:
            latanal.setNormalizer(new NormaMinkowski<RedDbl> (logDensity,
                  latanal.getDim()));
            break;
          case NORMA_GENERIC:
            latanal.setNormalizer(new Normalizer<RedDbl> (logDensity,
                  latanal.getDim(), "Norma_generic"));
            break;
          case PALPHA_N:
            latanal.setNormalizer(new NormaPalpha<NTL::ZZ, RedDbl> (
                  latanal.getReducer()->getIntLatticeBasis()->getModulo(),
                  alpha, latanal.getDim()));
            break;
          case L1:
            break;
          case L2:
            break;
          default:
            std::cout << "LatticeAnalysis::initNormalizer:   no such case";
            exit (2);
        }
      }
    };


  //===========================================================================
  // Utility functions

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

  //===========================================================================

  void eraseExtension (std::vector <std::string> & files)
  {
    for (unsigned int i = 0; i < files.size (); i++) {
      size_t j = files[i].rfind(".dat");
      if (j != std::string::npos)
        files[i].erase(j);
    }
  }

  //===========================================================================

  void printFileNames (std::vector <std::string> & files)
  {
    std::cout << "----------------------------------------------" << std::endl;
    for (unsigned int i = 0; i < files.size (); i++) {
      std::cout << files[i] << std::endl;
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
          Reducer<Int, BasInt, Dbl, RedDbl> 
        & reducer, CriterionType criterion, NormaType normaType,
        PreReductionType preRed, NormType norm, int alpha, std::int64_t maxNodesBB)
    {
      m_reducer = &reducer;
      m_criterion = criterion;
      m_normalizerType = normaType;
      initNormalizer(normaType, alpha);
      m_preRed = preRed;
      m_norm = norm;
      m_dim = m_reducer->getIntLatticeBasis()->getDim();
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
      spec.initNormalizer(*this, norma, alpha);
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
          m_reducer->preRedDieter(0);
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
              / m_normalizer->getPreComputedBound(m_dim);
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

      LatTestWriter<Int>* rw = createLatTestWriter (infile, config.outputType);

      // creating the Reducer object from input
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl> basis (config.basis,
          config.dim, config.norm);
      Reducer<Int, BasInt, Dbl, RedDbl> red (basis);
      // update parameters
      setReducer(red);
      setCriterion(config.test);
      setNormalizerType(config.normalizer);
      setDim(config.dim);
      initNormalizer (config.normalizer);
      setPreReduction(config.prereduction);
      setNorm(config.norm);
      setMaxNodesBB(config.maxNodesBB);

      if (!doTest(config.fact, config.precision, config.blocksize)) {
        MyExit(1, "error in LatticeAnalysis::doTestFromInputFile");
        exit(1);
      }

      // putting the results in the output stream
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
      LatTestWriter<Int>* LatticeAnalysis<Int, BasInt, Dbl, RedDbl>::
      createLatTestWriter (const char *infile, OutputType ot)
    {
      LatTestWriter<Int> *rw = 0;
      std::string fname;

      switch (ot) {
        case RES:
          fname = infile;
          fname += ".res";
          rw = new LatTestWriterRes<Int> (fname.c_str ());
          break;

        case TEX:
          fname = infile;
          fname += ".tex";
          //rw = new WriterTex(fname.c_str()); //EB Ne permet pas d'écrire en Tex
          std::cerr << "\n*** outputType:   TEX not implemented" << std::endl;
          return 0;
          break;

        case TERMINAL:
          rw = new LatTestWriterRes<Int> (&std::cout);
          break;

        default:
          std::cerr << "\n*** outputType:   no such case" << std::endl;
          return 0;
      }
      return rw;
    }

} // end namespace

#endif
