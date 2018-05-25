#ifndef LATTICETESTERCONFIG_H
#define LATTICETESTERCONFIG_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"


namespace LatticeTester {

  /**
   * This class is used to save the configuration of a lattice test. It is used
   * to keep all the parameters read in the data file and passed to different
   * methods for the spectral, Beyer or \f$P_{\alpha}\f$ tests.
   *
   */
  class LatticeTesterConfig {
    public:
      static const int MAX_WORD_SIZE = 64;

      /**
       * Constructor.
       */
      LatticeTesterConfig();

      /**
       * Destructor.
       */
      ~LatticeTesterConfig();

      /**
       * Frees the memory used by this object.
       */
      void kill();

      /**
       * For debugging: writes the configuration on the console.
       */
      void write();

      /**
       * This flag is set `true` if the generator is to be read from a file,
       * otherwise it is set `false`.
       */
      bool readGenFile;

      /**
       * If `readGenFile` is `true`, the name of the file from which to read
       * the generator.
       */
      std::string fileName;

      /*
       * The dimensin of the lattice
       */
      int dim;

      /**
       * The bound used for the normalization in the definition of \f$S_t\f$.
       * Only applicable for the spectral test.
       */
      NormaType normalizer;

      /**
       * The bound used for the normalization in the definition of \f$S_t\f$.
       * Only applicable for the spectral test.
       */
      PreReductionType prereduction;

      /**
       * This flag indicates to print more detailed results if `detailF`
       * \f$>0\f$. Default value: 0.
       */
      int detailF;

      /**
       * Norm used to measure the length of vectors. See module `Const` for a
       * definition of the possible norms.
       */
      NormType norm;

      /**
       * The criterion for which the test will be performed. See module
       * `Const` for the possible criterion types.
       */
      CriterionType test;

      /**
       * The precision in which the test will be performed. See module
       * `Const` for the possible precision types.
       */
      PrecisionType precision;

      /**
       * The maximum number of nodes to be examined in any given
       * branch-and-bound procedure when computing \f$d_t\f$ or \f$q_t\f$.
       */
      long maxNodesBB;

      /**
       * The basis of the Lattice where the test will be apply
       */
      BMat basis;

      /**
       * The modulo of the rank1 lattice for Palpha test
       */
      MScal modulo;

      /**
       * The alpha coefficient for Palpha test (for rank 1 lattices)
       */
      int alpha;

      /**
       * Used in LLL and BKZ. See Reducer. It must be smaller than
       * 1. If `fact` is closer to 1, the basis will be (typically) "more
       * reduced", but that will require more work.
       */
      double fact;

      /**
       * Used in BKZ. It stocks the number of blocks used for the BKZ reduction.
       * The more it is large, the more the basis will be reduced.
       */
      int blocksize;


      /**
       * File format used to store the results. See `Const` for a definition
       * of the possible output types.
       */
      OutputType outputType;
  };

}
#endif
