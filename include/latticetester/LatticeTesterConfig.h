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
    * Reinitializes this object and allocates enough memory for \f$j\f$
    * MRGs.
    */
   void setJ (int j);

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
   int dimension;

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
   LatticeTester::NormType norm;

   /**
    * This flag is set `true` if the test is applied for lacunary indices.
    * If it is `false`, the test is applied for successive indices.
    */
   bool lacunary;

   /**
    * Used for lacunary indices. If the respective values are \f$s\f$ and
    * \f$d\f$, then the program will analyze the lattice structure of vectors
    * formed by groups of \f$s\f$ successive values, taken \f$d\f$ values apart,
    * i.e. groups of the form \f$(u_{i+1}, …, u_{i+s}, u_{i+d+1}, …, u_{i+d+s},
    * u_{i+2d+1}, …, u_{i+2d+s}, …)\f$.
    */
   int lacGroupSize;

   /**
    * \copydoc lacGroupSize
    */
   NTL::ZZ lacSpacing;

   /**
    * The lacunary indices, either read explicitly or computed from
    * `lacGroupSize` and `lacSpacing`.
    */
   BVect Lac;

   /**
    * Is `true` when the modulus of congruence \f$m\f$ is a prime number,
    * is `false` otherwise.
    */
   bool primeM;

   /**
    * If `true`, the program will verify that the modulus \f$m\f$ is a
    * prime number. If `false`, will not verify it.
    */
   bool verifyM;

   /**
    * Is `true` when the generator has maximal period, is `false`
    * otherwise.
    */
   bool maxPeriod;

   /**
    * If `true`, the program will verify that the generator has maximal
    * period. If `false`, will not verify it.
    */
   bool verifyP;

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
