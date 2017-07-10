#ifndef LATTICEANALYSIS_H
#define LATTICEANALYSIS_H

#include <string>
#include <list>

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Weights.h"
#include "latticetester/Reducer.h"

//#include "latticetester/Writer.h"
#include "latticetester/LatticeTesterConfig.h"

namespace LatticeTester {

/**
 * This class gathers other classes of LatticeTester to create an object
 * performing tests on lattices. These tests are applied on lattices to
 * assess their structural properties and their qualities with respect to
 * different criteria. Included are well-known tests such as the *spectral*
 * test, the *Beyer* test, the \f$P_{\alpha}\f$ test. The corresponding
 * figures of merit for the lattice are the length of the shortest vector
 * in the *primal* or in the *dual* lattice computed with different norms,
 * the Beyer quotient, or the \f$P_{\alpha}\f$ criterion. For the standard
 * spectral test, the figure of merit is based on the length of the shortest
 * non-zero vector in the *dual* lattice, using the \f${\mathcal{L}}_2\f$ norm
 * to compute the length of vectors, and the inverse of this length gives the
 * maximal distance between successive hyperplanes covering all the points in
 * the *primal* lattice. If one computes the length of the shortest non-zero
 * vector in the *dual* lattice using the \f${\mathcal{L}}_1\f$ norm, one
 * obtains the minimal number of hyperplanes covering all the points of the
 * *primal* lattice.
 *
 */

 // PW_TODO update description

class LatticeAnalysis {

public:
  /**
    * Constructor. The test will be applied on `lattice`, with the selected
    * `normalizer`.
    */
   LatticeAnalysis (Reducer* reducer, NormType norm, 
      NormaType normaType, int alpha = 0);

   /**
    * Destructor.
    */
   ~LatticeAnalysis ();

   /**
    * Performs the test in dimension `dim`.
    * The method returns `false` if the test was interrupted for any reason
    * before completion, and it returns `true` upon success. The result of
    * the test is kept in <tt>m_merit</tt>.
    * PW_TODO : description à améliorer
    */
   bool performTest (PreReductionType PreRed, double fact, 
                                      long blockSize = 0);

   /**
    * Reads the parameters of the test in input text file `datafile`; then do
    * the test. The data file must always have the extension `".dat"`, but must
    * be given as argument here *without extension*. For example, if the data
    * file is named `myLattice.dat`, then the method must be called as
    * `doTest("myLattice")`. Returns 0 if the test completed successfully; returns
    * a negative integer if there was an error.
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
   double & getMerit () { return m_merit; }

   void setReducer (Reducer & red) {m_reducer = &red; }
   void setDim (int dim) { m_dim = dim; }
   void setNorm (NormType norm) { m_norm = norm; }
   void setNormalizerType (NormaType normalizerType) { m_normalizerType = normalizerType; }
    
   

private:

  Reducer* reducer;
  int m_dim;
  NormType m_norm;
  NormaType m_normalizerType;
  Normalizer* m_normalizer;
  

   /**
    * Contains the results of the test.
    */
   double m_merit;

   /**
    * Returns a `Writer` created from the input file `infile` and the given
    * `OutputType`.
    */
   Writer* createWriter (const char *infile, OutputType ot);

};

} // end namespace

#endif
