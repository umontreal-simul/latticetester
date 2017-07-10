#ifndef LATTICETEST_H
#define LATTICETEST_H
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Weights.h"
#include "latticetester/Reducer.h"
#include <string>
#include <list>


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


class LatticeTest {

public:

   /**
    * Constructor. The test will be applied on `lattice`, with the selected 
    * `normalizer`.
    */
   LatticeTest (Reducer & reducer, NormaType normaType, int alpha = 0);

   /**
    * Destructor.
    */
   ~LatticeTest ();

   /**
    * Gets the results of the applied test.
    */
   double & getMerit () { return m_merit; }

   /**
    * Performs the test in dimension `dim`.
    * The method returns `false` if the test was interrupted for any reason
    * before completion, and it returns `true` upon success. The result of 
    * the test is kept in <tt>m_merit</tt>.
    */
   bool performTest (double fact = 0.999999, long blockSize = 20);

   /**
    * Creates and returns the normalizer corresponding to criterion
    * `norma`. In the case of the \f$P_{\alpha}\f$ test, the argument
    * `alpha` = \f$\alpha\f$. In all other cases, it is unused.
    */
   void initNormalizer (NormaType norma, int alpha);

private:

   /**
    * The lattice on which the test is applied.
    */
   Reducer* m_reducer;

   /**
    * The type of normalizer used for the the test.
    */
   NormaType m_normaType;

   /**
    * The normalizer used for the the test.
    */

   Normalizer* m_normalizer;

   /**
    * Contains the results of the test.
    */
   double m_merit;

};

} // namespace LatticeTester

#endif
