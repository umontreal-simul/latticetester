#ifndef LATTICETEST_H
#define LATTICETEST_H
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Weights.h"
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
   LatticeTest (IntLatticeBasis * lattice, NormaType normaType);

   /**
    * Destructor.
    */
   ~LatticeTest ();

   /**
    * Gets the results of the applied test.
    */
   double & getMerit () { return m_merit; }

   /**
    * Starts the test in dimension `dim`.
    * Whenever the normalized value of the merit is smaller than `minVal`
    * for any dimension, the method returns `false` immediately. The
    * method returns `false` if the test was interrupted for any reason
    * before completion, and it returns `true` upon success. The results
    * of the test are kept in <tt>m_merit</tt>.
    */
   virtual bool test (int dim, double minVal[]) = 0;

   /**
    * Similar to `test` above, but with the weights `weights`.
    */
   virtual bool test (int dim, double minVal[], const double* weights);

private:

  /**
    * The lattice on which the test is applied.
    */
   IntLatticeBasis* m_lat;

  /**
    * The normalizer used for the the test.
    */
   normalizer* m_normalizer;

   /**
    * Contains the results of the test.
    */
   double m_merit;

};

} // namespace LatticeTester

#endif
