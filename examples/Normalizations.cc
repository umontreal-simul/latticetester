/** 
 * Suppose that we have a `m`-lattice and its `m`-dual. We can get the figure of
 * merit called the normalized spectral test by dividing the length of the
 * shortest vector in the dual lattice with an upper bound. LatticeTester
 * provides several modules which are different approximations of this
 * upper bound. This example introduces them and presents the syntax to use
 * them.
 * */

#define NTL_TYPES_CODE 2

#include <iostream>

// To read the matrix we will use in this example
#include "latticetester/ParamReader.h"
// Normalization modules for the spectral test
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
// Normalization for the spectral test with L1NORM
#include "latticetester/NormaMinkL1.h"

#include "latticetester/Types.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"

using namespace LatticeTester;
int main() {
  // Creating a Reader for file `44matrixEx`.
  ParamReader<MScal, BScal, RScal> reader("build/examples/44matrixEx.dat");
  // Storing all the lines in `44matrixEx`
  reader.getLines();
  // Reading a matrix and printing it
  BMat matrix(4,4); // The "recipient"
  unsigned int ln = 0; // The line counter
  reader.readBMat(matrix, ln, 0, 4);
  // matrix now contains a non-singular matrix. We will use it to do our test.
  // Basically we will assume that matrix is the basis of the dual of the lattice
  // we are studying. If we want to study a lattice, even though we only store
  // its rescalled counterpart, the dual we store is actually the dual of the
  // lattice (and the `m`-dual of the `m`-lattice). Therefore, in the
  // representation that is done of the lattice with the IntLattice class,
  // to do the spectral test, we simply have to use the dualize method before
  // passing the basis to a Reducer object and use `shortestVector()`. For a 
  // more detailed explanation, see IntLattice.h.
  IntLatticeBasis<MScal, BScal, NScal, RScal> lat_basis(matrix, 4);
  // This is to play nice with the Reducer.
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Give the lattice to the reducer.
  Reducer<MScal, BScal, NScal, RScal> red(lat_basis);
  // We do a pre-reduction even though it is not needed in this example.
  // That is because if we pass a real lattice in high dimension, this speeds up
  // the program a lot.
  red.redBKZ();
  // We find the shortest vector with the L2NORM (euclidian norm)
  red.shortestVector(L2NORM);
  double shortest = NTL::conv<double>(red.getMinLength());
  // We now instanciate the different Normalization modules and print the
  // "different" results (the results won't actually be different because the
  // bound is exact in dimension 4 and all the modules have the same).
  //NormaRogers<RScal>
  return 0;
}
