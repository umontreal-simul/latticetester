/**
 * This example uses the LatticeTester API to read and print from/to files.
 * */
#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBasis.h"

using namespace LatticeTester;

int main() {
  // Creating a Reader for file `44matrixEx`.
  ParamReader<MScal, BScal, RScal> reader("examples/44matrixEx.dat");
  // Storing all the lines in `44matrixEx`
  reader.getLines();
  /* To use a reader, you simply need, once you've called getLines(), to have
   * a variable in which to store the object you read, and to know where to get
   * it. In a few cases, it might be necessary to pass a pointer to the variable.
   * Also, in some functions, such as this one, it might not be possible to pass
   * the line as a temporary variable at the function call because it will be
   * incremented during execution.
   * */
  // Reading a matrix and printing it
  BMat matrix(4,4); // The "recipient"
  unsigned int ln = 0; // The line counter
  reader.readBMat(matrix, ln, 0, 4);
  // Printing the matrix to see if it reading it was successful.
  std::cout << "We read matrix:\n" <<matrix << std::endl;

  // Creating an IntLatticeBasis with the matrix we just read as a basis (it is
  // non-singular)
  // It is also possible to specify give a modulo and an dual basis when creating
  // an IntLattice, but we will not use them here.
  IntLatticeBasis<MScal, BScal, NScal, RScal> lat_basis(matrix, 4);
  // We can evaluate and print the lenght of the vectors of this basis
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  std::cout << "The basis and the vector norms before reductions:\n" <<
    lat_basis.toStringBasis();

  // Creating a Reducer
  Reducer<MScal, BScal, NScal, RScal> red(lat_basis);
  // Applying Dieter reduction and printing the new basis on standard output.
  red.redDieter(0);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  std::cout << "The basis and the vector norms after Dieter reduction:\n" <<
    lat_basis.toStringBasis();

  // This is not really satisfying. Let's do the reduction with LLL this time
  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, 4);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  red.redLLLNTL();
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  std::cout << "The basis and the vector norms after LLL reduction:\n" <<
    lat_basis.toStringBasis();

  // We got the same thing because we are on small dimensions. Still, let's
  // compare with BKZ.
  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, 4);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  red.redBKZ();
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  std::cout << "The basis and the vector norms after BKZ reduction:\n" <<
    lat_basis.toStringBasis();

  // Let's get the shortest vector to see how far we are with the reductions
  if (red.shortestVector(L2NORM)) {
    std::cout << "The basis and the vector norms with the shortest vector:\n" <<
      lat_basis.toStringBasis();
  }

  return 0;
}
