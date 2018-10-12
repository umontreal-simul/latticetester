/**
 * This example uses the LatticeTester API to read and print from/to files.
 * This does a really simple usage of the ParamReader and the WriterRes classes.
 * */
#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/WriterRes.h"

#include "NTL/LLL.h"

using namespace LatticeTester;

int main() {
  int size = 9;
  // Creating a Reader for file `44matrixEx`.
  ParamReader<MScal, BScal, RScal> reader("build/examples/44matrixEx.dat");
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
  BMat matrix(size, size); // The "recipient"
  unsigned int ln = 0; // The line counter
  reader.readBMat(matrix, ln, 0, size);

  // We want to save the results of the reductions to a file, rather than to the
  // standard output.
  WriterRes<MScal> writer("build/examples/IOExample.out");

  // Creating an IntLatticeBasis with the matrix we just read as a basis (it is
  // non-singular)
  // It is also possible to specify give a modulo and an dual basis when creating
  // an IntLattice, but we will not use them here.
  IntLatticeBasis<MScal, BScal, NScal, RScal> lat_basis(matrix, size);
  // We can evaluate and print the lenght of the vectors of this basis
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms before reductions:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  // Writing the same thing on standard output
  std::cout << "The basis and the vector norms before reductions:\n" <<
    lat_basis.toStringBasis();

  // Creating a Reducer
  Reducer<MScal, BScal, NScal, RScal> red(lat_basis);
  // Applying Dieter reduction and printing the new basis on standard output.
  red.redDieter(0);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms after Dieter reduction:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  // Writing the same thing on standard output
  std::cout << "The basis and the vector norms after Dieter reduction:\n" <<
    lat_basis.toStringBasis();

  // This is not really satisfying. Let's do the reduction with LLL this time
  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, size);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  red.redLLLNTL();
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms after LLL reduction:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  // Writing the same thing on standard output
  std::cout << "The basis and the vector norms after LLL reduction:\n" <<
    lat_basis.toStringBasis();

  // We got the same thing because we are on small dimensions. Still, let's
  // compare with BKZ.
  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, size);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  red.redBKZ();
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms after BKZ reduction:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  // Writing the same thing on standard output
  std::cout << "The basis and the vector norms after BKZ reduction:\n" <<
    lat_basis.toStringBasis();

  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, size);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Let's get the shortest vector to see how far we are with the reductions
  if (red.shortestVector(L2NORM)) {
    // Writing in the file
    writer.writeString("The basis and the vector norms with the shortest vector:");
    writer.newLine();
    writer.writeString(lat_basis.toStringBasis());
    // Writing the same thing on standard output
    std::cout << "The basis and the vector norms with the shortest vector:\n" <<
      lat_basis.toStringBasis();
  }

  std::cout << matrix << std::endl << lat_basis.getBasis() << std::endl;
  BVect x;
  x.SetLength(size);
  for (int i = 0; i < size; i++) {
    if (!NTL::LatticeSolve(x, lat_basis.getBasis(), matrix[i]))
      std::cout << "Not the same basis\n";
    else
      std::cout << "ok\n";
  }
  
  return 0;
}
