/**
 * This example uses the LatticeTester API to read and print from/to files.
 * This does a really simple usage of the ParamReader and the WriterRes classes.
 *
 * ParamReader is a class that can easily be used to read a file, split its lines
 * and "tokenize" them so that the program can read files word by word and
 * convert the words to different types. This class can also easily be extended
 * to read other types than those already present, as is done in LatMRG.
 *
 * WriterRes is anothers class that simplifies the creation of a stream to write
 * to a file. Just as ParamReader, this can also convert data types into strings
 * to write directly in the file.
 * */
#define NTL_TYPES_CODE 1

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
  // Those two lines create a reader and ready it to be used.
  ParamReader<MScal, BScal, RScal> reader("44matrixEx.dat");
  reader.getLines();
  /** To use a reader, you simply need, once you've called getLines(), to have
   * a variable in which to store the object you read, and to know where to get
   * it. In a few cases, it might be necessary to pass a pointer to the variable.
   * Also, in some functions, such as this one, it might not be possible to pass
   * the line as a temporary variable at the function call because it will be
   * incremented during execution.
   * */
  // Reading a matrix as an examples
  BMat matrix(size, size); // The "recipient"
  unsigned int ln = 0; // The line counter
  reader.readBMat(matrix, ln, 0, size);

  // Creation of a writer object for file IOExample.out
  WriterRes<MScal> writer("IOExample.out");

  //! Creating an IntLatticeBasis with the matrix we just read as a basis (it is
  //! non-singular)
  //! It is also possible to specify give a modulo and an dual basis when creating
  //! an IntLattice, but we will not use them here.
  IntLatticeBasis<MScal, BScal, NScal, RScal> lat_basis(matrix, size);
  //! We can evaluate and print the lenght of the vectors of this basis
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms before reductions:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  writer.newLine();
  //! Writing the same thing on standard output
  //! std::cout << "The basis and the vector norms before reductions:\n" <<
  //!   lat_basis.toStringBasis() << std::endl;

  //! Creating a Reducer
  Reducer<MScal, BScal, NScal, RScal> red(lat_basis);
  //! Applying Dieter reduction and printing the new basis on standard output.
  red.redDieter(0);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The basis and the vector norms after Dieter reduction:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  writer.newLine();
  //! Writing the same thing on standard output
  //! std::cout << "The basis and the vector norms after Dieter reduction:\n" <<
  //!   lat_basis.toStringBasis() << std::endl;

  //! This is not really satisfying. Let's do the reduction with LLL this time
  lat_basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, size);
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  red.redLLLNTL();
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  //! Writing in the file
  writer.writeString("The basis and the vector norms after LLL reduction:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  writer.newLine()
  //! Writing the same thing on standard output
  //! std::cout << "The basis and the vector norms after LLL reduction:\n" <<
  //!   lat_basis.toStringBasis() << std::endl;

  // Writing a specific type in the file
  writer.writeString("The original matrix");
  writer.newLine();
  writer.writeMMat(matrix);

  return 0;
}
