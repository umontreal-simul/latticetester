/**
 * This example uses the LatticeTester API to read and print from/to files.
 * */
#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"

using namespace LatticeTester;

int main() {
  // Creating a Reader for file `44matrixEx`.
  ParamReader<MScal, BScal, RScal> reader("44matrixEx.dat");
  // Storing all the lines in `44matrixEx`
  reader.getLines();
  // Reading a matrix and printing it
  BMat matrix(4,4);
  unsigned int ln = 0;
  reader.readBMat(matrix, ln, 0, 4);
  std::cout << matrix << std::endl;
  return 0;
}
