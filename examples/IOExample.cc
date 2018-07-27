/**
 * This example uses the LatticeTester API to read and print from/to files.
 * */
#define NTL_TYPES_CODE=2

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"

int main() {
  // Creating a Reader for file `44matrixEx`.
  ParamReader<MScal, BScal, RScal> reader("44matrixEx");
  // Storing all the lines in `44matrixEx`
  reader.getLines();
  return 0;
}
