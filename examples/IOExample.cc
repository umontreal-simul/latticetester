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
 *
 * The output of this program look like this : 
 * We can write messages in the file.
 * Here is one way to write a matrix.
 *    [[927 109 930 5 943 663 207 454 496]
 *     [334 1015 884 800 447 818 247 360 658]
 *     [830 391 721 95 1004 777 478 927 653]
 *     [900 565 378 327 598 114 662 407 483]
 *     [352 662 565 241 3 78 505 645 391]
 *     [725 990 80 500 931 301 107 717 972]
 *     [348 291 524 190 795 121 953 670 209]
 *     [958 367 927 803 755 872 144 329 57]
 *     [589 14 119 722 436 415 665 252 616]]
 * 
 * Here is another one.
 * [927 109 930 5 943 663 207 454 496]
 * [334 1015 884 800 447 818 247 360 658]
 * [830 391 721 95 1004 777 478 927 653]
 * [900 565 378 327 598 114 662 407 483]
 * [352 662 565 241 3 78 505 645 391]
 * [725 990 80 500 931 301 107 717 972]
 * [348 291 524 190 795 121 953 670 209]
 * [958 367 927 803 755 872 144 329 57]
 * [589 14 119 722 436 415 665 252 616]
 * 
 * It's also possible to use *ToString() methods to write stuff in the output
 * file.
 * The matrix and its vectors norms:
 * Primal Basis:
 *   Dim = 9 
 *     [ 352 662 565 241 3 78 505 645 391 ]
 *     [ 589 14 119 722 436 415 665 252 616 ]
 *     [ 348 291 524 190 795 121 953 670 209 ]
 *     [ 900 565 378 327 598 114 662 407 483 ]
 *     [ 927 109 930 5 943 663 207 454 496 ]
 *     [ 958 367 927 803 755 872 144 329 57 ]
 *     [ 334 1015 884 800 447 818 247 360 658 ]
 *     [ 725 990 80 500 931 301 107 717 972 ]
 *     [ 830 391 721 95 1004 777 478 927 653 ]
 *   Norms:
 *     [1769478 2130068 2563917 2586820 3559934 4019226 4055743 4189809 4496614 ]
 * 
 * Writers can be used as the program progresses, just likeusing a standard
 * output print, but also can be created at the end of theexecution to print the
 * result of the execution, making them quiteflexible.
 * */
#define NTL_TYPES_CODE 1

#include <iostream>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/WriterRes.h"

#include "NTL/LLL.h"

using namespace LatticeTester;

int main() {
  int size = 9;
  // Those two lines create a reader and ready it to be used.
  ParamReader<Int, RealRed> reader("44matrixEx.dat");
  reader.getLines();
  /** To use a reader, you simply need, once you've called getLines(), to have
   * a variable in which to store the object you read, and to know where to get
   * it. In a few cases, it might be necessary to pass a pointer to the variable.
   * Also, in some functions, such as this one, it might not be possible to pass
   * the line as a temporary variable at the function call because it will be
   * incremented during execution.
   * */
  // Reading a matrix as an examples
  IntMat matrix(size, size); // The "recipient"
  unsigned int ln = 0; // The line counter
  reader.readBMat(matrix, ln, 0, size);

  // Creation of a writer object for file IOExample.out
  WriterRes<Int> writer("IOExample.out");

  // Writing to the file
  writer.writeString("We can write messages in the file.");
  writer.newLine();
  writer.writeString("Here is one way to write a matrix.");
  writer.newLine();
  writer.writeMMat(matrix);

  // Adding indentation
  writer.beginTabbedSection();
  writer.newLine();
  writer.writeString("Here is another one.");
  writer.newLine();
  for (int i = 0; i < matrix.NumRows(); i++) {
    writer.writeString("[");
    for (int j = 0; j < matrix.NumCols(); j++) {
      writer.writeIntScal(matrix[i][j]);
      if (j == matrix.NumCols()-1) continue;
      writer.writeString(" ");
    }
    writer.writeString("]");
    writer.newLine();
  }
  // Printing \n and ending indented section.
  writer.newParagraph();

  writer.writeString("It's also possible to use *ToString() methods to write"
      " stuff in the output\nfile.");
  writer.newLine();
  //! Creating an IntLatticeBase with the matrix we just read as a basis (it is
  //! non-singular)
  //! It is also possible to specify give a modulo and an dual basis when creating
  //! an IntLattice, but we will not use them here.
  IntLatticeBase<Int, Real, RealRed> lat_basis(matrix, size);
  //! We can evaluate and print the lenght of the vectors of this basis
  lat_basis.updateVecNorm();
  lat_basis.sort(0);
  // Writing in the file
  writer.writeString("The matrix and its vectors norms:");
  writer.newLine();
  writer.writeString(lat_basis.toStringBasis());
  writer.newLine();

  writer.writeString("Writers can be used as the program progresses, just like"
      "using a standard\noutput print, but also can be created at the end of the"
      "execution to print the\nresult of the execution, making them quite"
      "flexible.");

  return 0;
}
