/**
 * This examples uses the different forms of reduction available in LatticeTester
 * to reduce the same basis and then prints the execution time and the reduced
 * basis.
 *
 * This example can also showcase the usage of IntLatticeBasis.
 * */
#define NTL_TYPES_CODE 2

#include <iostream>
#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;

int main() {
  BasisConstruction<BScal> constr;

  BMat matrix1;
  matrix1.SetDims(4,3);
  matrix1[0][0]=10; matrix1[0][1]=5; matrix1[0][2]=5;
  matrix1[1][0]=15; matrix1[1][1]=4; matrix1[1][2]=5;
  matrix1[2][0]=20; matrix1[2][1]=6; matrix1[2][2]=5;
  matrix1[3][0]=30; matrix1[3][1]=7; matrix1[3][2]=3;
  std::cout << "DualConstruction\n";
  // constr.GCDConstruction(matrix1);
  constr.DualConstruction(matrix1);

  matrix1.SetDims(4,5);
  matrix1[0][0]=10; matrix1[0][1]=5; matrix1[0][2]=5; matrix1[0][3]=2; matrix1[0][4]=3;
  matrix1[1][0]=15; matrix1[1][1]=4; matrix1[1][2]=5; matrix1[1][3]=2; matrix1[1][4]=3;
  matrix1[2][0]=20; matrix1[2][1]=6; matrix1[2][2]=5; matrix1[2][3]=2; matrix1[2][4]=3;
  matrix1[3][0]=20; matrix1[3][1]=6; matrix1[3][2]=5; matrix1[3][3]=2; matrix1[3][4]=3;
  constr.GCDConstruction(matrix1);
  return 0;
}
