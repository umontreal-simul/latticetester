
/**An example of program to compare the speed of two upper
 *triangularisation method in LatticeTester. The first is
 *Util::Triangularization @author Couture.
 +The second is BasisConstruction:: upperTriangular @author Lecuyer
 * We use 150 basis. We begin with 5x5 dimension
 *to 75x75 dimension. 10 different basis for each dimension
 *We compute triangularization time in time clocks machine.
 *The for each  method an two dimension array of 15 line and 10 column is used
  to save the triangularization time.
  *---The line 1 of of the array contain computation time of the 10 basis with dimension   5x5,
   file used /examples/bench/1073741827_5_0.dat to /examples/bench/1073741827_5_9.dat
  *--The line 2 of of the array contain computation time of the 10 basis with dimension 10x10,
   file used /examples/bench/1073741827_10_0.dat to /examples/bench/1073741827_10_9.dat
  *---The line 3 of of the array contain computation time of the 10 basis with dimension 15x15,
   file used /examples/bench/1073741827_15_0.dat to /examples/bench/1073741827_15_9.dat
  *---The line 4 of of the array contain computation time of the 10 basis with dimension 20x20,
   file used /examples/bench/1073741827_20_0.dat to /examples/bench/1073741827_20_9.dat
  *---The line 5 of of the array contain computation time of the 10 basis with dimension 25x25,
   file used /examples/bench/1073741827_25_0.dat to /examples/bench/1073741827_25_9.dat
  *---The line 6 of of the array contain computation time of the 10 basis with dimension 30x30,
   file used /examples/bench/1073741827_30_0.dat to /examples/bench/1073741827_30_9.dat
  *---The line 7 of of the array contain computation time of the 10 basis with dimension 35x35,
   file used /examples/bench/1073741827_35_0.dat to /examples/bench/1073741827_35_9.dat
  *---The line 8 of of the array contain computation time of the 10 basis with dimension 40x40,
   file used /examples/bench/1073741827_40_0.dat to /examples/bench/1073741827_40_9.dat
  *---The line 9 of of the array contain computation time of the 10 basis with dimension 45x45,
   file used /examples/bench/1073741827_45_0.dat to /examples/bench/1073741827_45_9.dat
  *---The line 10 of of the array contain computation time of the 10 basis with dimension 50x50,
   file used /examples/bench/1073741827_50_0.dat to /examples/bench/1073741827_50_9.dat
  *---The line 11 of of the array contain computation time of the 10 basis with dimension 55x55,
   file used /examples/bench/1073741827_55_0.dat to /examples/bench/1073741827_55_9.dat
  *---The line 12 of of the array contain computation time of the 10 basis with dimension 60x60,
   file used /examples/bench/1073741827_60_0.dat to /examples/bench/1073741827_60_9.dat
  *---The line 13 of of the array contain computation time of the 10 basis with dimension 65x65,
   file used /examples/bench/1073741827_65_0.dat to /examples/bench/1073741827_65_9.dat
  *---The line 14 of of the array contain computation time of the 10 basis with dimension 70x70,
   file used /examples/bench/1073741827_70_0.dat to /examples/bench/1073741827_70_9.dat
  *---The line 15 of of the array contain computation time of the 10 basis with dimension 75x75,
   file used /examples/bench/1073741827_75_0.dat to /examples/bench/1073741827_75_9.dat
 **/

#define NTL_TYPES_CODE 2
#include <iostream>
#include <ctime>
#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Reducer.h"
#include "latticetester/Const.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"
#include "latticetester/NTLWrap.h"
#include "Examples.h"

using namespace LatticeTester;
namespace
{
  const std::string prime = primes[2];
}

int main()
{
  // clock_t timer = clock();
  clock_t tmp;

  RealMat matTriangularTime;
  matTriangularTime.resize(15, 10);

  std::string prime = primes[2];
  double tps = 0;
  //

  for (int j = 0; j < 15; j++)
  {
    for (int k = 0; k < 10; k++)
    {
      //! Variables definition
      ParamReader<Int, RealRed> reader;
      BasisConstruction<Int> constr; // The basis constructor we will use
      std::string name;
      int numlines;
      IntMat matrix1, matrix2;
      unsigned int ln;
      Int m(1021);

      name = "bench/" + prime + "_" + std::to_string(5 * (j + 1)) + "_" + std::to_string(k);
      //std::cout << name << std::endl;
    
      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);

      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      tmp = clock();
      constr.upperTriangular(matrix1, matrix2, m);
      tps = (double)(clock() - tmp); //(CLOCKS_PER_SEC);
      matTriangularTime(j, k) = tps;
    }
  }

  std::cout << "Print the array that contain computation tine with the BasisConstruction::upperTriangular \n";
  for (int i = 0; i < 15; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      std::cout << matTriangularTime(i, j) << "   ";
    }
    std::cout << "" << std::endl;
  }

  // we repaeat the same code but this time we use Util::triangularization method

  for (int j = 0; j < 15; j++)
  {
    for (int k = 0; k < 10; k++)
    {

      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1, matrix2;
      unsigned int ln;
      // BasisConstruction<Int> constr;
      Int m(1021);
      name = "bench/" + prime + "_" + std::to_string(5 * (j + 1)) + "_" + std::to_string(k);
      //std::cout << name << std::endl;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);

      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      tmp = clock();

      Triangularization(matrix1, matrix2, numlines, numlines, m);

      double tps = (double)(clock() - tmp); //(CLOCKS_PER_SEC);
      matTriangularTime(j, k) = tps;
      // delete constr;
    }
  }

  std::cout << "Print the array that contain computation tine with the Util::Triangularization \n";
  for (int i = 0; i < 15; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      std::cout << matTriangularTime(i, j) << "   ";
    }
    std::cout << "" << std::endl;
  }

  return 0;
}
