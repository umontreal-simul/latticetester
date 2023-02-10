/*An example programm to compute a basis shortest vector using the
**Cholesky decomposition  and displays the computation time in minute. 
* we also compute the shortest vector using the triangular decomposition
* and displays the computation time in minute.
* In code, we use a basis have 4x4 dimension, the modulo 'm' is 101.
* The name of dat file is '1021_4_001.dat'. It is the 'examples' folder of LatticeTester.
* The absolue path is '/BuildPlace/examples/bench/1021_4_002'.
*/

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
#include "Examples.h"
#include "latticetester/WriterRes.h"

using namespace LatticeTester;
namespace
{
  const std::string prime = primes[0];
}

int main()
{
  std::string prime = primes[0];
  clock_t timer, tmp;
  IntLatticeBase<Int, Real, RealRed> *basis;
  Reducer<Int, Real, RealRed> *red;
  // BasisConstruction<Int> constr; // The basis constructor we will use
  //! Variables definition
  ParamReader<Int, RealRed> reader, reader2;
  std::string name;
  int numlines;
  // The primal basis
  IntMat matrix1;
  // The dual basis
  IntMat matrix2;
  unsigned int ln;

  std::string s1("cholesky");
  std::string s2("triangular");
  name = "bench/" + prime + "_4" + "_002";
  reader = ParamReader<Int, RealRed>(name + ".dat");
  reader.getLines();
  reader.readInt(numlines, 0, 0);
  //For the primal basis
  matrix1.SetDims(numlines, numlines);
  //for the dual basis
  matrix2.SetDims(numlines, numlines);
  ln = 1;
  reader.readBMat(matrix1, ln, 0, numlines);


  Int m(101);
 // basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, matrix2, m, numlines);
  basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, numlines);
  red = new Reducer<Int, Real, RealRed>(*basis);
  std::cout << " The base before reduction\n";
  printBase((red->getIntLatticeBase())->getBasis());
 
  // LLL reduction 
  red->redLLL(0.999, 1000000, numlines);
  // We can also the BKZ reduction by removing the comment of the next line
    //red->redBKZ();
 
  basis->updateVecNorm();

  std::cout << " The base after reduction\n";
  printBase((red->getIntLatticeBase())->getBasis());

 //  comput the shortest vector using Cholesky decompisition method
 //This is the recommanded method
  timer = clock();
  if (!red->shortestVector(L2NORM, s1))
  {
    std::cout << " The computation of shortest vector failed \n";
  }
  tmp=clock()-timer;
  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  std::cout << red->getMinLength() << std::endl;
  std::cout << "The time to compute shortest vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;


 //  comput the shortest vector using Cholesky decompisition method
 //The method can be very slow. It's not the recommanded one

  basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, matrix2, m, numlines);
  std::cout << " The base before reduction\n";
  printBase((red->getIntLatticeBase())->getBasis());
  red = new Reducer<Int, Real, RealRed>(*basis);
  red->redBKZ(0.9, 10, QUADRUPLE, 0);
  std::cout << " The base after reduction\n";
  printBase((red->getIntLatticeBase())->getBasis());
  timer = clock();
  if (!red->shortestVector(L2NORM, s2))
  {
    std::cout << " The computation of shortest vector failed \n";
  }
  tmp=clock()-timer;
  std::cout << " The shortest vector length wih triangular decomposition\n";
  std::cout << red->getMinLength() << std::endl;
  std::cout << " The time to compute shortest vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;


  return 0;
}
