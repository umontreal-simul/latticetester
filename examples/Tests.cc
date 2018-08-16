/**
 * This examples uses LatticeTester to perform tests on a basis.
 * There are two ways implemented in LatticeTester to perform tests on a 
 * lattice. You can either use the LatticeAnalysis class or the
 * LatticeTesterRoutines module. The prescribed usage for LatticeAnalysis is 
 * already presented in the lattester executable so this example will focus on
 * LatticeTesterRoutines.
 *
 * The LatticeTesterRoutines module is a high level API that that can be called
 * with a minimal number of parameters for each wanted test. Since those functions
 * do not keep track of the various parameters passed to them, they are best used
 * in small programs or scripts (as much as scripting is a thing in C++).
 *
 * */
#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/LatticeAnalysis.h"
#include "latticetester/Types.h"

using namespace LatticeTester;

int main() {
  // Creating a Reducer object. For details look at the other examples.
  ParamReader<MScal, BScal, RScal> reader("build/examples/44matrixEx.dat");
  reader.getLines();
  BMat matrix(4,4);
  unsigned int ln = 0;
  reader.readBMat(matrix, ln, 0, 4);
  IntLatticeBasis<MScal, BScal, NScal, RScal> basis(matrix, 4);
  Reducer<MScal, BScal, NScal, RScal> red(basis);

  // Creating an empty object. This is to showcase the setup from scratch.
  LatticeAnalysis<MScal, BScal, NScal, RScal> lattice;
  // Setting all the values needed to execute a test
  lattice.setReducer(red);
  lattice.setCriterion(SPECTRAL);
  lattice.setPreReduction(BKZ);
  lattice.setNorm(L2NORM);
  lattice.initNormalizer(ROGERS);
  lattice.doTest();
  lattice.printTestResults();
  // Here Reducer has been modified by doTest(). We can get the modified basis.
  std::cout << "The basis that we have:\n" << basis.toStringBasis();

  // Now, let's instanciate a whole new object with parameters in the constructor.
  basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, 4);
  // This instanciation lists all the possibilities. All the commented parameters
  // have somewhat sane defaults except normaType in the case of a spectral test.
  lattice = LatticeAnalysis<MScal, BScal, NScal, RScal>(red, BEYER
      /*, preRed, norm, maxNodesBB, normaType, alpha*/);
  lattice.doTest();
  lattice.printTestResults();
  // Here Reducer has been modified by doTest(). We can get the modified basis.
  std::cout << "The basis that we have:\n" << basis.toStringBasis();
  return 0;
}
