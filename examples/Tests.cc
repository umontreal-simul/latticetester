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

#include <iostream>

#include "latticetester/LatticeTesterRoutines.h"

using namespace LatticeTester;

int main() {
  return 0;
}
