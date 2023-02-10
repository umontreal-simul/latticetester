/**
 * This example showcases the usage of the BasisConstruction module. This reads
 * matrices from files and builds a basis and a dual for an `IntLatticeBase`
 * object. The files this is set to use are in the `bench.zip` archive. To
 * execute the program, the archive should be unziped and the `bench` folder
 * should be put in the same directory from which the executable is called.
 *
 * This example reads matrices from files and performs the different construction
 * algorithms in BasisConstruction on them. The program then prints the execution
 * time of the various algorithms. Note that the execution of the program is not
 * what you would expect in reality since bench contains random full matrices.
 * 
 * This is a sample output for NTL_TYPES_CODE 2:
 *                GCD    LLL     DUAL1    DUAL2
 * Dim     5     4418   3074       735     1002
 * Dim    10    13497   7900      2647     8151
 * Dim    15    38502  20984      9543    19052
 * Dim    20    94467  44949     88171    50834
 * Dim    25   152712  86751    154730   181654
 * Dim    30   594683 137168   2970433  1682890
 * Dim    35 21994254 221505 168412442 13860037
 * */

// This should always use Types 2 or 3, because we get too big numbers with GCD
// elimination.
#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/Types.h" 
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLatticeBase.h"

#include "Examples.h"

using namespace LatticeTester;

namespace {
  const std::string prime = primes[0];
}

int main() {
  clock_t timer = clock();
  // The different clocks we will use for benchmarking
  // We use ctime for implementation simplicity
  int max_dim = 6; // Actual max dim is 5*max_dim
  clock_t gcd_time[max_dim], lll_time[max_dim],
  dual1_time[max_dim], dual2_time[max_dim], totals[4];
  for (int i = 0; i < max_dim; i++) {
    gcd_time[i] = 0;
    lll_time[i] = 0;
    dual1_time[i] = 0;
    dual2_time[i] = 0;
  }

  // Defining constants for the execution of the algorithms
  BasisConstruction<Int> constr; // The basis constructor we will use
  IntMat bas_mat, dua_mat;
  Int mod(1021);

  clock_t tmp;
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      //! Reader shenanigans
      std::string name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      ParamReader<Int, RealRed> reader(name + ".dat");
      reader.getLines();
      int numlines;
      unsigned int ln;
      reader.readInt(numlines, 0, 0);
      IntMat bas_mat, dua_mat;
      bas_mat.SetDims(numlines, numlines);
      dua_mat.SetDims(numlines, numlines);
      ln = 1;
      //! Filling the matrix
      reader.readBMat(bas_mat, ln, 0, numlines);

      // Creating a lattice basis
      IntLatticeBase<Int, Real, RealRed> lattice(bas_mat, numlines);

      //! We want to avoid singular matrix because we can't compute the dual, and
      //! IntLatticeBase only really supports square matrices.
      if (NTL::determinant(bas_mat) == 0) {
        std::cout << name << " is singular\n";
        continue;
      }

      // Timing GCDConstruction first
      tmp = clock();
      constr.GCDTriangularBasis(bas_mat,mod);
      Int modulo(1);
      gcd_time[j] += clock() - tmp;

      // Timing mDualTriangular
      tmp = clock();
      constr.mDualTriangular(bas_mat, dua_mat, modulo);
      dual1_time[j] += clock() - tmp;

      // Timing LLLConstruction next
      tmp = clock();
      constr.LLLConstruction(lattice.getBasis());
      modulo = Int(1);
      lll_time[j] += clock() - tmp;

      // The following works, but does not set all the properties of lattice to
      // properly work with a dual.
      tmp = clock();
      constr.mDualTriangular(lattice.getBasis(), lattice.getDualBasis(), modulo);
      dual2_time[j] += clock() - tmp;
      // This sets the lattice to know it has a dual. Computing the norm of the
      // vectors in the lattice would also be wise.
      lattice.setDualFlag(true);
    }
  }

  std::cout << "         ";
  int width1 = getWidth(gcd_time, max_dim, "GCD", totals, 0);
  int width2 = getWidth(lll_time, max_dim, "LLL", totals, 1);
  int width3 = getWidth(dual1_time, max_dim, "DUAL1", totals, 2);
  int width4 = getWidth(dual2_time, max_dim, "DUAL2", totals, 3);
  std::cout << std::endl;

  std::cout << "Total time" << std::setw(width1) << totals[0]
    << std::setw(width2) << totals[1]
    << std::setw(width3) << totals[2]
    << std::setw(width4) << totals[3] << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim" << std::setw(6) << (i+1)*5
      << std::setw(width1) << gcd_time[i] << std::setw(width2) << lll_time[i]
      << std::setw(width3) << dual1_time[i] << std::setw(width4) << dual2_time[i];
    std::cout << std::endl;
  }
  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";

  return 0;
}
