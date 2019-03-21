/**
 * This example showcases the usage of the BasisConstruction module. This reads
 * matrices from files and builds a basis and a dual for an `IntLatticeBasis`
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
#include "latticetester/IntLattice.h"

#include "Examples.h"

using namespace LatticeTester;

namespace {
  const std::string prime = primes[0];
  // The function that tests one of the methods in BasisConstruction.
  // That note that this is a functionnal programming approach
  void testFunction(void (BasisConstruction<BScal>::*func)(BMat&),
      BasisConstruction<BScal>& constr, clock_t time1[], clock_t time2[],
      int max_dim){
    clock_t tmp;
    for (int j = 0; j < max_dim; j++) {
      for (int k = 0; k < 10; k++) {
        tmp = clock();
        // Reader shenanigans
        std::string name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
        ParamReader<MScal, BScal, RScal> reader(name + ".dat");
        reader.getLines();
        int numlines;
        unsigned int ln;
        reader.readInt(numlines, 0, 0);
        BMat bas_mat, dua_mat;
        bas_mat.SetDims(numlines, numlines);
        dua_mat.SetDims(numlines, numlines);
        ln = 1;
        reader.readBMat(bas_mat, ln, 0, numlines);

        // We want to avoid singular matrix because we can't compute the dual, and
        // IntLatticeBasis only really supports square matrices.
        if (NTL::determinant(bas_mat) == 0) {
          std::cout << name << " is singular\n";
          continue;
        }

        // Timing ma first
        (constr.*func)(bas_mat);
        // If you don't need the dual basis, the following line is sufficient
        // basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, numlines);
        MScal modulo(1);
        time1[j] += clock() - tmp;

        tmp = clock();
        constr.DualConstruction(bas_mat, dua_mat, modulo);
        time2[j] += clock() - tmp;
      }
    }
  }
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

  ParamReader<MScal, BScal, RScal> reader;

  // Defining constants for the execution of the algorithms
  BasisConstruction<BScal> constr; // The basis constructor we will use
  BMat bas_mat, dua_mat;

  testFunction(&BasisConstruction<BScal>::GCDConstruction, constr, gcd_time,
      dual1_time, max_dim);
  testFunction(&BasisConstruction<BScal>::LLLConstruction, constr, lll_time,
      dual2_time, max_dim);

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
