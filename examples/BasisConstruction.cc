/**
 * This example showcases the usage of the BasisConstruction module. This reads
 * matrices from files and builds a basis and a dual for an IntLatticeBasis
 * object. The files this is set to use are in the `bench.zip` archive. To
 * execute the program, the archive should be unziped and the `bench` folder
 * should be put in the same directory from which the executable is called.
 * 
 * This prints the time used by each strategy, but for the dimensions that we
 * use, we can conclude that there is not one method better than the other.
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

using namespace LatticeTester;

int main() {
  // The different clocks we will use for benchmarking
  // We use ctime for implementation simplicity
  int max_dim = 7; // Actual max dim is 5*max_dim
  clock_t tmp;
  clock_t gcd_time[max_dim], lll_time[max_dim];
  for (int i = 0; i < max_dim; i++) {
    gcd_time[i] = 0;
    lll_time[i] = 0;
  }

  std::string prime = "1021";
  ParamReader<MScal, BScal, RScal> reader;

  // The constructor we will use
  BasisConstruction<BScal> constr;
  BMat bas_mat, dua_mat;
  int numlines;
  unsigned int ln;
  std::string name;
  // This loop builds the basis with GCDConstruction straight up.
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      std::cout << name << std::endl;
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      bas_mat.kill();
      bas_mat.SetDims(numlines, numlines);
      dua_mat.kill();
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
      constr.GCDConstruction(bas_mat);
      // If you don't need the dual basis, the following line is sufficient
      // basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, numlines);
      MScal modulo(1);

      constr.DualConstruction(bas_mat, dua_mat, modulo);
      IntLatticeBasis <MScal, BScal, NScal, RScal> basis =
        IntLatticeBasis<MScal, BScal, NScal, RScal>(bas_mat, dua_mat, modulo, numlines);
      gcd_time[j] += clock() - tmp;
    }
  }
  // This loop builds the basis with LLLConstruction first.
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      bas_mat.kill();
      bas_mat.SetDims(numlines, numlines);
      dua_mat.kill();
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
      constr.GCDConstruction(bas_mat);
      // If you don't need the dual basis, the following line is sufficient
      // basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, numlines);
      MScal modulo(1);

      constr.DualConstruction(bas_mat, dua_mat, modulo);
      IntLatticeBasis <MScal, BScal, NScal, RScal> basis =
        IntLatticeBasis<MScal, BScal, NScal, RScal>(bas_mat, dua_mat, modulo, numlines);
      lll_time[j] += clock() - tmp;
    }
  }
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += gcd_time[i];
  }
  int width1 = log10(tmp)+2;
  std::cout << "           GCD";
  for (int i = 0; i<width1-3; i++) {
    std::cout << " ";
  }
  std::cout << "LLL\n";
  std::cout << "Total time " << tmp << " ";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += lll_time[i];
  }
  int width2 = log10(tmp)+2;
  std::cout << tmp << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 << std::setw(width1)
      << gcd_time[i] << std::setw(width2) << lll_time[i] << std::endl;
  }

  return 0;
}
