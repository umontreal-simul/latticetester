/**
 * For now, this example benchmarks the different implementations we have
 * available for basis reduction.
 *
 * Add an error counter to see if both methods work.
 * */
#define NTL_TYPES_CODE 1

#include <iostream>
#include <ctime>

#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"

using namespace LatticeTester;

int main() {
  clock_t ma_basis = 0, cou_basis = 0, ma_dual = 0, cou_dual = 0, tmp;
  std::string primes[6] = {"1021", "1048573", "1073741827", "1099511627791",
    "1125899906842597", "18446744073709551629"};
  ParamReader<MScal, BScal, RScal> reader;
  //long err_ma = 0, err_cou = 0;

  BasisConstruction<BScal> constr;
  BScal m;
  BMat matrix1, matrix2, matrix3;
  int numlines;
  unsigned int ln;
  std::string name;
  for (int i = 0; i < 6; i++) {
    for (int j = 5; j < 101; j+=5) {
      for (int k = 0; k < 10; k++) {
        // Reader shenanigans
        name = "bench/" + primes[i] + "_" + std::to_string(j) + "_" + std::to_string(k);
        reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
        reader.getLines();
        reader.readInt(numlines, 0, 0);
        matrix1.SetDims(numlines, numlines);
        matrix2.SetDims(numlines, numlines);
        matrix3.SetDims(numlines, numlines);
        ln = 1;
        reader.readBMat(matrix3, ln, 0, numlines);
        matrix1 = matrix2 = matrix3;

        // Timing ma first
        tmp = clock();
        constr.GCDConstruction(matrix1);
        ma_basis += clock() - tmp;
        tmp = clock();
        constr.DualConstruction(matrix1, BScal(1));
        ma_dual += clock() - tmp;

        // Timing cou
        tmp = clock();
        m = 1;
        Triangularization(matrix2, matrix1, numlines, numlines, m);
        cou_basis += clock() - tmp;
        tmp = clock();
        CalcDual(matrix1, matrix2, numlines, m);
        cou_dual += clock() - tmp;
      }
    }
  }
  std::cout << "ma_basis: " << ma_basis << std::endl;
  std::cout << "cou_basis: " << cou_basis << std::endl;
  std::cout << "ma_basis/cou_basis: " << (double)ma_basis/cou_basis << std::endl;
  std::cout << "ma_dual: " << ma_dual << std::endl;
  std::cout << "cou_dual: " << cou_dual << std::endl;
  std::cout << "ma_dual/cou_dual: " << (double)ma_dual/cou_dual << std::endl;

  return 0;
}
