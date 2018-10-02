/**
 * For now, this example benchmarks the different implementations we have
 * available for basis reduction.
 *
 * Add an error counter to see if both methods work.
 *
 * As of now the results look like this:
 * ma_basis: 21739192
 * cou_basis: 249924718
 * ma_basis/cou_basis: 0.086983
 * Dimension: 5
 * ma_basis: 300
 * cou_basis: 764
 * ma_basis/cou_basis: 0.39267
 * Dimension: 10
 * ma_basis: 1892
 * cou_basis: 3883
 * ma_basis/cou_basis: 0.487252
 * Dimension: 15
 * ma_basis: 7809
 * cou_basis: 18444
 * ma_basis/cou_basis: 0.42339
 * Dimension: 20
 * ma_basis: 48107
 * cou_basis: 186311
 * ma_basis/cou_basis: 0.258208
 * Dimension: 25
 * ma_basis: 70023
 * cou_basis: 329619
 * ma_basis/cou_basis: 0.212436
 * Dimension: 30
 * ma_basis: 449102
 * cou_basis: 6603505
 * ma_basis/cou_basis: 0.0680096
 * Dimension: 35
 * ma_basis: 21162057
 * cou_basis: 242782302
 * ma_basis/cou_basis: 0.0871647
 * ma_dual: 179732000
 * cou_dual: 5370
 * ma_dual/cou_dual: 33469.6
 *
 * This means that 1) Obviously the dual basis construction from solving a
 * linear system is crap. Doing the Euclid algorithm on the lines of the matrix
 * is about 10 times faster for the instances studied here but we should verify
 * if there is a difference in this ratio depending on the size of the instance.
 * Also, it is important to note that basis construction takes way less memory
 * in my method because it is done in place.
 *
 * Although it is not definitive, the triangularization seems faster in my
 * implementation, especially in larger dimensions.
 *
 * Aussi, clairement la construction de base par LLL est juste beaucoup trop
 * mieux que la triangularisation justement parce que les vecteurs n'escaladent
 * pas vers des tailles démoniaques.
 * */
#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"

using namespace LatticeTester;

int main() {
  clock_t ma_basis = 0, cou_basis = 0, ma_dual = 0, cou_dual = 0, LLL_basis = 0,
          LLL_dual = 0, tmp;
  clock_t basis_dim_ma[7], basis_dim_cou[7], basis_dim_LLL[7];
  for (int i = 0; i < 7; i++) {
    basis_dim_ma[i] = 0;
    basis_dim_cou[i] = 0;
    basis_dim_LLL[i] = 0;
  }
  std::string primes[6] = {"1021", "1048573", "1073741827", "1099511627791",
    "1125899906842597", "18446744073709551629"};
  ParamReader<MScal, BScal, RScal> reader;
  //long err_ma = 0, err_cou = 0;

  BasisConstruction<BScal> constr;
  BScal m;
  BMat matrix1, matrix2, matrix3, matrix4;
  int numlines;
  unsigned int ln;
  std::string name;
  for (int i = 0; i < 1; i++) {
    for (int j = 5; j < 36; j+=5) {
      for (int k = 0; k < 10; k++) {
        // Reader shenanigans
        name = "bench/" + primes[i] + "_" + std::to_string(j) + "_" + std::to_string(k);
        reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
        reader.getLines();
        reader.readInt(numlines, 0, 0);
        matrix1.kill();
        matrix2.kill();
        matrix3.kill();
        matrix4.kill();
        matrix1.SetDims(numlines, numlines);
        matrix2.SetDims(numlines, numlines);
        matrix3.SetDims(numlines, numlines);
        matrix4.SetDims(numlines, numlines);
        ln = 1;
        reader.readBMat(matrix4, ln, 0, numlines);
        matrix1 = matrix2 = matrix3 = matrix4;
        if (NTL::determinant(matrix1) == 0) continue;

        // Timing ma first
        tmp = clock();
        m = 1;
        constr.GCDConstruction(matrix1);
        ma_basis += clock() - tmp;
        basis_dim_ma[j/5-1] += clock()-tmp;
        tmp = clock();
        constr.DualConstruction(matrix1, matrix4, m);
        ma_dual += clock() - tmp;
        matrix1.kill();
        matrix1.SetDims(numlines, numlines);

        // Timing cou
        tmp = clock();
        m = 1;
        Triangularization(matrix2, matrix1, numlines, numlines, m);
        cou_basis += clock() - tmp;
        basis_dim_cou[j/5-1] += clock()-tmp;
        tmp = clock();
        constr.DualConstruction2(matrix1, matrix4, m);
        cou_dual += clock() - tmp;

        // Timing LLL
        tmp = clock();
        m = 1;
        constr.LLLConstruction(matrix3);
        LLL_basis += clock() - tmp;
        basis_dim_LLL[j/5-1] += clock()-tmp;
        tmp = clock();
        constr.DualConstruction2(matrix3, matrix4, m);
        LLL_dual += clock() - tmp;
      }
    }
  }
  std::cout << "ma_basis: " << ma_basis << std::endl;
  std::cout << "cou_basis: " << cou_basis << std::endl;
  std::cout << "ma_basis/cou_basis: " << (double)ma_basis/cou_basis << std::endl;
  for (int i = 0; i < 7; i++) {
    std::cout << "Dimension: " << (i+1)*5 << std::endl;
    std::cout << "ma_basis: " << basis_dim_ma[i] << std::endl;
    std::cout << "cou_basis: " << basis_dim_cou[i] << std::endl;
    std::cout << "ma_basis/cou_basis: " << (double)basis_dim_ma[i]/basis_dim_cou[i] << std::endl;
  }
  std::cout << "ma_dual: " << ma_dual << std::endl;
  std::cout << "cou_dual: " << cou_dual << std::endl;
  std::cout << "ma_dual/cou_dual: " << (double)ma_dual/cou_dual << std::endl;

  return 0;
}
