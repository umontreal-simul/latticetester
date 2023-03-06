
/**In this example, we show the use of most of the method in class BasisConstruction.
 * We begin by reading a file that contains the basis information and create an
 * 'IntLattice' object. 

 * We show a use of BasisContruction::upperTriangular, BasisContruction::lowerTriangular,
 * BasisContruction::LLLConstruction with two differrent parametter of delta,
 * BasisContruction::calcDualUpperTriangular, BasisContruction::calcDual
 *
 * In this example, we can compare the speed of  BasisConstruction::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and BasisConstruction::calcDualUpperTriangular method which compute an m-dual basis
 * with an upper triangular basis.
 *
 * We can also compare the speed of 'BasisConstruction::upperTriangular'
 * and the speed of 'BasisConstruction::LLLConstruction'
 * *
 **/

 /*
 *  The bases we used to test the functionnaly of LatticeTester are in a folder named 
 * 'examples/bench/'.
 *  Each file in 'examples/bench/' folder contain a basis, and the file is nameed as follows:
 * 'prime_dimBasis_exanpleNumber' where 'prime' is modulo value of the basis, 
 * 'dimBasis' is the dimension of the basis, and 'exampleNumber' is the number of the
 *  example for the bases of dimension 'dimBasis'.
 *   // This is a sample output for TYPES_CODE ZD:
 *  Total time  UPP    Low     LLL1   LLL2   Dual1 Dual2
 * Dim     5   2049   2257    466    361    878   223
 * Dim    10   7396   8419   1814   2089   2064   801
 * Dim    15  18459  21048   4122   5983   5256  2119
 * Dim    20  41071  47818   8670  13366  11358  4374
 * Dim    25  72438  80001  15633  32268  18913  7216
 * Dim    30 110689 134219  20085  40633  35397 11837
 * Dim    35 164154 186050  27150  58239  52095 18252
 * Dim    40 232165 267259  31969  70082  79393 25671
 * Total time: 0.0569002 minutes
*/



#define TYPES_CODE  ZD
//#define NTL_TYPES_CODE  2

#include <iostream>
#include <ctime>
#include <NTL/mat_GF2.h>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/EnumTypes.h"

using namespace LatticeTester;

int getWidth(clock_t time[], int dim, std::string message, clock_t totals[], int ind) {
  clock_t tmp = 0;
  for (int i = 0; i < dim; i++) {
    tmp += time[i];
  }
  int width = log10(tmp) + 2;
 // std::cout << std::setw(width) << message;
  totals[ind] = tmp;
  return width;
}

 //  The following array gives the possible modulo values for the basis examples.
 //  The nodulo are all prime values.
const int many_primes = 6;
const std::string primes[] = {"1021", "1048573", "1073741827", "1099511627791",
                  "1125899906842597", "18446744073709551629"};

 //use file use of basis values modulo 1021
 const std::string prime = primes[0];

 int main()
 {
  clock_t timer = clock();
  // The different clocks we will use for benchmarking
  // We use ctime for implementation simplicity
  int leng = 8; // Actual max dim is 5*leng
  clock_t upp_time[leng],low_time[leng], lll1_time[leng],lll2_time[leng],
  dual1_time[leng], dual2_time[leng], totals[6];
  for (int i = 0; i < leng; i++) {
    upp_time[i] = 0;
    low_time[i] = 0;
    lll1_time[i] = 0;
    lll2_time[i] = 0;
    dual1_time[i] = 0;
    dual2_time[i] = 0;
  }

  IntLattice<Int, Real> *lattice;
  BasisConstruction<Int> constr; // The basis constructor we will use
  IntMat bas_mat, dua_mat;
  NTL::Mat<NTL::ZZ> w_copie2, m_dual2;
  IntMat w_copie, m_v, m_v2,m_dual;
  Int m(1021);
  
  clock_t tmp;
  
 for (int j = 0; j < leng; j++) {
   for (int k = 0; k < 10; k++) {
    std::string name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
  
  ParamReader<Int, Real> reader(name + ".dat");

  reader.getLines();
  int numlines;
  unsigned int ln=1;
  reader.readInt(numlines, 0, 0);

  bas_mat.SetDims(numlines, numlines);
  dua_mat.SetDims(numlines, numlines);
  // for the upper triangular  basis
  m_v.SetDims(numlines, numlines);
  // for the lower triangular  basis
  m_v2.SetDims(numlines, numlines);
  // for the m-dual  basis
  m_dual.SetDims(numlines, numlines);
  // A copie of lattice basis
  w_copie.SetDims(numlines, numlines);
  m_dual2.SetDims(numlines, numlines);
  // A copie of lattice basis
  w_copie2.SetDims(numlines, numlines);
  //! Filling the matrix
  reader.readBMat(bas_mat, ln, 0, numlines);
  // Creating a lattice basis
  // lattice = new IntLattice<Int, Real>(bas_mat, bas_mat, m, numlines);
  lattice = new IntLattice<Int, Real>(bas_mat, numlines);
  copy(bas_mat, w_copie);
  tmp = clock();
  constr.upperTriangularBasis(w_copie, m_v, m);
  upp_time[j] += clock() - tmp;
  copy(bas_mat, w_copie);
  // The construction of the lower triangular basis
  tmp = clock();
  constr.lowerTriangularBasis(w_copie, m_v2, m);
  low_time[j] += clock() - tmp;

 // std::cout << " The LLL reduction basis with delta=0.8 \n";
  copy(bas_mat, w_copie);
  tmp = clock();
 constr.LLLConstruction(w_copie, 0.8);
  lll1_time[j] += clock() - tmp;

 // std::cout << " The LLL reduction basis with delta=0.99 \n";
  copy(bas_mat, w_copie);
  tmp = clock();
  constr.LLLConstruction(w_copie, 0.99999);
  lll2_time[j] += clock() - tmp;

  // NTL::ZZ mm(1021);
  // copyMatrixToMat(bas_mat, w_copie2); 
  copy(bas_mat, w_copie);
  tmp = clock();
  constr.mDualBasis(w_copie, m_dual, m);
  dual1_time[j] += clock() - tmp;
  
  copy(bas_mat, w_copie);
  constr.upperTriangularBasis(w_copie, m_v, m);
  tmp = clock();
  constr.mDualUpperTriangular(m_v, m_v2, numlines, m);
  dual2_time[j] += clock() - tmp;
    }
  }

  std::cout << "         ";
  int width1 = getWidth(upp_time, leng, "UPPTR", totals, 0);
  int width2 = getWidth(low_time, leng, "LOWTR", totals, 1);
  int width3 = getWidth(lll1_time, leng, "LLL1", totals, 2);
  int width4 = getWidth(lll2_time, leng, "LLL2", totals, 3);
  int width5 = getWidth(dual1_time, leng, "DUAL1", totals, 4);
  int width6 = getWidth(dual2_time, leng, "DUAL2", totals, 5);
  std::cout << std::endl;

  std::cout << "Total time" << std::setw(width1) << totals[0]
    << std::setw(width2) << totals[1]
    << std::setw(width3) << totals[2]
    << std::setw(width4) << totals[3]
    << std::setw(width5) << totals[4]
    << std::setw(width6) << totals[5] << std::endl;
  for (int i = 0; i < leng; i++) {
    std::cout << "Dim" << std::setw(6) << (i+1)*5
      << std::setw(width1) << upp_time[i] << std::setw(width2) << low_time[i]
      << std::setw(width3) << lll1_time[i] << std::setw(width4) << lll2_time[i]
      << std::setw(width5) << dual1_time[i] << std::setw(width6) << dual2_time[i];
    std::cout << std::endl;
  }
  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";


  return 0;
}