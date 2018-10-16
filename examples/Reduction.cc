/**
 * This is an example to the usage of the reducer class to perform reduction
 * of lattices. This also compares the average relative length reduction of each
 * method and the time taken by dimension.
 * */
/**
 * This example compares our implementation of the LLL algorithm with the one in
 * NTL. This comparison is done on a set of randomly generated matrices
 * available in the `bench.zip` archive in the example folder of the build tree.
 * Please note, before running this program that it takes quite a long time to
 * execute (in our test runs it took slightly more than an hour).
 *
 * Also, a lot of the information is hard coded in this example and we are aware
 * this is not a really good practice. This program is mostly a script we made
 * to have tangible proof that NTL is way faster than our algorithm.
 *
 * This program will do LLL reduction on the sets of matrix in the files in the
 * `bench.zip` archive and mesure the time it takes. The matrix where generated
 * with both different dimension and number sizes. For each of the combination
 * of dimension we consider, 10 different matrices are tested for the sake of 
 * consistency.
 *
 * Below is an example execution with NTL_TYPES_CODE 2 :
 *
 * ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS
 * Results depending on the size of the numbers
 *           Lower than            NTL             Us         NTL/Us
 *                 1021        3841541       73879949      0.0519971
 *              1048573        3892779       97965949       0.039736
 *           1073741827        7820678      256791162      0.0304554
 *        1099511627791        8234180      297252111       0.027701
 *     1125899906842597        9700107      348596562      0.0278262
 * 18446744073709551629       10460438      358252442      0.0291985
 * 
 * Results depending on the dimension
 *  Dimension            NTL             Us         NTL/Us
 *          5           4278          11313       0.378149
 *         10          21040         101709       0.206865
 *         15          42329         296306       0.142856
 *         20         100598         984737       0.102157
 *         25         196672        2429732      0.0809439
 *         30         302465        4725799      0.0640029
 *         35         460690        8367421      0.0550576
 *         40         658885       14334478      0.0459651
 *         45         829210       19395797       0.042752
 *         50        1134136       30368342       0.037346
 *         55        1394538       38832067       0.035912
 *         60        1869029       55709073      0.0335498
 *         65        2286715       68535529      0.0333654
 *         70        2632518       82066722      0.0320778
 *         75        3169892       98800234      0.0320839
 *         80        4040843      131201006      0.0307989
 *         85        4760273      158600524      0.0300142
 *         90        5791525      219756863      0.0263542
 *         95        6518928      238274269      0.0273589
 *        100        7734496      259945513      0.0297543
 * 
 * We took  1432735495 ticks
 * NTL took 43947377 ticks
 * We are globally 32.6012 times slower
 * */

// We define the numeric types.
// It is possible to use this example with TYPES 2 and 3. For now 1 calls the
// same function for both execution and we look forward to change that.
#define NTL_TYPES_CODE 1

#include <iostream>
#include <ctime>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/WriterRes.h"

using namespace LatticeTester;

int main() {
  // This is basically the C method of timing a program. We time globally, but
  // also for eache dimension and for each size of integers in the matrix.
  clock_t our_time = 0, ntl_time = 0, diff;
  clock_t our_time_num[6], ntl_time_num[6];
  for (int i = 0; i < 6; i++){
    our_time_num[i] = 0;
    ntl_time_num[i] = 0;
  }
  clock_t our_time_dim[20], ntl_time_dim[20];
  for (int i = 0; i < 20; i++){
    our_time_dim[i] = 0;
    ntl_time_dim[i] = 0;
  }

  // Variables declaration.
  std::string primes[6] = {"1021", "1048573", "1073741827", "1099511627791",
    "1125899906842597", "18446744073709551629"};
  ParamReader<MScal, BScal, RScal> reader;
  // We dynamically allocate memory to these two pointers every time we need to
  // create an object of their type. It is probably not the best way to do it,
  // but the time taken by that is negligible compared to the LLL execution.
  IntLatticeBasis<MScal, BScal, NScal, RScal>* basis;
  Reducer<MScal, BScal, NScal, RScal>* red;
  std::string name;
  int numlines;
  BMat matrix1, matrix2;
  unsigned int ln;

  for (int i = 0; i < 2; i++) {
    for (int j = 5; j < 101; j+=5) {
      for (int k = 0; k < 10; k++) {
        // Reader shenanigans
        name = "bench/" + primes[i] + "_" + std::to_string(j) + "_" + std::to_string(k);
        reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
        reader.getLines();
        reader.readInt(numlines, 0, 0);
        matrix1.SetDims(numlines, numlines);
        matrix2.SetDims(numlines, numlines);
        ln = 1;
        reader.readBMat(matrix1, ln, 0, numlines);
        matrix2 = matrix1;

        // We time NTL first
        basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix2, numlines);
        red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
        diff = clock();
        red->redLLLNTL();
        ntl_time += clock()-diff;
        ntl_time_dim[j/5-1] += clock()-diff;
        ntl_time_num[i] += clock()-diff;
        delete red;
        delete basis;

        // We time our implementation
        matrix2 = matrix1;
        basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix2, numlines);
        red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
        diff = clock();
        red->redLLL();
        our_time += clock()-diff;
        our_time_dim[j/5-1] += clock()-diff;
        our_time_num[i] += clock()-diff;
        delete red;
        delete basis;
      }
    }
  }

  // Printing the results in a somewhat formated way.
  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  std::cout << "Results depending on the size of the numbers\n";
  std::cout << std::setw(20) << "Lower than" << std::setw(15) << "NTL"
    << std::setw(15) << "Us" << std::setw(15) << "NTL/Us" << std::endl;
  for (int i = 0; i < 6; i++){
    std::cout << std::setw(20) << primes[i] << std::setw(15) << ntl_time_num[i]
      << std::setw(15) << our_time_num[i] << std::setw(15)
      << (double)ntl_time_num[i]/(double)our_time_num[i] << std::endl;
  }

  std::cout << "\nResults depending on the dimension\n";
  std::cout << std::setw(10) << "Dimension" << std::setw(15) << "NTL"
    << std::setw(15) << "Us" << std::setw(15) << "NTL/Us" << std::endl;
  for (int i = 0; i < 20; i++){
    std::cout << std::setw(10) << (i+1)*5 << std::setw(15) << ntl_time_dim[i]
      << std::setw(15) << our_time_dim[i] << std::setw(15)
      << (double)ntl_time_dim[i]/(double)our_time_dim[i] << std::endl;
  }

  std::cout << "\nWe took  " << our_time << " ticks\n";
  std::cout << "NTL took " << ntl_time << " ticks\n";
  std::cout << "We are globally " << (double)our_time/(double)ntl_time
    << " times slower\n";
  
  return 0;
}
