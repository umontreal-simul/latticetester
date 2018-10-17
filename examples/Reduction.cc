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
#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/WriterRes.h"

using namespace LatticeTester;

/*
 * On veut faire les réductions de Dieter, LLL, BKZ et de Minkowski pour
 * plusieurs bases, comparer le temps d'exécution et la réduction moyenne des
 * longueurs des vecteurs pour le fun.
 * */
int main() {
  int max_dim = 20; // Actual max dim is 5*max_dim
  // This is basically the C method of timing a program. We time globally, but
  // also for eache dimension and for each size of integers in the matrix.
  clock_t lll_time[max_dim], die_time[max_dim], bkz_time[max_dim],
  min_time[max_dim], tmp;
  clock_t total_times[4];
  for (int i = 0; i < max_dim; i++){
    lll_time[i] = 0;
    die_time[i] = 0;
    bkz_time[i] = 0;
    min_time[i] = 0;
  }

  // Variables declaration.
  std::string prime = "1021";
  ParamReader<MScal, BScal, RScal> reader;
  // We dynamically allocate memory to these two pointers every time we need to
  // create an object of their type. It is probably not the best way to do it,
  // but the time taken by that is negligible compared to the LLL execution.
  IntLatticeBasis<MScal, BScal, NScal, RScal>* basis;
  Reducer<MScal, BScal, NScal, RScal>* red;
  std::string name;
  int numlines;
  BMat matrix1;
  unsigned int ln;

  // Timing Dieter
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(j) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redDieter(0);
      delete red;
      delete basis;
      die_time[j] += clock() - tmp;
    }
  }
  // Timing LLL with default parameters
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(j) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redLLLNTL();
      delete red;
      delete basis;
      lll_time[j] += clock() - tmp;
    }
  }
  // Timing BKZ with default parameters
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(j) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redBKZ();
      delete red;
      delete basis;
      bkz_time[j] += clock() - tmp;
    }
  }
  // Timing Minkowski reduction
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(j) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->reductMinkowski(0);
      delete red;
      delete basis;
      min_time[j] += clock() - tmp;
    }
  }

  // Printing the results in a somewhat formated way.
  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += die_time[i];
  }
  int width1 = log10(tmp)+2;
  total_times[0] = tmp;
  std::cout << "          " << std::setw(width1) << "Dieter";
  for (int i = 0; i<width1-3; i++) {
    std::cout << " ";
  }
  std::cout << "Total time " << tmp << " ";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += lll_time[i];
  }
  int width2 = log10(tmp)+2;
  total_times[1] = tmp;
  std::cout << " " << std::setw(width2) << "LLL";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += bkz_time[i];
  }
  int width3 = log10(tmp)+2;
  total_times[2] = tmp;
  std::cout << " " << std::setw(width3) << "BKZ";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += min_time[i];
  }
  int width4 = log10(tmp)+2;
  total_times[3] = tmp;
  std::cout << " " << std::setw(width4) << "Minkowski";
  std::cout << "Total time " << std::setw(width1) << total_times[0]
    << std::setw(width2) << total_times[1]
    << std::setw(width3) << total_times[2]
    << std::setw(width4) << total_times[3];
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 << std::setw(width1)
      << die_time[i] << std::setw(width2) << lll_time[i] << std::setw(width3)
      << bkz_time[i] << std::setw(width4) << min_time[i];
  }
  
  return 0;
}
