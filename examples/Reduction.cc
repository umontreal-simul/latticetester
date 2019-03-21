/**
 * This is an example to the usage of the reducer class to perform reduction
 * of lattices. This serves as a comparison between the different reduction
 * methods when used before searching for the shortest vector. The output is
 * formated to include the number of clock ticks spent on each algorithm as well
 * as the number of times the program failed to find the shortest vector for
 * each pre-reduction.
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

#include "Examples.h"

using namespace LatticeTester;

int main() {
  clock_t timer = clock();
  int max_dim = 10; // Actual max dim is 5*max_dim
  // This is basically the C method of timing a program. We time globally, but
  // also for eache dimension and for each size of integers in the matrix.
  clock_t die_time[max_dim], lll_time[max_dim], bkz_time[max_dim],
  sho_die[max_dim], sho_lll[max_dim], sho_bkz[max_dim], tmp;
  clock_t total_times[6];
  for (int i = 0; i < max_dim; i++){
    lll_time[i] = 0;
    die_time[i] = 0;
    bkz_time[i] = 0;
    sho_die[i] = 0;
    sho_lll[i] = 0;
    sho_bkz[i] = 0;
  }
  int die_fails=0, lll_fails=0, bkz_fails=0;

  std::string prime = primes[0];

  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      // We dynamically allocate memory to these two pointers every time we need to
      // create an object of their type.
      IntLatticeBasis<MScal, BScal, NScal, RScal>* basis;
      Reducer<MScal, BScal, NScal, RScal>* red;
      // Variables definition
      ParamReader<MScal, BScal, RScal> reader;
      std::string name;
      int numlines;
      BMat matrix1;
      unsigned int ln;
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time Dieter first
      tmp = clock();
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redDieter(0);
      die_time[j] += clock() - tmp;
      tmp = clock();
      if (!red->shortestVector(L2NORM)) {
        die_fails++;
      }
      sho_die[j] += clock() - tmp;
      delete red;
      delete basis;

      // Next LLL
      tmp = clock();
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redLLLNTL();
      lll_time[j] += clock() - tmp;
      tmp = clock();
      if (!red->shortestVector(L2NORM)) {
        lll_fails++;
      }
      sho_lll[j] += clock() - tmp;
      delete red;
      delete basis;

      // Then BKZ
      tmp = clock();
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redBKZ();
      bkz_time[j] += clock() - tmp;
      tmp = clock();
      if (!red->shortestVector(L2NORM)) {
        bkz_fails++;
      }
      sho_bkz[j] += clock() - tmp;
      delete red;
      delete basis;
    }
  }

  // Printing the results in a somewhat formated way.
  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  std::cout << "          ";
  int width1 = getWidth(die_time, max_dim, "Dieter", total_times, 0);
  int width2 = getWidth(lll_time, max_dim, "LLL", total_times, 1);
  int width3 = getWidth(bkz_time, max_dim, "BKZ", total_times, 2);
  int width4 = getWidth(sho_die, max_dim, "SV Dieter", total_times, 3);
  int width5 = getWidth(sho_lll, max_dim, "SV LLL", total_times, 4);
  int width6 = getWidth(sho_bkz, max_dim, "SV BKZ", total_times, 5);
  std::cout << std::endl;

  std::cout << "Total time" << std::setw(width1) << total_times[0]
    << std::setw(width2) << total_times[1]
    << std::setw(width3) << total_times[2]
    << std::setw(width4) << total_times[3]
    << std::setw(width5) << total_times[4]
    << std::setw(width6) << total_times[5] << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 
      << std::setw(width1) << die_time[i] 
      << std::setw(width2) << lll_time[i] 
      << std::setw(width3) << bkz_time[i] 
      << std::setw(width4) << sho_die[i] 
      << std::setw(width5) << sho_lll[i] 
      << std::setw(width6) << sho_bkz[i] 
      << std::endl;
  }
  std::cout << "Fails     "
    << std::setw(width1) << die_fails 
    << std::setw(width2) << lll_fails 
    << std::setw(width3) << bkz_fails 
    << std::endl;

  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";
  
  return 0;
}
