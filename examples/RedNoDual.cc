/**
 * This compares the execution of our algorithms with and without dual.
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
 * This executes LLL reduction, BKZ reduction, Shortest vector search and
 * Minkowski reduction with and without dual and computes the time they take.
 * */
int main() {
  int max_dim = 20; // Actual max dim is 5*max_dim
  // This is basically the C method of timing a program. We time globally, but
  // also for each dimension and for each size of integers in the matrix.
  clock_t lll_dual[max_dim], lll_nodual[max_dim],
  bkz_dual[max_dim], bkz_nodual[max_dim],
  sho_dual[max_dim], sho_nodual[max_dim],
  sho_bkz_dual[max_dim], sho_bkz_nodual[max_dim],
  min_dual[max_dim], min_nodual[max_dim], tmp;
  clock_t total_times[10];
  for (int i = 0; i < max_dim; i++){
    lll_dual[i] = 0;
    bkz_dual[i] = 0;
    sho_dual[i] = 0;
    sho_bkz_dual[i] = 0;
    min_dual[i] = 0;
    lll_nodual[i] = 0;
    bkz_nodual[i] = 0;
    sho_nodual[i] = 0;
    sho_bkz_nodual[i] = 0;
    min_nodual[i] = 0;
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

  std::cerr << "LLL\n";
  // Timing LLL with default parameters
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
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
      lll_dual[j] += clock() - tmp;
    }
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redLLLNTLNoDual();
      delete red;
      delete basis;
      lll_nodual[j] += clock() - tmp;
    }
  }
  std::cerr << "BKZ\n";
  // Timing BKZ with default parameters
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
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
      bkz_dual[j] += clock() - tmp;
    }
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redBKZNoDual();
      delete red;
      delete basis;
      bkz_nodual[j] += clock() - tmp;
    }
  }
  std::cerr << "Shortest\n";
  // Timing shortest reduction
  for (int j = 0; j < max_dim; j++) {
    // The execution is too costly for great dimensions.
    if (j > 7) continue;
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->shortestVector(L2NORM);
      delete red;
      delete basis;
      sho_dual[j] += clock() - tmp;
    }
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->shortestVectorNoDual(L2NORM);
      delete red;
      delete basis;
      sho_nodual[j] += clock() - tmp;
    }
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
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
      red->shortestVector(L2NORM);
      delete red;
      delete basis;
      sho_bkz_dual[j] += clock() - tmp;
    }
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // We time NTL first
      basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
      red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
      red->redBKZNoDual();
      red->shortestVectorNoDual(L2NORM);
      delete red;
      delete basis;
      sho_bkz_nodual[j] += clock() - tmp;
    }
  }

  //std::cerr << "Minkowski\n";
  //// Timing Minkowski reduction
  //for (int j = 0; j < max_dim; j++) {
  //  // The execution is too costly for great dimensions.
  //  if (j > 7) continue;
  //  for (int k = 0; k < 10; k++) {
  //    tmp = clock();
  //    // Reader shenanigans
  //    name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
  //    reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
  //    reader.getLines();
  //    reader.readInt(numlines, 0, 0);
  //    matrix1.SetDims(numlines, numlines);
  //    ln = 1;
  //    reader.readBMat(matrix1, ln, 0, numlines);

  //    // We time NTL first
  //    basis = new IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix1, numlines);
  //    red = new Reducer<MScal, BScal, NScal, RScal>(*basis);
  //    red->redBKZ();
  //    red->reductMinkowski(0);
  //    delete red;
  //    delete basis;
  //    min_time[j] += clock() - tmp;
  //  }
  //}

  // Printing the results in a somewhat formated way.
  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += lll_dual[i];
  }
  int width1 = log10(tmp)+2;
  total_times[0] = tmp;
  std::cout << "          " << std::setw(width1) << "LLL 1";
  for (int i = 0; i<width1-3; i++) {
    std::cout << " ";
  }

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += lll_nodual[i];
  }
  int width2 = log10(tmp)+2;
  total_times[1] = tmp;
  std::cout << " " << std::setw(width2) << "LLL 2";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += bkz_dual[i];
  }
  int width3 = log10(tmp)+2;
  total_times[2] = tmp;
  std::cout << " " << std::setw(width3) << "BKZ 1";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += bkz_nodual[i];
  }
  int width4 = log10(tmp)+2;
  total_times[3] = tmp;
  std::cout << " " << std::setw(width4) << "BKZ 2";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += sho_dual[i];
  }
  int width5 = log10(tmp)+2;
  total_times[4] = tmp;
  std::cout << " " << std::setw(width4) << "Short 1";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += sho_nodual[i];
  }
  int width6 = log10(tmp)+2;
  total_times[5] = tmp;
  std::cout << " " << std::setw(width4) << "Short 2";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += sho_bkz_dual[i];
  }
  int width7 = log10(tmp)+2;
  total_times[6] = tmp;
  std::cout << " " << std::setw(width4) << "Short BKZ 1";

  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += sho_bkz_nodual[i];
  }
  int width8 = log10(tmp)+2;
  total_times[7] = tmp;
  std::cout << " " << std::setw(width4) << "Short BKZ 2\n";

  std::cout << "Total time" << std::setw(width1) << total_times[0]
    << std::setw(width2) << total_times[1]
    << std::setw(width3) << total_times[2]
    << std::setw(width4) << total_times[3]
    << std::setw(width5) << total_times[4]
    << std::setw(width6) << total_times[5]
    << std::setw(width7) << total_times[6]
    << std::setw(width8) << total_times[7] << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 << std::setw(width1)
      << lll_dual[i] << std::setw(width2) << lll_nodual[i] << std::setw(width3)
      << bkz_dual[i] << std::setw(width4) << bkz_nodual[i] << std::setw(width5)
      << sho_dual[i] << std::setw(width6) << sho_nodual[i] << std::setw(width7)
      << sho_bkz_dual[i] << std::setw(width8) << sho_bkz_nodual[i] << std::endl;
  }
  
  return 0;
}
