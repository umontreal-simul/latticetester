#ifndef EXAMPLES_H
#define EXAMPLES_H
/**
 * This file defines functions and constants used in all examples.
 * ***  Bad class name!
 * */

const int many_primes = 6;
const std::string primes[] = {"1021", "1048573", "1073741827", "1099511627791",
"1125899906842597", "18446744073709551629"};

/* This function does 3 things.
 * It computes the sum of the first `dim` elements in time and returns the width(+1)
 * (number of characters needed to write it) of the resulting integer in base 10.
 * This stores that sum in totals[ind].
 * Finally this prints message on std::cout in a column of the width computed
 * above.
 * We abuse this function to format the output of the examples.
 * */
int getWidth(clock_t time[], int dim, std::string message, clock_t totals[], int ind) {
  clock_t tmp = 0;
  for (int i = 0; i < dim; i++) {
    tmp += time[i];
  }
  int width = log10(tmp) + 2;
  std::cout << std::setw(width) << message;
  totals[ind] = tmp;
  return width;
}

#endif
