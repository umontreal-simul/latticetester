/*
 * This is an example that showcases the usage of CoordinateSets and benchmarks
 * the different implementations that we have in the
 * LatticeTester::CoordinateSets namespace.
 *
 * This should be expanded to include the AddCoordinate class, and to do a
 * proper benchmark.
 * */

#include "latticetester/CoordinateSets.h"

#include <ctime>

using namespace LatticeTester;

int main() {
  clock_t tt1 = 0, tt2 = 0, ranges_time, sets_time;
  long size = 27;
  CoordinateSets::FromRanges ranges(1, size, 1, size);
  long a[size];
  for (int i = 0; i < size; i++) {
    a[i] = i+1;
  }
  Coordinates coord (a, a+size);
  CoordinateSets::Subsets    sets(coord, 1, size);
  auto it1 = ranges.begin();
  auto it2 = sets.begin();
  for (int i = 0; i < 10; i++) {
    it1 = ranges.begin();
    it2 = sets.begin();
    ranges_time = clock();
    while (it1 != ranges.end()) {
      it1++;
    }
    tt1 += clock() - ranges_time;
    sets_time = clock();
    while (it2 != sets.end()) {
      it2++;
    }
    tt2 += clock() - sets_time;
  }
  std::cout << "FromRanges average time: " << (double)tt1/CLOCKS_PER_SEC/10 << std::endl;
  std::cout << "Subsets average time:    " << (double)tt2/CLOCKS_PER_SEC/10 << std::endl;
  return 0;
}
