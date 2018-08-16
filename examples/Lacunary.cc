/**
 * This examples showcases the usage of Lacunary indices sets construction to
 * perform tests with the same basis on different projection.
 * It then uses the Weights module to get a figure of merit based on those
 * projections.
 * */

#include <iostream>
#include "NTL/ZZ.h"

int main() {
  NTL::ZZ x(-10);
  NTL::ZZ y(3);
  std::cout << x/y << std::endl;
  std::cout << -10/3 << std::endl;
  std::cout << x%y << std::endl;
  std::cout << (-x)%(-y) << std::endl;
  std::cout << -10%3 << std::endl;
  return 0;
}
