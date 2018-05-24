#include "latticetester/IntFactor.h"

int main(){
    long x = 17179607041;
    std::cout << LatticeTester::IntFactor::isPrime(x, 3) << std::endl;
}