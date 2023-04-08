/**
 * This example showcases the usage of the example input creator.
 **/


#define TYPES_CODE  ZD

#include <iostream>
#include <cstdint>
#include <ctime>
#include <NTL/mat_GF2.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"

using namespace LatticeTester;


int main() {

	BasisConstruction<Int> constr; // The basis constructor we use.
	//Here we set up the matrix parameters
    int64_t d = 40; //dimension
    int64_t m = 1021; //modulus
    IntVec a; //vector for 
    int64_t k = 6; //length of reccurence in the LCG
    a.SetLength(k);
    
    //Define entries of a
    for (int64_t i =  0; i < k; i++) {
    	a[i] = (i+1)*(m-73);
    }
    
    IntMat V; //output variable
	constr.CreateExampleMatrix(m, d, a, V);
	std::cout << V;
    
    int64_t no = 0;//number of output file
	constr.CreateExampleMatrixToFile(m, d, a, no);
	

	return 0;
	}
