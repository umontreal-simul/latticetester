/**
 * A small test for LLL64
 **/

//#define TYPES_CODE  LD     // int64_t
#define TYPES_CODE  ZD     // ZZ

#include <iostream>
#include <cstdint>
#include <ctime>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
// #include "latticetester/IntLattice.h"
#include "latticetester/EnumTypes.h"

using namespace LatticeTester;

int64_t m = 101;      // Modulus
const PrecisionType prec = DOUBLE;  // For LLL construction.

int main() {

	BasisConstruction<Int> constr; // The basis constructor we use.
    long dim = 3;
	IntMat bas_mat, bas_copy;
	bas_mat.SetDims(dim, dim);
	bas_copy.SetDims(dim, dim);
	long korBase[][3] = {{1, 12, 43}, {0, 101, 0}, {0, 0, 101}};
	long i, j;
    for (i = 0; i < dim; i++) {
	    for (j = 0; j < dim; j++) {
            bas_mat[i][j] = korBase[i][j];
	    }
    }
	std::cout <<  bas_mat << "  \n";
    // copy(bas_mat, bas_copy);
    constr.LLLConstruction0(bas_mat, 0.9);

	std::cout <<  bas_mat << "  \n";
	std::cout << std::endl << std::endl;
	return 0;
	}

