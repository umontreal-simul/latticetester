/**
 * This example showcases the usage of the BasisConstruction module. This reads
 * matrices from files and builds a basis and a dual for an `IntLattice`
 * object. The files this is set to use are in the `bench.zip` archive. To
 * execute the program, the archive should be unziped and the `bench` folder
 * should be put in the same directory from which the executable is called.

 * The bases we used to showcase of BasisConstruction methods are in a folder named 
 * 'examples/bench/'. Each file in 'examples/bench/' folder contain a basis, and the 
 * file is nameed as follows: 'prime_dimBasis_exanpleNumber' where 'prime' is modulo 
 * value of the basis, 'dimBasis' is the dimension of the basis, and 'exampleNumber' 
 * is the number of the example for the bases of dimension 'dimBasis'.
 *
 * This example reads matrices from files and performs the different construction
 * algorithms in BasisConstruction on them. The program then prints the execution
 * time of the various algorithms. Note that the execution of the program is not
 * what you would expect in reality since bench contains random full matrices.
 *
 * We show a use of BasisContruction::upperTriangularBasis, 
 * BasisContruction::lowerTriangularBasis, BasisContruction::LLLConstruction with 
 * two differrent parametter of delta, BasisContruction::mDualUpperTriangular, 
 *  BasisContruction::mDualBasis
 *
 * In this example, we can compare the speed of BasisConstruction::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and BasisConstruction::mDualUpperTriangular method which compute an m-dual basis
 * with an upper triangular basis.
 *
 * We can also compare the speed of 'BasisConstruction::upperTriangularBasis'
 * and the speed of 'BasisConstruction::LLLConstruction'
 *
 * RESULTS with m = 1021:
 *
 *  dim:        10      20       30       40

UppTri       6590    40472   175472   561106
LowTri       7070    46714   194587   593768
TriGCD       5048    61863   438167  1762391
Tri96        5233    38514   166644   537749
LLL5          233      589     1272     2226
LLL8          263     1106     2801     3858
LLL9          313     1903     5016     9633
DualUT        909     3902    34482   347704
DualUT96     1443     6203    63702   502523
Dual         1261     6574    18648    45761

Total time: 7.55269 seconds

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
#include "latticetester/ParamReader.h"
// #include "latticetester/IntLattice.h"
#include "latticetester/EnumTypes.h"

using namespace LatticeTester;

//  The following array gives the possible modulo values for the basis examples.
const int many_primes = 6;
const std::string primes[] = { "1021", "1048573", "1073741827", "1099511627791",
		"1125899906842597", "18446744073709551629" };

//Use basis values modulo 1021
const std::string prime = primes[0];
Int m(1021);      // Modulus
const PrecisionType prec = DOUBLE;  // For LLL construction.
const int numSizes = 4; // Number of matrix sizes (choices of dimension).
const int numRep = 10;  // Number of replications for each case.
const int numMeth = 10;    // Number of methods, and their names.
std::string names[numMeth] = {"UppTri  ", "LowTri  ", "TriGCD  ", "Tri96   ",
		 "LLL5    ", "LLL8    ", "LLL9    ", "DualUT  ", "DualUT96", "Dual    "};
const int dimensions[numSizes] = {10, 20, 30, 40};

int main() {
	// We use ctime for implementation simplicity
	clock_t totalTime = clock();  // Global timer for total time.
	clock_t timer[numMeth][numSizes];
	clock_t tmp;

	BasisConstruction<Int> constr; // The basis constructor we use.
	IntMat bas_mat, bas_copy, m_v, m_v2;
    int d;

	for (d = 0; d < numSizes; d++) {  // Each matrix size
		unsigned int dim = dimensions[d]; // The corresponding dimension.
		bas_mat.SetDims(dim, dim);
		bas_copy.SetDims(dim, dim);
		m_v.SetDims(dim, dim);
		m_v2.SetDims(dim, dim);
		for (int meth = 0; meth < numMeth; meth++)
			timer[meth][d] = 0;
	    for (int r = 0; r < numRep; r++) {
	    	// We use a different file for each rep.
		    std::string name = "bench/" + prime + "_"
		       + std::to_string(dim) + "_" + std::to_string(r);
		    ParamReader<Int, Real> reader(name + ".dat");
	        reader.getLines();
	        int thisdim;
			reader.readInt(thisdim, 0, 0);
			unsigned int ln=1;
			reader.readBMat(bas_mat, ln, 0, dim);

			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.upperTriangularBasis(bas_copy, m_v, m);
			timer[0][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.lowerTriangularBasis(bas_copy, m_v, m);
			timer[1][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.GCDTriangularBasis(bas_copy, m);
			timer[2][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			// This one is in Util.h, it is the old method from 1996.
			Triangularization(bas_copy, m_v, dim, dim, m);
			timer[3][d] += clock() - tmp;

			// std::cout << " The LLL construction with delta=0.5 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.LLLConstruction0(bas_copy, 0.5, prec);
			timer[4][d] += clock() - tmp;

			// std::cout << " The LLL construction with delta=0.8 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.LLLConstruction0(bas_copy, 0.8, prec);
			timer[5][d] += clock() - tmp;

			// std::cout << " The LLL constructio with delta=0.99999 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			constr.LLLConstruction0(bas_copy, 0.99999, prec);
			timer[6][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			constr.upperTriangularBasis(bas_copy, m_v, m);
			tmp = clock();
			constr.mDualUpperTriangular(m_v, m_v2, m);
			timer[7][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			constr.upperTriangularBasis(bas_copy, m_v, m);
			tmp = clock();
			constr.mDualUpperTriangular96(m_v, m_v2, m);
			timer[8][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
            constr.mDualBasis(bas_copy, m_v, m);
			timer[9][d] += clock() - tmp;
		}
	}

	std::cout << " dim:  ";
	for (d = 0; d < numSizes; d++)
	    std::cout << std::setw(6) << dimensions[d] << " ";
	std::cout << std::endl << std::endl;
	for (int meth = 0; meth < numMeth; meth++) {
	    std::cout  << names[meth] << " ";
	    for (d = 0; d < numSizes; d++)
		   std::cout << std::setw(8) << timer[meth][d] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Total time: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
			<< " seconds\n";
	return 0;
	}

