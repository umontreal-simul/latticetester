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
 * OLD RESULTS on my laptop with m = 1021:   (outdated, with ZZ)
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
=========================================================

New results on my laptop with m = 1048573:

lecuyer@xubuntu-22042:~/git/latticetester/build/examples$ ./BasisManipulationLD  (int64_t)

 dim:          10       20       30       40

UppTri        385     1518     4779     9151
TriGCD        966     4632    14416    30021
Tri96         799     4587    14860    31337
LLL5          391     1508     4482     7286
LLL8          532     3402    10943    19267
LLL9          646     6316    25108    50430
DualUT        100      806     2029     3881
DualUT96      252     1860     6129    13911

Total time: 0.631648 seconds

lecuyer@xubuntu-22042:~/git/latticetester/build/examples$ ./BasisManipulation  (ZZ)

 dim:          10       20       30       40

UppTri        690     3075     8153    16477
TriGCD       2297    13100    38277    86614
Tri96        3919    25362    77597   178568
LLL5          286      604     1152     1882
LLL8          213     1074     2241     4184
LLL9          266     1929     4858     8752
DualUT        148     1018     3048     7040
DualUT96      328     1909     5818    13447
Dual         1263     7430    25956    65896

Total time: 1.07714 seconds

 **/


#define TYPES_CODE  ZD

#include <iostream>
#include <cstdint>
#include <ctime>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/ParamReader.h"

using namespace LatticeTester;

//  The following array gives the possible modulo values for the basis examples.
const int many_primes = 7;
const std::string primes[] = { "101", "1021", "1048573", "1073741827", "1099511627791",
		"1125899906842597", "18446744073709551629" };

//Use basis values modulo m
const std::string prime = primes[2];
// Int m = primes[2];      // Modulus near 2^{20}
Int m(1048573);      // Modulus
const PrecisionType prec = DOUBLE;  // For LLL construction.
const int numRep = 10;  // Number of replications for each case.
const int numMeth = 10;    // Number of methods, and their names.
std::string names[numMeth] = {"UppTri  ", "LowTri  ", "TriGCD  ", "Tri96   ",
		 "LLL5    ", "LLL8    ", "LLL9    ", "DualUT  ", "DualUT96", "Dual    "};
const int numSizes = 4; // Number of matrix sizes (choices of dimension).
// const int dimensions[numSizes] = {10, 20};
const int dimensions[numSizes] = {10, 20, 30, 40};

int main() {
	// We use ctime for implementation simplicity
	clock_t totalTime = clock();  // Global timer for total time.
	clock_t timer[numMeth][numSizes];
	clock_t tmp;

	// BasisConstruction<Int> constr; // The basis constructor we use.
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
			BasisConstruction<Int>::upperTriangularBasis(bas_copy, m_v, m);
			timer[0][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			//BasisConstruction<Int>::lowerTriangularBasis(bas_copy, m_v, m);
			timer[1][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::GCDTriangularBasis(bas_copy, m);
			timer[2][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			// This one is in Util.h, it is the old method from 1996.
			Triangularization(bas_copy, m_v, dim, dim, m);
			timer[3][d] += clock() - tmp;

			// std::cout << " The LLL construction with delta=0.5 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.5, prec);
			timer[4][d] += clock() - tmp;

			// std::cout << " The LLL construction with delta=0.8 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.8, prec);
			timer[5][d] += clock() - tmp;

			// std::cout << " The LLL constructio with delta=0.99999 \n";
			copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.99999, prec);
			timer[6][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			BasisConstruction<Int>::upperTriangularBasis(bas_copy, m_v, m);
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular(m_v, m_v2, m);
			timer[7][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			BasisConstruction<Int>::upperTriangularBasis(bas_copy, m_v, m);
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular96(m_v, m_v2, m);
			timer[8][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
            BasisConstruction<Int>::mDualBasis(bas_copy, m_v, m);
			timer[9][d] += clock() - tmp;
		}
	}

	std::cout << " dim:    ";
	for (d = 0; d < numSizes; d++)
	    std::cout << std::setw(8) << dimensions[d] << " ";
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

