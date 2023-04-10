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
 * In this example, we can compare the speed of BasisConstruction<Int>::calcDual method
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


 **/

//#define TYPES_CODE  LD     // int64_t
#define TYPES_CODE  ZD     // ZZ

#include <iostream>
#include <cstdint>
#include <ctime>
#include <type_traits>
#include <typeinfo>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"

using namespace LatticeTester;

//Int m(101);      // Modulus m = 101
//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Prime modulus near 2^{20}
//Int m(1073741827);  // Prime modulus near 2^{30}
//Int m(1099511627791);  // Prime modulus near 2^{40}
//Int m(1125899906842597);  // Prime modulus near 2^{50}
Int a;       // The LCG multiplier

const long numSizes = 5;    // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 10, 20, 30, 40, 50 };
//const long dimensions[numSizes] = {5};

const long numRep = 10;  // Number of replications for each case.
const long numMeth = 8;    // Number of methods, and their names.
std::string names[numMeth] = { "UppTri   ", "Tri96    ", "LLL5     ",
		"LLL8     ", "LLL9     ", "mDualUT  ", "mDualUT96", "mDual    " };

int main() {
	// We use ctime for the timings, for implementation simplicity
	clock_t totalTime = clock();  // Global timer for total time.
	clock_t timer[numMeth][numSizes];
	clock_t tmp;
	IntMat basis1, basis2, basis3, basisdual;
	Int sqlength;
	Rank1Lattice<Int, double> *korlat;    // Will be a Korobov lattice.

	long d;
	for (d = 0; d < numSizes; d++) {  // Each matrix size
		long dim = dimensions[d]; // The corresponding dimension.
		basis1.SetDims(dim, dim);
		basis2.SetDims(dim, dim);
		basis3.SetDims(dim, dim);
		basisdual.SetDims(dim, dim);
		for (int64_t meth = 0; meth < numMeth; meth++)
			timer[meth][d] = 0;
		for (int64_t r = 0; r < numRep; r++) {
			a = m / 5 + 10 * r;   // The multiplier we use for this rep.
			korlat = new Rank1Lattice<Int, Real>(m, a, dim);
			korlat->buildBasis(dim);
			copy(korlat->getBasis(), basis1); // This initial basis is triangular.
			// We apply LLL to change basis1.
			BasisConstruction<Int>::LLLConstruction0(basis1, 0.5);
			copy(basis1, basis2);

			tmp = clock();
			BasisConstruction<Int>::upperTriangularBasis(basis2, basis3, m);
			timer[0][d] += clock() - tmp;

			copy(basis1, basis2);
			// This one is in Util.h, it is the old method from 1996.
			tmp = clock();
			Triangularization(basis2, basis3, dim, dim, m);
			timer[1][d] += clock() - tmp;
			// This basis3 is upper triangular.

			copy(basis1, basis2);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.5);
			timer[2][d] += clock() - tmp;

			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.8);
			timer[3][d] += clock() - tmp;

			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.99999);
			timer[4][d] += clock() - tmp;

			copy(basis3, basis2);  // This one should be upper triangular.
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular(basis2, basisdual, m);
			timer[5][d] += clock() - tmp;

			copy(basis3, basis2);
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular96(basis2, basisdual, m);
			timer[6][d] += clock() - tmp;

#if TYPES_CODE  ==  ZD
			// copy(basis3, basis2);
			tmp = clock();
			BasisConstruction<Int>::mDualBasis(basis3, basisdual, m);
			timer[6][d] += clock() - tmp;
#endif
			}
	}
	std::cout << "Types: " << strFlexTypes << "\n";
	std::cout << "Timings for different methods, in basic clock units \n";
	std::cout << " dim:    ";
	for (d = 0; d < numSizes; d++)
		std::cout << std::setw(8) << dimensions[d] << " ";
	std::cout << std::endl << std::endl;
	for (int meth = 0; meth < numMeth; meth++) {
		std::cout << names[meth] << " ";
		for (d = 0; d < numSizes; d++)
			std::cout << std::setw(8) << timer[meth][d] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Total time: "
			<< (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
			<< " seconds\n";
	return 0;
}

