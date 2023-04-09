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

//  These are selected prime moduli to be used in examples.
const int64_t primes[] = { 101, 1021, 1048573, 1073741827, 1099511627791,
		                   1125899906842597 };

// int64_t m = primes[1];      // Modulus = 1021
Int m(primes[0]);           // Modulus near 2^{20}
const long numSizes = 1;    // Number of matrix sizes (choices of dimension).
//const long dimensions[numSizes] = {10, 20, 30, 40};
const long dimensions[numSizes] = {5};

const PrecisionType prec = DOUBLE;  // For LLL construction.
const long numRep = 1;  // Number of replications for each case.
const long numMeth = 10;    // Number of methods, and their names.
std::string names[numMeth] = {"UppTri  ", "LowTri  ", "TriGCD  ", "Tri96   ",
		 "LLL5    ", "LLL8    ", "LLL9    ", "DualUT  ", "DualUT96", "Dual    "};

int main() {
	// We use ctime for the timings, for implementation simplicity
	clock_t totalTime = clock();  // Global timer for total time.
	clock_t timer[numMeth][numSizes];
	clock_t tmp;

	IntMat bas_mat, bas_copy, m_v, m_v2;
	Rank1Lattice<Int, double> *korlat;    // Will be a Korobov lattice.
    Int a;  // The LCG multiplier.

    long d;
	for (d = 0; d < numSizes; d++) {  // Each matrix size
		long dim = dimensions[d]; // The corresponding dimension.
		bas_mat.SetDims(dim, dim);
		bas_copy.SetDims(dim, dim);
		m_v.SetDims(dim, dim);
		m_v2.SetDims(dim, dim);
		for (int64_t meth = 0; meth < numMeth; meth++)
			timer[meth][d] = 0;
	    for (int64_t r = 0; r < numRep; r++) {
            a = m / 5 + 10 * r;   // The multiplier we use for this rep.
            korlat = new Rank1Lattice<Int, Real>(m, a, dim);
			korlat->buildBasis (dim);
            copy(korlat->getBasis(), bas_mat);  // This basis is triangular.
    	    std::cout  << "korlat basis = \n" << bas_mat << "\n";

            // We apply LLL to change it.
			BasisConstruction<Int>::LLLConstruction0(bas_mat, 0.5);
    	    std::cout  << "LLL5 basis = \n" << bas_mat << "\n";

            copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::upperTriangularBasis(bas_copy, m_v, m);
			timer[0][d] += clock() - tmp;
    	    std::cout  << "UT basis = \n" << m_v << "\n";

			// copy(bas_mat, bas_copy);
			// tmp = clock();
			// BasisConstruction<Int>::lowerTriangularBasis(bas_copy, m_v, m);
			// timer[1][d] += clock() - tmp;

			copy(bas_mat, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::GCDTriangularBasis(bas_copy, m);
			timer[2][d] += clock() - tmp;
    	    std::cout  << "GCD-UT basis = \n" << bas_copy << "\n";
    	    std::cout  << "*** WRONG ***  This also changes the dimensions of bas_copy.\n" ;
    		bas_copy.SetDims(dim, dim);

			copy(bas_mat, bas_copy);
			tmp = clock();
			// This one is in Util.h, it is the old method from 1996.
			Triangularization(bas_copy, m_v, dim, dim, m);
			timer[3][d] += clock() - tmp;
    	    std::cout  << "UT96 basis = \n" << m_v << "\n";

			// The basis m_v is triangular.
			// std::cout << " The LLL construction with delta=0.5 \n";
			copy(m_v, bas_copy);
    	    std::cout  << "Before LLL5 m_v basis = \n" << m_v << "\n";
    	    std::cout  << "Before LLL5 basis = \n" << bas_copy << "\n";
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.5);
			timer[4][d] += clock() - tmp;
    	    std::cout  << "LLL5 basis = \n" << bas_copy << "\n";

			// std::cout << " The LLL construction with delta=0.5 \n";
			copy(m_v, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.8);
			timer[5][d] += clock() - tmp;

			// std::cout << " The LLL construction with delta=0.5 \n";
			copy(m_v, bas_copy);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(bas_copy, 0.99999);
			timer[6][d] += clock() - tmp;
    	    std::cout  << "LLL9 basis = \n" << bas_copy << "\n";

			copy(m_v, bas_copy);  // This one should be upper triangular.
    	    std::cout  << "UT basis = \n" << bas_copy << "\n";
			//if (!CheckTriangular(m_v, dim, m)) {
			//	std::cout << "Matrix not triangular! \n";
			//}
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular(bas_copy, m_v2, m);
			timer[7][d] += clock() - tmp;
    	    std::cout  << "dual LT basis = \n" << m_v2 << "\n";

			copy(m_v, bas_copy);
			tmp = clock();
			//  std::cout << "value of m = " << m << "\n";
            // This one was changing the value of m !!!!!
		    BasisConstruction<Int>::mDualUpperTriangular96(bas_copy, m_v2, m);
			timer[8][d] += clock() - tmp;
    	    std::cout  << "dual LT basis = \n" << m_v2 << "\n";

			copy(bas_mat, bas_copy);
			tmp = clock();
            BasisConstruction<Int>::mDualBasis(bas_copy, m_v, m);
			timer[9][d] += clock() - tmp;
    	    std::cout  << "dual basis = \n" << m_v << "\n";
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

