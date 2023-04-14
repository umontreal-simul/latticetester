/**
 * This example illustrates how to use the functions of the BasisConstruction class,
 * with a small five-dimensional lattice obtained from an LCG of modulus 101.
 * The basis is printed between function calls, to show what is going on.
 * An initial upper-triangular basis is constructed by Rank1Lattice in a standard way.
 * We then call LLL with `delta = 0.5` to obtain a basis with shorter vectors,
 * which is not triangular. Calling `upperTriangulatBasis` transforms this basis
 * to an upper triangular basis, which happens to be the same as the initial one.
 * Different triangularization methods are compared.  Then LLL with different
 * values of `delta`. Finally a dual basis is computed in different ways.
 *
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
#include "latticetester/Reducer.h"

using namespace LatticeTester;

Int m(101);      // Modulus m = 101
//Int m(1021);     // Modulus m = 1021
//Int m(1048573);  // Modulus m = 1048573 (prime number near 2^{20})
Int a(33);       // The LCG multiplier
const long dim = 5;  // Dimension

int main() {
	IntMat basis1, basis2, basis3, basisdual;
	basis1.SetDims(dim, dim);
	basis2.SetDims(dim, dim);
	basis3.SetDims(dim, dim);
	basisdual.SetDims(dim, dim);
    Int sqlength;

	// We construct a Korobov lattice.
	Rank1Lattice<Int, double> *korlat;
	korlat = new Rank1Lattice<Int, Real>(m, a, dim);
	korlat->buildBasis(dim);
	copy(korlat->getBasis(), basis1);  // This initial basis is triangular.
	std::cout << "Initial Korobov lattice basis = \n" << basis1 << "\n";
	ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
	std::cout << "Square length of shortest vector: " << sqlength << "\n\n";

	// We apply LLL to change basis1.
	BasisConstruction<Int>::LLLConstruction0(basis1, 0.5);
	std::cout << "Basis after applying LLL with delta=0.5: \n" << basis1 << "\n";
	ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
	std::cout << "Square length of shortest vector: " << sqlength << "\n\n";

	copy(basis1, basis2);
	BasisConstruction<Int>::upperTriangularBasis(basis2, basis3, m);
	std::cout << "Conversion with `upperTriangularBasis`: \n" << basis3 << "\n\n";

	copy(basis1, basis2);
	BasisConstruction<Int>::GCDTriangularBasis(basis2, m);
	std::cout << "Conversion with GCDTriangular: \n" << basis2 << "\n\n";
	if (basis2[0][0] != 1) 	std::cout << "The latter is WRONG *** \n\n";
	basis2.SetDims(dim, dim);

	copy(basis1, basis2);
	// This one is in Util.h, it is the old method from 1996.
	Triangularization(basis2, basis3, dim, dim, m);
	std::cout << "Conversion with upperTriangular96: \n" << basis3 << "\n\n";
	// This basis3 is upper triangular.

	copy(basis1, basis2);
	BasisConstruction<Int>::LLLConstruction0(basis2, 0.5);
	std::cout << "Basis after applying LLL with delta=0.5: \n" << basis2 << "\n";
	ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
	std::cout << "Square length of shortest vector: " << sqlength << "\n\n";

	BasisConstruction<Int>::LLLConstruction0(basis2, 0.8);
	std::cout << "Basis after applying LLL with delta=0.8: \n" << basis2 << "\n";
	ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
	std::cout << "Square length of shortest vector: " << sqlength << "\n\n";

	BasisConstruction<Int>::LLLConstruction0(basis2, 0.99999);
	std::cout << "Basis after applying LLL with delta=0.99999: \n" << basis2 << "\n";
	ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
	std::cout << "Square length of shortest vector: " << sqlength << "\n\n";

	copy(basis3, basis2);  // This one should be upper triangular.
	BasisConstruction<Int>::mDualUpperTriangular(basis2, basisdual, m);
	std::cout << "m-dual upperTriangular: \n" << basisdual << "\n\n";

	copy(basis3, basis2);
	//  std::cout << "value of m = " << m << "\n";
	// This one was changing the value of m !!!!!
	BasisConstruction<Int>::mDualUpperTriangular96(basis2, basisdual, m);
	std::cout << "m-dual upperTriangular96: \n" << basisdual << "\n\n";

	// copy(basis3, basis2);
	BasisConstruction<Int>::mDualBasis(basis3, basisdual, m);
	std::cout << "m-dual basis by general method: \n" << basisdual << "\n";
	std::cout << "The latter is WRONG *** \n\n";

	Reducer<Int, Real> *red = new Reducer<Int, Real>(*korlat);
	red->shortestVector();
	std::cout << " The shortest vector length:\n";
	std::cout << red->getMinLength() << std::endl;

	return 0;
}

