// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cmath>
#include "latticetester/LatticeTesterRoutines.h"

namespace LatticeTester {

double ShortestVector(BMat matrix, NormType norm, PreReductionType preRed, 
	PrecisionType doublePrecision, double fact, int blocsize)
{
	int dimension;
	if (matrix.size1() != matrix.size2()) {
		MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
		exit(1);
		// C'est pas un peu nul ça ?
	} else 
		dimension = matrix.size1();

	// creating objects needed to perform the test
	IntLatticeBasis basis (matrix, dimension, norm);
	Reducer red (basis);

	// performing pre-reduction
	switch (preRed) {
		case BKZ:
			red.redBKZ(doublePrecision, fact, blocksize);
			break;
		case LenstraLL:
			red.redLLLNTL(doublePrecision, fact);
			break;
		case PreRedDieter:
			red.preRedDieter(0);
			break;
		case NOPRERED:
			cout << "WARNING: no pre-reduction is performed before applying Branch-and-Bound";
			cout << " procedure. Running time could increase dramaticaly with the matrix dimension.";
			cout << endl;
			break;
		default:
			MyExit(1, "LatticeTesterRoutines::ShortestVector:   NO SUCH CASE FOR PreReductionType");
			exit(1);
	}

	// performing the Branch-and-Bound procedure to find the shortest non-zero vector.
	// foundShortest bool is set to true if the algorithm terminates without any error.
	bool foundShortest;
    foundShortest = red.shortestVector(norm);

    double length;
    if (foundShortest) {
	    if (norm == L2NORM) { //L2NORM is stored squared
	    	length = conv<double>(red.getMinLength());
	    	length = sqrt(length);
	    } else if (norm == L1NORM)
	    	length = conv<double>(red.getMinLength());
	} else 
		length = -1.0;

	return length;
}




} // end namespace LatticeTester




enum NormType { SUPNORM = 1, L1NORM = 2, L2NORM = 3, ZAREMBANORM = 4 };

PrecisionType { DOUBLE, QUADRUPLE, EXPONENT, ARBITRARY, EXACT };

PreReductionType {BKZ, PreRedDieter, LenstraLL, NOPRERED};



BKZ precision, fact, blocskize

LLL precision fact

PreRedDieter

NOpredred