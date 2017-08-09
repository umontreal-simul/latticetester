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

//*=============================================================================================	

double ShortestVector(BMat matrix, NormType norm, PreReductionType preRed, 
	PrecisionType doublePrecision, double fact, int blocksize)
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
			red.redBKZ(fact, blocksize, doublePrecision);
			break;
		case LenstraLL:
			red.redLLLNTL(fact, doublePrecision);
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
	// red.shortestVector(norm) is a bool set to *true* if the Branch-and-Bound algorithm
	// terminates without any error.
	if (red.shortestVector(norm))
		return conv<double>(red.getMinLength());
	else 
		return -1.0;
}

//*=============================================================================================

double ShortestVector(BMat matrix, NormType norm, long maxNodesBB, PreReductionType preRed, 
	PrecisionType doublePrecision, double fact, int blocksize)
{
	Reducer::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
	return ShortestVector(matrix, norm, preRed, doublePrecision, fact, blocksize);
}

//*=============================================================================================

double FigureOfMerit(BMat matrix, NormaType normalizerType, PreReductionType preRed,
	PrecisionType doublePrecision, double fact, int blocksize)
{
	double merit;

	int dimension;
	if (matrix.size1() != matrix.size2()) {
		MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
		exit(1);
		// C'est pas un peu nul ça ?
	} else 
		dimension = matrix.size1();

	// creation of the norm and normalizer objects
	NormType norm;
	Normalizer* normalizer;

	// calculation of the log-density of the matrix used to initialize the normalizer
	RScal logDensity;
#if NTL_TYPES_CODE > 1
	logDensity = - log( abs(NTL::determinant(matrix)) );
#else 
	// NTL library does not support matrix with double: we then use the boost library
	boost::numeric::ublas::matrix<long>  mat_tmps;
	mat_tmps.resize(dimension, dimension);
	for (unsigned int i = 0; i < dimension; i++) {
		for (unsigned int j = 0; j < dimension; j++) {
			mat_tmps(i,j) = matrix(i,j);
		}
	}
	logDensity = -log( abs(det_double(mat_tmps)) );
#endif

	// we initialize the norm and normalizer objects according to the normalizerType used
	switch (normalizerType) {
		case BESTLAT:
			norm = L2NORM;
			normalizer = new NormaBestLat (logDensity, dimension);
			break;
		case LAMINATED:
			norm = L2NORM;
			normalizer = new NormaLaminated (logDensity, dimension);
			break;
		case ROGERS:
			norm = L2NORM;
			normalizer = new NormaRogers (logDensity, dimension);
			break;
		case MINKOWSKI:
			norm = L2NORM;
			normalizer = new NormaMinkowski (logDensity, dimension);
			break;
		case MINKL1:
			norm = L1NORM;
			normalizer = new NormaMinkL1 (logDensity, dimension);
			break;
        default: //PALPHA_N, NORMA_GENERIC, L1, L2
        	MyExit(1, "LatticeTesterRoutines::FigureOfMerit:   NO SUCH CASE FOR *normalization type*");
			exit(1);
	}

	// compute the shortest non-zero vector
	merit = ShortestVector(matrix, norm, preRed, doublePrecision, fact, blocksize);
	
	if (merit == -1.0)
		return -1.0; // the BB procedure didn't terminated well

	// normalization step
	merit /= conv<double>(normalizer->getPreComputedBound(dimension));
	
	delete normalizer;
	return merit; 
}

//*=============================================================================================

double FigureOfMerit(BMat matrix, NormaType normalizerType, long maxNodesBB, 
	PreReductionType preRed, PrecisionType doublePrecision, 
	double fact, int blocksize)
{
	Reducer::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
	return FigureOfMerit(matrix, normalizerType, preRed, doublePrecision, fact, blocksize);
}

//*=============================================================================================

bool MinkowskiReduction(BMat & matrix, PreReductionType preRed, PrecisionType doublePrecision,
						double fact, int blocksize)
{
	int dimension;
	if (matrix.size1() != matrix.size2()) {
		MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
		exit(1);
		// C'est pas un peu nul ça ?
	} else 
		dimension = matrix.size1();

	// creating objects needed to perform the test
	IntLatticeBasis basis (matrix, dimension, L2NORM);
	Reducer red (basis);

	// performing pre-reduction
	red.preRedDieter(0);
	// PW_TODO: check if possible to use other (better) prereduction?

	// performing the Minkowski reduction. Returns *false* if the algorithm didn't terminated well, 
	// returns *true* if it a success.
	bool reductionSuccess = red.reductMinkowski (0);
	matrix = red.getIntLatticeBasis().getBasis();
	return reductionSuccess;
}

//*=============================================================================================

bool MinkowskiReduction(BMat & matrix, long maxNodesBB, PreReductionType preRed, 
						PrecisionType doublePrecision, double fact, int blocksize)
{
	Reducer::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
	return MinkowskiReduction(matrix, preRed, doublePrecision, fact, blocksize);
}

//*=============================================================================================

double FigureOfMeritBeyer(BMat matrix, PreReductionType preRed, PrecisionType doublePrecision, 
							double fact, int blocksize)
{
	int dimension;
	if (matrix.size1() != matrix.size2()) {
		MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
		exit(1);
		// C'est pas un peu nul ça ?
	} else 
		dimension = matrix.size1();

	bool reductionSuccess;
	//m_lat->dualize ();
   reductionSuccess = MinkowskiReduction(matrix, preRed, doublePrecision, fact, blocksize);
	//m_lat->dualize ();

	//if (m_dualF)
	//m_lat->dualize ();

	if (reductionSuccess) {
		IntLatticeBasis basis (matrix, dimension, L2NORM);
		basis.updateScalL2Norm (0);
		basis.updateScalL2Norm (dimension-1);
		double x1, x2; // maybe using something else than double (xdouble, RR?)
		conv (x1, basis.getVecNorm (0));
		conv (x2, basis.getVecNorm (dimension-1));
		return sqrt(x1 / x2);
	} else 
		return -1.0;
}

//*=============================================================================================

double FigureOfMeritBeyer(BMat matrix, long maxNodesBB, PreReductionType preRed,
							PrecisionType doublePrecision, double fact, int blocksize)
{
	Reducer::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
   return FigureOfMeritBeyer(matrix, preRed, doublePrecision, fact, blocksize);
}

//*=============================================================================================

} // end namespace LatticeTester

