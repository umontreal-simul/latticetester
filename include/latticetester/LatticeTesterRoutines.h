


#ifndef LATTICETESTERROUTINES_H
#define LATTICETESTERROUTINES_H

#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/LatticeAnalysis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaRogers.h"

#include "latticetester/Reducer.h"

namespace LatticeTester {


//PW_TODO: ou est le modulo ????????

/**
 * This function allows computation of the shortest non-zero vector in a lattice,
 * according to the selected norm. Many parameters can bet set by the user, otherwise
 * the function work with default values.
 * Returns -1.0 if there was an error in Branch-and-Bound procedure. Return the length
 * of the shortest non-zero vector otherwise.
 */
double ShortestVector(BMat matrix, NormType norm, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

/**
 * Same thing as before but with the possibility to set a different value for
 * the variable maxNodesBB.
 */
double ShortestVector(BMat matrix, NormType norm, long maxNodesBB,
	PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE,
	double fact = 0.999, int blocksize = 20);

/**
 * This function compute the Figure of Merit to a given matrix, according to a
 * normalization criteria. It first computes the shortest non-zero vector using the
 * above functions. It then normalizes this value.
 * Returns -1.0 if there was an error in Branch-and-Bound procedure while calculating
 * the length of shortest non-zero vector. Return the figure of merit otherwise.
 */
double FigureOfMerit(BMat matrix, NormaType normalizerType, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

/**
 * Same thing as before but with the possibility to set a different value for
 * the variable maxNodesBB.
 */
double FigureOfMerit(BMat matrix, NormaType normalizerType, long maxNodesBB,
	PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE,
	double fact = 0.999, int blocksize = 20);

/**
 * This function reduces a basis to a Minkowski-reduced basis. Such basis has strong
 * properties regarding the length of its vectors but will require a huge running time,
 * especially when the dimension of the basis increases. Such Minkowski-reduced basis
 * is usefull, for example, to calculate a Beyer quotient (as implemented in
 * FigureOfMerit()).
 */
bool MinkowskiReduction(BMat & matrix, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

/**
 * Same thing as before but with the possibility to set a different value for
 * the variable maxNodesBB.
 */
bool MinkowskiReduction(BMat & matrix, long maxNodesBB, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

/**
 * This function compute the Figure of Merit to a given matrix, according to the
 * Beyer criteria. It first computes the Minkowski-reduced basis of the lattice 
 * and then makes the quotient of shortest over longest vector.
 * Returns -1.0 if there was an error in Branch-and-Bound procedure while calculating
 * the Minkowski-reduced basis. Return the figure of merit otherwise.
 */
double FigureOfMeritBeyer(BMat matrix, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

/**
 * Same thing as before but with the possibility to set a different value for
 * the variable maxNodesBB.
 */
double FigureOfMeritBeyer(BMat matrix, long maxNodesBB, PreReductionType preRed = BKZ,
	PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20);

} // end namespace LatticeTester

#endif