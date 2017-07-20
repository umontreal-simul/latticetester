#ifndef LATTICETESTERROUTINES_H
#define LATTICETESTERROUTINES_H

#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/LatticeAnalysis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Reducer.h"

//PW_TODO vérifier que certains ne sont pas appelés inutilement

namespace LatticeTester {


// returns -1.0 if there was an error in redBB0.
double ShortestVector(BMat matrix, NormType norm, PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocsize = 20);
double ShortestVector(BMat matrix, NormType norm, int maxNodesBB, PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocsize = 20);

// Minkowki reduction of the basis


// Figure of Merit 	


} // end namespace LatticeTester

#endif