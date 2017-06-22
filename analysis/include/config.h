//
//  Tools.h
//  LatticeTester
//
//  Created by Erwan Bourceret on 21/06/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//

#ifndef CONFIG_H
#define CONFIG_H

#include "Tools.h"

// projection parameters definition
//----------------------------------------------------------------------------------------

/*
 * if true the output is written in a .txt file
 * if false the ouput is printed in the console
 * via the std::outFile command
 */
bool outFileRequested = false;


/*
 * Type of the basis input.
 * If true, all entry of the basis will be random. In that case,
 * the Dual can't be compute.
 * If false, the basis will be genereated from a MRG with random
 * parameters a = (a1, a2, ..., ak) and a modulus configurable below
 * and according to L'Ecuyer's publications.
 */
bool FullRandomMatrix = false;

/*
 * Select the interval of value for the full random basis.
 * The value will be chosen between minCoeff and maxCoeff.
 * FullRandomMatrix flag must be true.
 */
const int minCoeff = 40;
const int maxCoeff = 1000;

/*
 * All Reducer method, except redDieter, can be used
 * without the computation of the Dual. In that case,
 * we save memory and time on average, and the result is
 * the same. But, for a high period (about 2^40), the
 * Cholesky decomposition computation needs very large number.
 * Thus, the calcul without Dual can drive sometimes to a
 * longer commputing timing.
 * FullRandomMatrix flag must be false.
 */
bool WITH_DUAL = false;

/*
 * Order of the Basis according to L'Écuyer's paper.
 * FullRandomMatrix flag must be false.
 */
const int order = 3;

/*
 * Modulo of the Basis according to L'Écuyer's paper.
 * FullRandomMatrix flag must be false.
 */
const ZZ modulusRNG = power_ZZ(2, 19) - 1;

/*
 * The Dimension to be analysed.
 * Must be int value.
 */
const int dimension = 22;


/*
 * Iteration loop over basis of same dimension.
 * Each random basis is computed with a different seed.
 */
const int maxIteration = 10;

/*
 * a/b is the value of the delta in the LLL and BKZ
 * Reduction. NTL offers the possibility to compute
 * a LLL Reduction with the exact delta. We have noticed
 * only minor differences with this option.
 */
const double delta = 0.99999;
const double epsilon = 1.0 - delta;

/*
 * Reduction bound in the RedLLL algorithm.
 */
const int maxcpt = 10000000;

/*
 * Block Size in the BKZ algorithm. See NTL documention
 * for further information.
 */
const long blocksize = 20;

/*
 * Maximum number of Nodes in the Branch-and-Bound.
 */
const int maxNodesBB = 1000000;

/*
 * Selecting method of reducing.
 */
ReduceType Reduce_type[] ={
   Initial,                   // Calibration
   PairRed,                   // Performs Pairwise Reduction
   PairRedRandomized,         // Performs stochastic Pairwise Reduction
   RedLLL,                       // Performs LLL Reduction
   PairRed_LLL,               // Performs Pairwise Reduction
   PairRedRandomized_LLL,     // Perform Pairwise Reduction and then
   // LLL Reduction
   LLLNTL,                    // Performs LLL Reduction with NTL Library
   PairRed_LLLNTL,            // Perform Pairwise Reduction and then
   // LLL Reduction with NTL Library
   PairRedRandomized_LLLNTL,  // Perform stocastic Pairwise Reduction and
   // then LLL Reduction with NTL Library
   BKZNTL,                    // Performs BKZ Reduction with NTL Library
   PairRed_BKZNTL,            // Perform Pairwise Reduction and then
   // BKZ Reduction with NTL Library
   PairRedRandomized_BKZNTL,  // Perform stocastic Pairwise Reduction and
   // then BKZ Reduction with NTL Library
   //BB_Only,                   // Performs Branch-and-Bound Reduction
   BB_Classic,                // Perform Pairwise Reduction and then
   // Branch-and-Bound Reduction
   BB_BKZ                     // Performs BKZ Reduction with NTL Library
   // and then Branch-and-Bound Reduction
   //RedDieter,                    // Performs Dieter Reduction
   //WARNING USE DIETER ONLY FOR DIM < 6
   //RedMinkowski                  // Perform Minkowski Reduction with
   // Branch-and-Bound Algorithm.
};


#endif /* config_h */
