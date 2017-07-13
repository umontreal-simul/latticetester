//
//  Tools.h
//  LatticeTester
//
//  Created by Erwan Bourceret on 21/06/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//

#ifndef TOOLS_H
#define TOOLS_H

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Reducer.h"
#include "SimpleMRG.h"

enum ReduceType {
   Initial,                    // No reduction applied, use for calibration
   PairRed,                    // Performs Pairwise Reduction
   PairRedRandomized,          // Performs stochastic Pairwise Reduction
   RedLLL,                        // Performs LLL Reduction
   PairRed_LLL,                // Performs Pairwise Reduction
   PairRedRandomized_LLL,      // Perform Pairwise Reduction and then
                               // LLL Reduction
   LLLNTLProxyFP,              // Performs proxy LLL Reduction with NTL Library
                               // With double precision
   LLLNTLProxyRR,              // Performs proxy LLL Reduction with NTL Library
                               // With arbitrary precision
   LLLNTLExact,                // Performs exact LLL Reduction with NTL Library
   PairRed_LLLNTL,             // Perform Pairwise Reduction and then
                               // LLL Reduction with NTL Library
   PairRedRandomized_LLLNTL,   // Perform stocastic Pairwise Reduction and
                               // then LLL Reduction with NTL Library
   BKZNTL,                     // Performs BKZ Reduction with NTL Library
   PairRed_BKZNTL,             // Perform Pairwise Reduction and then
                               // BKZ Reduction with NTL Library
   PairRedRandomized_BKZNTL,   // Perform stocastic Pairwise Reduction and
                               // then BKZ Reduction with NTL Library
   BB_Only,                    // Performs Branch-and-Bound Reduction
   LLL_BB,                     // Perform LLL Reduction and then
                               // Branch-and-Bound Reduction
   PreRedDieter_BB,            // Perform Prereddieter Reduction and then
                               // Branch-and-Bound Reduction
   PreRedDieter_LLL_BB,        // Perform Dieter pre-Reduction, LLL reduction
                               // and then Branch-and-Bound Reduction
   BKZ_BB,                     // Performs BKZ Reduction with NTL Library
                               // and then Branch-and-Bound Reduction
   RedDieter,                     // Performs Dieter Reduction
                               //WARNING USE DIETER ONLY FOR DIM < 6
   RedMinkowski                   // Perform Minkowski Reduction with
                               // Branch-and-Bound Algorithm.
};

string toStringReduce (ReduceType reduce)
{
   switch(reduce) {
   case PairRed :
         return "PairRed";

   case PairRedRandomized :
         return "PairRedRandomized";

   case RedLLL :
         return "RedLLL";

   case PairRed_LLL :
         return "PairRed_LLL";

   case PairRedRandomized_LLL :
         return "PairRedRandomized_LLL";

   case LLLNTLProxyRR :
         return "LLLNTLProxy_RR";

   case LLLNTLProxyFP :
         return "LLLNTLProxy_FP";

   case LLLNTLExact :
         return "LLLNTLExact";

   case PairRed_LLLNTL :
         return "PairRed_LLLNTL";

   case PairRedRandomized_LLLNTL :
         return "PairRedRandomized_LLLNTL";

   case BKZNTL :
         return "BKZNTL";

   case PairRed_BKZNTL :
         return "PairRed_BKZNTL";

   case PairRedRandomized_BKZNTL :
         return "PairRedRandomized_BKZNTL";

   case LLL_BB :
         return "LLL_BB";

   case PreRedDieter_BB :
         return "PreRedDieter_BB";

   case PreRedDieter_LLL_BB :
         return "PreRedDieter_LLL_BB";

   case BB_Only :
         return "BB_Only";

   case BKZ_BB :
         return "BKZ_BB";

   case RedDieter :
         return "RedDieter";

   case RedMinkowski :
         return "RedMinkowski";

   case Initial :
         return "Initial";

   default :
         return "NOVALUE";
   }
}



BMat RandomMatrix (int dim, long min, long max, int seed)
{
   BMat basis;
   basis.SetDims(dim,dim);
   BScal det;

   srand(seed);

   do{
       for (int i = 0; i < dim; i++){
           for (int j = 0; j < dim; j++)
               basis[i][j] = min + (rand() % (max - min + 1));
       }

#if NTL_TYPES_CODE != 1
   det = determinant(basis);
#else
   // As NTL library does not support matrix with double
   // we compute the determinant with the boost library
   boost::numeric::ublas::matrix<long>  mat_tmps;

   mat_tmps.resize(dim, dim);
   for(unsigned int i = 0; i < dim; i++){
      for(unsigned int j = 0; j < dim; i++){
         mat_tmps(i,j) = basis[i][j];
      }
   }
   det = det_double(mat_tmps);
   //RScal logDensity(-log(abs(10000)));
#endif

   } while ( det == 0 );
   return basis;
}

/*
 * Compute the mean of the array of Type with the size Size.
 */
template<typename Type, unsigned long Size>
Type mean(const array <Type, Size> array) {
   Type sum (0);
   for(unsigned int i = 0; i<Size; i++)
       sum += array[i];
   return sum / Size;
}

/*
 * Compute the variance of the array of Type with the size Size.
 */
template<typename Type, unsigned long Size>
Type variance(const array <Type, Size> array) {
   Type sum (0);
   Type mean_tmp(mean(array));
   for(unsigned int i = 0; i<Size; i++)
       sum += (array[i] - mean_tmp) * (array[i] - mean_tmp);
   return sum / Size;
}


vec_ZZ canonicVector (int dimension, int i)
{
   vec_ZZ e;
   e.SetLength(dimension);
   e[i] = 1;
   return e;
}

vec_ZZ randomVector (int dimension, ZZ modulus, ZZ seed)
{
   vec_ZZ vector;
   vector.SetLength(dimension);
   SetSeed(seed);
   for (int i = 0; i < dimension; i++)
      vector[i] = RandomBnd(modulus);
   return vector;
}

BMat CreateRNGBasis (const ZZ modulus, const int order, const int dimension, vec_ZZ& a)
{
   BMat basis;
   basis.SetDims(dimension,dimension);

   for(unsigned int i = 0; i < dimension; i++){
      for(unsigned int j = 0; j < dimension; j++){
         basis[i][j] = 0;
      }
   }

   if (dimension < order+1) {
      // degenerate case: identity matrix only
      for (int i = 0; i < dimension; i++)
         basis[i][i] = 1;

   } else { //usual case

      // (a_i) coefficients

      for (int i = 0; i < order; i++) {
         // left upper block
         basis[i][i] = 1;

         //right upper block
         vec_ZZ initialState;
         initialState = canonicVector(order, i);
         SimpleMRG myMRG (modulus, order, a, initialState);

         for (int l = order; l < dimension; l++){
            ZZ_p nb_value(myMRG.getNextValue());
            basis[i][l] = conv<BScal>(myMRG.getNextValue());
         }
      }

      // right lower block
      for (int i = order; i < dimension; i++)
         conv(basis[i][i],modulus);

   } // end if

   return basis;
}

BMat Dualize (const BMat V, const ZZ modulus, const int k)
{
   BMat W;
   W.resize(V.NumRows(), V.NumRows());
#if NTL_TYPES_CODE != 1
   transpose(W,-V);
#else
   for(unsigned int i = 0; i < V.NumRows(); i++){
      for(unsigned int j = 0; j < V.NumRows(); j++){
         W[i][j] = -V[j][i];
      }
   }
#endif
   long rmax = k;
   if(k>V.NumRows()){ rmax = V.NumRows(); }
   for (int i = 0; i < rmax; i++)
      conv(W[i][i],modulus);
   for (int i = k; i < V.NumRows(); i++)
      W[i][i] = 1;
   return W;
}

#endif /* Tools_h */
