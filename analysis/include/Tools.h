//
//  Tools.h
//  LatticeTester
//
//  Created by Erwan Bourceret on 21/06/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//

#ifndef Tools_h
#define Tools_h

#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Util.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Reducer.h"

enum ReduceType {
   Initial,                    // No reduction applied, use for calibration
   PairRed,                    // Performs Pairwise Reduction
   PairRedRandomized,          // Performs stochastic Pairwise Reduction
   RedLLL,                        // Performs LLL Reduction
   PairRed_LLL,                // Performs Pairwise Reduction
   PairRedRandomized_LLL,      // Perform Pairwise Reduction and then
                               // LLL Reduction
   LLLNTLProxy,                // Performs proxy LLL Reduction with NTL Library
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

   case LLLNTLProxy :
         return "LLLNTLProxy";

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



mat_ZZ RandomMatrix (int dim, long min, long max, int seed)
{
   mat_ZZ basis;
   basis.SetDims(dim,dim);
   ZZ det;

   srand(seed);

   do{
       for (int i = 0; i < dim; i++){
           for (int j = 0; j < dim; j++)
               basis[i][j] = min + (rand() % (max - min + 1));
       }
       det = determinant(basis);

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

mat_ZZ CreateRNGBasis (const ZZ modulus, const int order, const int dimension, ZZ seed)
{
   mat_ZZ basis;
   basis.SetDims(dimension,dimension);

   if (dimension < order+1) {
      // degenerate case: identity matrix only
      for (int i = 0; i < dimension; i++)
         basis[i][i] = 1;

   } else { //usual case

      // (a_i) coefficients
      vec_ZZ a;
      a = randomVector(order, modulus, seed);

      for (int i = 0; i < order; i++) {
         // left upper block
         basis[i][i] = 1;

         //right upper block
         vec_ZZ initialState;
         initialState = canonicVector(order, i);
         SimpleMRG myMRG (modulus, order, a, initialState);

         for (int l = order; l < dimension; l++)
            basis[i][l] = conv<ZZ>(myMRG.getNextValue());
      }

      // right lower block
      for (int i = order; i < dimension; i++)
         basis[i][i] = modulus;

   } // end if

   return basis;
}

mat_ZZ Dualize (const mat_ZZ V, const ZZ modulus, const int k)
{
   mat_ZZ W;
   W.SetDims(V.NumRows(), V.NumRows());

   transpose(W,-V);
   long rmax = k;
   if(k>V.NumRows()){ rmax = V.NumRows(); }
   for (int i = 0; i < rmax; i++)
      W[i][i] = modulus;
   for (int i = k; i < V.NumRows(); i++)
      W[i][i] = 1;

   return W;
}


#endif /* Tools_h */
