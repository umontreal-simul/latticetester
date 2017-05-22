//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "latticetester/Util.h"
#include "latticetester/Basis.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"

#include <NTL/ctools.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include "NTL/vec_ZZ.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <RInside.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

#define PRINT_CONSOLE

using namespace std;
using namespace LatticeTester;

const int MaxDimension = 15;

#ifdef PRINT_CONSOLE
const int MinDimension = MaxDimension - 1;
#else
const int MinDimension = 5;
#endif


const int Interval_dim = MaxDimension - MinDimension;

template <typename type>
Rcpp::NumericMatrix toRcppMatrix(const type scal[][Interval_dim], const int & maxIteration)
{

   Rcpp::NumericMatrix mat(maxIteration, Interval_dim);
   for (int i = 0; i<maxIteration ; i++) {
      for (int j = 0; j<Interval_dim; j++) {
         conv(mat(i,j), scal[i][j]);
      }
   }
   return mat;
}



void RandomMatrix (mat_ZZ& A, ZZ& det, int min, int max, int seed){

   int dim = (int) A.NumRows() ;
   srand(seed);

   do{
       for (int i = 0; i < dim; i++){
           for (int j = 0; j < dim; j++)
               A[i][j] = min + (rand() % (int)(max - min + 1));
       }
       det = determinant(A);

   } while ( det == 0 );

}


template<typename Type, long Size>
void print(string name, Type const(& array)[Size], bool isIntegerOutput) {
   cout << name << " = ";
   for(int i=0; i<Size; i++){
       if (isIntegerOutput)
           cout << conv<ZZ>(array[i]) << " ";
       else
           cout << array[i] << " ";
   }
   //cout << endl;
}


template<typename Type, long Size>
Type Average(Type const(& array)[Size][Interval_dim]) {
   Type sum (0);
   for(int i=0; i<Size; i++)
       sum += array[i][0];
   return sum / Size;
}


//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main (int argc, char *argv[])
{

   // main parameters for the test
   int min = 30;
   int max = 100;


   long a = 999999;
   long b = 1000000;
   double delta = (double) a/b;
   double epsilon = 1.0 - delta;

   int maxcpt = 1000000;
   int d = 0;
   long blocksize = 20; // for BKZ insertions

   // iteration loop over matrices of same dimension
   const int maxIteration = 10;


   // print important information
   bool printMatricesDetails = false;
   cout << "epsilon = " << epsilon << endl;
   //cout << "dimension = " << dimension << endl;
   cout << "nombre de matrices testées = " << maxIteration << endl;
   cout << "dimension minimale : " << MinDimension << endl;
   cout << "dimension maximale : " << MaxDimension << endl;
   cout << endl;

   // to display progress bar
   boost::progress_display show_progress(5*maxIteration*Interval_dim);

   // arrays to store values

   double timing_PairRedPrimal [maxIteration][Interval_dim];
   double timing_PairRedPrimalRandomized [maxIteration][Interval_dim];

   double timing_LLL [maxIteration][Interval_dim];
   double timing_LLL_PairRedPrimal [maxIteration][Interval_dim];
   double timing_LLL_PostPairRedPrimal [maxIteration][Interval_dim];
   double timing_LLL_PairRedPrimalRandomized [maxIteration][Interval_dim];
   double timing_LLL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];

   double timing_LLLNTL [maxIteration][Interval_dim];
   double timing_LLLNTL_PairRedPrimal [maxIteration][Interval_dim];
   double timing_LLLNTL_PostPairRedPrimal [maxIteration][Interval_dim];
   double timing_LLLNTL_PairRedPrimalRandomized [maxIteration][Interval_dim];
   double timing_LLLNTL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];

   //double timing_LLL_NTL_Exact [maxIteration][Interval_dim];
   double timing_BB_Only [maxIteration][Interval_dim];
   double timing_BB_Classic1 [maxIteration][Interval_dim];
   double timing_BB_Classic2 [maxIteration][Interval_dim];
   double timing_BB_BKZ1 [maxIteration][Interval_dim];
   double timing_BB_BKZ2 [maxIteration][Interval_dim];

   double timing_BKZNTL [maxIteration][Interval_dim];
   double timing_BKZNTL_PairRedPrimal [maxIteration][Interval_dim];
   double timing_BKZNTL_PostPairRedPrimal [maxIteration][Interval_dim];
   double timing_BKZNTL_PairRedPrimalRandomized [maxIteration][Interval_dim];
   double timing_BKZNTL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];



   NScal length_Initial [maxIteration][Interval_dim];
   NScal length_PairRedPrimal [maxIteration][Interval_dim];
   NScal length_PairRedPrimalRandomized [maxIteration][Interval_dim];

   NScal length_LLL [maxIteration][Interval_dim];
   NScal length_LLL_PostPairRedPrimal [maxIteration][Interval_dim];
   NScal length_LLL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];

   NScal length_LLLNTL [maxIteration][Interval_dim];
   NScal length_LLLNTL_PostPairRedPrimal [maxIteration][Interval_dim];
   NScal length_LLLNTL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];

   NScal length_BB_Only [maxIteration][Interval_dim];
   NScal length_BB_Classic [maxIteration][Interval_dim];
   NScal length_BB_BKZ [maxIteration][Interval_dim];

   //NScal length_LLL_NTL_Exact [maxIteration][Interval_dim];

   NScal length_BKZNTL [maxIteration][Interval_dim];
   NScal length_BKZNTL_PostPairRedPrimal [maxIteration][Interval_dim];
   NScal length_BKZNTL_PostPairRedPrimalRandomized [maxIteration][Interval_dim];

   for (int dimension = MinDimension; dimension < MaxDimension; dimension++){

      for (int iteration = 0; iteration < maxIteration; iteration++){


         int seed = (iteration+1) * (iteration+1) * 123456789 * dimension;
         int seed_dieter = (iteration+1) * dimension * 12342;
         //int seed = (int) (iteration+1) * 12345 * time(NULL);

         int idx = dimension - MinDimension; // Stock the indices for table

         // We create copies of the same basis
         BMat basis_PairRedPrimal (dimension, dimension);
         ZZ det;
         RandomMatrix(basis_PairRedPrimal, det, min, max, seed);


         BMat basis_PairRedPrimalRandomized (basis_PairRedPrimal);
         BMat basis_LLL (basis_PairRedPrimal);
         BMat basis_PairRedPrimal_LLL (basis_PairRedPrimal);
         BMat basis_PairRedPrimalRandomized_LLL (basis_PairRedPrimal);
         BMat basis_LLLNTL (basis_PairRedPrimal);
         BMat basis_PairRedPrimal_LLLNTL (basis_PairRedPrimal);
         BMat basis_PairRedPrimalRandomized_LLLNTL (basis_PairRedPrimal);
         //BMat basis_LLLNTL_Exact (basis_PairRedPrimal);
         BMat basis_BKZNTL (basis_PairRedPrimal);
         BMat basis_PairRedPrimal_BKZNTL (basis_PairRedPrimal);
         BMat basis_PairRedPrimalRandomized_BKZNTL (basis_PairRedPrimal);
         BMat basis_BB_Only (basis_PairRedPrimal);
         BMat basis_BB_Classic (basis_PairRedPrimal);
         BMat basis_BB_BKZ (basis_PairRedPrimal);

         IntLatticeBasis lattice_PairRedPrimal (basis_PairRedPrimal,dimension);
         IntLatticeBasis lattice_PairRedPrimalRandomized (basis_PairRedPrimalRandomized,dimension);
         IntLatticeBasis lattice_LLL (basis_LLL, dimension);
         IntLatticeBasis lattice_PairRedPrimal_LLL (basis_PairRedPrimal_LLL, dimension);
         IntLatticeBasis lattice_PairRedPrimalRandomized_LLL (basis_PairRedPrimalRandomized_LLL, dimension);
         IntLatticeBasis lattice_LLLNTL (basis_LLLNTL,dimension);
         IntLatticeBasis lattice_PairRedPrimal_LLLNTL (basis_PairRedPrimal_LLLNTL,dimension);
         IntLatticeBasis lattice_PairRedPrimalRandomized_LLLNTL (basis_PairRedPrimalRandomized_LLLNTL,dimension);
         //IntLatticeBasis lattice_LLLNTL_Exact (basis_LLLNTL_Exact,dimension);
         IntLatticeBasis lattice_BKZNTL (basis_BKZNTL,dimension);
         IntLatticeBasis lattice_PairRedPrimal_BKZNTL (basis_PairRedPrimal_BKZNTL,dimension);
         IntLatticeBasis lattice_PairRedPrimalRandomized_BKZNTL (basis_PairRedPrimalRandomized_BKZNTL,dimension);
         IntLatticeBasis lattice_BB_Only (basis_BB_Only,dimension);
         IntLatticeBasis lattice_BB_Classic (basis_BB_Classic,dimension);
         IntLatticeBasis lattice_BB_BKZ (basis_BB_BKZ, dimension);


         lattice_PairRedPrimal.setNegativeNorm(true);
         lattice_PairRedPrimal.updateVecNorm();
         lattice_PairRedPrimal.sort(0);
         NScal initialShortestVectorLength = lattice_PairRedPrimal.getVecNorm(0);
         length_Initial [iteration][idx] = initialShortestVectorLength;

         if (printMatricesDetails) {
           cout << "\n*** Initial basis ***" << endl;
           cout << "det = " << det << endl;
           cout << "Shortest vector = " << initialShortestVectorLength << endl;
           lattice_PairRedPrimal.write();
         }

         Reducer reducer_PairRedPrimal (lattice_PairRedPrimal);
         Reducer reducer_PairRedPrimalRandomized (lattice_PairRedPrimalRandomized);
         Reducer reducer_LLL (lattice_LLL);
         Reducer reducer_PairRedPrimal_LLL (lattice_PairRedPrimal_LLL);
         Reducer reducer_PairRedPrimalRandomized_LLL (lattice_PairRedPrimalRandomized_LLL);
         Reducer reducer_LLLNTL (lattice_LLLNTL);
         Reducer reducer_PairRedPrimal_LLLNTL (lattice_PairRedPrimal_LLLNTL);
         Reducer reducer_PairRedPrimalRandomized_LLLNTL (lattice_PairRedPrimalRandomized_LLLNTL);
         //Reducer reducer_LLLNTL_Exact (lattice_LLLNTL_Exact);
         Reducer reducer_BKZNTL (lattice_BKZNTL);
         Reducer reducer_PairRedPrimal_BKZNTL (lattice_PairRedPrimal_BKZNTL);
         Reducer reducer_PairRedPrimalRandomized_BKZNTL (lattice_PairRedPrimalRandomized_BKZNTL);
         Reducer reducer_BB_Only (lattice_BB_Only);
         Reducer reducer_BB_Classic (lattice_BB_Classic);
         Reducer reducer_BB_BKZ (lattice_BB_BKZ);


         //------------------------------------------------------------------------------------
         // Pairwise reduction in primal basis only
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimal = clock();
         reducer_PairRedPrimal.preRedDieterPrimalOnly(d);
         clock_t end_PairRedPrimal = clock();

         lattice_PairRedPrimal.setNegativeNorm(true);
         lattice_PairRedPrimal.updateVecNorm();
         lattice_PairRedPrimal.sort(0);

         if (printMatricesDetails) {
           cout << "*** Pairwise reduction in primal basis only ***" << endl;
           cout << "Shortest vector = ";
           cout << lattice_PairRedPrimal.getVecNorm(0) << endl;
           lattice_PairRedPrimal.write();
         }
         ++show_progress;


         //------------------------------------------------------------------------------------
         // Randomized pairwise reduction in primal basis only
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimalRandomized = clock();
         reducer_PairRedPrimalRandomized.preRedDieterPrimalOnlyRandomized(d, seed_dieter);
         clock_t end_PairRedPrimalRandomized = clock();

         lattice_PairRedPrimalRandomized.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized.updateVecNorm();
         lattice_PairRedPrimalRandomized.sort(0);

         if (printMatricesDetails) {
           cout << "*** Randomized pairwise reduction in primal basis only ***" << endl;
           cout << "Shortest vector = ";
           cout << lattice_PairRedPrimalRandomized.getVecNorm(0) << endl;
           lattice_PairRedPrimalRandomized.write();
         }

         //------------------------------------------------------------------------------------
         // LLL Richard
         //------------------------------------------------------------------------------------

         clock_t begin_LLL = clock();
         reducer_LLL.redLLL(delta, maxcpt, dimension);
         clock_t end_LLL = clock();

         lattice_LLL.setNegativeNorm(true);
         lattice_LLL.updateVecNorm();
         lattice_LLL.sort(0);

         if (printMatricesDetails) {
           cout << "*** LLL only ***" << endl;
           cout << "Shortest vector = " << lattice_LLL.getVecNorm(0) << endl;
           lattice_LLL.write();
         }

         ++show_progress;

         //------------------------------------------------------------------------------------
         // Pairwise reduction (in primal basis only) and then LLL Richard
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimal_LLL1 = clock();
         reducer_PairRedPrimal_LLL.preRedDieterPrimalOnly(d);
         clock_t end_PairRedPrimal_LLL1 = clock();

         lattice_PairRedPrimal_LLL.setNegativeNorm(true);
         lattice_PairRedPrimal_LLL.updateVecNorm();
         lattice_PairRedPrimal_LLL.sort(0);

         clock_t begin_PairRedPrimal_LLL2 = clock();
         reducer_PairRedPrimal_LLL.redLLL(delta, maxcpt, dimension);
         clock_t end_PairRedPrimal_LLL2 = clock();

         lattice_PairRedPrimal_LLL.setNegativeNorm(true);
         lattice_PairRedPrimal_LLL.updateVecNorm();
         lattice_PairRedPrimal_LLL.sort(0);

         if (printMatricesDetails){
           cout << "*** Pairwise reduction in primal and LLL ***" << endl;
           cout << "Shortest vector = " << lattice_PairRedPrimal_LLL.getVecNorm(0) << endl;
           lattice_PairRedPrimal_LLL.write();
         }


         //------------------------------------------------------------------------------------
         // Randomized pairwise reduction (in primal basis only) and then LLL Richard
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimalRandomized_LLL1 = clock();
         reducer_PairRedPrimalRandomized_LLL.preRedDieterPrimalOnlyRandomized(d, seed_dieter);
         clock_t end_PairRedPrimalRandomized_LLL1 = clock();

         lattice_PairRedPrimalRandomized_LLL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_LLL.updateVecNorm();
         lattice_PairRedPrimalRandomized_LLL.sort(0);

         ++show_progress;

         clock_t begin_PairRedPrimalRandomized_LLL2 = clock();
         reducer_PairRedPrimalRandomized_LLL.redLLL(delta, maxcpt, dimension);
         clock_t end_PairRedPrimalRandomized_LLL2 = clock();

         lattice_PairRedPrimalRandomized_LLL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_LLL.updateVecNorm();
         lattice_PairRedPrimalRandomized_LLL.sort(0);


         if (printMatricesDetails){
           cout << "*** Randomized pairwise reduction in primal and LLL ***" << endl;
           cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_LLL.getVecNorm(0) << endl;
           lattice_PairRedPrimalRandomized_LLL.write();
         }

         //------------------------------------------------------------------------------------
         // LLL NTL reduction (floating point version = proxy)
         //------------------------------------------------------------------------------------

         clock_t begin_LLLNTL = clock();
         reducer_LLLNTL.redLLLNTLProxy(delta);
         clock_t end_LLLNTL = clock();

         lattice_LLLNTL.setNegativeNorm(true);
         lattice_LLLNTL.updateVecNorm();
         lattice_LLLNTL.sort(0);

         if (printMatricesDetails) {
           cout << "*** LLL NTL Proxy only ***" << endl;
           cout << "Shortest vector = " << lattice_LLLNTL.getVecNorm(0) << endl;
           lattice_LLLNTL.write();
         }


         //------------------------------------------------------------------------------------
         // Pairwise reduction (in primal basis only) and then LLL NTL proxy
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimal_LLLNTL1 = clock();
         reducer_PairRedPrimal_LLLNTL.preRedDieterPrimalOnly(d);
         clock_t end_PairRedPrimal_LLLNTL1 = clock();

         lattice_PairRedPrimal_LLLNTL.setNegativeNorm(true);
         lattice_PairRedPrimal_LLLNTL.updateVecNorm();
         lattice_PairRedPrimal_LLLNTL.sort(0);

         // useful ?
         //NScal intermediateLengthBis = lattice_PairRedPrimal_LLL_NTL.getVecNorm(0);

         clock_t begin_PairRedPrimal_LLLNTL2 = clock();
         reducer_PairRedPrimal_LLLNTL.redLLLNTLProxy(delta);
         clock_t end_PairRedPrimal_LLLNTL2 = clock();

         lattice_PairRedPrimal_LLLNTL.setNegativeNorm(true);
         lattice_PairRedPrimal_LLLNTL.updateVecNorm();
         lattice_PairRedPrimal_LLLNTL.sort(0);

         if (printMatricesDetails){
           cout << "*** Pairwise reduction in primal and LLL NTL ***" << endl;
           cout << "Shortest vector = " << lattice_PairRedPrimal_LLLNTL.getVecNorm(0) << endl;
           lattice_PairRedPrimal_LLLNTL.write();
         }


         //------------------------------------------------------------------------------------
         // Randomized pairwise reduction (in primal basis only) and then LLL NTL proxy
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimalRandomized_LLLNTL1 = clock();
         reducer_PairRedPrimalRandomized_LLLNTL.preRedDieterPrimalOnlyRandomized(d, seed_dieter);
         clock_t end_PairRedPrimalRandomized_LLLNTL1 = clock();


         lattice_PairRedPrimalRandomized_LLLNTL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_LLLNTL.updateVecNorm();
         lattice_PairRedPrimalRandomized_LLLNTL.sort(0);

         clock_t begin_PairRedPrimalRandomized_LLLNTL2 = clock();
         reducer_PairRedPrimalRandomized_LLLNTL.redLLLNTLProxy(delta);
         clock_t end_PairRedPrimalRandomized_LLLNTL2 = clock();

         lattice_PairRedPrimalRandomized_LLLNTL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_LLLNTL.updateVecNorm();
         lattice_PairRedPrimalRandomized_LLLNTL.sort(0);

         if (printMatricesDetails){
            cout << "*** Randomized pairwise reduction in primal and LLL NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_LLLNTL.getVecNorm(0) << endl;
            lattice_PairRedPrimalRandomized_LLLNTL.write();
         }

         ++show_progress;



         //------------------------------------------------------------------------------------
         // LLL NTL Exact reduction only
         //------------------------------------------------------------------------------------

         /*ZZ det2;
         clock_t begin_LLLNTL_Exact = clock();
         reducer_LLLNTL_Exact.redLLLNTLExact(det2, a, b);
         clock_t end_LLLNTL_Exact = clock();

         lattice_LLLNTL_Exact.setNegativeNorm(true);
         lattice_LLLNTL_Exact.updateVecNorm();
         lattice_LLLNTL_Exact.sort(0);

         if (printMatricesDetails) {
           cout << "*** LLL NTL Exact only ***" << endl;
           cout << "Shortest vector = " << lattice_LLLNTL_Exact.getVecNorm(0) << endl;
           lattice_LLLNTL_Exact.write();
         }*/


         //------------------------------------------------------------------------------------
         // BKZ NTL reduction
         //------------------------------------------------------------------------------------

         clock_t begin_BKZNTL = clock();
         reducer_BKZNTL.redBKZ(delta, blocksize);
         clock_t end_BKZNTL = clock();

         lattice_BKZNTL.setNegativeNorm(true);
         lattice_BKZNTL.updateVecNorm();
         lattice_BKZNTL.sort(0);

         if (printMatricesDetails) {
           cout << "*** BKZ NTL only ***" << endl;
           cout << "Shortest vector = " << lattice_BKZNTL.getVecNorm(0) << endl;
           lattice_BKZNTL.write();
         }

         //------------------------------------------------------------------------------------
         // Pairwise reduction (in primal basis only) and then BKZ NTL proxy
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimal_BKZNTL1 = clock();
         reducer_PairRedPrimal_BKZNTL.preRedDieterPrimalOnly(d);
         clock_t end_PairRedPrimal_BKZNTL1 = clock();

         lattice_PairRedPrimal_BKZNTL.setNegativeNorm(true);
         lattice_PairRedPrimal_BKZNTL.updateVecNorm();
         lattice_PairRedPrimal_BKZNTL.sort(0);

         // useful ?
         //NScal intermediateLengthBis = lattice_PairRedPrimal_LLL_NTL.getVecNorm(0);

         clock_t begin_PairRedPrimal_BKZNTL2 = clock();
         reducer_PairRedPrimal_BKZNTL.redBKZ(delta, blocksize);
         clock_t end_PairRedPrimal_BKZNTL2 = clock();

         lattice_PairRedPrimal_BKZNTL.setNegativeNorm(true);
         lattice_PairRedPrimal_BKZNTL.updateVecNorm();
         lattice_PairRedPrimal_BKZNTL.sort(0);

         if (printMatricesDetails){
           cout << "*** Pairwise reduction in primal and BKZ NTL ***" << endl;
           cout << "Shortest vector = " << lattice_PairRedPrimal_BKZNTL.getVecNorm(0) << endl;
           lattice_PairRedPrimal_BKZNTL.write();
         }

         lattice_PairRedPrimal_BKZNTL.setNegativeNorm(true);
         lattice_PairRedPrimal_BKZNTL.updateVecNorm();
         lattice_PairRedPrimal_BKZNTL.sort(0);


         //------------------------------------------------------------------------------------
         // Randomized pairwise reduction (in primal basis only) and then BKZ NTL proxy
         //------------------------------------------------------------------------------------

         clock_t begin_PairRedPrimalRandomized_BKZNTL1 = clock();
         reducer_PairRedPrimalRandomized_BKZNTL.preRedDieterPrimalOnlyRandomized(d, seed_dieter);
         clock_t end_PairRedPrimalRandomized_BKZNTL1 = clock();

         lattice_PairRedPrimalRandomized_BKZNTL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_BKZNTL.updateVecNorm();
         lattice_PairRedPrimalRandomized_BKZNTL.sort(0);

         NScal intermediateLengthRandomized_BKZNTL = lattice_PairRedPrimalRandomized_BKZNTL.getVecNorm(0);

         clock_t begin_PairRedPrimalRandomized_BKZNTL2 = clock();
         reducer_PairRedPrimalRandomized_BKZNTL.redBKZ(delta, blocksize);
         clock_t end_PairRedPrimalRandomized_BKZNTL2 = clock();
         if (printMatricesDetails){
            cout << "*** Randomized pairwise reduction in primal and BKZ NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_BKZNTL.getVecNorm(0) << endl;
            lattice_PairRedPrimalRandomized_BKZNTL.write();
         }


         lattice_PairRedPrimalRandomized_BKZNTL.setNegativeNorm(true);
         lattice_PairRedPrimalRandomized_BKZNTL.updateVecNorm();
         lattice_PairRedPrimalRandomized_BKZNTL.sort(0);

         //------------------------------------------------------------------------------------
         // Branch and Bound Only
         //------------------------------------------------------------------------------------
         clock_t begin_BB_Only = clock();
         reducer_BB_Only.shortestVector(L2NORM);
         clock_t end_BB_Only = clock();

         lattice_BB_Only.setNegativeNorm(true);
         lattice_BB_Only.updateVecNorm();
         lattice_BB_Only.sort(0);

         ++show_progress;

         //------------------------------------------------------------------------------------
         // Branch and Bound classic
         //------------------------------------------------------------------------------------

         // usual pre reduction heuristics
         clock_t begin_BB_Classic1 = clock();
         reducer_BB_Classic.preRedDieterPrimalOnly(d);
         reducer_BB_Classic.redLLL(delta, maxcpt, dimension);
         clock_t end_BB_Classic1 = clock();

         clock_t begin_BB_Classic2 = clock();
         reducer_BB_Classic.shortestVector(L2NORM);
         clock_t end_BB_Classic2 = clock();

         lattice_BB_Classic.setNegativeNorm(true);
         lattice_BB_Classic.updateVecNorm();
         lattice_BB_Classic.sort(0);


         //------------------------------------------------------------------------------------
         // Branch and Bound post BKZ
         //------------------------------------------------------------------------------------

         clock_t begin_BB_BKZ1 = clock();
         reducer_BB_BKZ.redBKZ(delta, blocksize);
         clock_t end_BB_BKZ1 = clock();

         clock_t begin_BB_BKZ2 = clock();
         reducer_BB_BKZ.shortestVector(L2NORM);
         clock_t end_BB_BKZ2 = clock();

         lattice_BB_BKZ.setNegativeNorm(true);
         lattice_BB_BKZ.updateVecNorm();
         lattice_BB_BKZ.sort(0);


         //------------------------------------------------------------------------------------
         // timing updating
         //------------------------------------------------------------------------------------


         double runningTime_PairRedPrimal = double (end_PairRedPrimal - begin_PairRedPrimal) / CLOCKS_PER_SEC;
         double runningTime_PairRedPrimalRandomized = double (end_PairRedPrimalRandomized - begin_PairRedPrimalRandomized) / CLOCKS_PER_SEC;

         double runningTime_LLL = double(end_LLL - begin_LLL) / CLOCKS_PER_SEC;
         double runningTime_LLL_PairRedPrimal = double(end_PairRedPrimal_LLL1 - begin_PairRedPrimal_LLL1) / CLOCKS_PER_SEC;
         double runningTime_LLL_PostPairRedPrimal = double(end_PairRedPrimal_LLL2 - begin_PairRedPrimal_LLL2) / CLOCKS_PER_SEC;
         double runningTime_LLL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_LLL1 - begin_PairRedPrimalRandomized_LLL1) / CLOCKS_PER_SEC;
         double runningTime_LLL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_LLL2 - begin_PairRedPrimalRandomized_LLL2) / CLOCKS_PER_SEC;

         double runningTime_LLLNTL = double(end_LLLNTL - begin_LLLNTL) / CLOCKS_PER_SEC;
         double runningTime_LLLNTL_PairRedPrimal = double(end_PairRedPrimal_LLLNTL1 - begin_PairRedPrimal_LLLNTL1) / CLOCKS_PER_SEC;
         double runningTime_LLLNTL_PostPairRedPrimal = double (end_PairRedPrimal_LLLNTL2 - begin_PairRedPrimal_LLLNTL2) / CLOCKS_PER_SEC;
         double runningTime_LLLNTL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_LLLNTL1 - begin_PairRedPrimalRandomized_LLLNTL1) / CLOCKS_PER_SEC;
         double runningTime_LLLNTL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_LLLNTL2 - begin_PairRedPrimalRandomized_LLLNTL2) / CLOCKS_PER_SEC;
         //double runningTime_LLLNTL_Exact = double(end_LLLNTL_Exact - begin_LLLNTL_Exact) / CLOCKS_PER_SEC;

         double runningTime_BKZNTL = double (end_BKZNTL - begin_BKZNTL) / CLOCKS_PER_SEC;
         double runningTime_BKZNTL_PairRedPrimal = double(end_PairRedPrimal_BKZNTL1 - begin_PairRedPrimal_BKZNTL1) / CLOCKS_PER_SEC;
         double runningTime_BKZNTL_PostPairRedPrimal = double (end_PairRedPrimal_BKZNTL2 - begin_PairRedPrimal_BKZNTL2) / CLOCKS_PER_SEC;
         double runningTime_BKZNTL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_BKZNTL1 - begin_PairRedPrimalRandomized_BKZNTL1) / CLOCKS_PER_SEC;
         double runningTime_BKZNTL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_BKZNTL2 - begin_PairRedPrimalRandomized_BKZNTL2) / CLOCKS_PER_SEC;

         double runningTime_BB_Only = double (end_BB_Only - begin_BB_Only) / CLOCKS_PER_SEC;
         double runningTime_BB_Classic1 = double (end_BB_Classic1 - begin_BB_Classic1) / CLOCKS_PER_SEC;
         double runningTime_BB_Classic2 = double (end_BB_Classic2 - begin_BB_Classic2) / CLOCKS_PER_SEC;
         double runningTime_BB_BKZ1 = double (end_BB_BKZ1 - begin_BB_BKZ1) / CLOCKS_PER_SEC;
         double runningTime_BB_BKZ2 = double (end_BB_BKZ2 - begin_BB_BKZ2) / CLOCKS_PER_SEC;

         //------------------------------------------------------------------------------------
         // timing and length arrays updating
         //------------------------------------------------------------------------------------
         timing_PairRedPrimal [iteration][idx] = runningTime_PairRedPrimal;
         timing_PairRedPrimalRandomized [iteration][idx] = runningTime_PairRedPrimalRandomized;

         timing_LLL [iteration][idx] = runningTime_LLL;
         timing_LLL_PairRedPrimal [iteration][idx] = runningTime_LLL_PairRedPrimal;
         timing_LLL_PostPairRedPrimal [iteration][idx] = runningTime_LLL_PostPairRedPrimal;
         timing_LLL_PairRedPrimalRandomized [iteration][idx] = runningTime_LLL_PairRedPrimalRandomized;
         timing_LLL_PostPairRedPrimalRandomized [iteration][idx] = runningTime_LLL_PostPairRedPrimalRandomized;

         timing_LLLNTL [iteration][idx] = runningTime_LLLNTL;
         timing_LLLNTL_PairRedPrimal [iteration][idx] = runningTime_LLLNTL_PairRedPrimal;
         timing_LLLNTL_PostPairRedPrimal [iteration][idx] = runningTime_LLLNTL_PostPairRedPrimal;
         timing_LLLNTL_PairRedPrimalRandomized [iteration][idx] = runningTime_LLLNTL_PairRedPrimalRandomized;
         timing_LLLNTL_PostPairRedPrimalRandomized [iteration][idx] = runningTime_LLLNTL_PostPairRedPrimalRandomized;

         //timing_LLLNTL_Exact [iteration][idx] = runningTime_LLLNTL_Exact;

         timing_BKZNTL [iteration][idx] = runningTime_BKZNTL;
         timing_BKZNTL_PairRedPrimal [iteration][idx] = runningTime_BKZNTL_PairRedPrimal;
         timing_BKZNTL_PostPairRedPrimal [iteration][idx] = runningTime_BKZNTL_PostPairRedPrimal;
         timing_BKZNTL_PairRedPrimalRandomized [iteration][idx] = runningTime_BKZNTL_PairRedPrimalRandomized;
         timing_BKZNTL_PostPairRedPrimalRandomized [iteration][idx] = runningTime_BKZNTL_PostPairRedPrimalRandomized;

         timing_BB_Only [iteration][idx] = runningTime_BB_Only;
         timing_BB_Classic1 [iteration][idx] = runningTime_BB_Classic1;
         timing_BB_Classic2 [iteration][idx] = runningTime_BB_Classic2;
         timing_BB_BKZ1 [iteration][idx] = runningTime_BB_BKZ1;
         timing_BB_BKZ2 [iteration][idx] = runningTime_BB_BKZ2;


         length_PairRedPrimal [iteration][idx] = lattice_PairRedPrimal.getVecNorm(0);
         length_PairRedPrimalRandomized [iteration][idx] = lattice_PairRedPrimalRandomized.getVecNorm(0);

         length_LLL [iteration][idx] = lattice_LLL.getVecNorm(0);
         length_LLL_PostPairRedPrimal [iteration][idx] = lattice_PairRedPrimal_LLL.getVecNorm(0);
         length_LLL_PostPairRedPrimalRandomized [iteration][idx] = lattice_PairRedPrimalRandomized_LLL.getVecNorm(0);

         length_LLLNTL [iteration][idx] = lattice_LLLNTL.getVecNorm(0);
         length_LLLNTL_PostPairRedPrimal [iteration][idx] = lattice_PairRedPrimal_LLLNTL.getVecNorm(0);
         length_LLLNTL_PostPairRedPrimalRandomized [iteration][idx] = lattice_PairRedPrimalRandomized_LLLNTL.getVecNorm(0);

         //length_LLLNTL_Exact [iteration][idx] = lattice_LLLNTL_Exact.getVecNorm(0);

         length_BKZNTL [iteration][idx] = lattice_BKZNTL.getVecNorm(0);
         length_BKZNTL_PostPairRedPrimal [iteration][idx] = lattice_PairRedPrimal_BKZNTL.getVecNorm(0);
         length_BKZNTL_PostPairRedPrimalRandomized [iteration][idx] = lattice_PairRedPrimalRandomized_BKZNTL.getVecNorm(0);


         length_BB_Only [iteration][idx] = lattice_BB_Only.getVecNorm(0);
         length_BB_Classic [iteration][idx] = lattice_BB_Classic.getVecNorm(0);
         length_BB_BKZ [iteration][idx] = lattice_BB_BKZ.getVecNorm(0);

      }

   } // end iteration loop over matrices of same dimension



#ifdef PRINT_CONSOLE

      //------------------------------------------------------------------------------------
      // Results printing in console
      //------------------------------------------------------------------------------------

      // print parameters used
      cout << "\n" << endl;
      cout << "epsilon = " << epsilon << endl;
      cout << "dimension = " << MinDimension << endl;
      cout << "nombre de matrices testées = " << maxIteration << endl;

      // print statistics
      cout << "\n---------------- TIMING AVG ----------------\n" << endl;

      cout << "       PairRedPrimal = " << Average(timing_PairRedPrimal) << endl;
      cout << " PairRedPrimalRandom = " << Average(timing_PairRedPrimalRandomized) << endl;
      cout << endl;

      cout << "                 LLL = " << Average(timing_LLL) << endl;
      cout << "         PairRed+LLL = " << Average(timing_LLL_PairRedPrimal) + Average(timing_LLL_PostPairRedPrimal);
      cout << " (" << Average(timing_LLL_PostPairRedPrimal) << ")" << endl;
      cout << "   PairRedRandom+LLL = " << Average(timing_LLL_PairRedPrimalRandomized) + Average(timing_LLL_PostPairRedPrimalRandomized);
      cout << " (" << Average(timing_LLL_PostPairRedPrimalRandomized) << ")" << endl;
      cout << endl;

      cout << "              LLLNTL = " << Average(timing_LLLNTL) << endl;
      cout << "      PairRed+LLLNTL = " << Average(timing_LLLNTL_PairRedPrimal) + Average(timing_LLLNTL_PostPairRedPrimal);
      cout << " (" << Average(timing_LLLNTL_PostPairRedPrimal) << ")" << endl;
      cout << "PairRedRandom+LLLNTL = " << Average(timing_LLLNTL_PairRedPrimalRandomized) + Average(timing_LLLNTL_PostPairRedPrimalRandomized);
      cout << " (" << Average(timing_LLLNTL_PostPairRedPrimalRandomized) << ")" << endl;
      cout << endl;

      //cout << "        LLLNTL_Exact = " << Average(timing_LLL_NTL_Exact) << endl;
      cout << endl;

      cout << "              BKZNTL = " << Average(timing_BKZNTL) << endl;
      cout << "      PairRed+BKZNTL = " << Average(timing_BKZNTL_PairRedPrimal) + Average(timing_BKZNTL_PostPairRedPrimal);
      cout << " (" << Average(timing_BKZNTL_PostPairRedPrimal) << ")" << endl;
      cout << "PairRedRandom+BKZNTL = " << Average(timing_BKZNTL_PairRedPrimalRandomized) + Average(timing_BKZNTL_PostPairRedPrimalRandomized);
      cout << " (" << Average(timing_BKZNTL_PostPairRedPrimalRandomized) << ")" << endl;
      cout << endl;

      cout << "             BB Only = " << Average(timing_BB_Only) << endl;
      cout << "          BB Classic = " << Average(timing_BB_Classic1) + Average(timing_BB_Classic2);
      cout << " (" << Average(timing_BB_Classic2) << ")" << endl,
      cout << "              BB BKZ = " << Average(timing_BB_BKZ1) + Average(timing_BB_BKZ2);
      cout << " (" << Average(timing_BB_BKZ2) << ")" << endl;

      cout << "\n--------------------------------------------" << endl;



      cout << "\n---------------- LENGTH AVG ----------------\n" << endl;

      cout << "             Initial = " << conv<ZZ>(Average(length_Initial)) << endl;
      cout << "       PairRedPrimal = " << conv<ZZ>(Average(length_PairRedPrimal)) << endl;
      cout << " PairRedPrimalRandom = " << conv<ZZ>(Average(length_PairRedPrimalRandomized)) << endl;
      cout << endl;

      cout << "                 LLL = " << conv<ZZ>(Average(length_LLL)) << endl;
      cout << "         PairRed+LLL = " << conv<ZZ>(Average(length_LLL_PostPairRedPrimal)) << endl;
      cout << "   PairRedRandom+LLL = " << conv<ZZ>(Average(length_LLL_PostPairRedPrimalRandomized)) << endl;
      cout << endl;

      cout << "              LLLNTL = " << conv<ZZ>(Average(length_LLLNTL)) << endl;
      cout << "      PairRed+LLLNTL = " << conv<ZZ>(Average(length_LLLNTL_PostPairRedPrimal)) << endl;
      cout << "PairRedRandom+LLLNTL = " << conv<ZZ>(Average(length_LLLNTL_PostPairRedPrimalRandomized)) << endl;
      cout << endl;

      //cout << "       LLL_NTL_Exact = " << conv<ZZ>(Average(length_LLL_NTL_Exact)) << endl;
      cout << endl;

      cout << "              BKZNTL = " << conv<ZZ>(Average(length_BKZNTL)) << endl;
      cout << "      PairRed+BKZNTL = " << conv<ZZ>(Average(length_BKZNTL_PostPairRedPrimal)) << endl;
      cout << "PairRedRandom+BKZNTL = " << conv<ZZ>(Average(length_BKZNTL_PostPairRedPrimalRandomized)) << endl;
      cout << endl;

      cout << "             BB Only = " << conv<ZZ>(Average(length_BB_Only)) << endl;
      cout << "          BB Classic = " << conv<ZZ>(Average(length_BB_Classic)) << endl;
      cout << "              BB BKZ = " << conv<ZZ>(Average(length_BB_BKZ)) << endl;


      cout << "\n--------------------------------------------" << endl;

#endif


#ifdef WITH_R
   /*----------------------------------------------*/
   // UTILISATION DE R
   /*----------------------------------------------*/


   RInside R(argc, argv);              // create an embedded R instance

   R["Mindimension"] = MinDimension;
   R["Maxdimension"] = MaxDimension;
   R["dimension"] = Interval_dim;
   //R["timing_Initial"] = toRcppMatrix(timing_Initial, maxIteration);
   R["timing_PairRedPrimal"] = toRcppMatrix<double>(timing_PairRedPrimal, maxIteration);
   R["timing_PairRedPrimalRandomized"] = toRcppMatrix<double>(timing_PairRedPrimalRandomized, maxIteration);
   R["timing_LLL"] = toRcppMatrix<double>(timing_LLL, maxIteration);
   R["timing_LLL_PostPairRedPrimal"] = toRcppMatrix<double>(timing_LLL_PostPairRedPrimal, maxIteration);
   R["timing_LLL_PostPairRedPrimalRandomized"] = toRcppMatrix<double>(timing_LLL_PostPairRedPrimalRandomized, maxIteration);
   R["timing_LLLNTL"] = toRcppMatrix<double>(timing_LLLNTL, maxIteration);
   R["timing_LLLNTL_PostPairRedPrimal"] = toRcppMatrix<double>(timing_LLLNTL_PostPairRedPrimal, maxIteration);
   R["timing_LLLNTL_PostPairRedPrimalRandomized"] = toRcppMatrix<double>(timing_LLLNTL_PostPairRedPrimalRandomized, maxIteration);
   R["timing_BKZNTL"] = toRcppMatrix<double>(timing_BKZNTL, maxIteration);
   R["timing_BKZNTL_PostPairRedPrimal"] = toRcppMatrix<double>(timing_BKZNTL_PostPairRedPrimal, maxIteration);
   R["timing_BKZNTL_PostPairRedPrimalRandomized"] = toRcppMatrix<double>(timing_BKZNTL_PostPairRedPrimalRandomized, maxIteration);
   R["timing_BB_Only"] = toRcppMatrix<double>(timing_BB_Only, maxIteration);
    //R["timing_LLLNTL"] = 4;

   /*std::string str =
      "cat('Running ls()\n'); print(ls()); "
      "cat('Showing M\n'); print(M); "
      "cat('Showing colSums()\n'); Z <- colSums(M); print(Z); "
      "Z";  */

    // by running parseEval, we get the last assignment back, here the filename

   std::string outPath = "~/Desktop";
   std::string outFile = "myPlot.png";
   R["outPath"] = outPath;
   R["outFile"] = outFile;

   // alternatively, by forcing a display we can plot to screen
   string library = "library(ggplot2); ";
   string build_data_frame =
   "df <- data.frame(indice = seq(1:dimension),"
     "PairRedPrimal=colMeans(timing_PairRedPrimal), "
     "PairRedPrimalRandomized = colMeans(timing_PairRedPrimalRandomized),"
     "LLL = colMeans(timing_LLL),"
     "LLL_PostPairRedPrimal = colMeans(timing_LLL_PostPairRedPrimal),"
     "LLL_PostPairRedPrimalRandomized = colMeans(timing_LLL_PostPairRedPrimalRandomized),"
     "LLLNTL = colMeans(timing_LLLNTL),"
     "LLLNTL_PostPairRedPrimal = colMeans(timing_LLLNTL_PostPairRedPrimal),"
     "LLLNTL_PostPairRedPrimalRandomized = colMeans(timing_LLLNTL_PostPairRedPrimalRandomized),"
     "BKZNTL = colMeans(timing_BKZNTL),"
     "BKZNTL_PostPairRedPrimal = colMeans(timing_BKZNTL_PostPairRedPrimal),"
     "BKZNTL_PostPairRedPrimalRandomized = colMeans(timing_BKZNTL_PostPairRedPrimalRandomized),"
     "BB_Only = colMeans(timing_BB_Only)"
     ");";

   string build_plot =
   "myPlot <- ggplot() + "
     "geom_line(data=df, aes(x=indice, y=PairRedPrimal, color ='PairRedPrimal')) + "
     "geom_line(data=df, aes(x=indice, y=PairRedPrimalRandomized, color ='PairRedPrimalRandomized')) +"
     "geom_line(data=df, aes(x=indice, y=LLL, color ='LLL')) +"
     "geom_line(data=df, aes(x=indice, y=LLL_PostPairRedPrimal, color ='LLL_PostPairRedPrimal')) +"
     "geom_line(data=df, aes(x=indice, y=LLL_PostPairRedPrimalRandomized, color ='LLL_PostPairRedPrimalRandomized')) +"
     "geom_line(data=df, aes(x=indice, y=LLLNTL, color ='LLLNTL')) +"
     "geom_line(data=df, aes(x=indice, y=LLLNTL_PostPairRedPrimal, color ='LLLNTL_PostPairRedPrimal')) +"
     "geom_line(data=df, aes(x=indice, y=LLLNTL_PostPairRedPrimalRandomized, color ='LLLNTL_PostPairRedPrimalRandomized')) +"
     "geom_line(data=df, aes(x=indice, y=BKZNTL, color ='BKZNTL')) +"
     "geom_line(data=df, aes(x=indice, y=BKZNTL_PostPairRedPrimal, color ='BKZNTL_PostPairRedPrimal')) +"
     "geom_line(data=df, aes(x=indice, y=BKZNTL_PostPairRedPrimalRandomized, color ='BKZNTL_PostPairRedPrimalRandomized')) +"
     "geom_line(data=df, aes(x=indice, y=BB_Only, color ='BB_Only')) +"
     "labs(color='Legend text'); ";


   string print_plot =
    "print(myPlot); "
     "ggsave(filename=outFile, path=outPath, plot=myPlot); ";
   // parseEvalQ evluates without assignment
    //R.parseEvalQ(cmd);
   R.parseEvalQ(library);

   R.parseEvalQ(build_data_frame);
   R.parseEvalQ(build_plot);
   R.parseEvalQ(print_plot);

#endif


    return 0;
}
