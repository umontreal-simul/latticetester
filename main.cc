//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//

#define PRINT_CONSOLE

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
#include "latticetester/Types.h"

#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>
#ifdef WITH_R
#include <RInside.h>
#endif
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

#include "SimpleMRG.h"


using namespace std;
using namespace NTL;
using namespace LatticeTester;

const int MaxDimension = 15;

#ifdef PRINT_CONSOLE
const int MinDimension = MaxDimension - 1;
#else
const int MinDimension = 5;
#endif


const int Interval_dim = MaxDimension - MinDimension;

#ifdef WITH_R
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
#endif


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


mat_ZZ CreateRNGBasis (const ZZ modulus, const int k, const int dimension, ZZ seed)
{
	mat_ZZ basis;
	basis.SetDims(dimension,dimension);

	if (dimension < k+1) {
		// degenerate case: identity matrix only
		for (int i = 0; i < dimension; i++)
			basis[i][i] = 1;

	} else { //usual case

		// (a_i) coefficients
		vec_ZZ a;
		a = randomVector(k, modulus, seed);

		for (int i = 0; i < k; i++) {
			// left upper block
			basis[i][i] = 1;

			//right upper block
			vec_ZZ initialState;
			initialState = canonicVector(k, i);
			SimpleMRG myMRG (modulus, k, a, initialState);

			for (int l = k; l < dimension; l++)
				basis[i][l] = conv<ZZ>(myMRG.getNextValue());
		}

		// right lower block
		for (int i = k; i < dimension; i++)
			basis[i][i] = modulus;

	} // end if

	return basis;
}


mat_ZZ Dualize (const mat_ZZ V, const ZZ modulus, const int k)
{
	mat_ZZ W;
	W.SetDims(V.NumRows(), V.NumRows());

	transpose(W,-V);

	for (int i = 0; i < k; i++)
		W[i][i] = modulus;
	for (int i = k; i < V.NumRows(); i++)
		W[i][i] = 1;

	return W;
}

void reduce(Reducer & red, const string & name, const int & d){

}

void reduce(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){

   //------------------------------------------------------------------------------------
   // Pairwise reduction in primal basis only
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimal" )
      red.preRedDieterPrimalOnly(d);

   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction in primal basis only
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimalRandomized")
      red.preRedDieterPrimalOnlyRandomized(d, seed_dieter);

   //------------------------------------------------------------------------------------
   // LLL Richard
   //------------------------------------------------------------------------------------
   if(name == "LLL")
      red.redLLL(delta, maxcpt, dimension);

   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then LLL Richard
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimal_LLL")
      red.preRedDieterPrimalOnly(d);


   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then LLL Richard
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimalRandomized_LLL")
      red.preRedDieterPrimalOnlyRandomized(d, seed_dieter);

   //------------------------------------------------------------------------------------
   // LLL NTL reduction (floating point version = proxy)
   //------------------------------------------------------------------------------------
   if(name == "LLLNTL")
      red.redLLLNTLProxy(delta);

   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then LLL NTL proxy
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimal_LLLNTL")
      red.preRedDieterPrimalOnly(d);

   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then LLL NTL proxy
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimalRandomized_LLLNTL")
      red.preRedDieterPrimalOnlyRandomized(d, seed_dieter);


   //------------------------------------------------------------------------------------
   // LLL NTL Exact reduction only
   //------------------------------------------------------------------------------------
   //if(name == "LLLNTL_Exact")
   //   red.redLLLNTLExact(det2,a,b);


   //------------------------------------------------------------------------------------
   // BKZ NTL reduction
   //------------------------------------------------------------------------------------
   if(name=="BKZNTL")
      red.redBKZ(delta, blocksize);

   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then BKZ NTL proxy
   //------------------------------------------------------------------------------------
   if(name =="PairRedPrimal_BKZNTL")
      red.preRedDieterPrimalOnly(d);

   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then BKZ NTL proxy
   //------------------------------------------------------------------------------------
   if(name =="PairRedPrimalRandomized_BKZNTL")
      red.preRedDieterPrimalOnlyRandomized(d, seed_dieter);

   //------------------------------------------------------------------------------------
   // Branch and Bound classic
   //------------------------------------------------------------------------------------
   if(name =="BB_Classic"){
      red.preRedDieterPrimalOnly(d);
      red.redLLL(delta, maxcpt, dimension);
   }


   //------------------------------------------------------------------------------------
   // Branch and Bound post BKZ
   //------------------------------------------------------------------------------------
   if(name =="BB_BKZ")
      red.redBKZ(delta, blocksize);
}


void reduce2(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){


   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then LLL Richard
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimal_LLL")
      red.redLLL(delta, maxcpt, dimension);


   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then LLL Richard
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimalRandomized_LLL")
      red.redLLL(delta, maxcpt, dimension);

   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then LLL NTL proxy
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimal_LLLNTL")
      red.redLLLNTLProxy(delta);

   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then LLL NTL proxy
   //------------------------------------------------------------------------------------
   if(name == "PairRedPrimalRandomized_LLLNTL")
      red.redLLLNTLProxy(delta);


   //------------------------------------------------------------------------------------
   // Pairwise reduction (in primal basis only) and then BKZ NTL proxy
   //------------------------------------------------------------------------------------
   if(name =="PairRedPrimal_BKZNTL")
      red.redBKZ(delta, blocksize);

   //------------------------------------------------------------------------------------
   // Randomized pairwise reduction (in primal basis only) and then BKZ NTL proxy
   //------------------------------------------------------------------------------------
   if(name =="PairRedPrimalRandomized_BKZNTL")
      red.redBKZ(delta, blocksize);

   //------------------------------------------------------------------------------------
   // Branch and Bound classic
   //------------------------------------------------------------------------------------
   if(name =="BB_Classic"){
      red.shortestVector(L2NORM);
   }


   //------------------------------------------------------------------------------------
   // Branch and Bound post BKZ
   //------------------------------------------------------------------------------------
   if(name =="BB_BKZ")
      red.shortestVector(L2NORM);
}



/* example new CreateRNGBasis

   int main ()
   {
      //ZZ modulus = conv<ZZ>("678956454545356865342357689098765324686546576787568");
      ZZ modulus;
      power(modulus, 2, 6);
      modulus-=1;

      int k = 3;
      int dimension = 10;
      ZZ seed = conv<ZZ>(123456);

      mat_ZZ V;
      V = CreateRNGBasis (modulus, k, dimension, seed);

      mat_ZZ W;
      W = Dualize (V, modulus, k);

      cout << "Test V*transpose(W) = m*ID : ";
      cout << IsIdent( V*transpose(W) - diag(dimension, modulus-1), dimension);
      cout << endl;

      return 0;
    }

*/


//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main (int argc, char *argv[])
{
   // printing total running time
   clock_t begin = clock();

   // main parameters for the test
   const int min = 30;
   const int max = 100;

   const long a = 999999;
   const long b = 1000000;
   const double delta = (double) a/b;
   const double epsilon = 1.0 - delta;

   const int maxcpt = 1000000;
   const int d = 0;
   const long blocksize = 20; // for BKZ insertions

   ZZ modulus;
   power(modulus, 2, 31);
   modulus--;
   int k = 3;

   string names[] = {
      "PairRedPrimal",
      "PairRedPrimalRandomized",
      "LLL",
      "PairRedPrimal_LLL",
      "PairRedPrimalRandomized_LLL",
      "LLLNTL",
      "PairRedPrimal_LLLNTL",
      "PairRedPrimalRandomized_LLLNTL",
      "BKZNTL",
      "PairRedPrimal_BKZNTL",
      "PairRedPrimalRandomized_BKZNTL",
      //"BB_Only",
      "BB_Classic",
      "BB_BKZ"};
   
   string names2[] = {
      "PairRedPrimal_LLL",
      "PairRedPrimalRandomized_LLL",
      "PairRedPrimal_LLLNTL",
      "PairRedPrimalRandomized_LLLNTL",
      "PairRedPrimal_BKZNTL",
      "PairRedPrimalRandomized_BKZNTL",
      "BB_Classic",
      "BB_BKZ"};

   // iteration loop over matrices of same dimension
   const int maxIteration = 10;


   map<string, map<int, NScal[maxIteration]> > length_results;
   map<string, map<int, double[maxIteration]> > timing_results;


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

         ZZ seed = conv<ZZ>((iteration+1) * (iteration+1) * 123456789 * dimension);
         //int seed = (iteration+1) * (iteration+1) * 123456789 * dimension;
         int seed_dieter = (iteration+1) * dimension * 12342;
         //int seed = (int) (iteration+1) * 12345 * time(NULL);

         int idx = dimension - MinDimension; // Stock the indices for table

         // We create copies of the same basis
         BMat basis_PairRedPrimal (dimension, dimension);
         ZZ det;
         //RandomMatrix(basis_PairRedPrimal, det, min, max, seed);


         mat_ZZ V;
         V = CreateRNGBasis (modulus, k, dimension, seed);

         mat_ZZ W;
         W = Dualize (V, modulus, k);

         map < string, BMat* > basis;
         map < string, BMat* > dualbasis;
         map < string, IntLatticeBasis* > lattices;
         map < string, Reducer* > reducers;

         for(const string &name : names){
            basis[name] = new BMat(V);
            dualbasis[name] = new BMat(W);
            lattices[name] = new IntLatticeBasis(*basis[name], dimension);
            reducers[name] = new Reducer(*lattices[name]);
         }

         map < string, clock_t > timing;






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

         lattice_PairRedPrimal.setNegativeNorm(true);
         lattice_PairRedPrimal.updateVecNorm();
         lattice_PairRedPrimal.sort(0);
         NScal initialShortestVectorLength = lattice_PairRedPrimal.getVecNorm(0);
         length_Initial [iteration][idx] = initialShortestVectorLength;


         clock_t begin = clock();
         clock_t end = clock();
         for(const string &name : names){
            begin = clock();
            reduce(*reducers[name], name, d, seed_dieter, blocksize, delta, maxcpt, dimension);
            end = clock();
            timing_results[name][dimension][iteration] = double (end - begin) / CLOCKS_PER_SEC;
            
            lattices[name]->setNegativeNorm();
            lattices[name]->updateVecNorm();
            lattices[name]->sort(0);
            ++show_progress;
         }

         for(const string &name : names2){
            begin = clock();
            reduce2(*reducers[name], name, d, seed_dieter, blocksize, delta, maxcpt, dimension);
            end = clock();
            timing_results[name+"2"][dimension][iteration] = double (end - begin) / CLOCKS_PER_SEC;
            lattices[name]->setNegativeNorm();
            lattices[name]->updateVecNorm();
            lattices[name]->sort(0);
            ++show_progress;
         }
         
         for(const string &name : names){
            length_results[name][dimension][iteration] = lattices[name]->getVecNorm(0);
         }
         


         for(const string &name : names){
            basis[name]->kill();
            dualbasis[name]->kill();
            delete lattices[name];
            delete reducers[name];
         }

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

#ifdef RESULTAT_BRUT



   map<string, map<int, NScal[iteration]> > length_results;
   map<string, map<int, double[iteration]> > timing_results;








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

   // printing total running time
   clock_t end = clock();
   cout << "\nTotal running time = " << (double) (end - begin) / CLOCKS_PER_SEC << endl;

   return 0;
}
