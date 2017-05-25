//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//


// select pre compiling options
//----------------------------------------------------------------------------------------

#define PRINT_CONSOLE

//----------------------------------------------------------------------------------------

#include <iostream>
#include <map>
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

#ifdef WITH_R
#include <RInside.h>
#endif

#include "SimpleMRG.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;

// projection parameters definition
//----------------------------------------------------------------------------------------

// Use of the Dual
bool WITH_DUAL = false;


// ireration loop over the dimension of lattices
const int MinDimension = 10;
#ifdef PRINT_CONSOLE
const int MaxDimension = MinDimension + 1;
#else
const int MaxDimension = 20;
#endif

// order
const int order = 5;

// iteration loop over matrices of same dimension
const int maxIteration = 10;

// Epsilon
const long a = 999999;
const long b = 1000000;
const double delta = (double) a/b;
const double epsilon = 1.0 - delta;

const int maxcpt = 1000000; // for redLLL
const int d = 0; // for preRedDieter
const long blocksize = 10; // for BKZ insertions

// modulus
const ZZ modulusRNG = power_ZZ(2, 31) - 1;

const int Interval_dim = MaxDimension - MinDimension;

// still usefull ?
const int minCoeff = 40;
const int maxCoeff = 1000;


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
   "BB_BKZ",
   //"DIETER", //WARNING USE DIETER ONLY FOR DIM < 6
   "MINKOWSKI"};

string names2[] = {
   "PairRedPrimal_LLL",
   "PairRedPrimalRandomized_LLL",
   "PairRedPrimal_LLLNTL",
   "PairRedPrimalRandomized_LLLNTL",
   "PairRedPrimal_BKZNTL",
   "PairRedPrimalRandomized_BKZNTL",
   "BB_Classic",
   "BB_BKZ"};



// functions used in main program
//----------------------------------------------------------------------------------------

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

template<typename Type, unsigned long Size>
Type mean(const array <Type, Size> array) {
   Type sum (0);
   for(int i = 0; i<Size; i++)
       sum += array[i];
   return sum / Size;
}

template<typename Type, unsigned long Size>
Type variance(const array <Type, Size> array) {
   Type sum (0);
   Type mean_tmp(mean(array));
   for(int i = 0; i<Size; i++)
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

	for (int i = 0; i < k; i++)
		W[i][i] = modulus;
	for (int i = k; i < V.NumRows(); i++)
		W[i][i] = 1;

	return W;
}

void reduce(Reducer & red, const string & name, const int & d){

}


void reduce(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){
   //cout << name << endl;

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

   //------------------------------------------------------------------------------------
   // Dieter Method
   //------------------------------------------------------------------------------------
   if(name == "DIETER" && WITH_DUAL)
      red.shortestVectorDieter(L2NORM);

   //------------------------------------------------------------------------------------
   // Minkowski reduction
   //------------------------------------------------------------------------------------
   if(name == "MINKOWSKI")
      red.reductMinkowski(d);
}


void reduce2(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){
   //cout << name << endl;


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
   if(name =="BB_BKZ" && !WITH_DUAL){
      red.shortestVector(L2NORM);
   }
}



//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main (int argc, char *argv[])


{
   // printing total running time
   clock_t begin = clock();

   // Print parameters
   cout << "epsilon = " << epsilon << endl;
   cout << "dimension = " << MinDimension << endl;
   cout << "nombre de matrices testées = " << maxIteration << endl;
   //cout << "dimension minimale : " << MinDimension << endl;
   //cout << "dimension maximale : " << MaxDimension << endl;
   cout << "ordre de la matrice : " << order << endl;
   cout << endl;

   // Stock Results
   map<string, map<int, array<NScal, maxIteration> > > length_results;
   map<string, map<int, array<double, maxIteration> > > timing_results;
   // old declaration using C-style arrays not accepted by some compilators
   //map<string, map<int, NScal[maxIteration]> > length_results;
   //map<string, map<int, double[maxIteration]> > timing_results;

   // to display progress bar
   boost::progress_display show_progress(maxIteration*Interval_dim);

   // arrays to store values

   for (int dimension = MinDimension; dimension < MaxDimension; dimension++){

      for (int iteration = 0; iteration < maxIteration; iteration++){

         ZZ seedZZ = conv<ZZ>((iteration+1) * (iteration+1) * 123456789 * dimension);
         int seed = (iteration+1) * (iteration+1) * 123456789 * dimension;
         int seed_dieter = (iteration+1) * dimension * 12342;
         //int seed = (int) (iteration+1) * 12345 * time(NULL);

         // We create copies of the same basis
         BMat basis_PairRedPrimal (dimension, dimension);
         ZZ det;
         RandomMatrix(basis_PairRedPrimal, det, minCoeff, maxCoeff, seed);

         mat_ZZ V;
         V = CreateRNGBasis (modulusRNG, order, dimension, seedZZ);

         mat_ZZ W;
         W = Dualize (V, modulusRNG, order);

         map < string, BMat* > basis;
         map < string, BMat* > dualbasis;
         map < string, IntLatticeBasis* > lattices;
         map < string, Reducer* > reducers;

         for(const string &name : names){
            basis[name] = new BMat(V);
            if(WITH_DUAL){
               dualbasis[name] = new BMat(W);
               lattices[name] = new IntLatticeBasis(*basis[name], *dualbasis[name], modulusRNG, dimension);
            }
            else{
               lattices[name] = new IntLatticeBasis(*basis[name], dimension);
            }
            // IF WE WANT FULL RANDOM MATRIX
            //basis[name] = new BMat(basis_PairRedPrimal);

            reducers[name] = new Reducer(*lattices[name]);
         }

         map < string, clock_t > timing;

         lattices["initial"] = new IntLatticeBasis(V, dimension);
         lattices["initial"]->setNegativeNorm(true);
         lattices["initial"]->updateVecNorm();
         lattices["initial"]->sort(0);
         length_results["initial"][dimension][iteration] = lattices["initial"]->getVecNorm(0);


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
            //cout << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;

         }

         for(const string &name : names2){
            //cout << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;
            begin = clock();
            reduce2(*reducers[name], name, d, seed_dieter, blocksize, delta, maxcpt, dimension);
            end = clock();
            timing_results[name+"2"][dimension][iteration] = double (end - begin) / CLOCKS_PER_SEC;
            lattices[name]->setNegativeNorm();
            lattices[name]->updateVecNorm();
            lattices[name]->sort(0);
            //cout << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;
         }

         for(const string &name : names){
            length_results[name][dimension][iteration] = lattices[name]->getVecNorm(0);
         }



         for(const string &name : names){
            basis[name]->BMat::clear();
            //dualbasis[name]->kill();
            delete lattices[name];
            delete reducers[name];
         }

         delete lattices["initial"];

         ++show_progress;
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
   if(WITH_DUAL)
      cout << "Dual utilisé" << endl;
   else
      cout << "Dual non utilisé" << endl;

   // print statistics
   cout << "\n---------------- TIMING AVG ----------------\n" << endl;

   cout << "       PairRedPrimal = " << mean(timing_results["PairRedPrimal"][MinDimension]) << "( +/- " << variance(timing_results["PairRedPrimal"][MinDimension]) << ")" << endl;
   cout << " PairRedPrimalRandom = " << mean(timing_results["PairRedPrimalRandomized"][MinDimension]) << "( +/- " << variance(timing_results["PairRedPrimalRandomized"][MinDimension]) << ")" << endl;
   cout << endl;

   cout << "                 LLL = " << mean(timing_results["LLL"][MinDimension]) << "( +/- " << variance(timing_results["LLL"][MinDimension]) << ")" << endl;
   cout << "         PairRed+LLL = " << mean(timing_results["PairRedPrimal_LLL"][MinDimension]) + mean(timing_results["PairRedPrimal_LLL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimal_LLL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimal_LLL2"][MinDimension]) << ")" << endl;
   cout << "   PairRedRandom+LLL = " << mean(timing_results["PairRedPrimalRandomized_LLL"][MinDimension]) + mean(timing_results["PairRedPrimalRandomized_LLL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimalRandomized_LLL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_LLL2"][MinDimension]) << ")" << endl;
   cout << endl;

   cout << "              LLLNTL = " << mean(timing_results["LLLNTL"][MinDimension]) << "( +/- " << variance(timing_results["LLLNTL"][MinDimension]) << ")" << endl;
   cout << "      PairRed+LLLNTL = " << mean(timing_results["PairRedPrimal_LLLNTL"][MinDimension]) + mean(timing_results["PairRedPrimal_LLLNTL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimal_LLLNTL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimal_LLLNTL"][MinDimension]) << ")" << endl;
   cout << "PairRedRandom+LLLNTL = " << mean(timing_results["PairRedPrimalRandomized_LLLNTL"][MinDimension]) + mean(timing_results["PairRedPrimalRandomized_LLLNTL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimalRandomized_LLLNTL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_LLLNTL2"][MinDimension]) << ")" << endl;
   cout << endl;

   //cout << "        LLLNTL_Exact = " << mean(timing_LLL_NTL_Exact) << endl;
   cout << endl;

   cout << "              BKZNTL = " << mean(timing_results["BKZNTL"][MinDimension]) << ") ( +/- " << variance(timing_results["BKZNTL"][MinDimension]) << ")" << endl;
   cout << "      PairRed+BKZNTL = " << mean(timing_results["PairRedPrimal_BKZNTL"][MinDimension]) + mean(timing_results["PairRedPrimal_BKZNTL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimal_BKZNTL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimal_BKZNTL2"][MinDimension]) << ")" << endl;
   cout << "PairRedRandom+BKZNTL = " << mean(timing_results["PairRedPrimalRandomized_BKZNTL"][MinDimension]) + mean(timing_results["PairRedPrimalRandomized_BKZNTL2"][MinDimension]);
   cout << " (" << mean(timing_results["PairRedPrimalRandomized_BKZNTL2"][MinDimension]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_BKZNTL2"][MinDimension]) << ")" << endl;
   cout << endl;

   cout << "          BB Classic = " << mean(timing_results["BB_Classic"][MinDimension]) + mean(timing_results["BB_Classic2"][MinDimension]);
   cout << " (" << mean(timing_results["BB_Classic2"][MinDimension]) << ") ( +/- " << variance(timing_results["BB_Classic2"][MinDimension]) << ")" ")" << endl,
   cout << "              BB BKZ = " << mean(timing_results["BB_BKZ"][MinDimension]) + mean(timing_results["BB_BKZ2"][MinDimension]);
   cout << " (" << mean(timing_results["BB_BKZ2"][MinDimension]) << ") ( +/- " << variance(timing_results["BB_BKZ2"][MinDimension]) << ")";
   if (WITH_DUAL)
      cout << " BB NON EFFECTUER CAR UTILISATION DU DUAL" << endl;
   else
      cout << endl;

   cout << endl;

   //cout << "    Dieter Reduction = " << mean(timing_results["DIETER"][MinDimension]);
   //if (!WITH_DUAL)
   //   cout << " DIETER NON EFFECTUER CAR DUAL NECESSAIRE" << endl;
   //else
   //   cout << endl;
   cout << " Minkowski Reduction = " << mean(timing_results["MINKOWSKI"][MinDimension]) << endl;

   cout << "\n--------------------------------------------" << endl;



   cout << "\n---------------- LENGTH AVG ----------------\n" << endl;

   cout << "             Initial = " << conv<ZZ>(mean(length_results["initial"][MinDimension])) << endl;

   cout << "       PairRedPrimal = " << conv<ZZ>(mean(length_results["PairRedPrimal"][MinDimension])) << endl;
   cout << " PairRedPrimalRandom = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized"][MinDimension])) << endl;
   cout << endl;

   cout << "                 LLL = " << conv<ZZ>(mean(length_results["LLL"][MinDimension])) << endl;
   cout << "         PairRed+LLL = " << conv<ZZ>(mean(length_results["PairRedPrimal_LLL"][MinDimension])) << endl;
   cout << "   PairRedRandom+LLL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_LLL"][MinDimension])) << endl;
   cout << endl;

   cout << "              LLLNTL = " << conv<ZZ>(mean(length_results["LLLNTL"][MinDimension])) << endl;
   cout << "      PairRed+LLLNTL = " << conv<ZZ>(mean(length_results["PairRedPrimal_LLLNTL"][MinDimension])) << endl;
   cout << "PairRedRandom+LLLNTL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_LLLNTL"][MinDimension])) << endl;
   cout << endl;
   cout << endl;

   cout << "              BKZNTL = " << conv<ZZ>(mean(length_results["BKZNTL"][MinDimension])) << endl;
   cout << "      PairRed+BKZNTL = " << conv<ZZ>(mean(length_results["PairRedPrimal_BKZNTL"][MinDimension])) << endl;
   cout << "PairRedRandom+BKZNTL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_BKZNTL"][MinDimension])) << endl;
   cout << endl;

   //cout << "             BB Only = " << conv<ZZ>(mean(length_results["PairRedPrimal_BKZNTL"][MinDimension])) << endl;
   cout << "          BB Classic = " << conv<ZZ>(mean(length_results["BB_Classic"][MinDimension])) << endl,
   cout << "              BB BKZ = " << conv<ZZ>(mean(length_results["BB_BKZ"][MinDimension])) << endl;
   cout << endl;

   //cout << "    Dieter Reduction = " << conv<ZZ>(mean(length_results["DIETER"][MinDimension])) << endl,
   cout << " Minkowski Reduction = " << conv<ZZ>(mean(length_results["MINKOWSKI"][MinDimension])) << endl;

   cout << "\n--------------------------------------------" << endl;

#endif



#ifdef WITH_R
   /*----------------------------------------------*/
   // UTILISATION DE R
   /*----------------------------------------------*/


   RInside R(argc, argv);              // create an embedded R instance

   R["MinDimension"] = MinDimension;
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
