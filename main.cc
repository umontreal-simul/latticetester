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

// if true the output is written in a .txt file
// if false the ouput is printed in the console
// via the std::outFile command
bool outFileRequested = true;

// Use of the Dual
bool WITH_DUAL = true;

// ireration loop over the dimension of lattices
const int MinDimension = 20;
#ifdef PRINT_CONSOLE
const int MaxDimension = MinDimension + 1;
#else
const int MaxDimension = 12;
#endif


// order
const int order = 5;

// iteration loop over matrices of same dimension
const int maxIteration = 1;

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
   "BB_BKZ"
   //"DIETER", //WARNING USE DIETER ONLY FOR DIM < 6
   //"MINKOWSKI"
   };

string names2[] = {
   "PairRedPrimal_LLL",
   "PairRedPrimalRandomized_LLL",
   "PairRedPrimal_LLLNTL",
   "PairRedPrimalRandomized_LLLNTL",
   "PairRedPrimal_BKZNTL",
   "PairRedPrimalRandomized_BKZNTL",
   "BB_Classic",
   "BB_BKZ"
   };



// functions used in main program
//----------------------------------------------------------------------------------------

#ifdef WITH_R
template <typename Type, unsigned long Size, unsigned long Size2>
Rcpp::NumericMatrix toRcppMatrix( const array<array<Type, Size>, Size2> array)
{

   Rcpp::NumericMatrix mat(maxIteration, Interval_dim);
   for (int i = 0; i<Size ; i++) {
      for (int j = 0; j<Size2; j++) {
         conv(mat(i,j), array[j][i]);
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

template<typename Type, unsigned long Size>
Type mean(const array <Type, Size> array) {
   Type sum (0);
   for(unsigned int i = 0; i<Size; i++)
       sum += array[i];
   return sum / Size;
}

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

	for (int i = 0; i < k; i++)
		W[i][i] = modulus;
	for (int i = k; i < V.NumRows(); i++)
		W[i][i] = 1;

	return W;
}

bool reduce(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){

   bool ok(true);
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
   //if(name == "DIETER" && WITH_DUAL)
      //red.shortestVectorDieter(L2NORM);

   //------------------------------------------------------------------------------------
   // Minkowski reduction
   //------------------------------------------------------------------------------------
   if(name == "MINKOWSKI")
      ok = red.reductMinkowski(d);

   return ok;
}


bool reduce2(Reducer & red, const string & name, const int & d, int & seed_dieter, const int & blocksize, const double & delta, const int maxcpt, int dimension){
   //outFile << name << endl;

   bool ok(true);
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
      ok = red.shortestVector(L2NORM);
   }

   //------------------------------------------------------------------------------------
   // Branch and Bound post BKZ
   //------------------------------------------------------------------------------------
   if(name =="BB_BKZ" && !WITH_DUAL){
      ok = red.shortestVector(L2NORM);
   }
   return ok;
}



//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main (int argc, char *argv[])
{
   BMat dualbasis(4,4);
   BMat basis;
   basis.resize(4, 4);
   basis(0,0) = 2;
   basis(0,1) = 3;
   basis(0,2) = 4;
   basis(0,3) = 5;
   for (int i = 1; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
         if (i == j) {
            basis (i, j) = 10;
         } else {
            basis (i, j) = 0;
         }
      }
   }
   MScal m(10);
   Triangularization < BMat > (basis, dualbasis, 4, 4, m);
   CalcDual < BMat > (dualbasis, basis, 4, m);
   
   cout << " Base Primal \n" << basis << endl;
   cout << " Base Dual \n" << dualbasis << endl;
   
   
   
   
#if 0
   ofstream realOutFile;
   string fileName;
   if (outFileRequested) {
      cout << "Enter file name (without .txt extension): ";
      cin >> fileName;
   }
   ostream & outFile = outFileRequested ? realOutFile.open(fileName+".txt", std::ios::out), realOutFile : std::cout;

   //ofstream realOutFile;
   //ostream & outFile = outFileRequested ? realOutFile.open("nik.txt", std::ios::out), realOutFile : std::cout;

   // printing total running time
   clock_t begin = clock();

   outFile << endl;

   // Stock Results
   map<string, array< array<NScal, maxIteration>, Interval_dim > > length_results;
   map<string, array< array<double, maxIteration>, Interval_dim > > timing_results;
   map<string, int> nb_diff;

   // old declaration using C-style arrays not accepted by some compilators
   //map<string, map<int, NScal[maxIteration]> > length_results;
   //map<string, map<int, double[maxIteration]> > timing_results;

   // to display progress bar
   boost::progress_display show_progress(maxIteration*Interval_dim);

   // arrays to store values
   int id_dimension = 0;
   bool all_BB_over = true;
   int nb_error = 0;

   for (int dimension = MinDimension; dimension < MaxDimension; dimension++){

      id_dimension = dimension - MinDimension;

      for (int iteration = 0; iteration < maxIteration; iteration++){
         do{
            if(!all_BB_over){
               outFile << "/";
               nb_error++; // Pour la boucle
            }
            all_BB_over = true;
            ZZ seedZZ = conv<ZZ>((iteration+1) * (iteration+1) * 123456789 * dimension * (nb_error+1));
            int seed = (iteration+1) * (iteration+1) * 123456789 * dimension * (nb_error+1);
            int seed_dieter = (iteration+1) * dimension * 12342 * (nb_error+1) ;
            //int seed = (int) (iteration+1) * 12345 * time(NULL);

            // We create copies of the same basis
            BMat basis_PairRedPrimal (dimension, dimension);
            ZZ det;
            RandomMatrix(basis_PairRedPrimal, det, minCoeff, maxCoeff, seed);

            mat_ZZ V;
            V = CreateRNGBasis (modulusRNG, order, dimension, seedZZ);

            

            mat_ZZ W;
            W = Dualize (V, modulusRNG, order);
            
            mat_ZZ Wtmp;
            transpose(Wtmp, W);

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
            length_results["initial"][id_dimension][iteration] = lattices["initial"]->getVecNorm(0);


            clock_t begin = clock();
            clock_t end = clock();

            bool ok = true;
            for(const string &name : names){
               begin = clock();
               ok = reduce(*reducers[name], name, d, seed_dieter, blocksize, delta, maxcpt, dimension);
               end = clock();
               all_BB_over = all_BB_over && ok;
               timing_results[name][id_dimension][iteration] = double (end - begin) / CLOCKS_PER_SEC;
               lattices[name]->setNegativeNorm();
               lattices[name]->updateVecNorm();
               lattices[name]->sort(0);
               //outFile << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;

            }

            for(const string &name : names2){
               //outFile << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;
               begin = clock();
               ok = reduce2(*reducers[name], name, d, seed_dieter, blocksize, delta, maxcpt, dimension);
               end = clock();
               all_BB_over = all_BB_over && ok;
               timing_results[name+"2"][id_dimension][iteration] = double (end - begin) / CLOCKS_PER_SEC;
               lattices[name]->setNegativeNorm();
               lattices[name]->updateVecNorm();
               lattices[name]->sort(0);
               //outFile << "Norm of " << name << " : " << lattices[name]->getVecNorm(0) << endl;
            }



            for(const string &name : names){
               if((lattices[name]->getVecNorm(0) - lattices["BB_Classic"]->getVecNorm(0)) > 0.1)
                  nb_diff[name]++;
               length_results[name][id_dimension][iteration] = lattices[name]->getVecNorm(0);
            }

            if((lattices["initial"]->getVecNorm(0) - lattices["BB_Classic"]->getVecNorm(0)) > 0.1)
                  nb_diff["initial"]++;

            for(const string &name : names){
               basis[name]->BMat::clear();
               //dualbasis[name]->kill();
               delete lattices[name];
               delete reducers[name];
            }

            delete lattices["initial"];

         } while(!all_BB_over);
         ++show_progress;
      }

   } // end iteration loop over matrices of same dimension


#ifdef PRINT_CONSOLE

   //------------------------------------------------------------------------------------
   // Results printing in console
   //------------------------------------------------------------------------------------

   // print parameters used
   outFile << "\n" << endl;
   outFile << "epsilon = " << epsilon << endl;
   outFile << "dimension = " << MinDimension << endl;
   outFile << "nombre de matrices testées = " << maxIteration << endl;
   outFile << "ordre de la matrice : " << order << endl;
   outFile << "Nombre de Branch and Bound non terminés : " << nb_error << endl;
   if(WITH_DUAL)
      outFile << "Dual utilisé" << endl;
   else
      outFile << "Dual non utilisé" << endl;

   // print statistics
   outFile << "\n---------------- TIMING AVG ----------------\n" << endl;

   outFile << "       PairRedPrimal = " << mean(timing_results["PairRedPrimal"][0]) << "( +/- " << variance(timing_results["PairRedPrimal"][0]) << ")" << endl;
   outFile << " PairRedPrimalRandom = " << mean(timing_results["PairRedPrimalRandomized"][0]) << "( +/- " << variance(timing_results["PairRedPrimalRandomized"][0]) << ")" << endl;
   outFile << endl;

   outFile << "                 LLL = " << mean(timing_results["LLL"][0]) << "( +/- " << variance(timing_results["LLL"][0]) << ")" << endl;
   outFile << "         PairRed+LLL = " << mean(timing_results["PairRedPrimal_LLL"][0]) + mean(timing_results["PairRedPrimal_LLL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimal_LLL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimal_LLL2"][0]) << ")" << endl;
   outFile << "   PairRedRandom+LLL = " << mean(timing_results["PairRedPrimalRandomized_LLL"][0]) + mean(timing_results["PairRedPrimalRandomized_LLL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimalRandomized_LLL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_LLL2"][0]) << ")" << endl;
   outFile << endl;

   outFile << "              LLLNTL = " << mean(timing_results["LLLNTL"][0]) << "( +/- " << variance(timing_results["LLLNTL"][0]) << ")" << endl;
   outFile << "      PairRed+LLLNTL = " << mean(timing_results["PairRedPrimal_LLLNTL"][0]) + mean(timing_results["PairRedPrimal_LLLNTL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimal_LLLNTL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimal_LLLNTL"][0]) << ")" << endl;
   outFile << "PairRedRandom+LLLNTL = " << mean(timing_results["PairRedPrimalRandomized_LLLNTL"][0]) + mean(timing_results["PairRedPrimalRandomized_LLLNTL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimalRandomized_LLLNTL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_LLLNTL2"][0]) << ")" << endl;
   outFile << endl;

   //outFile << "        LLLNTL_Exact = " << mean(timing_LLL_NTL_Exact) << endl;
   outFile << endl;

   outFile << "              BKZNTL = " << mean(timing_results["BKZNTL"][0]) << ") ( +/- " << variance(timing_results["BKZNTL"][0]) << ")" << endl;
   outFile << "      PairRed+BKZNTL = " << mean(timing_results["PairRedPrimal_BKZNTL"][0]) + mean(timing_results["PairRedPrimal_BKZNTL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimal_BKZNTL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimal_BKZNTL2"][0]) << ")" << endl;
   outFile << "PairRedRandom+BKZNTL = " << mean(timing_results["PairRedPrimalRandomized_BKZNTL"][0]) + mean(timing_results["PairRedPrimalRandomized_BKZNTL2"][0]);
   outFile << " (" << mean(timing_results["PairRedPrimalRandomized_BKZNTL2"][0]) << ") ( +/- " << variance(timing_results["PairRedPrimalRandomized_BKZNTL2"][0]) << ")" << endl;
   outFile << endl;

   outFile << "          BB Classic = " << mean(timing_results["BB_Classic"][0]) + mean(timing_results["BB_Classic2"][0]);
   outFile << " (" << mean(timing_results["BB_Classic2"][0]) << ") ( +/- " << variance(timing_results["BB_Classic2"][0]) << ")" ")" << endl,
   outFile << "              BB BKZ = " << mean(timing_results["BB_BKZ"][0]) + mean(timing_results["BB_BKZ2"][0]);
   outFile << " (" << mean(timing_results["BB_BKZ2"][0]) << ") ( +/- " << variance(timing_results["BB_BKZ2"][0]) << ")";
   if (WITH_DUAL)
      outFile << " BB NON EFFECTUER CAR UTILISATION DU DUAL" << endl;
   else
      outFile << endl;

   outFile << endl;

   //outFile << "    Dieter Reduction = " << mean(timing_results["DIETER"][0]);
   //if (!WITH_DUAL)
   //   outFile << " DIETER NON EFFECTUER CAR DUAL NECESSAIRE" << endl;
   //else
   //   outFile << endl;
   //outFile << " Minkowski Reduction = " << mean(timing_results["MINKOWSKI"][0]) << endl;

   outFile << "\n--------------------------------------------" << endl;



   outFile << "\n---------------- LENGTH AVG ----------------\n" << endl;

   outFile << "             Initial = " << conv<ZZ>(mean(length_results["initial"][0])) << " Error Rate : " << (double) nb_diff["initial"]/maxIteration << endl;

   outFile << "       PairRedPrimal = " << conv<ZZ>(mean(length_results["PairRedPrimal"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimal"]/maxIteration << endl;
   outFile << " PairRedPrimalRandom = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimalRandomized"]/maxIteration << endl;
   outFile << endl;

   outFile << "                 LLL = " << conv<ZZ>(mean(length_results["LLL"][0])) << " Error Rate : " << (double) nb_diff["LLL"]/maxIteration << endl;
   outFile << "         PairRed+LLL = " << conv<ZZ>(mean(length_results["PairRedPrimal_LLL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimal_LLL"]/maxIteration << endl;
   outFile << "   PairRedRandom+LLL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_LLL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimalRandomized_LLL"]/maxIteration << endl;
   outFile << endl;

   outFile << "              LLLNTL = " << conv<ZZ>(mean(length_results["LLLNTL"][0])) << " Error Rate : " << (double) nb_diff["LLLNTL"]/maxIteration << endl;
   outFile << "      PairRed+LLLNTL = " << conv<ZZ>(mean(length_results["PairRedPrimal_LLLNTL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimal_LLLNTL"]/maxIteration << endl;
   outFile << "PairRedRandom+LLLNTL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_LLLNTL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimalRandomized_LLLNTL"]/maxIteration << endl;
   outFile << endl;
   outFile << endl;

   outFile << "              BKZNTL = " << conv<ZZ>(mean(length_results["BKZNTL"][0])) << " Error Rate : " << (double) nb_diff["BKZNTL"]/maxIteration << endl;
   outFile << "      PairRed+BKZNTL = " << conv<ZZ>(mean(length_results["PairRedPrimal_BKZNTL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimal_BKZNTL"]/maxIteration << endl;
   outFile << "PairRedRandom+BKZNTL = " << conv<ZZ>(mean(length_results["PairRedPrimalRandomized_BKZNTL"][0])) << " Error Rate : " << (double) nb_diff["PairRedPrimalRandomized_BKZNTL"]/maxIteration << endl;
   outFile << endl;

   //outFile << "             BB Only = " << conv<ZZ>(mean(length_results["PairRedPrimal_BKZNTL"][0])) << " Error Rate : " << (double) nb_diff[name]/maxIteration << endl;
   outFile << "          BB Classic = " << conv<ZZ>(mean(length_results["BB_Classic"][0])) << " Error Rate : " << (double) nb_diff["BB_Classic"]/maxIteration << endl,
   outFile << "              BB BKZ = " << conv<ZZ>(mean(length_results["BB_BKZ"][0])) << " Error Rate : " << (double) nb_diff["BB_BKZ"]/maxIteration << endl;
   outFile << endl;

   //outFile << "    Dieter Reduction = " << conv<ZZ>(mean(length_results["DIETER"][0])) << " Error Rate : " << (double) nb_diff[name]/maxIteration << endl,
   //outFile << " Minkowski Reduction = " << conv<ZZ>(mean(length_results["MINKOWSKI"][0])) << " Error Rate : " << (double) nb_diff["MINKOWSKI"]/maxIteration << endl;

   outFile << "\n--------------------------------------------\n" << endl;


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
   for(const string &name : names){
      R[name] = toRcppMatrix(timing_results[name]);
   }

    // by running parseEval, we get the last assignment back, here the filename

   std::string outPath = "~/Desktop";
   std::string outFile = "myPlot.png";
   R["outPath"] = outPath;
   R["outFile"] = outFile;

   // alternatively, by forcing a display we can plot to screen
   string library = "library(ggplot2); ";
   string build_data_frame = "df <- data.frame(indice = seq(1:dimension)";
   for(const string &name : names){
      build_data_frame += ", col_" + name + " =colMeans(" + name + ")";
   }
   build_data_frame += ");";


   string build_plot = "myPlot <- ggplot() + ";
   for(const string &name : names){
      build_plot += "geom_line(data=df, aes(x=indice, y=col_" + name + ", color ='col_" + name + "')) + ";
   }
   build_plot += "labs(color='Legend text'); ";

   string print_plot =
    "print(myPlot); "
     "ggsave(filename=outFile, path=outPath, plot=myPlot); ";
   // parseEvalQ evluates without assignment
   R.parseEvalQ(library);
   R.parseEvalQ(build_data_frame);
   R.parseEvalQ(build_plot);
   R.parseEvalQ(print_plot);

#endif

   // printing total running time
   clock_t end = clock();
   outFile << "\nTotal running time = " << (double) (end - begin) / CLOCKS_PER_SEC << endl;


   if (outFileRequested)
      realOutFile.close();
   
#endif

   return 0;
}

