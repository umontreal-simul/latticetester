//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//


#include "latticetester/Util.h"
#include "latticetester/Basis.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/IntLatticeBasis.h"
#include <NTL/ctools.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <RInside.h>
#include <iomanip>



#include "latticetester/Reducer.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <time.h>
#include <boost/progress.hpp>

// for LLL test
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






using namespace std;
using namespace LatticeTester;

void RandomMatrix (mat_ZZ& A, ZZ& det, int min, int max, int seed){
    
    int dim = (int) A.NumRows() ;
    srand(seed);
    
    do{
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++)
                A[i][j] = min + (rand() % (int)(max - min + 1));
        }
        // Richard implementation of basis include an empty first column
        // and an empty first lign, so we had a 1 on top left position
        // to preserve determinant calculation
        //A[0][0] = 1;
        
        det = determinant(A);
        
    } while ( det == 0 );
}

void RealLatticePrint (IntLatticeBasis& lattice)
{
    int dim = lattice.getDim();
    
    for (int i = 0; i < dim+1; i++){
        cout << "[";
        for (int j = 0; j < dim+1; j++)
            cout << lattice.getBasis()[i][j] << " ";
        cout << "]" << endl;
    }
    cout << endl;
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
Type Average(Type const(& array)[Size]) {
    Type sum (0);
    for(int i=0; i<Size; i++)
        sum += array[i];
    return sum / Size;
}




//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main()
{
    bool printMatricesDetails = false;
    bool printLists = false;
    
    // main parameters for the test
    int dimension = 10;
    int min = 10;
    int max = 100;
    
    long a = 999999;
    long b = 1000000;
    double delta = (double) a/b;
    double epsilon = 1.0 - delta;
    
    int maxcpt = 1000000;
    int d = 0;
    long blocksize = 30; // for BKZ insertions
    
    // iteration loop over matrices of same dimension
    const int maxIteration = 1;
    
    // print important information
    cout << "epsilon = " << epsilon << endl;
    cout << "dimension = " << dimension << endl;
    cout << "nombre de matrices testées = " << maxIteration << endl;
    cout << endl;
    
    // to display progress bar
    boost::progress_display show_progress(maxIteration);
    
    // arrays to store values
    double timing_PairRedPrimal [maxIteration];
    double timing_PairRedPrimalRandomized [maxIteration];
    double timing_LLL [maxIteration];
    double timing_LLL_PostPairRedPrimal [maxIteration];
    double timing_LLL_NTL_Proxy [maxIteration];
    double timing_PairRedPrimal_NTL [maxIteration];
    double timing_LLL_NTL_PostPairRedPrimal [maxIteration];
    double timing_LLL_NTL_Exact [maxIteration];
    double timing_BKZ_NTL [maxIteration];
    
    NScal length_PairRedPrimalRandomized [maxIteration];
    NScal length_Initial [maxIteration];
    NScal length_LLL [maxIteration];
    NScal length_PairRedPrimal [maxIteration];
    NScal length_LLL_PostPairRedPrimal [maxIteration];
    NScal length_LLL_NTL_Proxy [maxIteration];
    NScal length_PairRedPrimal_NTL [maxIteration];
    NScal length_LLL_NTL_PostPairRedPrimal [maxIteration];
    NScal length_LLL_NTL_Exact [maxIteration];
    NScal length_BKZ_NTL [maxIteration];
    
    
    for (int iteration = 0; iteration < maxIteration; iteration++){
        ++show_progress;
        
        int seed = (iteration+1) * (iteration+1) * 123456789;
        //int seed = (int) (iteration+1) * 12345 * time(NULL);
        
        // We create copies of the same basis for: LLL,
        // pairwisePrimalRed + LLL, LLL NTL floating point,
        // LLL NTL exact version, pairwiseRed + LLL, BKZ NTL,
        BMat basis_LLL (dimension, dimension);
        ZZ det;
        RandomMatrix(basis_LLL, det, min, max,seed);
        
        BMat basis_PairRedPrimalRandomized (basis_LLL);
        BMat basis_PairRedPrimal_LLL (basis_LLL);
        BMat basis_LLL_NTL_Proxy (basis_LLL);
        BMat basis_PairRedPrimal_LLL_NTL (basis_LLL);
        BMat basis_LLL_NTL_Exact (basis_LLL);
        BMat basis_BKZ_NTL (basis_LLL);
        
        IntLatticeBasis lattice_LLL (basis_LLL, dimension);
        IntLatticeBasis lattice_PairRedPrimalRandomized (basis_PairRedPrimalRandomized, dimension);
        IntLatticeBasis lattice_PairRedPrimal_LLL (basis_PairRedPrimal_LLL, dimension);
        IntLatticeBasis lattice_LLL_NTL_Proxy (basis_LLL_NTL_Proxy, dimension);
        IntLatticeBasis lattice_PairRedPrimal_LLL_NTL (basis_PairRedPrimal_LLL_NTL, dimension);
        IntLatticeBasis lattice_LLL_NTL_Exact (basis_LLL_NTL_Exact, dimension);
        IntLatticeBasis lattice_BKZ_NTL (basis_BKZ_NTL, dimension);
        
        lattice_LLL.setNegativeNorm();
        lattice_LLL.updateVecNorm();
        lattice_LLL.sort(0);
        NScal initialShortestVectorLength = lattice_LLL.getVecNorm(0);
        length_Initial [iteration] = initialShortestVectorLength;
        
        if (printMatricesDetails) {
            cout << "\n*** Initial basis ***" << endl;
            cout << "det = " << det << endl;
            cout << "Shortest vector = " << initialShortestVectorLength << endl;
            lattice_LLL.write();
        }
        
        Reducer reducer_LLL (lattice_LLL);
        Reducer reducer_PairRedPrimalRandomized (lattice_PairRedPrimalRandomized);
        Reducer reducer_PairRedPrimal_LLL (lattice_PairRedPrimal_LLL);
        Reducer reducer_LLL_NTL_Proxy (lattice_LLL_NTL_Proxy);
        Reducer reducer_PairRedPrimal_LLL_NTL (lattice_PairRedPrimal_LLL_NTL);
        Reducer reducer_LLL_NTL_Exact (lattice_LLL_NTL_Exact);
        Reducer reducer_BKZ_NTL (lattice_BKZ_NTL);
        
        
        //------------------------------------------------------------------------------------
        // Randomized pairwise reduction in primal basis only
        //------------------------------------------------------------------------------------
        
        clock_t begin_PairRedPrimalRandomized = clock();
        reducer_PairRedPrimalRandomized.preRedDieterPrimalOnlyRandomized(d);
        clock_t end_PairRedPrimalRandomized = clock();
        
        lattice_PairRedPrimalRandomized.setNegativeNorm();
        lattice_PairRedPrimalRandomized.updateVecNorm();
        lattice_PairRedPrimalRandomized.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** Randomized pairwise reduction only ***" << endl;
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
        
        
        //------------------------------------------------------------------------------------
        // Pairwise reduction (in primal basis only) and then LLL Richard
        //------------------------------------------------------------------------------------
        
        clock_t begin_PairRedPrimal_LLL1 = clock();
        reducer_PairRedPrimal_LLL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal_LLL1 = clock();
        
        lattice_PairRedPrimal_LLL.setNegativeNorm();
        lattice_PairRedPrimal_LLL.updateVecNorm();
        lattice_PairRedPrimal_LLL.sort(0);
        
        NScal intermediateShortestVectorLength = lattice_PairRedPrimal_LLL.getVecNorm(0);
        
        clock_t begin_PairRedPrimal_LLL2 = clock();
        reducer_PairRedPrimal_LLL.redLLL(delta, maxcpt, dimension);
        clock_t end_PairRedPrimal_LLL2 = clock();
        
        lattice_PairRedPrimal_LLL.setNegativeNorm();
        lattice_PairRedPrimal_LLL.updateVecNorm();
        lattice_PairRedPrimal_LLL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Pairwise reduction in primal and LLL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimal_LLL.getVecNorm(0) << endl;
            lattice_PairRedPrimal_LLL.write();
        }
        
        
        //------------------------------------------------------------------------------------
        // LLL NTL reduction (floating point version = proxy)
        //------------------------------------------------------------------------------------
        
        clock_t begin_LLL_NTL_Proxy = clock();
        reducer_LLL_NTL_Proxy.redLLLNTLProxy(delta);
        clock_t end_LLL_NTL_Proxy = clock();
        
        lattice_LLL_NTL_Proxy.setNegativeNorm();
        lattice_LLL_NTL_Proxy.updateVecNorm();
        lattice_LLL_NTL_Proxy.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** LLL NTL Proxy only ***" << endl;
            cout << "Shortest vector = " << lattice_LLL_NTL_Proxy.getVecNorm(1) << endl;
            lattice_LLL_NTL_Proxy.write();
        }
        
        
        //------------------------------------------------------------------------------------
        // Pairwise reduction (in primal basis only) and then LLL NTL proxy
        //------------------------------------------------------------------------------------
        
        clock_t begin_PairRedPrimal_LLL_NTL1 = clock();
        reducer_PairRedPrimal_LLL_NTL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal_LLL_NTL1 = clock();
        
        lattice_PairRedPrimal_LLL_NTL.setNegativeNorm(true);
        lattice_PairRedPrimal_LLL_NTL.updateVecNorm();
        lattice_PairRedPrimal_LLL_NTL.sort(0);
        
        NScal intermediateShortestVectorLengthBis = lattice_PairRedPrimal_LLL_NTL.getVecNorm(0);
        
        clock_t begin_PairRedPrimal_LLL_NTL2 = clock();
        reducer_PairRedPrimal_LLL_NTL.redLLLNTLProxy(delta);
        clock_t end_PairRedPrimal_LLL_NTL2 = clock();
        
        lattice_PairRedPrimal_LLL_NTL.setNegativeNorm();
        lattice_PairRedPrimal_LLL_NTL.updateVecNorm();
        lattice_PairRedPrimal_LLL_NTL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Pairwise reduction in primal and LLL NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimal_LLL_NTL.getVecNorm(0) << endl;
            lattice_PairRedPrimal_LLL_NTL.write();
        }
        
        
        //------------------------------------------------------------------------------------
        // LLL NTL Exact reduction only
        //------------------------------------------------------------------------------------
        
         ZZ det2;
         clock_t begin_LLL_NTL_Exact = clock();
         reducer_LLL_NTL_Exact.redLLLNTLExact(det2, a, b);
         clock_t end_LLL_NTL_Exact = clock();
         
         lattice_LLL_NTL_Exact.setNegativeNorm();
         lattice_LLL_NTL_Exact.updateVecNorm();
         lattice_LLL_NTL_Exact.sort(0);
         
         if (printMatricesDetails) {
         cout << "*** LLL NTL Exact only ***" << endl;
         cout << "Shortest vector = " << lattice_LLL_NTL_Exact.getVecNorm(0) << endl;
         lattice_LLL_NTL_Exact.write();
         }
        
        
        //------------------------------------------------------------------------------------
        // BKZ NTL reduction
        //------------------------------------------------------------------------------------
        
        clock_t begin_BKZ_NTL = clock();
        reducer_BKZ_NTL.redBKZ(delta, blocksize);
        clock_t end_BKZ_NTL = clock();
        
        lattice_BKZ_NTL.setNegativeNorm();
        lattice_BKZ_NTL.updateVecNorm();
        lattice_BKZ_NTL.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** BKZ NTL only ***" << endl;
            cout << "Shortest vector = " << lattice_BKZ_NTL.getVecNorm(0) << endl;
            lattice_BKZ_NTL.write();
        }
        
        //------------------------------------------------------------------------------------
        // Branch and Bound reduction
        //------------------------------------------------------------------------------------
        
        clock_t begin_BKZ_NTL = clock();
        reducer_PairRedPrimalRandomized.redBB0(L2NORM);
        
        
        
        //------------------------------------------------------------------------------------
        // timing updating
        //------------------------------------------------------------------------------------
        
        double runningTime_LLL = double(end_LLL - begin_LLL) / CLOCKS_PER_SEC;
        double runningTime_PairRedPrimalRandomized = double (end_PairRedPrimalRandomized - begin_PairRedPrimalRandomized) / CLOCKS_PER_SEC;
        double runningTime_PairRedPrimal = double(end_PairRedPrimal_LLL1 - begin_PairRedPrimal_LLL1) / CLOCKS_PER_SEC;
        double runningTime_LLL_PostPairRedPrimal = double(end_PairRedPrimal_LLL2 - begin_PairRedPrimal_LLL2) / CLOCKS_PER_SEC;
        double runningTime_LLL_NTL_Proxy = double(end_LLL_NTL_Proxy - begin_LLL_NTL_Proxy) / CLOCKS_PER_SEC;
        double runningTime_PairRedPrimal_NTL = double(end_PairRedPrimal_LLL_NTL1 - begin_PairRedPrimal_LLL_NTL1) / CLOCKS_PER_SEC;
        double runningTime_LLL_NTL_PostPairRedPrimal = double (end_PairRedPrimal_LLL_NTL2 - begin_PairRedPrimal_LLL_NTL2) / CLOCKS_PER_SEC;
        double runningTime_LLL_NTL_Exact = double(end_LLL_NTL_Exact - begin_LLL_NTL_Exact) / CLOCKS_PER_SEC;
        double runningTime_BKZ_NTL = double (end_BKZ_NTL - begin_BKZ_NTL) / CLOCKS_PER_SEC;
        
        
        //------------------------------------------------------------------------------------
        // timing and length arrays updating
        //------------------------------------------------------------------------------------
        
        timing_PairRedPrimalRandomized [iteration] = runningTime_PairRedPrimalRandomized;
        timing_LLL [iteration] = runningTime_LLL;
        timing_PairRedPrimal [iteration] = runningTime_PairRedPrimal;
        timing_LLL_PostPairRedPrimal [iteration] = runningTime_LLL_PostPairRedPrimal;
        timing_LLL_NTL_Proxy [iteration] = runningTime_LLL_NTL_Proxy;
        timing_PairRedPrimal_NTL [iteration] = runningTime_PairRedPrimal_NTL;
        timing_LLL_NTL_PostPairRedPrimal [iteration] = runningTime_LLL_NTL_PostPairRedPrimal;
        timing_LLL_NTL_Exact [iteration] = runningTime_LLL_NTL_Exact;
        timing_BKZ_NTL [iteration] = runningTime_BKZ_NTL;
        
        length_PairRedPrimalRandomized [iteration] = lattice_PairRedPrimalRandomized.getVecNorm(0);
        length_LLL [iteration] = lattice_LLL.getVecNorm(0);
        length_PairRedPrimal [iteration] = intermediateShortestVectorLength;
        length_LLL_PostPairRedPrimal [iteration] = lattice_PairRedPrimal_LLL.getVecNorm(0);
        length_LLL_NTL_Proxy [iteration] = lattice_LLL_NTL_Proxy.getVecNorm(0);
        length_PairRedPrimal_NTL [iteration] = intermediateShortestVectorLengthBis;
        length_LLL_NTL_PostPairRedPrimal [iteration] = lattice_PairRedPrimal_LLL_NTL.getVecNorm(0);
        length_LLL_NTL_Exact [iteration] = lattice_LLL_NTL_Exact.getVecNorm(0);
        length_BKZ_NTL [iteration] = lattice_BKZ_NTL.getVecNorm(0);
        
        
        
    } // end iteration loop over matrices of same dimension
    
    
    
    // print arrays
    
    if (printLists) {
        
        cout << "\nTIMING LIST ---------" << endl;
        print("PairRedPrimalRandomized", timing_PairRedPrimalRandomized, false);
        print("                    LLL", timing_LLL, false);
        print("          PairRedPrimal", timing_PairRedPrimal, false);
        print("      PostPairRedPrimal", timing_LLL_PostPairRedPrimal, false);
        print("          LLL_NTL_Proxy", timing_LLL_NTL_Proxy, false);
        print("      PairRedPrimal_NTL", timing_PairRedPrimal_NTL, false);
        print("  PostPairRedPrimal_NTL", timing_LLL_NTL_PostPairRedPrimal, false);
        print("          LLL_NTL_Exact", timing_LLL_NTL_Exact, false);
        print("                BKZ_NTL", timing_BKZ_NTL, false);
        
        cout << "\nLENGTH LIST ---------" << endl;
        print("Initial",length_Initial,true);
        print("PairRedPrimalRandomized", length_PairRedPrimalRandomized, true);
        print("                    LLL", length_LLL, true);
        print("          PairRedPrimal", length_PairRedPrimal, true);
        print("      PostPairRedPrimal", length_LLL_PostPairRedPrimal, true);
        print("          LLL_NTL_Proxy", length_LLL_NTL_Proxy, true);
        print("      PairRedPrimal_NTL", length_PairRedPrimal_NTL, true);
        print("  PostPairRedPrimal_NTL", length_LLL_NTL_PostPairRedPrimal, true);
        print("          LLL_NTL_Exact", length_LLL_NTL_Exact, true);
        print("                BKZ_NTL", length_BKZ_NTL, true);
    }
    
    
    //------------------------------------------------------------------------------------
    // Results printing in console
    //------------------------------------------------------------------------------------
    
    // print parameters used
    cout << "\n" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "dimension = " << dimension << endl;
    cout << "nombre de matrices testées = " << maxIteration << endl;
    
    // print statistics
    
    
    cout << "\n---------------- TIMING AVG ----------------\n" << endl;
    
    cout << " PairwisePrimal = " << Average(timing_PairRedPrimal) << endl;
    cout << " PairwiseRandom = " << Average(timing_PairRedPrimalRandomized) << endl;
    cout << endl;
    
    cout << "            LLL = " << Average(timing_LLL) << endl;
    cout << "    PairRed+LLL = " << Average(timing_PairRedPrimal) + Average(timing_LLL_PostPairRedPrimal);
    cout << " (" << Average(timing_LLL_PostPairRedPrimal) << ")" << endl;
    cout << endl;
    
    cout << "        LLL_NTL = " << Average(timing_LLL_NTL_Proxy) << endl;
    cout << "PairRed+LLL_NTL = " << Average(timing_PairRedPrimal_NTL) + Average(timing_LLL_NTL_PostPairRedPrimal);
    cout << " (" << Average(timing_LLL_NTL_PostPairRedPrimal) << ")" << endl;
    cout << "  LLL_NTL_Exact = " << Average(timing_LLL_NTL_Exact) << endl;
    cout << endl;
    
    cout << "        BKZ NTL = " << Average(timing_BKZ_NTL) << endl;
    
    cout << "\n--------------------------------------------" << endl;
    
    
    
    cout << "\n---------------- LENGTH AVG ----------------\n" << endl;
    
    cout << "        Initial = " << conv<ZZ>(Average(length_Initial)) << endl;
    cout << " PairwiseRandom = " << conv<ZZ>(Average(length_PairRedPrimalRandomized)) << endl;
    cout << "  PairRedPrimal = " << conv<ZZ>(Average(length_PairRedPrimal)) << endl;
    cout << endl;
    
    cout << "            LLL = " << conv<ZZ>(Average(length_LLL)) << endl;
    cout << "    PairRed+LLL = " << conv<ZZ>(Average(length_LLL_PostPairRedPrimal)) << endl;
    cout << endl;
    
    cout << "        LLL_NTL = " << conv<ZZ>(Average(length_LLL_NTL_Proxy)) << endl;
    cout << "PairRed+LLL_NTL = " << conv<ZZ>(Average(length_LLL_NTL_PostPairRedPrimal)) << endl;
    //cout << "  LLL_NTL_Exact = " << conv<ZZ>(Average(length_LLL_NTL_Exact)) << endl;
    cout << endl;
    
    cout << "        BKZ NTL = " << conv<ZZ>(Average(length_BKZ_NTL)) << endl;
    
    cout << "\n--------------------------------------------" << endl;
    
    
    
    
    
    
    
    
    
    
    
    //ZZ ab = conv<ZZ>("10888937647076736567376563809027278474092774049093877637265396456737653652456136387647828336378903893780");
    //for (int k = 0; k < 50; k++)
    //    cout << RandomBnd(ab) << endl;
    
    
    //CalcDual (MyLattice1.getPrimalBasis(), MyLattice1.getDualBasis(), dimension, det);
    //CalcDual (MyLattice2.getPrimalBasis(), MyLattice2.getDualBasis(), dimension, det);
    
    
    
    /*
     
     // printing matrices
     bool printMatrices = true;
     
     // Loop over dimension
     const int min_dimension = 6;
     const int max_dimension = 6;
     
     double TimerLLL [max_dimension - min_dimension + 1];
     double TimerLLLNTL [max_dimension - min_dimension + 1];
     bool EqualReducedBasis [max_dimension - min_dimension + 1];
     
     boost::progress_display show_progress(max_dimension);
     
     for (int dimension = min_dimension; dimension <= max_dimension; dimension++){
     
     ++show_progress;
     
     int modulus = 2;
     double epsilon = 0.000001;
     double delta = 1 - epsilon;
     
     int min = 1;
     int max = 500;
     
     Basis MyPrimalBasis (dimension, L2NORM);
     srand (54321);
     for (int i = 0; i < dimension; i++){
     for (int j = i; j < dimension; j++)
     //MyPrimalBasis[i+1][j+1] = power_ZZ(i+1,j);
     MyPrimalBasis[i+1][j+1] = min + (rand() % (int)(max - min + 1));
     }
     if (printMatrices){
     cout << "\nold Basis = " << endl;
     MyPrimalBasis.write();
     }
     
     IntLattice MyLattice (MyPrimalBasis, modulus);
     Reducer MyReducer (MyLattice);
     
     Basis MyPrimalBasisNTL (MyPrimalBasis);
     IntLattice MyLatticeNTL (MyPrimalBasisNTL,modulus);
     Reducer MyReducerNTL (MyLatticeNTL);
     if (printMatrices){
     cout << "\nold NTLBasis = " << endl;
     MyPrimalBasis.write();
     }
     
     clock_t begin = clock();
     
     MyReducer.redLLL(delta,100000000,dimension);
     if (printMatrices){
     cout << "\nnew Basis = " << endl;
     MyLattice.getPrimalBasis().write();
     }
     
     clock_t end = clock();
     
     MyReducerNTL.redLLLNTL(delta);
     if (printMatrices){
     cout << "\nnew NTLBasis = " << endl;
     MyLatticeNTL.getPrimalBasis().write();
     }
     
     clock_t end2 = clock();
     
     double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
     double elapsed_secs2 = double(end2 - end) / CLOCKS_PER_SEC;
     bool test = EqualityTest(MyLatticeNTL.getPrimalBasis(), MyLattice.getPrimalBasis(), dimension);
     
     if (printMatrices)
     cout << "Egalite reduced basis = " << test << endl;
     
     TimerLLL[dimension-1] = elapsed_secs;
     TimerLLLNTL[dimension-1] = elapsed_secs2;
     EqualReducedBasis[dimension-1] = test;
     
     //cout << "Durée pour LatMRG = " << elapsed_secs << endl;
     //cout << "Durée NTl = " << elapsed_secs2 << endl;
     //cout << "Égalité des new Basis = " << test << endl;
     
     } // end dimension loop
     
     
     cout << "\nTimerLLL = [";
     for (int j = 0; j < max_dimension; j++)
     cout << TimerLLL[j] << " ";
     cout << "]" << endl;
     
     cout << "\nTimerLLLNTL = [";
     for (int j = 0; j < max_dimension; j++)
     cout << TimerLLLNTL[j] << " ";
     cout << "]" << endl;
     
     cout << "\nEqualReducedBasis = [";
     for (int j = 0; j < max_dimension; j++)
     if (EqualReducedBasis[j] == 0)
     cout << EqualReducedBasis[j] << "_[" << j+1 << "] ";
     else
     cout << EqualReducedBasis[j] << " ";
     cout << "]" << endl;
     */
    
    /*
     
     // Printing matrices for specific dimensions
     
     int modulus = 2;
     int dim = 11;
     double epsilon = 0.000001;
     double delta = 1 - epsilon;
     
     int min = 1;
     int max = 500;
     
     Basis MyPrimalBasis (dim, L2NORM);
     srand (12345);
     for (int i = 0; i < dim; i++){
     for (int j = i; j < dim; j++)
     //MyPrimalBasis[i+1][j+1] = power_ZZ(i+1,j);
     MyPrimalBasis[i+1][j+1] = min + (rand() % (int)(max - min + 1));
     }
     cout << "old Basis = " << endl;
     MyPrimalBasis.write();
     
     IntLattice MyLattice (MyPrimalBasis, modulus);
     Reducer MyReducer (MyLattice);
     
     Basis MyPrimalBasisNTL (MyPrimalBasis);
     IntLattice MyLatticeNTL (MyPrimalBasisNTL,modulus);
     Reducer MyReducerNTL (MyLatticeNTL);
     cout << "old NTLBasis = " << endl;
     MyPrimalBasis.write();
     
     bool test;
     //test = MyLatticeNTL.getPrimalBasis() == MyLattice.getPrimalBasis();
     test = EqualityTest(MyLatticeNTL.getPrimalBasis(), MyLattice.getPrimalBasis(), dim);
     cout << "*** Equality of input matrices = " << test << " ***" << endl;
     
     MyReducer.redLLL(delta,10000,dim);
     cout << "\nnew Basis = " << endl;
     MyLattice.getPrimalBasis().write();
     
     MyReducerNTL.redLLLNTL(delta);
     cout << "new NTLBasis = " << endl;
     MyLatticeNTL.getPrimalBasis().write();
     
     test = EqualityTest(MyLatticeNTL.getPrimalBasis(), MyLattice.getPrimalBasis(), dim);
     cout << "*** Equality of the output matrices = " << test << " ***" << endl;
     
     */
    
    
    
    
    
    
    
    
    
    /*
     
     bool EqualityTest (mat_ZZ& A, mat_ZZ& B, int dimension){
     bool result = true;
     for (int i = 0; i < dimension; i++){
     for (int j = 0; j < dimension; j++){
     if (A[i][j] != B[i][j])
     result = false;
     if (!result)
     break;
     }
     if (!result)
     break;
     }
     return result;
     }
     
     ZZ TriangularMaxtrixDeterminant (mat_ZZ& A, int dimension){
     ZZ det (1);
     for (int i = 1; i < dimension+1; i++)
     det *= A[i][i];
     return det;
     }
     
     */
    
    
/*
    // UTILISATION DE R
    RInside R(argc, argv);              // create an embedded R instance
    
    R["M"] = lat.toRccpMatrix();                  // eval command, no return
    std::string str =
        "cat('Running ls()\n'); print(ls()); "
        "cat('Showing M\n'); print(M); "
        "cat('Showing colSums()\n'); Z <- colSums(M); print(Z); "
        "Z";
    
    Rcpp::NumericVector v = R.parseEval(str);
*/
    
    
/*
    // Slide Bar of progression
    boost::progress_display show_progress(max_dimension);
    ++show_progress;
 
 
*/
    
    
    
    return 0;
}
