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
#include "latticetester/ReduceFct.cpp"
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ_p.h>






using namespace std;
using namespace LatticeTester;

int main(int argc, const char * argv[]) {
    
    

    // printing matrices
    bool printMatrices = true;
    
    // Loop over dimension
    const int min_dimension = 7;
    const int max_dimension = 7;
    
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
        int max = 10;
        
        IntLatticeBasis MyPrimalLattice (dimension, L2NORM);
        srand (1);
        for (int i = 0; i < dimension; i++){
            for (int j = i; j < dimension; j++)
                //MyPrimalBasis[i+1][j+1] = power_ZZ(i+1,j);
                MyPrimalLattice.getBasis()(j,i) = min + (rand() % (int)(max - min + 1));
        }
        MyPrimalLattice.getBasis()(0,0) = 1;
        if (printMatrices){
            cout << "\nold Basis = " << endl;
            MyPrimalLattice.updateVecNorm();
            MyPrimalLattice.write();
        }
        RMat gram = calculGram( MyPrimalLattice );
        RMat cho = calculCholeski(MyPrimalLattice, gram);
        calculCholeskiuntiln(cho, gram, dimension, 1);
        
        RMat D(dimension,dimension);
        for(int i = 0; i < dimension; i++){
            for(int j = 0; j < dimension; j++){
                D(i,j) = 0;
                for(int k = 0; k < dimension; k++){
                    D(i,j) += cho(i,k)*cho(j,k);
                }
            }
        }
        //cho(0,0) = 1;
        //RMat cho2 = prod(cho, trans(cho));
        //cout << gram << "\n \n" << endl;
        //cout << cho  << "\n \n" << endl;
        //cout << D << "\n \n" << endl;
        IntLatticeBasis lat = pairwiseRed (MyPrimalLattice, 0, 0);
        lat.write();
        bool test(true);
        for(int i = 0; i < dimension; i++){
            for(int j = 0; j < dimension; j++){
                if((D(i,j) - gram(i,j)) > 0.0000001){
                    test = false;
                    cout << "D(i,j) = " << D(i,j) << " et gram(i,j) = " << gram(i,j) << endl;
                }
            }
        }

        
        //cout << test << endl;
    }/*
      
      
    
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
        
        cout << "Egalite reduced basis = " << test << endl;
        
        TimerLLL[dimension-1] = elapsed_secs;
        TimerLLLNTL[dimension-1] = elapsed_secs2;
        EqualReducedBasis[dimension-1] = test;
        
        //cout << "Durée pour LatMRG = " << elapsed_secs << endl;
        //cout << "Durée NTl = " << elapsed_secs2 << endl;
        //cout << "Égalité des new Basis = " << test << endl;
        
    } // end dimension loop
    
    /*
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
    
    
    // Création de la première base.
    /*
    Basis W(4);
    
    W[0][0] = 1;  W[0][1] = 0;  W[0][2] = 3;  W[0][3] = 0;
    W[1][0] = 0;  W[1][1] = 1;  W[1][2] = 1;  W[1][3] = 0;
    W[2][0] = 0;  W[2][1] = 0;  W[2][2] = 1;  W[2][3] = 0;
    W[3][0] = 0;  W[3][1] = 0;  W[3][2] = 0;  W[3][3] = 1;
    
    
    Basis W1(W);
    
    
    Basis V(4);

    V[0][0] = 19; V[0][1] = -40;  V[0][2] = 23;  V[0][3] = 4;
    V[1][0] = 0;  V[1][1] = -1;  V[1][2] = 2;  V[1][3] = 2;
    V[2][0] = 0;  V[2][1] = 0;  V[2][2] = 19; V[2][3] = 0;
    V[3][0] = 0;  V[3][1] = 10;  V[3][2] = -20;  V[3][3] = -1;

    Basis V1(V);

    IntLattice reseau1(W, 19);
    reseau1.getPrimalBasis().updateVecNorm();
    reseau1.getPrimalBasis().write();

    
    Reducer teston(reseau1);
    
    teston.redLLL(0.99, 1000000, 2);
    
    reseau1.getPrimalBasis().write();
    
    BKZ_FP(W1);
    W1.write();
     */
    
    //boost::numeric::ublas::identity_matrix<long> m(3);
    
    
    IntLatticeBasis A(10);
    
    A.getBasis()(0,2) = 2;
    
    A.updateVecNorm();
    
    
    
    
    //A.write();
    
    
    
    //cout << A.getBasis()[1][2] << endl;
    
    //BMat hey(boost::numeric::ublas::identity_matrix<long>(10));
    
    //cout << hey << endl;



#ifdef WITH_NTL

    cout << "yeah" << endl;

#endif


    return 0;
}
