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
#include "latticetester/ReduceFct.cpp"
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ_p.h>






using namespace std;
using namespace LatticeTester;

int main(int argc, const char * argv[]) {




    //RInside R(argc, argv);              // create an embedded R instance

    //R["txt"] = "Hello, world!\n";      // assign a char* (string) to 'txt'

    //R.parseEvalQ("cat(txt)");           // eval the init string, ignoring any returns




    // printing matrices
    bool printMatrices = true;

    // Loop over dimension
    const int min_dimension = 40;
    const int max_dimension = 4;

    
    int dimension = min_dimension;

    //for (int dimension = min_dimension; dimension <= max_dimension; dimension++){


    int modulus = 2;
    double epsilon = 0.000001;
    double delta = 1 - epsilon;

    int min = 10;
    int max = 50;

    IntLatticeBasis MyPrimalLattice (dimension, L2NORM);
    
    // Remplissage aléatoire de la lattice
    srand (2);
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++)
            //MyPrimalBasis[i+1][j+1] = power_ZZ(i+1,j);
            MyPrimalLattice.getBasis()(i,j) = min + (rand() % (int)(max - min + 1));
        
    }
    
    
    if (printMatrices){
        cout << "\nold Basis = " << endl;
        MyPrimalLattice.updateVecNorm();
        MyPrimalLattice.write();
    }
    
    //RMat gram = calculGram( MyPrimalLattice );
    //RMat cho = calculCholeski(MyPrimalLattice, gram);
    //calculCholeskiuntiln(cho, gram, dimension, 1);
    
    long maxcpt = 1000000;
    
    IntLatticeBasis MySecondLattice(MyPrimalLattice);
    
    
    
    redLLL(MyPrimalLattice, delta, maxcpt, dimension);
    
    LLL_XD(MySecondLattice.getBasis(), delta, 0, 0, 0);
    MySecondLattice.setNegativeNorm();
    MySecondLattice.updateVecNorm();
    
    cout << "\new Basis with ReduceFonction = " << endl;
    MyPrimalLattice.write();
    
    cout << "\new Basis with NTL = " << endl;
    MySecondLattice.write();
    
    
    
    
    
    
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
