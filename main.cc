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
        //RMat gram = calculGram( MyPrimalLattice );
        //RMat cho = calculCholeski(MyPrimalLattice, gram);
        //calculCholeskiuntiln(cho, gram, dimension, 1);
        
        
        IntLatticeBasis lat = preRedDieter(MyPrimalLattice, 0);
        lat.write();
        bool test(true);
        
    }

    return 0;
}
