//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//


#include <iostream>
#include "latticetester/Util.h"
#include "latticetester/Basis.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"


using namespace std;
using namespace LatticeTester;

int main(int argc, const char * argv[]) {
 //   Basis base_test(2,7);
//    for(NormType k = SUPNORM; k <=  ZAREMBANORM; k++){
//        cout << k << endl;
//    }

 //   enum NormType { SUPNORM = 1, L1NORM = 3, L2NORM = 3, ZAREMBANORM = 4 };

//    cout << norme_test << "  " << norme_test2 << "   " << norme_test3 << "   " << norme_test4 << endl;
//    cout << base_test(1, 0) << endl;
    MScal n = 20;
    //long t = 11;

    MVect a_vector(5,2);
    a_vector[2] = 1;
    a_vector[1] = 1;

    //cout << a_vector[0] << endl;
    //a = { 1,3,9 };
    //PrimeType tes = nombre.getStatus();
    //nombre.setStatus(tes);
    //cout << nombre.isPrime(s, t) << endl;
    //Rank1Lattice reseau(n, a_vector, 10);
    //reseau.buildBasis(4);
    
    Basis base_lattice(2);
    int s = 0;
    /* 
     for(int i = 0; i<2; i++){
        for(int j = 0; j<2; j++){
            base_lattice(i, j) = s;
            ++s;
        }
    }
     */
    base_lattice(0,0) = 1;
    base_lattice(1,0) = 0;
    base_lattice(0,1) = 0;
    base_lattice(1,1) = 1;
    
    
    IntLattice reseau(base_lattice, 2);
    
    //cout << reseau.checkDuality() << endl;
    
    int myints[]= {1};
    Coordinates second (myints,myints+5);
    
    /*while (!second.empty()) {
        std::cout << ' ' << *second.begin();
        second.erase(second.begin());
    } */
    cout << second.size() << endl;
    reseau.buildProjection(second);
    
    
    
    
    
    
    
    //char *mess("hey");
    //reseau.trace(mess);
    return 0;
}
