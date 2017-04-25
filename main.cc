//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//


#include <iostream>
#include "latcommon/Util.h"
#include "latcommon/Base.h"
#include "latcommon/Const.h"
#include "latcommon/Types.h"
#include "latcommon/IntFactor.h"
#include "latcommon/IntLattice.h"
#include "latcommon/Rank1Lattice.h"


using namespace std;
using namespace LatCommon;

int main(int argc, const char * argv[]) {
    Base base_test(2,7);
//    for(NormType k = SUPNORM; k <=  ZAREMBANORM; k++){
//        cout << k << endl;
//    }
    
    enum NormType { SUPNORM = 1, L1NORM = 3, L2NORM = 3, ZAREMBANORM = 4 };
//    NormType norme_test;
//    NormType norme_test2(L1NORM);
//    NormType norme_test3;
//    NormType norme_test4;
//    cout << norme_test << "  " << norme_test2 << "   " << norme_test3 << "   " << norme_test4 << endl;
//    cout << base_test(1, 0) << endl;
    MScal n = 20;
    //long t = 11;

    MVect a_vector(5,2);
    a_vector[2] = 1;
    a_vector[1] = 1;
    cout << a_vector[0] << endl;
    
    //cout << a_vector[0] << endl;
    //a = { 1,3,9 };
    //PrimeType tes = nombre.getStatus();
    //nombre.setStatus(tes);
    //cout << nombre.isPrime(s, t) << endl;
    Rank1Lattice reseau(n, a_vector, 10);
    reseau.buildBasis(4);
    Base base_lattice = reseau.getPrimalBasis();
    for(int i = 0; i<5; i++){
        for(int j = 0; j<5; j++){
            cout << base_lattice(i, j);
        }
        cout << endl;
        
    }
    //char *mess("hey");
    //reseau.trace(mess);
    return 0;
}
