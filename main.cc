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
#include <NTL/ctools.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include "latticetester/Reducer.h"


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
    //MScal n = 20;
    //long t = 11;

   // MVect a_vector(5,2);
   // a_vector[2] = 1;
   // a_vector[1] = 1;

    //cout << a_vector[0] << endl;
    //a = { 1,3,9 };
    //PrimeType tes = nombre.getStatus();
    //nombre.setStatus(tes);
    //cout << nombre.isPrime(s, t) << endl;
    //Rank1Lattice reseau(n, a_vector, 10);
    //reseau.buildBasis(4);


    int s = 0;
    /*
     for(int i = 0; i<2; i++){
        for(int j = 0; j<2; j++){
            base_lattice(i, j) = s;
            ++s;
        }
    }
     */

    Basis W(4);
    /*

    W[0][0] = 1;  W[1][2] = 0;  W[1][3] = 10;
    W[2][1] = 2;  W[2][2] = 1;  W[2][3] = 0;
    W[3][1] = 2;  W[3][2] = 0;  W[3][3] = 1;
    */
    
    ZZ scal;
    scal = 0;
    
    W[0][0] = 1;  W[0][1] = 0;  W[0][2] = 3;  W[0][3] = 0;
    W[1][0] = 0;  W[1][1] = 1;  W[1][2] = 1;  W[1][3] = 0;
    W[2][0] = 0;  W[2][1] = 0;  W[2][2] = 1;  W[2][3] = 0;
    W[3][0] = 0;  W[3][1] = 0;  W[3][2] = 0;  W[3][3] = 1;
    
    Basis W1(W);
    
  
    
   // W.write();
    
    
    
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
    

    /*
    reseau1.trace("hey");

    Triangularization(W1, V1, 4, 4, 19);
    CalcDual(V1, W1, 4, 19);

    V1.write();
    W1.write();
    //cout << V1.toString() << endl;
    //cout << W1.toString() << endl;

    //IntLattice reseau1(V, 11);

    IntLattice reseau2(V1, W1, 19);

    cout << reseau1.baseEquivalence(reseau2) << endl;

    */
    //cout <<trace << endl;




    //BMat C(3,3);

    //C = prod(V,W);
    //cout << C << endl;

    //cout << C.toString() << endl;

#ifdef WITH_NTL

    cout << "yeah" << endl;

#endif






    //IntLattice reseau(primal_basis, 2);

    //cout << reseau.checkDuality() << endl;

    //int myints[]= {0, 1};
    //Coordinates second (myints,myints+2);

    /*while (!second.empty()) {
        std::cout << ' ' << *second.begin();
        second.erase(second.begin());
    } */
    /*for(Coordinates::const_iterator iter = second.begin(); iter != second.end(); ++iter)
    {
        cout << *iter << endl;
    }*/
    //reseau.buildProjection(second);

    //int lin = 2;
    //int col = 2;
    //int m = 5;

    //MScal T1, T2, T3, T4, T5, T6, T7, T8;

   /*for (int j = 0; j < col; j++) {
      for (int i = 0; i < lin; i++)
         Modulo (W(i,j), m, W(i,j));
      int r = j;
      while (r < lin-1) {
         while (IsZero (W(r,j)) && r < lin-1)
            ++r;
         if (r < lin-1) {
            int s = r + 1;
            while (IsZero (W(s,j)) && s < lin-1)
               ++s;
            if (!IsZero (W(s,j))) {
               Euclide (W(r,j), W(s,j), T1, T2, T3, T4, W(s,j)); //pivot de gausse?
               clear (W(r,j));
               for (int j1 = j + 1; j1 < col; j1++) {
                  T5 = T1 * W(r,j1);
                  T6 = T2 * W(s,j1);
                  T7 = T3 * W(r,j1);
                  T8 = T4 * W(s,j1);
                  W(s,j1) = T5 + T6;
                  Modulo (W(s,j1), m, W(s,j1));
                  W(r,j1) = T7 + T8;
                  Modulo (W(r,j1), m, W(r,j1));
               }
            } else {
               for (int j1 = j; j1 < col; j1++) {
                  std::swap (W(r,j1), W(s,j1));
               }
            }
            r = s;
         }
      }
      cout << W.toString() << endl;
      if (IsZero (W(lin-1,j))) {
         for (int j1 = 0; j1 < col; j1++) {
            if (j1 != j)
               clear (V(j,j1));
            else
               V(j,j1) = m;
         }
      } else {
         Euclide (W(lin-1,j), m, T1, T2, T3, T4, V(j,j));
         for (int j1 = 0; j1 < j; j1++)
            clear (V(j,j1));
         for (int j1 = j + 1; j1 < col; j1++) {
            T2 = W(lin-1,j1) * T1;
            Modulo (T2, m, V(j,j1));
         }
          T1 = m/V(j,j);
         //Quotient (m, V(j,j), T1);
         for (int j1 = j + 1; j1 < col; j1++) {
            W(lin-1,j1) *= T1;
            Modulo (W(lin-1,j1), m, W(lin-1,j1));
         }
      }
   }*/









    //char *mess("hey");
    //reseau.trace(mess);
    return 0;
}
