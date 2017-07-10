//
//  normalizationTester.cc
//  Lattice Tester
//
//

/**
 * This main program aims to test the normalization
 * as performed in LatticeTester
 */

// include headers
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

// include LatticeTester headers
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"
#include "latticetester/Types.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/LatticeTest.h"

// include NTL headers
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

// include Boost headers
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

// include local files
#include "SimpleMRG.h"
#include "Tools.h"


using namespace std;
using namespace NTL;
using namespace LatticeTester;


// parameters 
//----------------------------------------------------------------------------------------

/*
 * Modulo of the Basis according to L'Écuyer's paper.
 * FullRandomMatrix flag must be false.
 */
//const ZZ modulusRNG = conv<ZZ>("32749");
const ZZ modulusRNG = power_ZZ(2,63) - 2247;
//const ZZ modulusRNG = power_ZZ(2, 2) - 1;
//const ZZ modulusRNG = power_ZZ(2, 3) - 1;
//const ZZ modulusRNG = power_ZZ(2, 5) - 1;
//const ZZ modulusRNG = power_ZZ(2, 7) - 1;
//const ZZ modulusRNG = power_ZZ(2, 13) - 1;
//const ZZ modulusRNG = power_ZZ(2, 17) - 1;
//const ZZ modulusRNG = power_ZZ(2, 19) - 1;
//const ZZ modulusRNG = power_ZZ(2, 61) - 1;

/*
 * Order of the Basis according to L'Écuyer's paper.
 * FullRandomMatrix flag must be false.
 */
const int order = 3;

/*
 * The Dimension to be analysed.
 * Must be int value.
 */
const int dimension = 6;

/*
 * a/b is the value of the delta in the LLL and BKZ
 * Reduction. NTL offers the possibility to compute
 * a LLL Reduction with the exact delta. We have noticed
 * only minor differences with this option.
 */
const double delta = 0.99999;
const double epsilon = 1.0 - delta;

/*
 * Block Size in the BKZ algorithm. See NTL documention
 * for further information.
 */
const long blocksize = 20;

/*
 * Maximum number of Nodes in the Branch-and-Bound.
 */
const int maxNodesBB = 1000000;



//----------------------------------------------------------------------------------------

int main ()
{
    // Seed initialization
    ZZ seedZZ = conv<ZZ>(123456789 * dimension);
    int seed = (123456789 * dimension);

    //int seed_dieter = (iteration+1) * dimension * 12342 * (nb_error+1) ;

    // coefficients (a_i)
    vec_ZZ a;
    a.SetLength(dimension);
    a[0] = conv<ZZ>("1145902849652723");
    a[1] = conv<ZZ>("0");
    a[2] = conv<ZZ>("-1184153554609676");

    // parameters printing
    cout << "\n------------------------" << endl;
    cout << "  dimension = " << dimension << endl;
    cout << "  epsilon = " << epsilon << endl;
    cout << "  blocksize = " << blocksize << endl;
    cout << endl;
    cout << "  modulo = " << modulusRNG << endl;
    cout << "  ordre = " << order << endl;
    cout << "  a = [";
    for (int i = 0; i < order-1; i++)
        cout << a[i] << ", ";
    cout << a[order-1] << "]" << endl;
    cout << "------------------------" << endl;

    // lattice basis creation
    BMat V;
    BMat W;

    V = CreateRNGBasis (modulusRNG, order, dimension, a, seedZZ);
    W = Dualize (V, modulusRNG, order);
    cout << "V =\n" << V << endl;
    cout << "W =\n" << W << endl;

    // attention
    IntLatticeBasis basis (W, V, modulusRNG, dimension);
    Reducer red (basis);



    // Ancienne façon de faire ---------------------------
    // ---------------------------------------------------

    /*
    // BKZ NTL
    red.redBKZ(delta, blocksize);
    basis.setNegativeNorm();
    basis.updateVecNorm();
    basis.sort(0);
    cout << "reduced =\n" << basis.getBasis() << endl;

    // Branch and Bound
    red.redBB0(L2NORM);
    basis.setNegativeNorm();
    basis.updateVecNorm();
    basis.sort(0);
    cout << "reduced bis =\n" << basis.getBasis() << endl;

    // normalizer
    RScal logDensity;
   
    logDensity = - log(determinant(W));
    //LogDensity = log( conv<double>(determinant(V)) / conv<double>(power(modulusRNG,dimension)) );
    cout << "Density = " << exp(logDensity) << endl;
   
    NormaBestLat normalizer (logDensity, 48);
   
    double shortestLength = conv<double>(red.getMinLength());
    shortestLength *= shortestLength; // squared

    // Results printing : a comparer table 7
    cout << "\nShortest Length real = " << shortestLength << endl;
    cout << "Bound on Length = " << normalizer.getBound(dimension) << "(" << normalizer.getPreComputedBound(dimension) << ")" << endl;

    cout << "FoM = " << sqrt( shortestLength / normalizer.getBound(dimension) ) << endl;
    // sqrt because L2 norm
    */



    // Nouvelle façon de faire ---------------------------
    // ---------------------------------------------------

    NormaType normaType = BESTLAT;
    LatticeTest latTest (red, normaType);
    latTest.performTest();

    cout << "FoM = " << latTest.getMerit() << endl;

    return 0;
}

