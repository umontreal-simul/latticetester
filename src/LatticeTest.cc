#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <typeinfo>
#include <list>
#include "NTL/ZZ.h"
#include "NTL/LLL.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"

#include "latticetester/Util.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"


using namespace std;
using namespace NTL;

///**************************************************************************

namespace LatticeTester
{

LatticeTest::LatticeTest (LatMRG::IntLattice * lat): m_merit(lat->getDim())
{
   m_lat = lat;
   m_dualF = true;
   m_invertF = false;
   m_maxAllDimFlag = true;
   m_detailF = 0;
   Reducer::maxNodesBB = m_maxNodesBB = 10000000;
   timer.init ();
}


//===========================================================================

LatticeTest::~LatticeTest ()
{
};


//===========================================================================

bool LatticeTest::test (int minDim, int maxDim, double minVal[],
                        const double* weights)
{
   throw "Not implemented with weights";
   // compiler warnings
   minDim = maxDim = -1;
   minVal[0] =  weights[0];
}


//===========================================================================

void LatticeTest::resetFromDim(int order, int & fromDim)
{
   if (fromDim <= order) {
      // fromDim = order + 1;
     // cout << "********* WARNING LatticeTest:  fromDim is <= order.\n";
      //  cout << " It is now reset to: fromDim = order + 1\n" << endl;
   }
}


//===========================================================================

void LatticeTest::setDualFlag (bool dualF)
{
   m_dualF = dualF;
   m_lat->fixLatticeNormalization (dualF);
}


//===========================================================================

void LatticeTest::setMaxAllDimFlag (bool maxAllDimF)
{
   m_maxAllDimFlag = maxAllDimF;
}


//===========================================================================

void LatticeTest::setInvertFlag (bool invertF)
{
   m_invertF = invertF;
}


//===========================================================================

void LatticeTest::setDetailFlag (int d)
{
   m_detailF = d;
}


//===========================================================================

void LatticeTest::setMaxNodesBB (long maxNodesBB)
{
   Reducer::maxNodesBB = m_maxNodesBB = maxNodesBB;
}


//===========================================================================


} // end namespace LatticeTester
