//
//  inputFileTest.cpp
//  main program to test LatticeTester execution 
//  on an external input file.
//


// Include Header
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

// Include LatticeTester Header
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"
#include "latticetester/Types.h"
#include "latticetester/ParamReader.h"
#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/LatticeAnalysis.h"

// Include NTL Header
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

// Include Boost Header
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

// Include Random Generator of MRG Matrix an tools
#include "SimpleMRG.h"
#include "Tools.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;


//==================================================================================

int main (int argc, char *argv[])
{
   struct stat buf; // properties of a file or directory
   LatticeAnalysis latAnalysis;
   int status = 0;

   stat("latticeAnalysis_test1", &buf);
   if (0 != S_ISDIR(buf.st_mode)) // directory
      status |= latAnalysis.doTestFromDirectory ("latticeAnalysis_test1");
   else {
      string dataname("latticeAnalysis_test1");
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);

      status |= latAnalysis.doTestFromInputFile ("/Users/paulwambergue/UdeM/latticetester/latticeAnalysis_test1");
   }

   return 0;
}

