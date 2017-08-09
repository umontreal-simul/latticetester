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

#if 0
int main (int argc, char *argv[])
{
   // Erwan
   //string testLocation = "ton_path_vers_dossier_input_files";
   //string testLocation = "ton_path_vers_dossier_input_files/input_file";
   
   // Paul
   //string testLocation = "/Users/paulwambergue/UdeM/latticetester/inputTestFiles";
   string testLocation = "/Users/paulwambergue/UdeM/latticetester/inputTestFiles/latticeAnalysis_test3";
   
   //PW_TODO ça marche pas
   
   struct stat buf; // properties of a file or directory
   LatticeAnalysis latAnalysis;
   int status = 0;
   
   stat(testLocation.c_str(), &buf);
   
   if (0 != S_ISDIR(buf.st_mode)) //directory
      status |= latAnalysis.doTestFromDirectory (testLocation.c_str());
   else { //file
      string dataname(testLocation.c_str());
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      status |= latAnalysis.doTestFromInputFile (testLocation.c_str());
   }
   
   return 0;
}
#endif

