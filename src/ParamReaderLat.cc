/**
ParamReaderLat.cc for ISO C++
version
modified
authors: Hicham Wahbi
      Frederik Rozon
      Richard Simardr
*/

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/ParamReaderLat.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstring>



using namespace std;
using namespace LatticeTester;



namespace LatticeTester
{

//===========================================================================

ParamReaderLatticeTester::ParamReaderLatticeTester (string fileName): ParamReader (fileName)
{}


//===========================================================================

ParamReaderLatticeTester::~ParamReaderLatticeTester ()
{}


//===========================================================================

void ParamReaderLatticeTester::read (LatticeTesterConfig & config)
{
   getLines ();
   unsigned int ln = 1;

   readNormType (config.norm, ++ln, 0);
   readNormaType (config.normalizer, ++ln, 0);
   readPreRed (config.prereduction, ++ln, 0);

   if(config.prereduction == BKZ){
      readDouble (config.fact, ++ln, 0);
      readInt (config.blocksize, ++ln, 0);
   }
   else if(config.prereduction == PreRedDieter){

   else if(config.prereduction == LLL){
      readDouble (config.fact, ++ln, 0);
   }

   readInt (config.dimension, ++ln, 0);

   readBMat(config.basis, ++ln, 0, config.dimension);

   readLong (config.maxNodesBB, ++ln, 1);

   readOutputType(config.outputType, ++ln, 0);




   //config.write();


}


//===========================================================================

}
