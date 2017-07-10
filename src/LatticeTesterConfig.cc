#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/Util.h"

using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

LatticeTesterConfig::LatticeTesterConfig()
{
   fileName.reserve(MAX_WORD_SIZE);
   //a.SetLength(1+maxDim);
   comp = 0;
   genType = 0;
   td = 0;
   invertF = false;
   detailF = 0;
}

LatticeTesterConfig::~LatticeTesterConfig()
{
   kill();
}


void LatticeTesterConfig::kill()
{
   if (comp != 0) {
      for (int i = 0; i < J; i++) {
         delete comp[i];
      }
      delete[] comp;
   }

   if (genType != 0)
      delete[] genType;
   delete[] td;
}

void LatticeTesterConfig::write()
{
   cout << "readGenFile: " << boolalpha << readGenFile << endl;
   if (readGenFile)
      cout << "fileName: " << fileName << endl;

   cout << "Dimension:   " << dimension << endl;
   cout << "normalizer: " << toStringNorma (normalizer) << endl;
   cout << "Prereduction: " << toStringPreReduction (prereduction) << endl;
   cout << "Lattice Basis" << basis << endl;
   cout << "maxNodesBB: " << maxNodesBB << endl;
   cout << "outputType: " << toStringOutput (outputType) << endl;

}

} // namespace
