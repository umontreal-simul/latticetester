#include "latmrg/LatConfig.h"
#include "latticetester/Util.h"

using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

LatConfig::LatConfig()
{
   fileName.reserve(MAX_WORD_SIZE);
   //a.SetLength(1+maxDim);
   comp = 0;
   genType = 0;
   td = 0;
   invertF = false;
   detailF = 0;
}

LatConfig::~LatConfig()
{
   kill();
}

void LatConfig::setJ(int j)
{
   kill();
   genType = new GenType[j];
   comp = new MRGComponent * [j];
}

void LatConfig::kill()
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

void LatConfig::write()
{
   cout << "readGenFile: " << boolalpha << readGenFile << endl;
   if (readGenFile)
      cout << "fileName: " << fileName << endl;
   cout << "J: " << J << endl << endl;

   for (int i = 0; i < J; i++) {
      if (J > 1)
         cout << "================ Component " << i+1 << " =================\n";
      cout << "   genType: " << toStringGen (genType[i]) << endl;
      cout << "   m: " << comp[i]->getM() << endl;
      cout << "   verifyM: " << boolalpha << verifyM << endl;
      cout << "   k: " << comp[i]->k << endl;
      cout << "   a: " << toString<MVect>(comp[i]->a, comp[i]->k) << endl << endl;
   }

//   cout << "fromDim: " << fromDim << endl;
//   cout << "toDim:   " << toDim << endl;
   cout << "td:   " << toString (td, 0, d) << endl;
   cout << "criter: " << toStringCriterion(criter) << endl;
   cout << "norma: " << toStringNorma (norma) << endl;
   cout << "latType: " << toStringLattice (latType) << endl;
   if (dualF)
      cout << "lattice: DUAL" << endl;
   else
      cout << "lattice: PRIMAL" << endl;
   cout << "lacGroupSize: " << lacGroupSize << endl;
   cout << "lacSpacing: " << lacSpacing << endl;
   cout << "maxPeriod: " << boolalpha << maxPeriod << endl;
   cout << "verifyP: " << boolalpha << verifyP << endl;
   cout << "maxNodesBB: " << maxNodesBB << endl;
   cout << "outputType: " << toStringOutput (outputType) << endl;

}

} // namespace
