#include "latticetester/LatConfig.h"
#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/Util.h"

using namespace std;
using namespace LatticeTester;

namespace LatticeTester
{

LatticeTesterConfig::LatticeTesterConfig()
{
   fileName.reserve(MAX_WORD_SIZE);
   //a.SetLength(1+maxDim);
   detailF = 0;
}

LatticeTesterConfig::~LatticeTesterConfig()
{
   kill();
}

void LatticeTesterConfig::kill()
{
}

void LatticeTesterConfig::write()
{
   cout << "Dimension:   " << dimension << endl;
   cout << "normalizer: " << toStringNorma (normalizer) << endl;
   cout << "Prereduction: " << toStringPreRed (prereduction) << endl;
   cout << "Lattice Basis" << basis << endl;
   cout << "maxNodesBB: " << maxNodesBB << endl;
   cout << "outputType: " << toStringOutput (outputType) << endl;

}

} // namespace
