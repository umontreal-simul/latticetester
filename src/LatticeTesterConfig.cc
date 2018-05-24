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

   switch(test) {

      case SPECTRAL:
         cout << "\n----- LatticeTesterConfig::write -----" << endl;
         cout << "Criterion: " << toStringCriterion(test) << endl;
         cout << "Norm: " << toStringNorm(norm) << endl;
         cout << "Normalizer: " << toStringNorma (normalizer) << endl;

         cout << "Prereduction: " << toStringPreRed (prereduction);
         cout << " (with " << toStringPrecision(precision) << " precision";
         cout << " and fact = " << fact << ")" << endl;

         cout << "Dimension of the basis: " << dim << endl;
         cout << "Lattice Basis = \n" << basis << endl;
         cout << "maxNodesBB: " << maxNodesBB << endl;
         break;

      case PALPHA:
         cout << "\n----- LatticeTesterConfig::write -----" << endl;
         cout << "Criterion: " << toStringCriterion(test) << endl;
         cout << "alpha: " << alpha << endl;
         cout << "modulo: " << modulo << endl;

         cout << "Dimension of the basis: " << dim << endl;
         cout << "Lattice Basis = \n" << basis << endl;
         cout << "maxNodesBB: " << maxNodesBB << endl;
         break;

      case BEYER:
         cout << "\n----- LatticeTesterConfig::write -----" << endl;
         cout << "Criterion: " << toStringCriterion(test) << endl;
         cout << "Dimension of the basis: " << dim << endl;
         cout << "Lattice Basis = \n" << basis << endl;
         cout << "maxNodesBB: " << maxNodesBB << endl;
         break;

      default:
         MyExit(1, "LatticeTesterConfig::write:  NO SUCH CASE");
         break;
   }
}

} // namespace LatticeTester
