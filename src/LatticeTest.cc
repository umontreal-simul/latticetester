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

#include "latticetester/LatticeTest.h"
#include "latticetester/Util.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaRogers.h"

using namespace std;
using namespace NTL;

namespace LatticeTester
{


//===========================================================================

LatticeTest::LatticeTest (Reducer & reducer, NormaType normaType, int alpha) 
{
   m_reducer = &reducer;
   m_normaType = normaType;
   initNormalizer(normaType, alpha);

}

//===========================================================================

LatticeTest::~LatticeTest ()
{
   delete m_normalizer;
}

//===========================================================================

bool LatticeTest::performTest (double fact, long blockSize)
{
   bool result;
   NormType norm = m_reducer->getIntLatticeBasis().getNorm();
   int dim = m_reducer->getIntLatticeBasis().getDim();

   // finding shortest non-zero vector in the lattice
   m_reducer->shortestVectorWithBKZ(norm, fact, blockSize);
   //m_reducer.redBKZ(fact,blockSize); // BKZ reduction is performed
   //result = m_reducer.redBB0(norm) // Then Brand-and-Bpund procedure is performed

   // calculating the Figure of Merit
   double length = conv<double>(m_reducer->getMinLength());
   if (norm == L2NORM)
      length *= length;

   double maxLength = m_normalizer->getPreComputedBound(dim);

   if (norm == L2NORM)
      m_merit = sqrt(length / maxLength);
   else 
      m_merit = length / maxLength;

   return result;
  
}

//===========================================================================

void LatticeTest::initNormalizer (NormaType norma, int alpha)
{
   int dim = m_reducer->getIntLatticeBasis().getDim(); 

   // PW_TODO
   // ok si la matrice est directement construite comme m-dual mais probleme si 
   // travail direct sur primale re-scaled ?
   // version avec logDensity en parametre (pou m^k sans besoin de calcul det)

   RScal logDensity;
   logDensity = - log( determinant(m_reducer->getIntLatticeBasis().getBasis()) );

   switch (norma) {
      case BESTLAT:
         m_normalizer = new NormaBestLat (logDensity, dim);
         break;
      case LAMINATED:
         m_normalizer = new NormaLaminated (logDensity, dim);
         break;
      case ROGERS:
         m_normalizer = new NormaRogers (logDensity, dim);
         break;
      case MINKL1:
         m_normalizer = new NormaMinkL1 (logDensity, dim);
         break;
      case MINKOWSKI:
         m_normalizer = new NormaMinkowski (logDensity, dim);
         break;
      case NORMA_GENERIC:
         m_normalizer = new Normalizer (logDensity, dim, "Norma_generic");
         break;
      case PALPHA_N:
         m_normalizer = new NormaPalpha (m_reducer->getIntLatticeBasis().getModulo(), alpha, dim);
         //PW_TODO : c'est bien ça ?
         break;
      default:
         cout << "normalizer:   no such case";
         exit (2);
   }
}

//=========================================================================







} // end namespace LatticeTester
