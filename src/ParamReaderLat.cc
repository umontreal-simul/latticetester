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

void ParamReaderLatticeTester::read (LatConfig & config)
{
   getLines ();
   unsigned int ln = 1;

   readCriterionType (config.criter, ln, 0);
   if (config.criter == SPECTRAL)
      readNormaType (config.norma, ln, 1);

   if (config.criter == PALPHA) {
      readCalcType (config.calcPalpha, ++ln, 0);
      config.J = 1;
   } else {
      readNormType (config.norm, ++ln, 0);
      readBool (config.readGenFile, ++ln, 0);
      if (config.readGenFile)
         readString (config.fileName, ln, 1);
      readInt (config.J, ++ln, 0); // J
   }

   config.setJ (config.J);

   MScal m;
   long m1(0), m2(0), m3(0);
   int k(0);
   MVect a;
   int maxOrder = 0;
   string coefkind (" ");  // = NOCOND, EQUAL, NONZERO

   for (int i = 0; i < config.J; i++) {
      readGenType (config.genType[i], ++ln, 0);
      readNumber3 (m, m1, m2, m3, ++ln, 0);
      if (config.criter == PALPHA)
         k = 1;
      else
         readInt (k, ++ln, 0);
      if (maxOrder < k)
         maxOrder = k;
      a.SetLength (k);
      SetZero (a, k);

      if (config.criter != PALPHA) {
         readString (coefkind, ++ln, 0);
      }
      if (0 == strcasecmp(coefkind.c_str(), "NONZERO")) {
         int s;
         readInt (s, ln, 1);
         int* T = new int[s+3];
         T[0] = 0;
         readIntVect (T, ln++, 2, s+3, 0);
         readMVect (a, ln, 1, s, 0);
         //   ln++;
         for (int j = s; j >= 1; j--) {
            int r = T[j];
            a[r] = a[j-1];
            a[j-1] = 0;
         }
         delete[] T;
      } else if (0 == strcasecmp(coefkind.c_str(), "EQUAL")) {
         int s;
         readInt (s, ln, 2);
         int* T = new int[s+3];
         readIntVect (T, ln++, 3, s+3, 1);
         int p0 = 1;
         for (int j = 1; j <= s; j++) {
            int r = T[j];
            readMScal (a[r], ln, j);
            for (int p = p0; p < r; p++) {
               a[p] = a[r];
            }
            p0 = r + 1;
         }
         delete[] T;
      } else   // NoCond
         readMVect (a, ++ln, 0, k, 0);
      checkBound (m, a, k);

      if (config.genType[i] == KOROBOV) {
         if (1 != k)
             MyExit(1, "KOROBOV must have k = 1");
         config.comp[i] = new MRGComponent (m, a, 1);
      } else if (config.genType[i] == RANK1) {
         config.comp[i] = new MRGComponent (m, a, k);
      } else if (config.genType[i] == MRG || config.genType[i] == LCG) {
         config.comp[i] = new MRGComponent (m1, m2, m3, a, k);
         if (1 == k)
            config.comp[i]->module.reduceM (config.comp[i]->a[0]);
      } else if (config.genType[i] == MWC) {
         config.comp[i] = MRGComponentFactory::fromMWC (m, a, k);
      } else
         assert(0);
   }

   readInt (config.d, ++ln, 0);
   if (config.d < 1)
      MyExit (1, "ParamReaderLat:   config.d < 1");
   config.td = new int[config.d];
   readIntVect (config.td, ++ln, 0, 1 + config.d, 0);
   int fromDim = config.td[0];
   int toDim = config.td[1];

   if (config.criter == PALPHA) {
      readBool (config.primeM, ++ln, 0);
      readBool (config.verifyM, ln, 1);
      readBool (config.maxPeriod, ++ln, 0);
      readBool (config.verifyP, ln, 1);
      readInt (config.alpha, ++ln, 0);
      readInt (config.seed, ++ln, 0);
      config.Beta.reserve(1 + toDim);
      ++ln;
      for (int i = 0; i <= toDim; i++)
         readDouble (config.Beta[i], ln, 1 + i);
      config.lacunary = false;

   } else {
      readBool (config.dualF, ++ln, 0);
      readLatticeType (config.latType, ++ln, 0);
      if ((config.genType[0] == RANK1 || config.genType[0] == KOROBOV) &&
           config.latType != FULL) {
         MyExit (1,
         "ParamReaderLat:   latType must be FULL for KOROBOV or RANK1 lattices");
         }
      checkPrimePower (config.latType, m2, m3, k);
      if (config.latType == ORBIT) {
         MyExit (1, "case ORBIT is not finished");
         readOrbit (config.J, config.comp, ++ln);
         --ln;
      }
      readLacunary (k, fromDim, toDim, ln, config.lacunary,
                    config.lacGroupSize, config.lacSpacing, config.Lac,
                    config.genType[0]);
      if (config.genType[0] != RANK1 && config.td[0] <= maxOrder &&
           !config.lacunary)
         config.td[0] = maxOrder + 1;

      readLong (config.maxNodesBB, ++ln, 0);
      readBool (config.invertF, ++ln, 0);
      readInt (config.detailF, ++ln, 0);
   }
   readOutputType (config.outputType, ++ln, 0);
}


//===========================================================================

}
