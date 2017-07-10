#include <fnmatch.h>
#include <dirent.h>
#include <sys/types.h>
#include <cerrno>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <typeinfo>
#include <list>

#ifdef WITH_NTL
#include "NTL/ZZ.h"
#include "NTL/LLL.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#endif

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/LatticeAnalysis.h"

#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/ParamReaderLat.h"

using namespace std;
using namespace NTL;
//using namespace NTL;


namespace LatticeTester
{

//===========================================================================

LatticeAnalysis::LatticeAnalysis (Reducer & reducer, NormType norm, 
      NormaType normaType, int alpha) 
{
   m_reducer = &reducer;
   m_dim = m_reducer->getIntLatticeBasis().getDim(); 
   m_norm = norm;
   m_normalizerType = normaType;
   initNormalizer(normaType, alpha);

}

//===========================================================================

LatticeAnalysis::~LatticeAnalysis ()
{
   delete m_normalizer;
}

//===========================================================================

void LatticeAnalysis::initNormalizer (NormaType norma, int alpha)
{
   // PW_TODO
   // ok si la matrice est directement construite comme m-dual mais probleme si
   // travail direct sur primale re-scaled ?
   // version avec logDensity en parametre (pou m^k sans besoin de calcul det)

   RScal logDensity;
   logDensity = - log( determinant(m_reducer->getIntLatticeBasis().getBasis()) );

   switch (norma) {
      case BESTLAT:
         m_normalizer = new NormaBestLat (logDensity, m_dim);
         break;
      case LAMINATED:
         m_normalizer = new NormaLaminated (logDensity, m_dim);
         break;
      case ROGERS:
         m_normalizer = new NormaRogers (logDensity, m_dim);
         break;
      case MINKL1:
         m_normalizer = new NormaMinkL1 (logDensity, m_dim);
         break;
      case MINKOWSKI:
         m_normalizer = new NormaMinkowski (logDensity, m_dim);
         break;
      case NORMA_GENERIC:
         m_normalizer = new Normalizer (logDensity, m_dim, "Norma_generic");
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

//===========================================================================

bool LatticeAnalysis::performTest (double fact, long blockSize)
{
   bool result;

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

//==========================================================================

/*
 * Reads the test parameters in infile; then do the test.
 * infile is the data file name without extension: if the data file is named
 * "poil.dat", then infile is "poil".
 * Data files must always have the extension "dat".
 */
int LatticeAnalysis::doTestFromInputFile (const char *infile)
{
   // parameters reading
   string fname (infile);
   fname += ".dat";
   ParamReaderLat paramRdr (fname.c_str ());
   fname.clear ();

   LatConfig config;
   paramRdr.read (config);
   //config.write();

   Writer* rw = createWriter (infile, config.outputType);



   // creating the Reducer object
   IntLatticeBasis basis (config.matrix, config.dim, config.norm);
   Reducer red (basis);

   // finding shortest non-zero vector in the lattice
   // PW_TODO : disjonction de cas selon la pre-reduction choisie
   double fact = 1.0 - config.epsilon;
   red.shortestVectorWithBKZ(config.norm, fact, config.blockSize);

   // calculating the density
   // PW_TODO
   // ok si la matrice est directement construite comme m-dual mais probleme si 
   // travail direct sur primale re-scaled ?
   // version avec logDensity en parametre (pou m^k sans besoin de calcul det)
   RScal logDensity;
   logDensity = - log( determinant(red.getIntLatticeBasis().getBasis()) );


   // initializing the normalizer

   Normalizer* normalizer;

   switch (config.normalizer) {
      case BESTLAT:
         m_normalizer = new NormaBestLat (logDensity, config.dim);
         break;
      case LAMINATED:
         m_normalizer = new NormaLaminated (logDensity, config.dim);
         break;
      case ROGERS:
         m_normalizer = new NormaRogers (logDensity, config.dim);
         break;
      case MINKL1:
         m_normalizer = new NormaMinkL1 (logDensity, config.dim);
         break;
      case MINKOWSKI:
         m_normalizer = new NormaMinkowski (logDensity, config.dim);
         break;
      case NORMA_GENERIC:
         m_normalizer = new Normalizer (logDensity, config.dim, "Norma_generic");
         break;

      /*
      case PALPHA_N:
         m_normalizer = new NormaPalpha (m_reducer->getIntLatticeBasis().getModulo(), alpha, dim);
         //PW_TODO : c'est bien ça ?
         break;
      */

      default:
         cout << "normalizer:   no such case";
         exit (2);
   }


   // Calculating the Figure of Merit using the selected normalization
   m_merit = conv<double>(red.getMinLength()) / normalizer->getPreComputedBound(dim);

   // deleting the normalizer pointer
   delete normalizer;




   ReportHeaderLat header (rw, &config, lattice);
   ReportFooterLat footer (rw);
   ReportLat report (rw, &config, &header, &footer);

   double minVal[1 + toDim];
   SetZero (minVal, toDim);

   Normalizer *normal = 0;

   if (config.criter == SPECTRAL) {
      normal = lattice->getNormalizer (config.norma, 0);
      // creates and returns the normalizer corresponding to config.norma
      normal->setNorm (config.norm);
   } else if (config.criter == PALPHA &&
              (config.calcPalpha == NORMPAL || config.calcPalpha == BAL)) {
      normal = new NormaPalpha (lattice->getModulo(), config.alpha, toDim);
   }

   if (!memLacF && config.lacunary) {
      plac = new Lacunary (config.Lac, toDim);
      lattice->setLac (*plac);
   }


   switch (config.criter) {
   case SPECTRAL: {
         LatTestSpectral spectralTest (normal, lattice);
         lattice->buildBasis (fromDim - 1);

         spectralTest.attach (&report);
         report.printHeader ();
         spectralTest.setDualFlag (config.dualF);
         spectralTest.setInvertFlag (config.invertF);
         spectralTest.setDetailFlag (config.detailF);
         spectralTest.setMaxAllDimFlag (true);
         spectralTest.setMaxNodesBB (config.maxNodesBB);

         if (1 == config.d) {
            spectralTest.test (fromDim, toDim, minVal);
            // lattice->write();
            footer.setLatticeTest (&spectralTest);
            report.printTable ();
            report.printFooter ();

         } else {
            if (config.genType[0] == MRG || config.genType[0] == LCG)
               master = new MRGLattice (*(MRGLattice *) lattice);
            else if (config.genType[0] == KOROBOV)
               master = new KorobovLattice (*(KorobovLattice *) lattice);
            else if (config.genType[0] == RANK1)
               master = new Rank1Lattice (*(Rank1Lattice *) lattice);

            master->buildBasis (toDim);
            TestProjections proj (master, lattice, &spectralTest, config.td,
                                  config.d);
            proj. setOutput (rw);
            // proj.setDualFlag (config.dualF);
            proj.setPrintF (true);
            double merit = proj.run (stationary, false, minVal);
            int nbProj = proj.getNumProjections ();
            rw->writeString ("\nMin merit:   ");
            rw->writeDouble (sqrt (merit));
            rw->newLine ();
            rw->writeString ("Num projections:   ");
            rw->writeInt (nbProj);
            rw->newLine ();
            // nbProj = proj.calcNumProjections(stationary, false);
            // cout << "Num projections2:  " << nbProj << endl << endl;
            delete master;
         }
      }
      break;

   case BEYER: {
         LatTestBeyer beyerTest (lattice);
         lattice->buildBasis (fromDim - 1);
         beyerTest.attach (&report);
         report.printHeader ();
         beyerTest.setDualFlag (config.dualF);
         beyerTest.setMaxAllDimFlag (true);
         beyerTest.setMaxNodesBB (config.maxNodesBB);
         beyerTest.setDetailFlag (config.detailF);
         beyerTest.test (fromDim, toDim, minVal);
         footer.setLatticeTest (&beyerTest);
         report.printTable ();
         report.printFooter ();
         //rw->writeString (lattice->toStringDualBasis ());
      }
      break;

   case PALPHA: {
         LatTestPalpha palphaTest (normal, lattice);
         palphaTest.setConfig (&config);
         palphaTest.attach (&report);
         report.printHeader ();
         if (1 == config.d) {
            palphaTest.test (fromDim, toDim, minVal);
            footer.setLatticeTest (&palphaTest);
            report.printTable ();
            report.printFooter ();
         } else {
            MRGLattice master = MRGLattice (*(MRGLattice *) lattice);
            master.buildBasis (toDim);
            TestProjections proj (&master, lattice, &palphaTest, config.td,
                                  config.d);
            proj. setOutput (rw);
            double merit = proj.run (true, false, minVal);
            int nbProj = proj.getNumProjections ();
            rw->writeString ("\n\nMin merit:   ");
            rw->writeDouble (sqrt (merit));
            rw->newLine ();
            rw->writeString ("Num projections:   ");
            rw->writeInt (nbProj);
            rw->newLine ();
            rw->newLine ();
         }
      }
      break;

   default:
      cerr << "Default case for config.criter" << endl;
      return -1;
   }

   if (normal != 0)
      delete normal;
   if (!memLacF && config.lacunary)
      delete plac;
   delete lattice;
   delete rw;
   return 0;
}


//==========================================================================

int LatticeAnalysis::doTestFromDirectory (const char *dirname)
{
   string dir = string (dirname);
   std::vector <string> files = std::vector <string> ();

   getDir (dir, files);
   printFileNames (files);
   eraseExtension (files);

   int flag = 0;
   for (unsigned int i = 0; i < files.size (); i++)
      flag |= doTest (files[i].c_str());

   return flag;
}

//==========================================================================

Writer* LatticeAnalysis::createWriter (const char *infile, OutputType ot)
{
   Writer *rw = 0;
   string fname;

   switch (ot) {
   case RES:
      fname = infile;
      fname += ".res";
      rw = new WriterRes (fname.c_str ());
      break;

   case TEX:
      fname = infile;
      fname += ".tex";
      // rw = new WriterTex(fname.c_str()); //EB Ne permet pas d'écrire en Tex
      break;

   case TERMINAL:
      rw = new WriterRes (&cout);
      break;

   default:
      cerr << "\n*** outputType:   no such case" << endl;
      return 0;
   }
   return rw;
}

//==========================================================================

} // end namespace LatticeTester




namespace
{

//==========================================================================

int getDir (string dir, std::vector <string> & files)
{
   DIR *dp;
   struct dirent *dirp;
   if ((dp = opendir (dir.c_str())) == NULL) {
      cerr << "Directory: " << dir << endl;
      perror ("Couldn't open the directory");
      return errno;
   }

   // Does directory name ends with /
   size_t j = dir.rfind('/');
   string SEP("");
   // if not, add one /
   if (dir.size() != (1 + j))
      SEP += "/";

   while ((dirp = readdir (dp)) != NULL) {
      if (0 == fnmatch("*.dat", dirp->d_name, 0))
         // keeps full name including directory name
         files.push_back (string (dir + SEP + dirp->d_name));
   }
   closedir (dp);
   return 0;
}

//==========================================================================

void eraseExtension (std::vector <string> & files)
{
   for (unsigned int i = 0; i < files.size (); i++) {
      size_t j = files[i].rfind(".dat");
      if (j != string::npos)
         files[i].erase(j);
   }
}

//==========================================================================

void printFileNames (std::vector <string> & files)
{
   cout << "----------------------------------------------" << endl;
   for (unsigned int i = 0; i < files.size (); i++) {
      cout << files[i] << endl;
   }
}

//==========================================================================

} // end namespace
