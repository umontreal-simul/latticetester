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
#include "latticetester/ParamReader.h"

using namespace std;
using namespace NTL;
//using namespace NTL;



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



namespace LatticeTester
{

//===========================================================================

LatticeAnalysis::LatticeAnalysis ()
{
}

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
   logDensity = - log( abs( determinant(m_reducer->getIntLatticeBasis().getBasis()) ) );

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
         m_normalizer = new NormaPalpha (m_reducer->getIntLatticeBasis().getModulo(), alpha, m_dim);
         //PW_TODO : c'est bien ça ?
         break;
      default:
         cout << "LatticeAnalysis::initNormalizer:   no such case";
         exit (2);
   }
}

//===========================================================================

/*

bool LatticeAnalysis::performTest (CriterionType criterion, PreReductionType preRed, 
      NormType norm, int dim, double fact, PrecisionType precision, int blocksize)
{

   //PW_TODO utiliser les member parameters

   bool result;

   // performing pre-reduction
   switch (preRed) {
      case BKZ:
         m_reducer->redBKZ(fact, blocksize, precision);
         break;
      case LenstraLL:
         m_reducer->redLLLNTL(fact, precision);
         break;
      case PreRedDieter:
         m_reducer->preRedDieter(0);
         break;
      default:
         MyExit(1, "LatticeLatticeAnalysis::doTestFromInputFile:   no such case");
         exit(1);
   }

   // performing the Branch-and-Bound procedure to find the shortest non-zero vector
   switch (criterion) {
      case SPECTRAL:
         // performing the Branch-and-Bound procedure to find the shortest non-zero vector
         result = m_reducer->shortestVector(norm);

         // calculating the Figure of Merit
         m_merit = conv<double>(m_reducer->getMinLength())
                  / m_normalizer->getPreComputedBound(dim);
         break;

      case BEYER:
         //performing the Branch-and-Bound procedure to find the Minkowski-reduced matrix
         result = m_reducer->reductMinkowski(0);

         // calculating the Figure of Merit
         m_merit = conv<double>(m_reducer->getMinLength())
                  /conv<double>(m_reducer->getMaxLength());
         //cout << "plus court vecteur : " << conv<double>(m_reducer->getMinLength()) << endl;
         //cout << "plus long vecteur : " << conv<double>(m_reducer->getMaxLength()) << endl;
         break;

      case PALPHA:
         MyExit(1, "PALPHA:   NOT YET");
         break;

      case BOUND_JS:
         MyExit(1, "BOUND_JS:   NOT YET");
         break;

      default:
         MyExit(1, "BOUND_JS:   NOT YET");
         exit(1);
   }

   return result;

}

//==========================================================================

void LatticeAnalysis::printTestResults () 
{
   cout << "\n----------------------------------------------------------" << endl;
   cout << "Prereduction : " << toStringPreRed(config.prereduction) << endl;
   cout << "Length of shortest non-zero vector = " << conv<double>(m_reducer->getMinLength());
   cout << " (" << toStringNorm(config.norm) << ")" << endl;
   cout << "Criterion : " << toStringCriterion(config.test) << endl;
   cout << "Figure of Merit = " << m_merit;
   cout << " (" << toStringNorma(config.normalizer) << " normalization)" << endl;
   cout << "----------------------------------------------------------" << endl;

}

*/

//==========================================================================

/*
 * Reads the test parameters in infile; then do the test.
 * infile is the data file name without extension: if the data file is named
 * "poil.dat", then infile is "poil".
 * Data files must always have the extension "dat".
 */
int LatticeAnalysis::doTestFromInputFile (const char *infile)
{

   // VERSION DEGUEUE AVEC WRITER
   // PW_TODO à supprimer ou non plus tard
   /*
   // parameters reading
   string fname (infile);
   fname += ".dat";
   ParamReader paramRdr (fname.c_str ());
   fname.clear ();

   LatConfig config;
   paramRdr.read (config);
   //config.write();

   Writer* rw = createWriter (infile, config.outputType);

   // creating the Reducer object from input
   IntLatticeBasis basis (config.matrix, config.dim, config.norm);
   Reducer red (basis);

   // update parameters
   setReducer(red);
   setDim(config.dim);
   setNorm(config.norm);
   setNormalizerType(config.normalizer);
   initNormalizer (config.normalizer);

   double fact = 1.0 - config.epsilon;

   // performing pre-reduction
   switch (config.prered) {
      case BKZ:
         m_reducer->redBKZ(fact, config.blockSize);
         break;
      case LenstraLL:
         m_reducer->redLLLNTLProxyFP(fact);
         break;
      case PreRedDieter:
         m_reducer->preRedDieter(0);
         break;
      default:
         cout << "LatticeLatticeAnalysis::doTestFromInputFile:   no such case";
         exit (2);
   }

   // performing the Branch-and-Bound procedure to find the shortest non-zero vector
   m_reducer->shortestVector(config.norm);

   // calculating the Figure of Merit
   m_merit = conv<double>(m_reducer->getMinLength())
            / m_normalizer->getPreComputedBound(config.dim);





   ReportHeader header (rw, &config, lattice);
   ReportFooter footer (rw);
   Report report (rw, &config, &header, &footer);

   double minVal[1 + toDim];
   SetZero (minVal, toDim);


   report.printHeader ();
   //footer.setLatticeTest (&spectralTest);
   report.printTable ();
   report.printFooter ();

   delete rw;
   return 0;
   */



   // VERSION SANS WRITER

   // parameters reading
   string fname (infile);
   fname += ".dat";
   ParamReader paramRdr (fname.c_str ());
   fname.clear ();

   LatticeTesterConfig config;
   paramRdr.read (config);

   //cout << "Writing config = " << endl;
   //config.write();

   // creating the Reducer object from input
   IntLatticeBasis basis (config.basis, config.dim, config.norm);
   Reducer red (basis);

   // update parameters
   setReducer(red);
   setDim(config.dim);
   setNorm(config.norm);
   setNormalizerType(config.normalizer);
   initNormalizer (config.normalizer);


   // performing pre-reduction
   switch (config.prereduction) {
      case BKZ:
         m_reducer->redBKZ(config.fact, config.blocksize, config.precision);
         break;
      case LenstraLL:
         m_reducer->redLLLNTL(config.fact, config.precision);
         break;
      case PreRedDieter:
         m_reducer->preRedDieter(0);
         break;
      default:
         MyExit(1, "LatticeLatticeAnalysis::doTestFromInputFile:   no such case");
         exit(1);
   }

   // performing the Branch-and-Bound procedure to find the shortest non-zero vector
   switch (config.test) {
      case SPECTRAL:
         // performing the Branch-and-Bound procedure to find the shortest non-zero vector
         m_reducer->shortestVector(config.norm);

         // calculating the Figure of Merit
         m_merit = conv<double>(m_reducer->getMinLength())
                  / m_normalizer->getPreComputedBound(config.dim);
         break;

      case BEYER:
         //performing the Branch-and-Bound procedure to find the Minkowski-reduced matrix
         m_reducer->reductMinkowski(0);

         // calculating the Figure of Merit
         m_merit = conv<double>(m_reducer->getMinLength())
                  /conv<double>(m_reducer->getMaxLength());
         cout << "plus court vecteur : " << conv<double>(m_reducer->getMinLength()) << endl;
         cout << "plus long vecteur : " << conv<double>(m_reducer->getMaxLength()) << endl;
         break;

      case PALPHA:
         MyExit(1, "PALPHA:   NOT YET");
         break;

      case BOUND_JS:
         MyExit(1, "BOUND_JS:   NOT YET");
         break;

      default:
         MyExit(1, "BOUND_JS:   NOT YET");
         exit(1);
   }

   cout << "\n----------------------------------------------------------" << endl;
   cout << "Prereduction : " << toStringPreRed(config.prereduction) << endl;
   cout << "Length of shortest non-zero vector = " << conv<double>(m_reducer->getMinLength());
   cout << " (" << toStringNorm(config.norm) << ")" << endl;
   cout << "Criterion : " << toStringCriterion(config.test) << endl;
   cout << "Figure of Merit = " << m_merit;
   cout << " (" << toStringNorma(config.normalizer) << " normalization)" << endl;
   cout << "----------------------------------------------------------" << endl;

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
      flag |= doTestFromInputFile (files[i].c_str());

   return flag;
}

//==========================================================================

/*
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
*/

//==========================================================================

} // end namespace LatticeTester

