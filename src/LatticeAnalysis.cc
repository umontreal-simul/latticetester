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

#include "NTL/ZZ.h"
#include "NTL/LLL.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"

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
   m_normalizer = 0;
   m_reducer = 0;
}

//===========================================================================

LatticeAnalysis::LatticeAnalysis (Reducer & reducer, CriterionType criterion, 
   NormaType normaType, PreReductionType preRed, NormType norm, int alpha,
   long maxNodesBB)
{
   m_reducer = &reducer;
   m_criterion = criterion;
   m_normalizerType = normaType;
   initNormalizer(normaType, alpha);
   m_preRed = preRed;
   m_norm = norm;
   m_dim = m_reducer->getIntLatticeBasis().getDim();
   Reducer::maxNodesBB = m_maxNodesBB = maxNodesBB;
}

//===========================================================================

LatticeAnalysis::~LatticeAnalysis ()
{
   if (m_normalizer != 0)
      delete m_normalizer;
}

//===========================================================================

void LatticeAnalysis::initNormalizer (NormaType norma, int alpha)
{
 
#if NTL_TYPES_CODE > 1
   RScal logDensity;
   logDensity = - log( abs( NTL::determinant(m_reducer->getIntLatticeBasis().getBasis()) ) );
#else
   // As NTL library does not support matrix with double
   // we compute the determinant with the boost library
   boost::numeric::ublas::matrix<long>  mat_tmps;
   mat_tmps.resize(m_dim, m_dim);
   for(unsigned int i = 0; i < m_dim; i++){
      for(unsigned int j = 0; j < m_dim; j++){
         mat_tmps(i,j) = m_reducer->getIntLatticeBasis().getBasis()(i,j);
      }
   }
   RScal logDensity(-log( abs( det_double(mat_tmps) ) ) );
#endif

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
         break;
      case L1:
         break;
      case L2:
         break;
      default:
         cout << "LatticeAnalysis::initNormalizer:   no such case";
         exit (2);
   }
}

//===========================================================================

bool LatticeAnalysis::doTest (double fact, PrecisionType precision, int blocksize)
{
   bool result = false;

   // performing pre-reduction
   switch (m_preRed) {
      case BKZ:
         m_reducer->redBKZ(fact, blocksize, precision);
         break;
      case LenstraLL:
         m_reducer->redLLLNTL(fact, precision);
         break;
      case PreRedDieter:
         m_reducer->preRedDieter(0);
         break;
      case NOPRERED:
         break;
      default:
         MyExit(1, "LatticeLatticeAnalysis::doTest:   no such case");
         exit(1);
   }

   // performing the Branch-and-Bound procedure to find the shortest non-zero vector
   switch (m_criterion) {
      case SPECTRAL:
         // performing the Branch-and-Bound procedure to find the shortest non-zero vector
         result = m_reducer->shortestVector(m_norm);
         // calculating the Figure of Merit
         if(m_normalizerType == L1 || m_normalizerType == L2) {
            m_merit = conv<double>(m_reducer->getMinLength());
         } else {
            m_merit = conv<double>(m_reducer->getMinLength())
                  / m_normalizer->getPreComputedBound(m_dim);
         }
         break;
      case BEYER:
         //performing the Branch-and-Bound procedure to find the Minkowski-reduced matrix
         result = m_reducer->reductMinkowski(0);
         // calculating the Figure of Merit
         m_merit = conv<double>(m_reducer->getMinLength())
                  /conv<double>(m_reducer->getMaxLength());
         break;
      case PALPHA:
         MyExit(1, "PALPHA:   to be implemented");
         break;
      case BOUND_JS:
         MyExit(1, "BOUND_JS:   NOT YET");
         break;
      default:
         MyExit(1, "LatticeAnalysis::doTest:   NO SUCH CASE");
         exit(1);
   }

   return result;
}

//==========================================================================

void LatticeAnalysis::printTestResults ()
{
   cout << "\n----------------------------------------------------------" << endl;
   cout << "Criterion: " << toStringCriterion(m_criterion) << endl;
   cout << "Prereduction used: " << toStringPreRed(m_preRed) << endl;
   cout << "Length of shortest non-zero vector = " << conv<double>(m_reducer->getMinLength());
   cout << " (" << toStringNorm(m_norm) << ")" << endl;
   cout << "Figure of Merit = " << m_merit;
   cout << " (" << toStringNorma(m_normalizerType) << " normalization)" << endl;
   cout << "----------------------------------------------------------\n" << endl;
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
   ParamReader paramRdr (fname.c_str ());
   fname.clear ();

   LatticeTesterConfig config;
   paramRdr.read (config);
   //config.write();

   LatTestWriter* rw = createLatTestWriter (infile, config.outputType);

   // creating the Reducer object from input
   IntLatticeBasis basis (config.basis, config.dim, config.norm);
   Reducer red (basis);
   // update parameters
   setReducer(red);
   setCriterion(config.test);
   setNormalizerType(config.normalizer);
   setDim(config.dim);
   initNormalizer (config.normalizer);
   setPreReduction(config.prereduction);
   setNorm(config.norm);
   setMaxNodesBB(config.maxNodesBB);
   
   if (!doTest(config.fact, config.precision, config.blocksize)) {
      MyExit(1, "error in LatticeAnalysis::doTestFromInputFile");
      exit(1);
   }

   // putting the results in the output stream
   rw->writeString("\n----------------------------------------------------------");
   rw->newLine();
   rw->writeString("Criterion: "); 
   rw->writeString(toStringCriterion(m_criterion));
   rw->newLine();
   rw->writeString("Prereduction used: ");
   rw->writeString(toStringPreRed(m_preRed));
   rw->newLine();
   rw->writeString("Length of shortest non-zero vector = ");
   rw->writeDouble(conv<double>(m_reducer->getMinLength()));
   rw->writeString(" (");
   rw->writeString(toStringNorm(m_norm));
   rw->writeString(")");
   rw->newLine();
   rw->writeString("Figure of Merit = ");
   rw->writeDouble(m_merit);
   rw->writeString(" (");
   rw->writeString(toStringNorma(m_normalizerType));
   rw->writeString(" normalization)");
   rw->newLine();
   rw->writeString("----------------------------------------------------------\n");
   rw->newLine();

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
      flag |= doTestFromInputFile (files[i].c_str());

   return flag;
}

//==========================================================================

LatTestWriter* LatticeAnalysis::createLatTestWriter (const char *infile, OutputType ot)
{
   LatTestWriter *rw = 0;
   string fname;

   switch (ot) {
   case RES:
      fname = infile;
      fname += ".res";
      rw = new LatTestWriterRes (fname.c_str ());
      break;

   case TEX:
      fname = infile;
      fname += ".tex";
      //rw = new WriterTex(fname.c_str()); //EB Ne permet pas d'écrire en Tex
      cerr << "\n*** outputType:   TEX not implemented" << endl;
      return 0;
      break;

   case TERMINAL:
      rw = new LatTestWriterRes (&cout);
      break;

   default:
      cerr << "\n*** outputType:   no such case" << endl;
      return 0;
   }
   return rw;
}

//==========================================================================

} // end namespace LatticeTester

