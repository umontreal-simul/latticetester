/**
ParamReader.cc for ISO C++
version
modified 28/05/06 10:36
authors: Hicham Wahbi
       Frederik Rozon
       Richard Simard
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <cstdlib>

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
<<<<<<< HEAD
=======
#include "latticetester/LatticeTesterConfig.h"
>>>>>>> new_class_LatticeTest

using namespace std;
using namespace LatticeTester;


namespace LatticeTester
{

//===========================================================================

ParamReader::ParamReader()
{
   m_fileName.reserve(MAX_WORD_SIZE);
}


//===========================================================================

ParamReader::ParamReader(string fileName)
{
   m_fileName.reserve(fileName.length());
   m_fileName = fileName;
}


//===========================================================================

ParamReader::~ParamReader()
{
   for (int i = 0; i < (int) m_lines.size(); i++)
      m_lines[i].clear();
   m_lines.clear();
}


//===========================================================================

void ParamReader::getLines()
{
   ifstream inFile(m_fileName.c_str());
   if (inFile.fail()) {
      cerr << "An error occurred. Unable to read input file:" <<
               m_fileName << endl;
      exit(1);
   }
   m_lines.clear();
   const string blancs = " \t";
   string line;
   string::size_type i;

   while (getline(inFile, line, '\n')) {
      i = line.find_first_not_of(blancs);  // find first non blank char
      if (i != string::npos && line[i] != '#')  // comment line begins with a #
         m_lines.push_back(line);
   }

   inFile.close();
}

//===========================================================================

void ParamReader::getToken(string& field, unsigned int ln, unsigned int pos)
{
   std::vector<string> tokens;
   tokens.reserve(20);
   tokens.clear();

   tokenize(tokens, ln - 1);
   if (pos > tokens.size()) {
      cerr << "Warning: position " << pos << " exceeds number of params at line "
      << ln << endl;
      field.clear();
   } else
      field = tokens[pos];
   return;
}


//===========================================================================

int ParamReader::tokenize(std::vector<string>& tokens, unsigned int ln)
{
   if (ln >= m_lines.size()) {
      cerr << "Warning: line " << ln << " doesn't exist in param file\n";
      exit(0);
   }
   string str = m_lines[ln];
   string word;
   int wnum = 0;
   bool firstFound = false;
   char c = 0;
   for (unsigned int i = 0; i < str.length(); i++) {
      c = str[i];
      if (!isDelim(c)) {
         string tmp(1, c);
         word.append(tmp);
         //tokens[wnum] += c;
         firstFound = true;

      } else {
         if (firstFound) {
            tokens.push_back(word);
            word.clear();
            wnum++;
         }
         firstFound = false;
      }
   }

   if (firstFound) {
      tokens.push_back(word);
      word.clear();
      wnum++;
   }
   if (isDelim(c))
      wnum--;
   return wnum;
}


//===========================================================================

bool ParamReader::isDelim(char c)
{
   const string delim = " ?!,\t\n";
   bool bRetVal = 0;
   for (int i = 0; delim[i] != 0; ++i) {
      if (delim[i] == c) {
         bRetVal = 1;
         break;
      }
   }

   return bRetVal;
}


//===========================================================================

void ParamReader::readBool(bool & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (strcasecmp(val.c_str(), "true") == 0)
      field = true;
   else if (strcasecmp(val.c_str(), "false") == 0)
      field = false;
   else
      MyExit(1, "readBool:   NO SUCH CASE");
}


//===========================================================================

void ParamReader::readString(string & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = string(val);
}


//===========================================================================

void ParamReader::readChar(char & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = val[0];
}


//===========================================================================

void ParamReader::readNumber3 (MScal & m, long & m1, long & m2, long & m3,
                               unsigned int ln, unsigned int pos)
{
   m2 = m3 = 0;
   string val;
   getToken(val, ln, pos++);
   istringstream in (val);
   in >> m;
   getToken(val, ln, pos++);
   if (val.empty())
      return;
   istringstream in2 (val);
   conv (m1, m);
   in2 >> m2;
   if (! in2)
      return;
   assert (m2 > 0);
   getToken(val, ln, pos++);
   if (val.empty())
      return;
   istringstream in3 (val);
   in3 >> m3;
   if (! in3)
      return;
   int sign;
   if (m < 0) {
      sign = -1;
      m = -m;
   } else
      sign = 1;
   m = sign*(power (m, m2) + m3);
}


//===========================================================================

void ParamReader::readInt(int& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = (int) strtol(val.c_str(), (char **)NULL, 10);
//   field = atoi(val.c_str());
}

//===========================================================================

void ParamReader::readLong(long& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = strtol(val.c_str(), (char **)NULL, 10);
//   field = atol(val.c_str());
}

//===========================================================================

void ParamReader::readZZ(NTL::ZZ & field, unsigned int ln, int pos)
{
   string val;
   getToken(val, ln, pos);
   field = to_ZZ(val.c_str());
}

//===========================================================================

void ParamReader::readDouble(double& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = strtod(val.c_str(), (char **)NULL);
//   field = atof(val.c_str());
}

//===========================================================================

<<<<<<< HEAD
void ParamReader::readGenType(GenType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (strcasecmp(val.c_str(), "LCG") == 0)
      field = LCG;
   else if (strcasecmp(val.c_str(), "MRG") == 0)
      field = MRG;
   else if (strcasecmp(val.c_str(), "MWC") == 0)
      field = MWC;
   else if (strcasecmp(val.c_str(), "KOROBOV") == 0)
      field = KOROBOV;
   else if (strcasecmp(val.c_str(), "RANK1") == 0)
      field = RANK1;
   else
      MyExit(1, "readGenType:   NO SUCH CASE");
}


//===========================================================================

=======
>>>>>>> new_class_LatticeTest
void ParamReader::readMScal(MScal & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   conv(field, val.c_str());
}

//===========================================================================

void ParamReader::readBScal(BScal& field, unsigned int ln, int pos)
{
   string val;
   getToken(val, ln, pos);
   conv(field, val.c_str());
}

//===========================================================================

void ParamReader::readMVect(MVect & fields, unsigned int & ln, unsigned int pos,
   unsigned int numPos, int j)
{
   for (unsigned int i = pos; i < numPos; i++) {
      readMScal(fields[j], ln, i);
      j++;
   }
}

//===========================================================================

<<<<<<< HEAD
=======
void ParamReader::readBMat(BMat & fields, unsigned int & ln, unsigned int pos,
   unsigned int numPos)
{
   for (unsigned int j = pos; j < numPos; j++){
      for (unsigned int i = pos; i < numPos; i++) {
         readBScal(fields[j][i], ln, i);
      }
      ++ln;
   }
}

//===========================================================================

>>>>>>> new_class_LatticeTest
void ParamReader::readIntVect (int* fields, unsigned int ln, unsigned int pos,
                               unsigned int num, int j)
{
   for (unsigned int i = pos; i < num; i++) {
      readInt(fields[j], ln, i);
      j++;
   }
}

//===========================================================================

void ParamReader::readDoubleVect(double* fields, unsigned int ln,
   unsigned int pos, unsigned int numPos, int j)
{
   for (unsigned int i = pos; i <= numPos; i++) {
      readDouble(fields[j], ln, i);
      j++;
   }
}


//===========================================================================

void ParamReader::readInterval (MVect & B, MVect & C, unsigned int & ln, int k)
{
   long m1, m2, m3;
   for (int i = 1; i <= k; i++) {
      readNumber3 (B[i], m1, m2, m3, ++ln, 1);
      readNumber3 (C[i], m1, m2, m3, ++ln, 1);
      assert (C[i] >= B[i]);
   }
}

<<<<<<< HEAD

//===========================================================================

void ParamReader::readCriterionType(CriterionType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "SPECTRAL"))
      field = SPECTRAL;
   else if (0 == strcasecmp(val.c_str(), "BEYER"))
      field = BEYER;
   else if (0 == strcasecmp(val.c_str(), "PALPHA"))
      field = PALPHA;
   else
      MyExit(1, "readCriterionType:   NO SUCH CASE");
}


=======
>>>>>>> new_class_LatticeTest
//===========================================================================

void ParamReader::readNormType (NormType & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "SUPNORM"))
      field = SUPNORM;
   else if (0 == strcasecmp(val.c_str(), "L1NORM"))
      field = L1NORM;
   else if (0 == strcasecmp(val.c_str(), "L2NORM"))
      field = L2NORM;
   else if (0 == strcasecmp(val.c_str(), "ZAREMBANORM"))
      field = ZAREMBANORM;
    else
      MyExit(1, "readNormType:   NO SUCH CASE");
}


//===========================================================================

<<<<<<< HEAD
void ParamReader::readCalcType (CalcType & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "PAL"))
      field = PAL;
   else if (0 == strcasecmp(val.c_str(), "NORMPAL"))
      field = NORMPAL;
   else if (0 == strcasecmp(val.c_str(), "BAL"))
      field = BAL;
   else if (0 == strcasecmp(val.c_str(), "SEEKPAL"))
      field = SEEKPAL;
   else
      MyExit(1, "readCalcType:   NO SUCH CASE");
}


//===========================================================================

void ParamReader::readDecompType (DecompType & field, unsigned int line,
        unsigned int pos)
{
   string val;
   getToken(val, line, pos);

   if (0 == strcasecmp(val.c_str(), "decomp"))
      field = DECOMP;
   else if (0 == strcasecmp(val.c_str(), "write"))
      field = DECOMP_WRITE;
   else if (0 == strcasecmp(val.c_str(), "read"))
      field = DECOMP_READ;
   else if (0 == strcasecmp(val.c_str(), "prime"))
      field = DECOMP_PRIME;
   else
      MyExit(1, "readDecompType:   NO SUCH CASE");
}


//===========================================================================

=======
>>>>>>> new_class_LatticeTest
void ParamReader::readNormaType(NormaType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "BESTLAT"))
      field = BESTLAT;
   else if (0 == strcasecmp(val.c_str(), "LAMINATED"))
      field = LAMINATED;
   else if (0 == strcasecmp(val.c_str(), "ROGERS"))
      field = ROGERS;
   else if (0 == strcasecmp(val.c_str(), "MINKL1"))
      field = MINKL1;
   else if (0 == strcasecmp(val.c_str(), "MINKOWSKI"))
      field = MINKOWSKI;
   else if (0 == strcasecmp(val.c_str(), "NORMA_GENERIC"))
      field = NORMA_GENERIC;
   else
      MyExit(1, "readNormaType:   NO SUCH CASE");
}


//===========================================================================

<<<<<<< HEAD
void ParamReader::readLatticeType(LatticeType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "FULL"))
      field = FULL;
   else if (0 == strcasecmp(val.c_str(), "RECURRENT"))
      field = RECURRENT;
   else if (0 == strcasecmp(val.c_str(), "ORBIT"))
      field = ORBIT;
   else if (0 == strcasecmp(val.c_str(), "PRIMEPOWER"))
      field = PRIMEPOWER;
   else
      MyExit(1, "readLatticeType:   NO SUCH CASE");
}


//===========================================================================

void ParamReader::readOrbit (int J, MRGComponent **comp, unsigned int & ln)
{
   for (int j = 0; j < J; j++) {
      unsigned int k = comp[j]->k;
      readMVect (comp[j]->orbitSeed, ln, 1U, k, 1);
      ln++;
   }
}


//===========================================================================

=======
>>>>>>> new_class_LatticeTest
void ParamReader::readOutputType(OutputType & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "TERMINAL"))
      field = TERMINAL;
   else if (0 == strcasecmp(val.c_str(), "RES"))
      field = RES;
   else if (0 == strcasecmp(val.c_str(), "GEN")) {
      field = GEN;
      MyExit(1, "readOutputType:   GEN case not ready");
   } else if (0 == strcasecmp(val.c_str(), "TEX")) {
      field = TEX;
      MyExit(1, "readOutputType:   TEX case not ready");
   } else
      MyExit(1, "readOutputType:   NO SUCH CASE");
}

<<<<<<< HEAD
//===========================================================================

void ParamReader::readImplemCond(ImplemCond& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "NoCond"))
      field = NO_COND;
   else if (0 == strcasecmp(val.c_str(), "AppFact"))
      field = APP_FACT;
   else if (0 == strcasecmp(val.c_str(), "NonZero")) {
      field = ZERO_COEF;
   } else if (0 == strcasecmp(val.c_str(), "Equal"))
      field = EQUAL_COEF;
   else if (0 == strcasecmp(val.c_str(), "PowerTwo"))
      field = POWER_TWO;
   else
      MyExit(1, "readImplemCond:   NO SUCH CASE");
}
=======
>>>>>>> new_class_LatticeTest

//===========================================================================

void ParamReader::readPreRed(PreReductionType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "BKZ"))
      field = BKZ;
   else if (0 == strcasecmp(val.c_str(), "PreRedDieter"))
<<<<<<< HEAD
      field = APP_FACT;
   else if (0 == strcasecmp(val.c_str(), "LLL"))
      field = ZERO_COEF;
=======
      field = PreRedDieter;
   else if (0 == strcasecmp(val.c_str(), "LLL"))
      field = LLL;
>>>>>>> new_class_LatticeTest
   else
      MyExit(1, "readPreRed:   NO SUCH CASE");
}

<<<<<<< HEAD
//===========================================================================

void ParamReader::readSearchMethod(SearchMethod& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "Exhaust"))
      field = EXHAUST;
   else if (0 == strcasecmp(val.c_str(), "Random"))
      field = RANDOM;
   else
      MyExit(1, "readSearchMethod:   NO SUCH CASE");
}

//===========================================================================

void ParamReader::readLacunary(int ordre, int fromDim, int toDim,
   unsigned int & ln, bool & lacunary, int & lacGroupSize, NTL::ZZ & lacSpacing,
   BVect & Lac, GenType genType)
{
   readInt (lacGroupSize, ++ln, 0);
   const int t = lacGroupSize;
   if (t > 0)
      readZZ(lacSpacing, ln, 1);

   if (((t == 1) && (lacSpacing == 1)) || (t > toDim)) {
      lacunary = false;
      if (RANK1 == genType)
         return;
      if (toDim <= ordre)
         MyExit(2, "ParamReader::ReadLacunary:   toDim <= k");
      return;
   }

   lacunary = true;
   CreateVect (Lac, toDim);
   int i;
   if (t < 0) {
      for (i = 0; i < toDim; i++)
         readBScal (Lac[i], ++ln, 0);
      return;
   }

   NTL::ZZ Q1, Q;
   if (t > 0)
      Q1 = lacSpacing;

   Q = 0;
   i = 1;
   while (true) {
      for (int j = 0; j < t; j++) {
         if (i < toDim) {
            conv (Lac[i], Q + j);
            i++;
         } else
            return;
      }
      Q += Q1;
   }
}


//===========================================================================

bool ParamReader::checkBound (const MScal & m, const MVect & A, int k)
{
   for (int i = 0; i < k; i++) {
      assert (A[i] < m);
      assert (A[i] > -m);
   }
   return true;
}


//===========================================================================

bool ParamReader::checkPrimePower(LatticeType lat, long m2, long m3, int k)
{
   if (lat != PRIMEPOWER)
      return true;
   if (m2 <= 0)
      MyExit(1, "LatticeType = PRIMEPOWER:   e <= 0");
   if (m3 != 0)
      MyExit(1, "LatticeType = PRIMEPOWER:   c != 0");
   if (k > 1)
      MyExit(1, "LatticeType = PRIMEPOWER:   k > 1");
   return true;
=======


//===========================================================================

bool ParamReader::checkBound (const MScal & m, const MVect & A, int k)
{
   for (int i = 0; i < k; i++) {
      assert (A[i] < m);
      assert (A[i] > -m);
   }
   return true;
}


void ParamReader::read (LatticeTesterConfig & config)
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

>>>>>>> new_class_LatticeTest
}


//===========================================================================

} // End namespace LatMRG
