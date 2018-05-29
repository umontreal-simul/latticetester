#ifndef PARAMREADER_H
#define PARAMREADER_H
#include "NTL/ZZ.h"

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/LatticeTesterConfig.h"

#include <string>
#include <vector>
#include <cassert>
#include <fstream>


namespace LatticeTester {

  /**
   * Utility class used to read basic parameter fields in a configuration file.
   * Lines whose first non-blank character is a <tt>#</tt> are considered as
   * comments and discarded.
   *
   */
  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    class ParamReader { 

      public:
        static const int MAX_WORD_SIZE = 64;

        /**
         * Constructor.
         */
        ParamReader();

        /**
         * Constructor. Opens the file `fileName`.
         */
        ParamReader (std::string fileName);

        /**
         * Destructor.
         */
        ~ParamReader();

        /**
         * Reads all the lines from the file and stores them into this object’s
         * buffer. Lines whose first non-blank character is a <tt>#</tt> are
         * considered as comments and discarded. Empty lines are also
         * discarded.
         */
        void getLines();

        /**
         * Puts into `field` the <tt>pos</tt>-th string token from line `ln`.
         */
        void getToken (std::string & field, unsigned int ln, unsigned int pos);

        /**
         * Splits line `ln` from the file into several string tokens. Separator
         * characters are defined in function `IsDelim`. Tokens are stored in
         * vector `tokens`.
         */
        int tokenize (std::vector<std::string> & tokens, unsigned int ln);

        /**
         * Reads a string from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readString (std::string & field, unsigned int ln, unsigned int pos);

        /**
         * Reads a boolean from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readBool (bool & field, unsigned int ln, unsigned int pos);

        /**
         * Reads a character from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readChar (char & field, unsigned int ln, unsigned int pos);

        /**
         * Reads an integer from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readInt (int & field, unsigned int ln, unsigned int pos);

        /**
         * Reads a long from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readLong (long & field, unsigned int ln, unsigned int pos);

        /**
         * Reads a large integer from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readZZ (NTL::ZZ & field, unsigned int ln, int pos);

        /**
         * Reads a double from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readDouble (double & field, unsigned int ln, unsigned int pos);


        /**
         * Reads \f$b\f$, \f$e\f$ and \f$c\f$, starting at the <tt>pos</tt>-th
         * token of the <tt>ln</tt>-th line and uses them to define \f$r\f$.
         * The numbers in the data file may be given in one of the two
         * following formats:
         *
         * \f$\bullet\f$ A single integer giving the value of \f$r=b\f$
         * directly on a line. In that case, one sets \f$e=c=0\f$. <br>
         * \f$\bullet\f$ Three integers \f$b\f$, \f$e\f$, \f$c\f$ on the same
         * line, separated by at least one blank. The \f$r\f$ value will be set
         * as \f$r=b^e+c\f$ if \f$b>0\f$, and \f$r= -(|b|^e+c)\f$ if \f$b<0\f$.
         * One must have \f$e\ge0\f$. For example, \f$(b, e, c) = (2,
         * 5, -1)\f$ will give \f$r=31\f$, while \f$(b, e, c) = (-2, 5, -1)\f$
         * will give \f$r=-31\f$.
         */
        void readNumber3 (Int & r, long & b, long & e, long & c, unsigned int ln, unsigned int pos);

        /**
         * Reads a `BScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readBScal (BasInt & field, unsigned int ln, int pos);

        /**
         * Reads a `BMat` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readBMat (BasIntMat & fields, unsigned int & ln, unsigned int pos,
            unsigned int numPos);

        /**
         * Reads a `Int` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readMScal (Int & field, unsigned int ln, unsigned int pos);

        /**
         * Reads a `RedDbl` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readRedDbl (RedDbl & field, unsigned int ln, unsigned int pos);

        /**
         * Reads `num` tokens (from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
         * `field`.
         */
        void readMVect (IntVec & field, unsigned int & ln, unsigned int pos,
            unsigned int num, int j);

        /**
         * Reads `num` tokens (from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
         * `field`.
         */
        void readIntVect (int* field, unsigned int ln, unsigned int pos,
            unsigned int num, int j);

        /**
         * Reads `num` tokens (from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
         * `field`.
         */
        void readDoubleVect (double* field, unsigned int ln, unsigned int pos,
            unsigned int num, int j);

        /**
         * Reads `2k MScal` tokens into vectors `B` and `C`, starting at the
         * <tt>ln</tt>-th line. These represent a box \f$[B_i, C_i]\f$, \f$i =
         * 1, 2, …, k\f$. The \f$B_i, C_i\f$ must be given in the order \f$B_1,
         * C_1, B_2, C_2, …, B_k, C_k\f$, each on a line of its own. Each
         * coefficient may be given in the form described in `readNumber3`
         * above.
         */
        //void readInterval (MVect & B, MVect & C, unsigned int & ln, int k);

        /**
         * Reads a norm from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readNormType (NormType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a test type from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readCriterionType (CriterionType& field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a type of normalization from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readNormaType (NormaType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a type of precision from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readPrecisionType (PrecisionType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a type of PreReduction from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readPreRed (PreReductionType & field, unsigned int ln,
            unsigned int pos);


        /**
         * Reads an output form from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readOutputType (OutputType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads the configuration file into `config` for the Beyer and the
         * spectral tests.
         */
        void read (LatticeTesterConfig<Int, BasIntMat> & config);
      private:

        /**
         * Internal line buffer.
         */
        std::vector<std::string> m_lines;

        /**
         * The path of the opened file.
         */
        std::string m_fileName;

        /**
         * Checks if the character `c` is to be considered as a token separator
         * or not.
         */
        bool isDelim (char c);

        /**
         * Does nothing for now.
         */
        void init() {}
    }; // End class ParamReader

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::ParamReader()
    {
      m_fileName.reserve(MAX_WORD_SIZE);
    }


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::ParamReader(string fileName)
    {
      m_fileName.reserve(fileName.length());
      m_fileName = fileName;
    }


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::~ParamReader()
    {
      for (int i = 0; i < (int) m_lines.size(); i++)
        m_lines[i].clear();
      m_lines.clear();
    }


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::getLines()
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::getToken(string& field, unsigned int ln, unsigned int pos)
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    int ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::tokenize(std::vector<string>& tokens, unsigned int ln)
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    bool ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::isDelim(char c)
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readBool(bool & field, unsigned int ln, unsigned int pos)
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readString(string & field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = string(val);
    }


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readChar(char & field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = val[0];
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readNumber3 (Int & m, long & m1, long & m2, long & m3,
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

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readInt(int& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = (int) strtol(val.c_str(), (char **)NULL, 10);
      //   field = atoi(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readLong(long& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = strtol(val.c_str(), (char **)NULL, 10);
      //   field = atol(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readZZ(NTL::ZZ & field, unsigned int ln, int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = to_ZZ(val.c_str());
    }


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readDouble(double& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      field = strtod(val.c_str(), (char **)NULL);
      //   field = atof(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readMScal(Int & field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      conv(field, val.c_str());
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readBScal(BasInt& field, unsigned int ln, int pos)
    {
      string val;
      getToken(val, ln, pos);
      conv(field, val.c_str());
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readMVect(IntVec & fields, unsigned int & ln, unsigned int pos,
        unsigned int numPos, int j)
    {
      for (unsigned int i = pos; i < numPos; i++) {
        readMScal(fields[j], ln, i);
        j++;
      }
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readBMat(BasIntMat & fields, unsigned int & ln, unsigned int pos,
        unsigned int numPos)
    {
      for (unsigned int i = pos; i < numPos; i++){
        for (unsigned int j = pos; j < numPos; j++){
          readBScal(fields(i,j), ln, j);
        }
        ln++;
      }

    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readIntVect (int* fields, unsigned int ln, unsigned int pos,
        unsigned int num, int j)
    {
      for (unsigned int i = pos; i < num; i++) {
        readInt(fields[j], ln, i);
        j++;
      }
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readDoubleVect(double* fields, unsigned int ln,
        unsigned int pos, unsigned int numPos, int j)
    {
      for (unsigned int i = pos; i <= numPos; i++) {
        readDouble(fields[j], ln, i);
        j++;
      }
    }


  //===========================================================================

  /*
     void ParamReader::readInterval (MVect & B, MVect & C, unsigned int & ln, int k)
     {
     long m1, m2, m3;
     for (int i = 1; i <= k; i++) {
     readNumber3 (B[i], m1, m2, m3, ++ln, 1);
     readNumber3 (C[i], m1, m2, m3, ++ln, 1);
     assert (C[i] >= B[i]);
     }
     }
     */

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readCriterionType(CriterionType& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);

      if (0 == strcasecmp(val.c_str(), "SPECTRAL"))
        field = SPECTRAL;
      else if (0 == strcasecmp(val.c_str(), "BEYER"))
        field = BEYER;
      else if (0 == strcasecmp(val.c_str(), "PALPHA")){
        field = PALPHA;
      }
      else
        MyExit(1, "readCriterionType:   NO SUCH CASE");
    }

  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readNormType (NormType & field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);
      if (0 == strcasecmp(val.c_str(), "SUPNORM")){
        field = SUPNORM;
        MyExit(1, "readNormType:   SUPNORM case not ready");
      }
      else if (0 == strcasecmp(val.c_str(), "L1NORM"))
        field = L1NORM;
      else if (0 == strcasecmp(val.c_str(), "L2NORM"))
        field = L2NORM;
      else if (0 == strcasecmp(val.c_str(), "ZAREMBANORM")){
        field = ZAREMBANORM;
        MyExit(1, "readNormType:   ZAREMBANORM case not ready");
      }
      else
        MyExit(1, "readNormType:   NO SUCH CASE");
    }


  //===========================================================================


  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readNormaType(NormaType& field, unsigned int ln, unsigned int pos)
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
      else if (0 == strcasecmp(val.c_str(), "L1"))
        field = L1;
      else if (0 == strcasecmp(val.c_str(), "L2"))
        field = L2;
      else
        MyExit(1, "readNormaType:   NO SUCH CASE");
    }

  //===========================================================================


  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readPrecisionType(PrecisionType& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);

      if (0 == strcasecmp(val.c_str(), "DOUBLE"))
        field = DOUBLE;
      else if (0 == strcasecmp(val.c_str(), "QUADRUPLE"))
        field = QUADRUPLE;
      else if (0 == strcasecmp(val.c_str(), "EXPONENT"))
        field = EXPONENT;
      else if (0 == strcasecmp(val.c_str(), "ARBITRARY"))
        field = ARBITRARY;
      else if (0 == strcasecmp(val.c_str(), "EXACT"))
        field = EXACT;
      else
        MyExit(1, "readPrecisionType:   NO SUCH CASE");
    }



  //===========================================================================


  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readOutputType(OutputType & field, unsigned int ln, unsigned int pos)
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


  //===========================================================================

  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::readPreRed(PreReductionType& field, unsigned int ln, unsigned int pos)
    {
      string val;
      getToken(val, ln, pos);

      if (0 == strcasecmp(val.c_str(), "BKZ"))
        field = BKZ;
      else if (0 == strcasecmp(val.c_str(), "PreRedDieter"))
        field = PreRedDieter;
      else if (0 == strcasecmp(val.c_str(), "LenstraLL"))
        field = LenstraLL;
      else if (0 == strcasecmp(val.c_str(), "NOPRERED"))
        field = NOPRERED;
      else
        MyExit(1, "readPreRed:   NO SUCH CASE");
    }

  //===========================================================================

  template<typename Int, typename BasIntMat>
    void initNorm(LatticeTesterConfig<Int, BasIntMat> & config)
    {
      switch(config.normalizer){
        case BESTLAT: case LAMINATED: case ROGERS: case NORMA_GENERIC: case MINKOWSKI: case L2:
          config.norm = L2NORM;
          break;
        case MINKL1: case L1:
          config.norm = L1NORM;
          break;
        default:
          MyExit(1, "initNorm NORMA:   NO SUCH CASE");
          break;
      }
    }

  //===========================================================================


  template<typename Int, typename IntVec, typename BasInt, typename BasIntMat, typename RedDbl>
    void ParamReader<Int, IntVec, BasInt, BasIntMat, RedDbl>::read (LatticeTesterConfig<Int, BasIntMat> & config)
    {
      getLines ();
      unsigned int ln = 1;

      readCriterionType (config.test, ln, 0);

      switch(config.test) {

        case SPECTRAL:
          readNormaType (config.normalizer, ln, 1);
          initNorm(config);
          readPreRed (config.prereduction, ++ln, 0);
          if (config.prereduction == BKZ) {
            readPrecisionType (config.precision, ln, 1);
            readDouble (config.fact, ln, 2);
            readInt (config.blocksize, ln, 3);
          } else if(config.prereduction == LenstraLL) {
            readPrecisionType (config.precision, ln, 1);
            readDouble (config.fact, ln, 2);
          }
          readInt (config.dim, ++ln, 0);
          config.basis.resize(config.dim, config.dim);
          readBMat(config.basis, ++ln, 0, config.dim);
          readLong (config.maxNodesBB, ln, 0);
          readOutputType(config.outputType, ++ln, 0);
          break;

        case BEYER: 
          MyExit(1, "BEYER not implemented yet. Need to change pre-red in reductMinkowski (either preRedDieter with dual basis calculation or redBKZ pre-redduction directly"); 
          break;

        case PALPHA:
          long m1, m2, m3;
          readNumber3(config.modulo, m1, m2, m3, ++ln, 0);
          readInt (config.alpha, ++ln, 0);
          readInt (config.dim, ++ln, 0);
          config.basis.resize(config.dim, config.dim);
          readBMat(config.basis, ++ln, 0, config.dim);
          readLong (config.maxNodesBB, ln, 0);
          readOutputType(config.outputType, ++ln, 0);
          break;

        default:
          MyExit(1, "Criterion type:   NO SUCH CASE");
          break;
      }
    }

} // End namespace LatticeTester
#endif
