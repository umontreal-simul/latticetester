// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PARAMREADER_H
#define PARAMREADER_H

#include "NTL/ZZ.h"

#include "latticetester/ntlwrap.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/LatticeTesterConfig.h"

#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdint>

namespace LatticeTester {

  /**
   * Utility class that can be used to read different kind of data from a file.
   * This class has to be initialized with a `fileName` to be used. After the 
   * object is initialized, it is then necessary to read the entire file with
   * `getLines`. This will make a vector containing all the lines of the file,
   * removing those with a `#` at the start. There then are many methods 
   * allowing the conversion of string tokens to various data types by 
   * specifying a line and a token number. In every line, a token is a string 
   * delimited by any of these caracters: ` ` (blank), `?`, `!`, `,`, `\t`, and
   * `\n`. For example, `ReadString(string, 2, 3)` will put the 4th token 
   * (numerotation starts at 0) of the 3rd line in string.  The user of this 
   * class has to be aware of the format of the file that will be read.
   *
   * \remark This class should have proper error management in the `GetToken` 
   * method because the program curently simply crashes without explanation if 
   * this class is misued.
   */
  template<typename Int, typename BasInt, typename RedDbl>
    class ParamReader { 
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<BasInt> BasIntMat;
      public:
        static const int MAX_WORD_SIZE = 64;

        /**
         * Empty constructor. Allocates space for a word.
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
         * Reads a string from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readString (std::string & field, unsigned int ln,
            unsigned int pos);

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
         * Reads a std::int64_t from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readLong (std::int64_t & field, unsigned int ln, unsigned int pos);

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
        void readNumber3 (Int & r, std::int64_t & b, std::int64_t & e, std::int64_t & c,
            unsigned int ln, unsigned int pos);

        /**
         * Reads a `BScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readBScal (BasInt & field, unsigned int ln, int pos);

        /**
         * Reads a square `BMat` of size `numPos*numPos` from the 
         * <tt>pos</tt>-th token of the <tt>ln</tt>-th line into `field`. The 
         * lines in the matrix will be read from the subsequent lines in the 
         * file, starting from token at the position `pos`. `fields` has to be
         * initialized to the right size.
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
         * Reads the configuration file into `config` for. This is a preset 
         * reading order that is used in `LatticeAnalysis` to get data from a
         * file with the standard format described in `Using the executable`.
         */
        void read (LatticeTesterConfig<Int, BasIntMat> & config);

      protected:

        /**
         * Puts into `field` the <tt>pos</tt>-th string token from line `ln`.
         * This is intended to be used by other methods the get a specific 
         * string token that can then be converted. This method uses `tokenize`
         * to split the line it has to operate on.
         */
        void getToken (std::string & field, unsigned int ln, unsigned int pos);

      private:

        /**
         * Internal line buffer that stores every line not starting with a `#`
         * after `getLines()` has been called.
         */
        std::vector<std::string> m_lines;

        /**
         * The path of the file that this object can read. Note that there are
         * no methods allowing to modify it in this class so it has to be 
         * initialized with the constructor.
         */
        std::string m_fileName;

        /**
         * Checks if the character `c` is to be considered as a token separator
         * or not. This is a method used internally by `tokenize`.
         */
        bool isDelim (char c);

        /**
         * Splits line `ln` in `m_lines` into several string tokens. This
         * function calls `isDelim` on every character to split the line at any
         * of ` ?!,\t\n`. Tokens are stored in vector `tokens`.
         */
        int tokenize (std::vector<std::string> & tokens, unsigned int ln);

        /**
         * Does nothing for now.
         */
        void init() {}
    }; // End class ParamReader

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    ParamReader<Int, BasInt, RedDbl>::ParamReader()
    {
      m_fileName.reserve(MAX_WORD_SIZE);
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    ParamReader<Int, BasInt, RedDbl>::ParamReader(
        std::string fileName)
    {
      m_fileName.reserve(fileName.length());
      m_fileName = fileName;
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    ParamReader<Int, BasInt, RedDbl>::~ParamReader()
    {
      for (int i = 0; i < (int) m_lines.size(); i++)
        m_lines[i].clear();
      m_lines.clear();
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::getLines()
    {
      std::ifstream inFile(m_fileName.c_str());
      if (inFile.fail()) {
        std::cerr << "An error occurred. Unable to read input file:" <<
          m_fileName << std::endl;
        exit(1);
      }
      m_lines.clear();
      const std::string blancs = " \t";
      std::string line;
      std::string::size_type i;

      while (getline(inFile, line, '\n')) {
        i = line.find_first_not_of(blancs);  // find first non blank char
        if (i != std::string::npos && line[i] != '#')  // comment line begins with a #
          m_lines.push_back(line);
      }

      inFile.close();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::getToken(
        std::string& field, unsigned int ln, unsigned int pos)
    {
      std::vector<std::string> tokens;
      tokens.reserve(20);
      tokens.clear();

      tokenize(tokens, ln);
      if (pos > tokens.size()) {
        std::cerr << "Warning: position " << pos << " exceeds number of params at line "
          << ln << std::endl;
        field.clear();
      } else
        field = tokens[pos];
      return;
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    int ParamReader<Int, BasInt, RedDbl>::tokenize(
        std::vector<std::string>& tokens, unsigned int ln)
    {
      if (ln >= m_lines.size()) {
        std::cerr << "Warning: line " << ln << " doesn't exist in param file\n";
        exit(0);
      }
      std::string str = m_lines[ln];
      std::string word;
      int wnum = 0;
      bool firstFound = false;
      char c = 0;
      for (unsigned int i = 0; i < str.length(); i++) {
        c = str[i];
        if (!isDelim(c)) {
          std::string tmp(1, c);
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

  template<typename Int, typename BasInt, typename RedDbl>
    bool ParamReader<Int, BasInt, RedDbl>::isDelim(char c)
    {
      const std::string delim = " ?!,\t\n";
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

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readBool(
        bool & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);

      if (strcasecmp(val.c_str(), "true") == 0)
        field = true;
      else if (strcasecmp(val.c_str(), "false") == 0)
        field = false;
      else
        MyExit(1, "readBool:   NO SUCH CASE");
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readString(
        std::string & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = std::string(val);
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readChar(
        char & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = val[0];
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readNumber3 (
        Int & m, std::int64_t & m1, std::int64_t & m2, std::int64_t & m3, unsigned int ln, 
        unsigned int pos)
    {
      m2 = m3 = 0;
      std::string val;
      getToken(val, ln, pos++);
      std::istringstream in (val);
      in >> m;
      getToken(val, ln, pos++);
      if (val.empty())
        return;
      std::istringstream in2 (val);
      NTL::conv (m1, m);
      in2 >> m2;
      if (! in2)
        return;
      assert (m2 > 0);
      getToken(val, ln, pos++);
      if (val.empty())
        return;
      std::istringstream in3 (val);
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

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readInt(
        int& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = (int) strtol(val.c_str(), (char **)NULL, 10);
      //   field = atoi(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readLong(
        std::int64_t& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = strtol(val.c_str(), (char **)NULL, 10);
      //   field = atol(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readZZ(
        NTL::ZZ & field, unsigned int ln, int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = NTL::to_ZZ(val.c_str());
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readDouble(
        double& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      field = strtod(val.c_str(), (char **)NULL);
      //   field = atof(val.c_str());
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readMScal(
        Int & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      NTL::conv(field, val.c_str());
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readBScal(
        BasInt& field, unsigned int ln, int pos)
    {
      std::string val;
      getToken(val, ln, pos);
      NTL::conv(field, val.c_str());
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readMVect(
        IntVec & fields, unsigned int & ln, unsigned int pos,
        unsigned int numPos, int j)
    {
      for (unsigned int i = pos; i < numPos; i++) {
        readMScal(fields[j], ln, i);
        j++;
      }
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readBMat(
        BasIntMat & fields, unsigned int & ln, unsigned int pos, 
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

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readIntVect (
        int* fields, unsigned int ln, unsigned int pos, unsigned int num, int j)
    {
      for (unsigned int i = pos; i < num; i++) {
        readInt(fields[j], ln, i);
        j++;
      }
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readDoubleVect(
        double* fields, unsigned int ln, unsigned int pos, unsigned int numPos,
        int j)
    {
      for (unsigned int i = pos; i <= numPos; i++) {
        readDouble(fields[j], ln, i);
        j++;
      }
    }


  //===========================================================================

  /*
   * void ParamReader::readInterval (MVect & B, MVect & C, unsigned int & ln, int k)
   * {
   * std::int64_t m1, m2, m3;
   * for (int i = 1; i <= k; i++) {
   * readNumber3 (B[i], m1, m2, m3, ++ln, 1);
   * readNumber3 (C[i], m1, m2, m3, ++ln, 1);
   * assert (C[i] >= B[i]);
   * }
   * }
   * */

  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readCriterionType(
        CriterionType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
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

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readNormType (
        NormType & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
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


  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readNormaType(
        NormaType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
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


  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readPrecisionType(
        PrecisionType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
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


  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readOutputType(
        OutputType & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      getToken(val, ln, pos);

      if (0 == strcasecmp(val.c_str(), "TERMINAL"))
        field = TERMINAL;
      else if (0 == strcasecmp(val.c_str(), "RES"))
        field = RES;
      else if (0 == strcasecmp(val.c_str(), "GEN")) {
        field = GEN;
        // MyExit(1, "readOutputType:   GEN case not ready");
      } else if (0 == strcasecmp(val.c_str(), "TEX")) {
        field = TEX;
        MyExit(1, "readOutputType:   TEX case not ready");
      } else
        MyExit(1, "readOutputType:   NO SUCH CASE");
    }


  //===========================================================================

  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::readPreRed(
        PreReductionType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
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

  /**
   * Sets the norm in config depending on it's normalizer.
   * */
  template<typename Int, typename BasIntMat>
    void initNorm(LatticeTesterConfig<Int, BasIntMat> & config)
    {
      switch(config.normalizer){
        case BESTLAT: 
        case LAMINATED: 
        case ROGERS: 
        case NORMA_GENERIC: 
        case MINKOWSKI: 
        case L2:
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


  template<typename Int, typename BasInt, typename RedDbl>
    void ParamReader<Int, BasInt, RedDbl>::read (
        LatticeTesterConfig<Int, BasIntMat> & config)
    {
      getLines ();
      unsigned int ln = 0;

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
          std::int64_t m1, m2, m3;
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
