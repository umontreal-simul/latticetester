#ifndef PARAMREADER_H
#define PARAMREADER_H
#include "NTL/ZZ.h"
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/LatticeTesterConfig.h"

#include <string>
#include <vector>


namespace LatticeTester {

/**
 * Utility class used to read basic parameter fields in a configuration file.
 * Lines whose first non-blank character is a <tt>#</tt> are considered as
 * comments and discarded.
 *
 */
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
   void readNumber3 (MScal & r, long & b, long & e, long & c,
                     unsigned int ln, unsigned int pos);

   /**
    * Reads a `BScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readBScal (BScal & field, unsigned int ln, int pos);

   /**
    * Reads a `BMat` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readBMat (BMat & fields, unsigned int & ln, unsigned int pos,
                 unsigned int numPos);

   /**
    * Reads a `MScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readMScal (MScal & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a `RScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readRScal (RScal & field, unsigned int ln, unsigned int pos);

   /**
    * Reads `num` tokens (from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
    * `field`.
    */
   void readMVect (MVect & field, unsigned int & ln, unsigned int pos,
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
   void readInterval (MVect & B, MVect & C, unsigned int & ln, int k);

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
   void read (LatticeTesterConfig & config);
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
};

}
#endif
