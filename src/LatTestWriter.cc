/*
 Writer.cc for ISO C++
 version 1.00
 authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard

*/

#include <iostream>
#include <fstream>
#include <sstream>

#include "latticetester/LatTestWriter.h"
#include "latticetester/Types.h"
#include "latticetester/Const.h"

using namespace std;

namespace LatticeTester
{

LatTestWriter::LatTestWriter(const char* fileName)
{
   m_stream = new ofstream(fileName);
   //_stream = dynamic_cast<ostream*>(new ofstream(fileName));
   if (m_stream->fail()) {
      throw 1;
   }
   m_clean = true;
}

LatTestWriter::LatTestWriter(ostream* stream)
{
   m_stream = stream;
   m_clean = false;
}

LatTestWriter::~LatTestWriter()
{
   if (m_clean && m_stream) {
      delete m_stream;
   }
}

void LatTestWriter::writeInt(const int & value)
{
   *m_stream << value;
}

void LatTestWriter::writeBool(const bool & value)
{
   *m_stream << (value ? "true" : "false");
}

void LatTestWriter::writeString(const string & value)
{
   *m_stream << value;
}

void LatTestWriter::writeMScal(const MScal & value)
{
   *m_stream << value;
}

void LatTestWriter::writeMMat(const MMat & A)
{
   int sizeA = A.size1();
   *m_stream << "   [";
   for (int i = 0; i < sizeA; i++) {
      if (i == 0) { *m_stream << "["; }
      else { *m_stream << "    ["; }

      for (int j = 0; j < (sizeA-1); j++)
        *m_stream << A[i][j] << " ";

      if (i == (sizeA-1)) { *m_stream << A[i][sizeA-1] << "]"; }
      else { *m_stream << A[i][sizeA-1] << "]" << endl; }
   }
   *m_stream << "]" << endl;
}

void LatTestWriter::writeDouble(const double & value)
{
   *m_stream << value;
}

}
