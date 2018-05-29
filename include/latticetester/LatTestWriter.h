#ifndef LATTESTWRITER_H
#define LATTESTWRITER_H
#include "latticetester/Types.h"

#include <iostream>
#include <string>
#include <fstream>


namespace LatticeTester {

  /**
   * This is the abstract class that does the writing of basic elements
   * (<tt>string</tt>’s, <tt>int</tt>’s, <tt>double</tt>’s, etc.) into a file
   * or into an `ostream`. Derived classes must be implemented to write in
   * different formats, for instance text or LaTeX.
   *
   */
  template<typename Int, typename IntMat>
    class LatTestWriter {
      public:

        /**
         * Constructor. Opens a `Writer` to write in the file `filename`.
         */
        LatTestWriter (const char* fileName);

        /**
         * Constructor. Opens a `Writer` to write directly in an `ostream`.
         */
        LatTestWriter (std::ostream* stream);

        /**
         * Destructor.
         */
        virtual ~LatTestWriter();

        /**
         * Begins a tabbed section. In a tabbed section, every element is
         * aligned at a tab start.
         */
        virtual void beginTabbedSection() = 0;

        /**
         * Ends a tabbed section.
         */
        virtual void endTabbedSection() = 0;

        /**
         * Advances the tab start once.
         */
        virtual void addTab() = 0;

        /**
         * Moves back the tab start once.
         */
        virtual void removeTab() = 0;

        /**
         * Resets the tab start to the beginning of the line.
         */
        virtual void clearTab() = 0;

        /**
         * Starts a new line. If in a tabbed section, the new line starts at
         * the current tab start position.
         */
        virtual void newLine() = 0;

        /**
         * Starts a new paragraph. Ends automaticaly a tabbed section.
         */
        virtual void newParagraph() = 0;

        /**
         * Writes a `bool`.
         */
        virtual void writeBool (const bool & value);

        /**
         * Writes an `int`.
         */
        virtual void writeInt (const int & value);

        /**
         * Writes a `string`.
         */
        virtual void writeString (const std::string & value);

        /**
         * Writes a `double`.
         */
        virtual void writeDouble (const double & value);

        /**
         * Writes a `IntScal`.
         */
        virtual void writeIntScal (const Int & value);

        /**
         * Writes a `IntMat`.
         */
        virtual void writeMMat(const IntMat & A);

        /**
         * Writes a `string` in a mathematical format using LaTeX notation.
         */
        virtual void writeMathString (const std::string) = 0;

        /**
         * Writes a `string` just as in an `equation` environment in LaTeX.
         */
        virtual void writeStandOutMathString (const std::string) = 0;

        /**
         * Writes a formatted table. A column horizontal alignment can be
         * modified using the string `alignment`. The first character of the
         * string is the alignment of the first column, the second character is
         * the alignment of the second column and so on. Three alignments are
         * possible for a column: `r` for right alignment, `c` for center
         * alignment, `l` for left alignment. If nothing is given, or if the
         * string length is smaller than the number of columns of the table,
         * then left alignment is used by default for the rest of the columns.
         */
        //virtual void writeTable (Table &, const std::string alignment) = 0;

        virtual std::ostream& getStream() { return *m_stream; }

      protected:

        /**
         * The stream in which the writings are done.
         */
        std::ostream* m_stream;

      private:

        /**
         * Set to `true` if this object needs to clean <tt>_stream</tt> upon
         * destruction.
         */
        bool m_clean;
    }; // class LatTestWriter

  //===========================================================================

  template<typename Int, typename IntMat>
    LatTestWriter<Int, IntMat>::LatTestWriter(const char* fileName)
    {
      m_stream = new ofstream(fileName);
      //_stream = dynamic_cast<ostream*>(new ofstream(fileName));
      if (m_stream->fail()) {
        throw 1;
      }
      m_clean = true;
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    LatTestWriter<Int, IntMat>::LatTestWriter(ostream* stream)
    {
      m_stream = stream;
      m_clean = false;
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    LatTestWriter<Int, IntMat>::~LatTestWriter()
    {
      if (m_clean && m_stream) {
        delete m_stream;
      }
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeInt(const int & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeBool(const bool & value)
    {
      *m_stream << (value ? "true" : "false");
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeString(const string & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeIntScal(const Int & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeMMat(const IntMat & A)
    {
      long sizeA = A.size1();
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

  //===========================================================================

  template<typename Int, typename IntMat>
    void LatTestWriter<Int, IntMat>::writeDouble(const double & value)
    {
      *m_stream << value;
    }

} // End namespace LatticeTester
#endif
