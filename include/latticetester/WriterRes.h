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

#ifndef LATTESTWRITERRES_H
#define LATTESTWRITERRES_H
#include "latticetester/Writer.h"

#include <iostream>
#include <fstream>
#include <string>


namespace LatticeTester {

  /**
   * This class is a simple implementation of the `Writer` abstract class
   * to write in plain text format on the stream.
   */
  template<typename Int>
    class WriterRes : public Writer<Int> {
      public:

        /**
         * Constructor. Opens the writer to write in file `fileName`. `margins` 
         * is the number of white spaces that will be used as a margin inside a 
         * table’s cell when printing a table.
         */
        WriterRes (const char* fileName, unsigned int margins = 5);

        /**
         * Same as above, except that the writer is opened to write directly
         * into `stream`.
         */
        WriterRes (std::ostream* stream, unsigned int margins = 5);

        /**
         * Defined in `Writer`.
         */
        void beginTabbedSection() {}

        /**
         * Defined in `Writer`.
         */
        void endTabbedSection();

        /**
         * Defined in `Writer`.
         */
        void addTab();

        /**
         * Defined in `Writer`.
         */
        void removeTab();

        /**
         * Defined in `Writer`.
         */
        void clearTab();

        /**
         * Defined in `Writer`.
         */
        void newLine();

        /**
         * Defined in `Writer`.
         */
        void newParagraph();

        /**
         * Defined in `Writer`.
         */
        void writeMathString (const std::string);

        /**
         * Defined in `Writer`.
         */
        void writeStandOutMathString (const std::string);

      private:

        /**
         * The number of white space used as a margin in a table’s cell when
         * printing.
         */
        unsigned int m_margins;

        /**
         * Used to remember the state of the tabs and white spaces needed to
         * align correctly a line, for instance when in a tabbed section.
         */
        std::string m_prefix;

    }; // End class WriterRes

  //===========================================================================

  template<typename Int>
    WriterRes<Int>::WriterRes (const char * fileName, unsigned int margins):
      Writer<Int> (fileName)
  {
    m_margins = margins;
  }


  //===========================================================================

  template<typename Int>
    WriterRes<Int>::WriterRes (std::ostream * stream, unsigned int margins):
      Writer<Int> (stream)
  {
    m_margins = margins;
  }


  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::endTabbedSection ()
    {
      m_prefix = "";
    }


  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::addTab ()
    {
      m_prefix += "\t";
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::removeTab ()
    {
      if (m_prefix.length () > 0) {
        m_prefix.erase (m_prefix.length () - 1, 1);
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::clearTab ()
    {
      m_prefix = "";
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::newLine ()
    {
      *this->m_stream << m_prefix << std::endl;
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::newParagraph ()
    {
      *this->m_stream << std::endl << std::endl;
      endTabbedSection ();
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::writeMathString (const std::string s)
    {
      *this->m_stream << m_prefix << s;
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::writeStandOutMathString (const std::string s)
    {
      *this->m_stream << std::endl << "\t" << s << std::endl;
    }

} // End namespace LatticeTester
#endif
