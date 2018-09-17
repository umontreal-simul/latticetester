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

#include <iostream>
#include <fstream>
#include <string>

#include "latticetester/Writer.h"

namespace LatticeTester {

  /**
   * This class is a simple implementation of the `Writer` abstract class
   * to write in plain text format on the stream. This class defines no other
   * method than the interface and can be viewed as an example to what the
   * `Writer` interface is intended to be.
   */
  template<typename Int>
    class WriterRes : public Writer<Int> {
      public:

        /**
         * Constructor that opens the writer to write at the beginning of file
         * `fileName`.
         */
        WriterRes (const char* fileName);

        /**
         * Constructor that will use `stream` as an output stream. This will
         * make the writer start writing at the end of the stream.
         */
        WriterRes (std::ostream* stream);

        /**
         * Implementation of the method defined in `Writer`. In this class, by
         * default, a tabbed section has 2 spaces at the beginning of lines.
         * Calling this methods resets the prefix of tabbed sections to the
         * default.
         */
        void beginTabbedSection();

        /**
         * Implementation of the method defined in `Writer`.
         */
        void endTabbedSection();

        /**
         * Implementation of the method defined in `Writer`. This adds 2 spaces
         * to the prefix of tabbed sections.
         */
        void addTab();

        /**
         * Implementation of the method defined in `Writer`. This removes 2
         * spaces to the prefix of tabbed sections.
         */
        void removeTab();

        /**
         * Implementation of the method defined in `Writer`. This sets the
         * prefix of tabbed sections to an empty string.
         */
        void clearTab();

        /**
         * Implementation of the method defined in `Writer`. This sets the
         * prefix of tabbed sections to a string with 2 spaces.
         */
        void defaultTab();

        /**
         * Implementation of the method defined in `Writer`. Prints a newline
         * character and the prefix of a tabbed section if a the writer is in
         * a tabbed section.
         */
        void newLine();

        /**
         * Implementation of the method defined in `Writer`. Prints 2 newline
         * characters and ends the tabbed section if the program is in a tabbed
         * section.
         */
        void newParagraph();

        /**
         * Implementation of the method defined in `Writer`. This writes the
         * string passed as an argument to the output string.
         */
        void writeMathString (const std::string);

        /**
         * Implementation of the method defined in `Writer`. This writes the 
         * string passed as argument in a tabbed section with a newline
         * character before and after.
         */
        void writeStandOutMathString (const std::string);

      protected:

        /**
         * Indicates if we currently are in a tabbed section. The methods of
         * this class will print m_prefix at each newline if m_isTabbed is `true`.
         * */
        bool m_isTabbed = false;

        /**
         * Used to remember the state of the tabs and white spaces needed to
         * align correctly a line, for instance when in a tabbed section.
         */
        std::string m_prefix;

    }; // End class WriterRes

  //===========================================================================

  template<typename Int>
    WriterRes<Int>::WriterRes (const char * fileName):
      Writer<Int> (fileName)
  { }


  //===========================================================================

  template<typename Int>
    WriterRes<Int>::WriterRes (std::ostream * stream):
      Writer<Int> (stream)
  { }


  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::beginTabbedSection ()
    {
      m_prefix = "  ";
      m_isTabbed = true;
    }


  //===========================================================================
  template<typename Int>
    void WriterRes<Int>::endTabbedSection ()
    {
      m_prefix = "";
      m_isTabbed = false;
    }


  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::addTab ()
    {
      if (m_isTabbed) {
        m_prefix += "  ";
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::removeTab ()
    {
      if (m_isTabbed) {
        if (m_prefix.length () > 0) {
          m_prefix.erase (m_prefix.length () - 1, 2);
        }
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::clearTab ()
    {
      if (m_isTabbed) {
        m_prefix = "";
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::defaultTab ()
    {
      if (m_isTabbed) {
        m_prefix = "  ";
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::newLine ()
    {
      *this->m_stream << std::endl;
      if (m_isTabbed) {
        *this->m_stream << m_prefix;
      }
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::newParagraph ()
    {
      endTabbedSection ();
      newLine();
      newLine();
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::writeMathString (const std::string s)
    {
      *this->m_stream << s;
    }

  //===========================================================================

  template<typename Int>
    void WriterRes<Int>::writeStandOutMathString (const std::string s)
    {
      beginTabbedSection();
      newLine();
      *this->m_stream << s;
      endTabbedSection();
      newLine();
    }

} // End namespace LatticeTester
#endif
