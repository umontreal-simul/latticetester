// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
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

#ifndef LATTESTWRITER_H
#define LATTESTWRITER_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>

#include "latticetester/NTLWrap.h"

namespace LatticeTester {

  /**
   * This is an abstract class that represents an interface to `Writer` classes.
   * That is, this class methods are intended to be used as a way to format an
   * output on an `ostream` or a file. It already implements functions to write
   * basic data types. A subclass of this class **must** be implemented if you
   * want to use it. `WriterRes` offers utilities to print in a file or on the
   * command line.
   * Since this class must be overridden, it is possible to implement subclasses
   * to format text for different kind of output such as plain text, \f$\LaTeX\f$
   * or HTML.
   *
   * If this class is used on a file, it will automatically print at the start
   * of the file and overwrite it. But if a stream is passed to this object, it
   * will not be modified. That means that it is possible to write on a single
   * stream with multiple writers.
   */
  template<typename Int>
    class Writer {
      private:
        typedef NTL::matrix<Int> IntMat;
      public:

        /**
         * Constructor that opens a `Writer` to write in the file `filename`.
         * This will overwrite the file and start printing at the start. The
         * program will make/ovewrite a file in the working directory from which
         * the program was called.
         */
        Writer (const char* fileName);

        /**
         * Constructor that opens a `Writer` to write directly in an `ostream`.
         * This will not change the `stream` object. That means that it is
         * possible to create a stream and pass it to different implementations
         * of this interface via this constructor to use different printing
         * formats.
         */
        Writer (std::ostream* stream);

        /**
         * Destructor.
         */
        virtual ~Writer();

        /**
         * Begins a tabbed section. In a tabbed section, every newline character
         * is followed by a series of whitespace or tab characters such that the
         * first non-blank characters of the lines in this section are aligned.
         * The amount of tabs in the begginning of the section can be modified
         * with the `addTab()`, `removeTab()` and `clearTab()` methods.
         *
         * This is a pure virtual method that has no implementation.
         */
        virtual void beginTabbedSection() = 0;

        /**
         * Ends a tabbed section.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void endTabbedSection() = 0;

        /**
         * After this is called, newlines in tabbed sections will have an
         * additional tab.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void addTab() = 0;

        /**
         * After this is called, newlines in tabbed sections will have one less
         * tab.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void removeTab() = 0;

        /**
         * After this is called, newlines in tabbed sections will have no tab.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void clearTab() = 0;

        /**
         * After this is called, newlines in tabbed sections will use the
         * default amount of tabs and spaces of the class.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void defaultTab() = 0;

        /**
         * Starts a new line and adds the right amount of tabs if in a tabbed 
         * section.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void newLine() = 0;

        /**
         * Starts a new paragraph. Automatically ends the tabbed section if the
         * writer was in one.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void newParagraph() = 0;

        /**
         * Writes a `bool` on the stream.
         */
        virtual void writeBool (const bool & value);

        /**
         * Writes an `int` on the stream.
         */
        virtual void writeInt (const int & value);

        /**
         * Writes a `string` on the stream.
         */
        virtual void writeString (const std::string & value);

        /**
         * Writes a `double` on the stream.
         */
        virtual void writeDouble (const double & value);

        /**
         * Writes a `IntScal` on the stream.
         */
        virtual void writeIntScal (const Int & value);

        /**
         * Writes a `IntMat` on the stream.
         */
        virtual void writeMMat(const IntMat & A);

        /**
         * Writes a `string` that is formated and in LaTeX math mode. This can
         * be used to print mathemetical formulas in a standardised format.
         *
         * This method can be implemented in various way. For example, it is
         * possible to have it print with or without the LaTeX $ signs.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void writeMathString (const std::string) = 0;

        /**
         * Writes a `string` that is formated and in LaTeX math mode with a
         * newline after and before the string. This can be used to print
         * mathemetical formulas in a standardised format.
         * 
         * This method can be implemented in various way. For example, it is
         * possible to have it print with or without the LaTeX \\equation
         * delimiters.
         * 
         * This is a pure virtual method that has no implementation.
         */
        virtual void writeStandOutMathString (const std::string) = 0;

        /**
         * Returns the stream on which this object writes. If the stream is a
         * simple `ostream`, this can be used to then print on standard output
         * for example. If this is a file. It can be possible to used different
         * writers in succession.
         * */
        virtual std::ostream& getStream() { return *m_stream; }

      protected:

        /**
         * The stream in which the writings are done.
         */
        std::ostream* m_stream;

      private:

        /**
         * Set to `true` if this object needs to clean <tt>m_stream</tt> upon
         * destruction.
         */
        bool m_clean;
    }; // class Writer

  //===========================================================================

  template<typename Int>
    Writer<Int>::Writer(const char* fileName)
    {
      m_stream = new std::ofstream(fileName);
      //_stream = dynamic_cast<ostream*>(new ofstream(fileName));
      if (m_stream->fail()) {
        throw 1;
      }
      m_clean = true;
    }

  //===========================================================================

  template<typename Int>
    Writer<Int>::Writer(std::ostream* stream)
    {
      m_stream = stream;
      m_clean = false;
    }

  //===========================================================================

  template<typename Int>
    Writer<Int>::~Writer()
    {
      if (m_clean && m_stream) {
        delete m_stream;
      }
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeInt(const int & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeBool(const bool & value)
    {
      *m_stream << (value ? "true" : "false");
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeString(const std::string & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeIntScal(const Int & value)
    {
      *m_stream << value;
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeMMat(const IntMat & A)
    {
      std::int64_t sizeA = A.size1();
      *m_stream << "   [";
      for (int i = 0; i < sizeA; i++) {
        if (i == 0) { *m_stream << "["; }
        else { *m_stream << "    ["; }

        for (int j = 0; j < (sizeA-1); j++)
          *m_stream << A[i][j] << " ";

        if (i == (sizeA-1)) { *m_stream << A[i][sizeA-1] << "]"; }
        else { *m_stream << A[i][sizeA-1] << "]" << std::endl; }
      }
      *m_stream << "]" << std::endl;
    }

  //===========================================================================

  template<typename Int>
    void Writer<Int>::writeDouble(const double & value)
    {
      *m_stream << value;
    }

} // End namespace LatticeTester
#endif
