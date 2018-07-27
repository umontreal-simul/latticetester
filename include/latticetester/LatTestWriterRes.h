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
#include "latticetester/LatTestWriter.h"

#include <iostream>
#include <fstream>
#include <string>


namespace LatticeTester {

  /**
   * This class implements the `Writer` abstract class to write basic elements
   * in plain text format.
   *
   */
  template<typename Int>
    class LatTestWriterRes : public LatTestWriter<Int> {
      public:

        /**
         * Constructor. Opens the writer to write in file `fileName`. `margins` 
         * is the number of white spaces that will be used as a margin inside a 
         * table’s cell when printing a table.
         */
        LatTestWriterRes (const char* fileName, unsigned int margins = 5);

        /**
         * Same as above, except that the writer is opened to write directly
         * into `stream`.
         */
        LatTestWriterRes (std::ostream* stream, unsigned int margins = 5);

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

        /**
         * Defined in `Writer`.
         */
        //void writeTable (Table & data, const std::string pos);

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

        /**
         * Set to `true` if the writer is in a tabbed section, otherwise it is
         * `false`.
         */
        //bool m_inTabbed;
        // remark: unused variable

    }; // End class LatTestWriterRes

  //===========================================================================

  template<typename Int>
    LatTestWriterRes<Int>::LatTestWriterRes (const char * fileName,
        unsigned int margins):
      LatTestWriter<Int> (fileName)
  {
    m_margins = margins;
  }


  //===========================================================================

  template<typename Int>
    LatTestWriterRes<Int>::LatTestWriterRes (std::ostream * stream, 
        unsigned int margins): LatTestWriter<Int> (stream)
  {
    m_margins = margins;
  }


  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::endTabbedSection ()
    {
      m_prefix = "";
    }


  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::addTab ()
    {
      m_prefix += "\t";
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::removeTab ()
    {
      if (m_prefix.length () > 0) {
        m_prefix.erase (m_prefix.length () - 1, 1);
      }
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::clearTab ()
    {
      m_prefix = "";
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::newLine ()
    {
      *this->m_stream << m_prefix << std::endl;
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::newParagraph ()
    {
      *this->m_stream << std::endl << std::endl;
      endTabbedSection ();
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::writeMathString (const std::string s)
    {
      *this->m_stream << m_prefix << s;
    }

  //===========================================================================

  template<typename Int>
    void LatTestWriterRes<Int>::writeStandOutMathString (const std::string s)
    {
      *this->m_stream << std::endl << "\t" << s << std::endl;
    }

  //===========================================================================

  /*

     void WriterRes::writeTable (Table & data, const string pos)
     {
     int t = data.size ();       // Nb de colonne
     int n = data.getHeight () + 1; // Nb de ligne

     int max_length[t];
     int nb_lig_cell[t][n];
     int nb_lig[n];
     string::size_type curpos;
     int curpos_b;
     int lig;
     int total_length = 0;

     memset (max_length, 0, t * sizeof (int));
     memset (nb_lig, 0, n * sizeof (int));

     string temp;

  // On calcule le nombre de sous ligne dans chaque case et la longueur
  // maximale de chaque colonne.

  for (int i = 0; i < t; i++) {
  for (int j = -1; j < data[i]->size (); j++) {
  lig = 0;
  curpos = 0;
  curpos_b = -1;

  temp = data[i]->getFormattedValue (j);

  while ((curpos = temp.find ("\n", curpos)) != string::npos) {
  int tmp = curpos - curpos_b - 1;
  max_length[i] = std::max (max_length[i], tmp);
  lig++;
  curpos_b = curpos;
  curpos++;
  }
  if (curpos_b == -1) {
  curpos_b = 0;
  }

  max_length[i] = std::max (max_length[i],
  static_cast <int>(temp.length () - curpos_b));
  nb_lig[j + 1] = std::max (nb_lig[j + 1], ++lig);
  nb_lig_cell[i][j + 1] = lig;
  }

  total_length += max_length[i] + 2 * m_margins;
  }

  int nb_sub_lig = 0;

  for (int i = 0; i < n; i++) {
  nb_sub_lig += nb_lig[i];
  }

  int k, offset;



  //string sub_data[t][nb_sub_lig];
  // this is not compiling with old versions of clang

  string **sub_data = new string* [t];
  for (int i = 0; i < t; ++i)
  sub_data[i] = new string [nb_sub_lig];
  // better could be : string *sub_data = new string [t * nb_sub_lig]



  // On décompose les chaines en sous chaines

  for (int i = 0; i < t; i++) {
    offset = 0;
    for (int j = -1; j < data[i]->size (); j++) {
      curpos = 0;
      curpos_b = 0;

      // padding with empty strings
      for (k = 0; k < nb_lig[j + 1] - nb_lig_cell[i][j + 1]; k++) {
        sub_data[i][j + offset + k + 1] = "";
      }

      temp = data[i]->getFormattedValue (j);

      while ((curpos = temp.find ("\n", curpos)) != string::npos) {
        sub_data[i][j + offset + k + 1] =
          temp.substr (curpos_b, curpos);
        k++;
        curpos_b = curpos;
        curpos++;
      }
      if (curpos_b > 0) {
        curpos_b++;
      }
      sub_data[i][j + offset + k + 1] =
        temp.substr (curpos_b, temp.length () - curpos_b);
      sub_data[i][j + offset + k + 1] = temp.substr (curpos_b, curpos);
      offset += k;
    }
  }

  // On dessine le tableau, finalement ...

  int length;
  char align;
  for (int i = 0; i < nb_sub_lig; i++) {
    for (int j = 0; j < t; j++) {

      align = pos[j];
      length = sub_data[j][i].length ();

      switch (align) {
        case 'c':
          offset = (max_length[j] >> 1) - (length >> 1);
          break;
        case 'r':
          offset = max_length[j] - length;
          break;
        default:               // case l
          offset = 0;

      }

      for (unsigned int k = 0; k < m_margins + offset; k++) {
        *this->m_stream << " ";
      }

      *this->m_stream << sub_data[j][i];

      for (unsigned int k = 0;
          k < max_length[j] - length - offset + m_margins; k++) {
        *this->m_stream << " ";
      }
    }
    if (i == nb_lig[0] - 1) {
      *this->m_stream << std::endl;
      for (int k = 0; k < total_length; k++) {
        *this->m_stream << "-";
      }
    }
    *this->m_stream << std::endl;
  }

  // deleting pointers
  for (int i = 0; i < t; ++i)
    delete [] sub_data[i];
  delete [] sub_data;

}

*/

} // End namespace LatticeTester
#endif
