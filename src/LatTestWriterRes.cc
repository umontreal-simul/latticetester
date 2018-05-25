/**
 WriterRes.cc for ISO C++
 version 1.00
 modified 02/06/06 15:03
 authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard
*/

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "latticetester/LatTestWriterRes.h"

using namespace std;

namespace LatticeTester
{

//=========================================================================

LatTestWriterRes::LatTestWriterRes (const char * fileName, unsigned int margins):
      LatTestWriter (fileName)
{
   m_margins = margins;
}


//=========================================================================

LatTestWriterRes::LatTestWriterRes (ostream * stream, unsigned int margins): LatTestWriter (stream)
{
   m_margins = margins;
}


//=========================================================================

void LatTestWriterRes::endTabbedSection ()
{
   m_prefix = "";
}


//=========================================================================

void LatTestWriterRes::addTab ()
{
   m_prefix += "\t";
}


//=========================================================================
void LatTestWriterRes::removeTab ()
{
   if (m_prefix.length () > 0) {
      m_prefix.erase (m_prefix.length () - 1, 1);
   }
}


//=========================================================================

void LatTestWriterRes::clearTab ()
{
   m_prefix = "";
}


//=========================================================================

void LatTestWriterRes::newLine ()
{
   *m_stream << m_prefix << endl;
}


//=========================================================================

void LatTestWriterRes::newParagraph ()
{
   *m_stream << endl << endl;
   endTabbedSection ();
}


//=========================================================================

void LatTestWriterRes::writeMathString (const string s)
{
   *m_stream << m_prefix << s;
}


//=========================================================================

void LatTestWriterRes::writeStandOutMathString (const string s)
{
   *m_stream << endl << "\t" << s << endl;
}


//=========================================================================

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
            *m_stream << " ";
         }

         *m_stream << sub_data[j][i];

         for (unsigned int k = 0;
               k < max_length[j] - length - offset + m_margins; k++) {
            *m_stream << " ";
         }
      }
      if (i == nb_lig[0] - 1) {
         *m_stream << endl;
         for (int k = 0; k < total_length; k++) {
            *m_stream << "-";
         }
      }
      *m_stream << endl;
   }

   // deleting pointers
   for (int i = 0; i < t; ++i)
      delete [] sub_data[i];
   delete [] sub_data;

}

*/

//=========================================================================

}
