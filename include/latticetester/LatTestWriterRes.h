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
class LatTestWriterRes : public LatTestWriter {
public:

/**
 * Constructor. Opens the writer to write in file `fileName`. `margins` is
 * the number of white spaces that will be used as a margin inside a table’s
 * cell when printing a table.
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
   
};

}   
#endif
