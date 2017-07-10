#ifndef PARAMREADERLAT_H
#define PARAMREADERLAT_H
#include "latticetester/ParamReader.h"
#include "latticetester/LatticeTesterConfig.h"

#include <string>


namespace LatticeTester {

/**
 * This class is used to read a configuration file for the executable
 * programs `lat*`, created from the `LatMain` main program. The format of
 * the configuration file is described in this guide for the program
 * `LatMain` on page (FIXME: page#).
 *
 */
class ParamReaderLat : public ParamReader {
public:

/**
 * Constructor. Opens configuration file `fileName`.
 */
ParamReaderLat (std::string fileName);

   /**
    * Destructor.
    */
   ~ParamReaderLat();

   /**
    * Reads the configuration file into `config` for the Beyer and the
    * spectral tests.
    */
   void read (LatConfig & config);
};

}
#endif
