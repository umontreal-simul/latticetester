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

#ifndef LATTICETESTERCONFIG_H
#define LATTICETESTERCONFIG_H

#include "latticetester/Const.h"
#include "latticetester/Util.h"

#include <cstdint>

namespace LatticeTester {

  /**
   * This class contains a configuration that can be passed to a
   * `LatticeAnalysis` object to perform a computation. This class is a little
   * bare bones to avoid the confusion of having too many parameters when
   * working with higher level classes.
   * */
  template<typename Int, typename BasIntMat>
    class Config {
      public:
        class BasisConfig {
          public:
            bool method;
        };

        class DualConfig {
        };

        class ReductionConfig {
          public:
            PreReductionType method;
        };

        class ShortestConfig {
          public:
            bool reduction;
            PreReductionType method;
        };

        class MeritConfig {
          public:
            CriterionType figure;
            NormaType norma;
            bool reduction;
            PreReductionType method;
        };

        /**
         * The name of the file that populated this object.
         * */
        std::string filename;

        /**
         * The problem this object stores configuration for. All the fields that
         * are not relevant to this problem might contain meaningless
         * information.
         * */
        ProblemType prob;

        /**
         * File format used to store the results. See `Const` for a definition
         * of the possible output types.
         */
        OutputType outputType;

        /**
         * */
        int NumCols;

        /**
         * */
        int NumRows;

        /**
         * */
        BasIntMat basis;
        BasIntMat dual_basis;
        Int m;

        union Configuration {
          BasisConfig basis;
          DualConfig dual;
          ReductionConfig reduct;
          ShortestConfig shortest;
          MeritConfig merit;
        } config;

    };

} // End namespace LatticeTester
#endif
