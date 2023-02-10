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

#ifndef LATTICETESTER_CONFIG_H
#define LATTICETESTER_CONFIG_H

#include "latticetester/Const.h"
#include "latticetester/Util.h"

#include <cstdint>

namespace LatticeTester {

  /**
   * This class contains a configuration that can be passed to a
   * `LatticeAnalysis` object to perform a computation. This class is a little
   * bare bones to avoid the confusion of having too many parameters when
   * working with higher level classes.
   *
   * Note: A "configuration" contains a basis and dual basis (but not an IntLattice object),
   * info on what we want to compute with this lattice (construct a basis, the m-dual,
   * compute a shortest vector, compute an FOM, etc.), and which method we want to use for that.
   *
   * */
  template<typename Int, typename IntMat>
    class Config {
      public:
        /*
         * Each `Config` object stores exactly one of the following classes to store the
         * information that is specific to its problem.
         * */
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
         * File format used to print the results. See `Const` for a definition
         * of the possible output types.
         */
        OutputType outputType;

        /**
         * The number of columns of the matrix stored in this object.
         * */
        int NumCols;

        /**
         * The number of rows of the matrix stored in this object.
         * */
        int NumRows;

        /**
         * The basis matrix read from the file.
         * */
        IntMat basis;

        /**
         * The dual matrix to the basis matrix read from the file. This stores
         * the dual so that the results can be reprinted.
         * */
        IntMat dual_basis;

        /**
         * The scaling factor m used (if it was needed by the test).
         * */
        Int m;

        /**
         * This will store the information specific to the problem this `Config`
         * object is for.
         * */
        union Configuration {
          BasisConfig basis;
          DualConfig dual;
          ReductionConfig reduct;
          ShortestConfig shortest;
          MeritConfig merit;
        } config;

    };

  extern template class Config<std::int64_t, NTL::matrix<std::int64_t>>;
  extern template class Config<NTL::ZZ, NTL::matrix<NTL::ZZ>>;

} // End namespace LatticeTester

#endif
