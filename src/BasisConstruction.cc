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

#include "NTL/ZZ.h"

#include "latticetester/BasisConstruction.h"

namespace LatticeTester {
  // Instantiation of supported types for the class
  template class BasisConstruction<std::int64_t>;
  template class BasisConstruction<NTL::ZZ>;
  
  //============================================================================
  // Specialization of LLLConstr
  //============================================================================

  template <>
    void LLLConstr<NTL::matrix<std::int64_t>>::LLLConstruction(
        NTL::matrix<std::int64_t>& matrix) {
      std::cerr << "LLL Construction can only be done with NTL::ZZ integers.\n";
      std::cerr << "Aborting.\n";
      exit(1);
    }

  //============================================================================

  template <>
    void LLLConstr<NTL::matrix<NTL::ZZ>>::LLLConstruction(
        NTL::matrix<NTL::ZZ>& matrix) {
      long rank = NTL::LLL_XD(matrix);
      long num = matrix.NumRows();
      for (long i = 0; i < rank; i++) {
        NTL::swap(matrix[i], matrix[num-rank+i]);
      }
      matrix.SetDims(rank, matrix.NumCols());
    }

}
