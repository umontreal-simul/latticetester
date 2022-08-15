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

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/Reducer.h"
#include <iostream>

namespace LatticeTester{

void negativeCholeski()
  {
    std::cout << "\n***** Negative diagonal element in Choleski Decomposition\n"
      << std::endl;
  }

  template class Reducer<std::int64_t, double, double>;
  template class Reducer<NTL::ZZ, double, double>;
  template class Reducer<NTL::ZZ, NTL::RR, NTL::RR>;

}
