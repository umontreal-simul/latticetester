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

#ifndef LATTICETESTER__FLEXTYPES_H
#define LATTICETESTER__FLEXTYPES_H

  /**
   * The following macros define shorter names for vectors and matrices of Int or Real.
   * This file is used only to avoid repeating these four lines in several class files.
   */
// Simple macros for the preprocessor:
//
// #define IntVec NTL::vector<Int>
// #define IntMat NTL::matrix<Int>
// #define RealVec NTL::vector<Real>
// #define RealMat NTL::matrix<Real>

// Other option (typedef):
//
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;
	typedef NTL::matrix<Real> RealMat;

#endif
