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

/*
 * \file latticetester/EnumTypes.h
 *
 * This class contains enumeration types and global constants used in LatticeTester.
 *
 * TO DO: I am not sure if it is a good idea to collect all enumeration types in one file
 * like this. It may be better to define them where they belong.    ***************
 */

#ifndef LATTICETESTER__ENUMTYPES_H
#define LATTICETESTER__ENUMTYPES_H
#include <string>
#include <array>

namespace LatticeTester {

  /**
   * The available norm types to measure the length of vectors.
   * For \f$X = (x_1,…,x_t)\f$:<br>
   * `SUPNORM` corresponds to \f$\Vert X\Vert= \max(|x_1|,…,|x_t|)\f$.<br>
   * `L1NORM` corresponds to \f$\Vert X\Vert= |x_1|+\cdots+|x_t|\f$.<br>
   * `L2NORM` corresponds to \f$\Vert X\Vert= (x_1^2+\cdots+x_t^2)^{1/2}\f$.<br>
   * `ZAREMBANORM` corresponds to \f$\Vert X\Vert= \max(1, |x_1|)\cdots\max(1, |x_t|)\f$.
   * Note: Use the names, not the numbers!   *******
   */
  enum NormType { SUPNORM = 1, L1NORM = 2, L2NORM = 3, ZAREMBANORM = 4 };

  /**
   * Different choices of output formats.
   *
   * `TERM`: the results will appear only on the terminal screen.<br>
   * `RES`: the results will be in plain text and sent to a `.res` file.<br>
   * `TEX`: the results will be in a LaTeX file with extension `.tex`.<br>
   * `GEN`: the results will be sent to a file with extension `.gen`.
   */
  enum OutputType { TERM, RES, TEX, GEN };

  /**
   * Types of problems that LatticeTester can handle.
   */
  enum ProblemType {BASIS, DUAL, REDUCTION, SHORTEST, MERIT};

  /**
   * Types of precision that the NTL can use for real numbers:
   * `DOUBLE` -- double
   * `QUADRUPLE` -- quad_float (quasi quadruple precision)
   *         this is useful when roundoff errors can cause problems
   * `XDOUBLE` -- xdouble (extended exponent doubles)
   *         this is useful when numbers get too big
   * `RR` -- RR (arbitrary precision floating point).
   * The choice `DOUBLE` is usually the fastest,
   * but may be prone to roundoff errors and/or overflow.
   */
 // enum PrecisionType { DOUBLE, QUADRUPLE, EXPONENT, ARBITRARY, EXACT };
  enum PrecisionType { DOUBLE, QUADRUPLE, XDOUBLE, RR};

  /**
   * Indicates whether an integer is prime, probably prime, composite or its
   * status is unknown (or we do not care).
   */
  enum PrimeType { PRIME, PROB_PRIME, COMPOSITE, UNKNOWN };

  /**
   * Merit criteria to measure the quality of generators or lattices.
   * TO DO: this list is not very clear.    ****************
   *
   * `LENGTH`: Only using the length of the shortest vector as a criterion.
   * `SPECTRAL`: figure of merit \f$S_T\f$ based on the spectral test.<br>
   * `BEYER`: figure of merit is the Beyer quotient \f$Q_T\f$.<br>
   * `PALPHA`: figure of merit based on \f$P_{\alpha}\f$.<br>
   * <tt>BOUND_JS</tt>: figure of merit based on
   *     the Joe-Sinescu bound \cite rSIN08a.<br>   ???
   */
  enum CriterionType { LENGTH, SPECTRAL, BEYER, PALPHA, BOUND_JS };

  /**
   * Different types of normalizations that can be used for shortest-vector lengths.
   * Corresponds to different ways of approximating the Hermite constants `gamma_t`.
   *
   * `BESTLAT`: the value used for \f$d_t^*\f$ corresponds to the best
   * lattice.<br>
   * `BESTBOUND`: the value used for \f$d_t^*\f$ corresponds to the best
   * bound known to us.<br>
   * `LAMINATED`: the value used for \f$d_t^*\f$ corresponds to the best
   * *laminated* lattice.<br>
   * `ROGERS`: the value for \f$d_t^*\f$ is obtained from *Rogers’* bound on the
   * density of sphere packing.<br>
   * `MINKL1`: the value for \f$d_t^*\f$ is obtained from the theoretical bounds
   * on the length of the shortest nonzero vector in the lattice using the
   * \f${\mathcal{L}}_1\f$ norm.<br>
   * `MINKL2`: the value for \f$d_t^*\f$ is obtained from *Minkowski’*
   * theoretical bounds on the length of the shortest nonzero vector in the
   * lattice using the \f${\mathcal{L}}_2\f$ norm.<br>
   * `NONE`: no normalization will be used.<br>
   */
 // enum NormaType { BESTLAT, BESTBOUND, LAMINATED, ROGERS, MINKOWSKI, MINKL1, MINK,L1,L2, NONE };
  enum NormaType { BESTLAT, BESTBOUND, LAMINATED, ROGERS, MINKL1, MINKL2, NONE };

  /**
   * Indicates which type of calculation is considered for the
   * \f$P_{\alpha}\f$ test. \anchor REF__Const_CalcType_def
   *
   * `PAL` is for the \f$P_{\alpha}\f$ test. <br>
   * `BAL` is for the bound on the \f$P_{\alpha}\f$ test. <br>
   * `NORMPAL` is for the \f$P_{\alpha}\f$ test `PAL`, with the result normalized
   *      over the `BAL` bound. <br>
   * `SEEKPAL` is for the \f$P_{\alpha}\f$ seek, which searches
   *      for good values of the multiplier.
   */
  enum CalcType { PAL, NORMPAL, BAL, SEEKPAL };

  /**
   * A list of all the possible lattice reductions implemented in `LatticeTester`.
   *
   * `NORPRERED`: no reduction
   * `DIETER`: Pairwise reduction
   * `LLL`: LLL reduction
   * `BKZ`: block Korkine-Zolotarev reduction
   * `FULL`: shortest vector search.
   */
  enum PreReductionType { NOPRERED, DIETER, LLL, BKZ, FULL };

  /**
   * Two possible ways of obtaining a triangular matrix to compute the bounds
   * in the BB algorithm.
   *
   * `CHOLESKY`: use a lower-triangular matrix obtained as the Cholesky decomposition
   *             of the matrix of scalar products.
   * `TRIANGULAR`: use a lower-triangular basis
   */
  enum DecompType { CHOLESKY, TRIANGULAR };

  /**
   * \name toString functions
   *
   * Useful functions for printing the `enum` constants in this module.
   *
   * @{
   *
   * Returns the value of the `enum` variable given as input as a string.
   */
  std::string toStringNorm (NormType);
  std::string toStringPrime (PrimeType);
  std::string toStringCriterion (CriterionType);
  std::string toStringProblem (ProblemType);
  std::string toStringNorma (NormaType);
  std::string toStringCalc (CalcType);
  std::string toStringPreRed (PreReductionType);
  std::string toStringOutput (OutputType);
  std::string toStringPrecision (PrecisionType);

  /**
   * @}
   */

  static constexpr uint64_t NB_PRIMES = 6543;
  extern const std::array<uint64_t, NB_PRIMES> PRIMES_ARRAY;

}
#endif
