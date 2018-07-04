// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
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

/**
 * \file latticetester/Const.h
 *
 * This module contains global constants used in LatticeTester.
 */

#ifndef LATTICETESTER__CONST_H
#define LATTICETESTER__CONST_H
#include <string>


namespace LatticeTester {

  /**
   * Indicates which norm is used to measure the length of vectors. For \f$X =
   * (x_1,…,x_t)\f$,
   *
   * `SUPNORM` corresponds to
   * \f$\Vert X\Vert= \max(|x_1|,…,|x_t|)\f$.<br> `L1NORM` corresponds to
   * \f$\Vert X\Vert= |x_1|+\cdots+|x_t|\f$.<br> `L2NORM` corresponds to
   * \f$\Vert X\Vert= (x_1^2+\cdots+x_t^2)^{1/2}\f$.<br> `ZAREMBANORM`
   * corresponds to \f$\Vert X\Vert= \max(1, |x_1|)\cdots\max(1, |x_t|)\f$.
   */
  enum NormType { SUPNORM = 1, L1NORM = 2, L2NORM = 3, ZAREMBANORM = 4 };

  /**
   * Indicates in which form and where the results will be sent.
   * \anchor REF__Const_co_output
   *
   * `TERMINAL`: the results will appear only on the terminal screen.<br>
   * `RES`: the results will be in plain text format and sent to a file with
   * extension `.res`.<br>
   * `TEX`: the results will be in LaTeX format and sent to a file with extension
   * `.tex`.\n
   * `GEN`: the results will be sent to a file with extension `.gen`.
   */
  enum OutputType { TERMINAL, RES, TEX, GEN };

  /**
   * Indicates in which precision the NTL algorithms will be perfoms :
   * `FP` -- double
   * `QP` -- quad_float (quasi quadruple precision)
   *         this is useful when roundoff errors can cause problems
   * `XD` -- xdouble (extended exponent doubles)
   *         this is useful when numbers get too big
   * `RR` -- RR (arbitrary precision floating point)
   this is useful for large precision and magnitudes
   * Generally speaking, the choice FP will be the fastest,
   * but may be prone to roundoff errors and/or overflow.
   */
  enum PrecisionType { DOUBLE, QUADRUPLE, EXPONENT, ARBITRARY, EXACT };

  /**
   * Indicates whether an integer is prime, probably prime, composite or its
   * status is unknown (or don’t care).
   */
  enum PrimeType { UNKNOWN, PRIME, PROB_PRIME, COMPOSITE };

  /**
   * Gives the merit criterion for ranking generators or lattices.
   *
   * `BEYER`: the figure of merit is the Beyer quotient \f$Q_T\f$.<br>
   * `SPECTRAL`: the figure of merit \f$S_T\f$ is based on the spectral
   * test.<br>
   * `PALPHA`: the figure of merit is based on \f$P_{\alpha}\f$.<br>
   * <tt>BOUND_JS</tt>: the figure of merit is based on
   * the Joe-Sinescu bound \cite rSIN08a&thinsp;.
   */
  enum CriterionType { SPECTRAL, BEYER, PALPHA, BOUND_JS };

  /**
   * Indicates which normalization is used to compute \f$S_t\f$ in the spectral
   * test, for each dimension \f$t\f$.
   *
   * `BESTLAT`: the value used for \f$d_t^*\f$ corresponds to the best
   * lattice.<br>
   * `LAMINATED`: the value used for \f$d_t^*\f$ corresponds to the best
   * *laminated* lattice.<br>
   * `ROGERS`: the value for \f$d_t^*\f$ is obtained from *Rogers’* bound on the
   * density of sphere packing.<br>
   * `MINKOWSKI`: the value for \f$d_t^*\f$ is obtained from *Minkowski’*
   * theoretical bounds on the length of the shortest nonzero vector in the
   * lattice using the \f${\mathcal{L}}_2\f$ norm.<br>
   * `MINKL1`: the value for \f$d_t^*\f$ is obtained from the theoretical bounds
   * on the length of the shortest nonzero vector in the lattice using the
   * \f${\mathcal{L}}_1\f$ norm.<br>
   * <tt>PALPHA_N</tt>: the case of the \f$P_{\alpha}\f$ test.<br>
   * <tt>NORMA_GENERIC</tt>: the trivial normalization (= 1) used for the generic
   * case when no useful normalization constant is known.
   */
  enum NormaType { BESTLAT, LAMINATED, ROGERS, MINKOWSKI, MINKL1,
    PALPHA_N, NORMA_GENERIC, L1, L2 };

  /**
   * Indicates which type of calculation is considered for the
   * \f$P_{\alpha}\f$ test. \anchor REF__Const_CalcType_def
   *
   * `PAL` is for the \f$P_{\alpha}\f$ test. <br>
   * `BAL` is for the bound on the \f$P_{\alpha}\f$ test. <br>
   * `NORMPAL` is for the \f$P_{\alpha}\f$ test `PAL`, with the result normalized
   * over the `BAL` bound. <br>
   * `SEEKPAL` is for the \f$P_{\alpha}\f$ seek, which searches
   * for good values of the multiplier.
   */
  enum CalcType { PAL, NORMPAL, BAL, SEEKPAL };


  /**
   * Indicates the Prereduction Type (BKZ, LenstraLL, ...) used before applying the 
   * Branch and Bound procedure.
   */

  enum PreReductionType {BKZ, PreRedDieter, LenstraLL, NOPRERED};

  /**
   * \name toString functions
   *
   * Useful functions for printing the `enum` constants in this module.
   *
   * @{
   *
   * Returns the `enum` constants in this module as strings.
   */
  std::string toStringNorm (NormType);
  std::string toStringPrime (PrimeType);
  std::string toStringCriterion (CriterionType);
  std::string toStringNorma (NormaType);
  std::string toStringCalc (CalcType);
  std::string toStringPreRed (PreReductionType);
  std::string toStringOutput (OutputType);
  std::string toStringPrecision (PrecisionType);


  /**
   * @}
   */
}
#endif
