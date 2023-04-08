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
 *
 * \file latticetester/Num.h
 *
 * This module implements a few mathematical functions of obscure relevance.
 */

#ifndef LATTICETESTER__NUM_H
#define LATTICETESTER__NUM_H

#include <cstdint>

namespace LatticeTester {

  /**
   * Calculates \f$t!\f$, the factorial of \f$t\f$. Might throw if `t` is too
   * large or if std::int64_t can't contain the factorial asked for.
   */
  std::int64_t lFactorial (int64_t t);

  /**
   * Returns the value of the logarithmic derivative of the Gamma function
   * \f$\psi(x) = \Gamma'(x) / \Gamma(x)\f$.
   */
  double Digamma (double x);

  /**
   * Evaluates the Bernoulli polynomial \f$B_n(x)\f$ of degree \f$n\f$
   * at \f$x\f$.
   * The first Bernoulli polynomials are:
   * \f{align*}{
   * B_0(x) &= 1  \\
   * B_1(x) &= x - 1/2  \\
   * B_2(x) &= x^2-x+1/6  \\
   * B_3(x) &= x^3 - 3x^2/2 + x/2  \\
   * B_4(x) &= x^4-2x^3+x^2-1/30 \\
   * B_5(x) &= x^5 - 5x^4/2 + 5x^3/3 - x/6  \\
   * B_6(x) &= x^6-3x^5+5x^4/2-x^2/2+1/42  \\
   * B_7(x) &= x^7 - 7x^6/2 +  7x^5/2 - 7x^3/6 + x/6  \\
   * B_8(x) &= x^8-4x^7+14x^6/3 - 7x^4/3 +2x^2/3-1/30.
   * \f}
   * Only degrees \f$n \le 8\f$ are programmed for now.
   */
  double BernoulliPoly (int64_t n, double x);

  /**
   * Computes the \f$n\f$-th harmonic number \f$H_n  = \sum_{j=1}^n 1/j\f$.
   */
  double Harmonic (std::int64_t n);

  /**
   * Computes the sum
   * \f[
   * \sideset{}{'}\sum_{-n/2<j\le n/2}\; \frac 1{|j|},
   * \f]
   * where the symbol \f$\sum^\prime\f$ means that the term with \f$j=0\f$ is excluded
   * from the sum.
   */
  double Harmonic2 (std::int64_t n);

  /**
   * Computes and returns the value of the series (see \cite vJOE92b)
   * \f[
   * S(x, n) = \sum_{j=1}^{n} \frac{\cos(2\pi j x)}{j}.
   * \f]
   * Restrictions: \f$n>0\f$ and \f$0 \le x \le 1\f$.
   */
  double FourierC1 (double x, std::int64_t n);

  /**
   * Computes and returns the value of the series
   * \f[
   * G(x, n) = \sideset{}{'}\sum_{-n/2<h\le n/2}\;  \frac{e^{2\pi i h x}}{|h|},
   * \f]
   * where the symbol \f$\sum^\prime\f$ means that the term with \f$h=0\f$ is excluded
   * from the sum, and assuming that the imaginary part of \f$G(x, n)\f$ vanishes.
   * Restrictions: \f$n>0\f$ and \f$0 \le x \le 1\f$.
   */
  double FourierE1 (double x, std::int64_t n);

}

#endif

