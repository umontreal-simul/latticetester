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

/**
 * \file latticetester/Util.h
 *
 * This module describes various useful functions as well as functions
 * implementing various interactions with NTL.
 */


#ifndef LATTICETESTER__UTIL_H
#define LATTICETESTER__UTIL_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdint>

#include <type_traits>

#include "latticetester/Const.h"

#include "NTL/tools.h"
#include "NTL/ZZ.h"

namespace NTL {

  /**
   * \name Some other compatibility utilities
   * @{
   * These functions perform convertions between different types. Most of them
   * do not really need explainations, but sometimes a specific logic is used
   * when doing the convertion.
   * 
   * \todo Maybe some of this can be rewritten as templates. This could be
   * usefull because it is easy to see how this bit of code could not work with
   * certain compilators.
   */

  /**
   * Converts the array of characters (string) `c` into an `std::int64_t` `l`
   * using the strtol() function of cstdlib.h.
   */
  inline void conv (std::int64_t & l, const char* c) {
    l = strtol(c, (char **) NULL, 10);
  }

  /**
   * Converts the array of characters (string) `c` into a `double` `r` using the
   * strtod() function of cstdlib.h.
   */
  inline void conv (double & r, const char* c) {
    r = strtod(c, (char **) NULL);
  }

  /**
   * Converts a `long long` to a `double.
   * */
  inline void  conv(double &x, long long a)
  {
    x = static_cast<double>(a);
  }

  /**
   * Converts a `double` to a `long long`. This truncates the decimals of a.
   * */
  inline void conv(long long &x, double a)
  {
    x = static_cast<long long>(a);
  }

  /**
   * Converts a `long long` to a `NTL::ZZ`.
   * */
  inline void conv(ZZ &x, long long a)
  {
    x = NTL::conv<ZZ>(a);
  }

  /**
   * Converts a `NTL::ZZ` to a `long long`. This will truncate a if it has to
   * many digits.
   * */
  inline void conv(long long &x, ZZ a)
  {
    x = NTL::conv<long long>(a);
  }

  /**
   * Converts a `long long` to a `long`. This will truncate a if it has to
   * many digits.
   * */
  inline void conv(long &x, long long a)
  {
    x = static_cast<long>(a);
  }

  /**
   * Converts a `long` to a `long long`.
   * */
  inline void conv(long long &x, long a)
  {
    x = static_cast<long long>(a);
  }

  /**
   * Since both are of the same type, this assigns a to x.
   * */
  inline void conv(long long &x, long long a)
  {
    x = a;
  }

  /**
   * @}
   */

}     // namespace NTL


namespace LatticeTester {

#define SEPAR "===============================================================\n"

  /**
   * Maximum integer that can be represented exactly as a `double`:
   * \f$2^{53}\f$.
   */
  const double MAX_LONG_DOUBLE = 9007199254740992.0;

  /**
   * Table of powers of 2: <tt>TWO_EXP[</tt>\f$i\f$<tt>]</tt> \f$= 2^i\f$,
   * \f$i=0, 1, …, 63\f$.
   */
  extern const std::int64_t TWO_EXP[];

  /**
   * Takes references to two variables of a generic type and swaps their
   * content. This uses the assignment operator, so it might not always work
   * well if this operator's implementation is not thorough.
   */
  template <typename T>
    inline void swap9 (T & x, T & y) { T t = x; x = y; y = t; } // WARNING: can we rename this?

  /**
   * \name Random numbers
   *
   * @{
   * All the functions of this module use LFSR258 as an underlying source for
   * pseudo-random numbers. A free (as in freedom) implementation of this
   * generator can be found at http://simul.iro.umontreal.ca/ as well as the
   * article presenting it. All the functions generating some sort of random
   * number will advance an integer version of LFSR258 by one state and output
   * a transformation of the state to give a double, an int or bits.
   */
  /**
   * Returns a random number in \f$[0, 1)\f$. The number will have 53
   * pseudo-random bits.
   */
  double RandU01();

  /**
   * Return a uniform pseudo-random integer in \f$[i, j]\f$. Note that the
   * numbers \f$i\f$ and \f$j\f$ are part of the possible output. It is
   * important that \f$i < j\f$ because the underlying arithmetic uses unsigned
   * integers to store j-i+1 and that will be undefined behavior.
   *
   * \todo For some reason, this only uses the first 62 bits of the generator. Why?
   */
  int RandInt (int i, int j);

  /**
   * Returns the first s pseudo-random bits of the underlying RNG in the form of
   * a s-bit integer. It is imperative that \f$1 \leq s \leq 64\f$ because
   * the RNG is 64 bits wide.
   */
  std::uint64_t RandBits (int s);

  /**
   * Sets the seed of the generator. Because of the constraints on the state,
   * `seed` has to be \f$ > 2\f$. If this is never called, a default seed will be
   * used.
   */
  void SetSeed (std::uint64_t seed);
  /**
   * @}
   */

  /**
   * \name Mathematical functions
   * @{
   * These are basic wrapper and overloads to common mathematical functions
   * such as logarithm, exponentiation, and square-root.
   */

  /**
   * Returns \f$p^i\f$.
   */
  inline std::int64_t power (std::int64_t p, std::int64_t i)
  {
    return NTL::power_long (p, i);
  }

  /**
   * Sets \f$z = 2^i\f$.
   */
  inline void power2 (std::int64_t & z, std::int64_t i)
  {
    z = NTL::power_long (2, i);
  }
  /**
   * Sets \f$z = 2^i\f$.
   */
  inline void power2 (NTL::ZZ& z, std::int64_t i)
  {
    z = NTL::power_ZZ (2, i);
  }

  /**
   * Returns \f$\sqrt{x}\f$ for \f$x\ge0\f$, and \f$-1\f$ for \f$x < 0\f$.
   */
  inline double mysqrt (double x)
  {
    if (x < 0.0)
      return -1.0;
    return sqrt (x);
  }

  /**
   * Returns the logarithm of \f$x\f$ in base 2.
   */
  template <typename T>
    inline double Lg (const T & x)
    {
      return log (x) / 0.6931471805599453094;
    }

  /**
   * Returns the logarithm of \f$x\f$ in base 2.
   */
  inline double Lg (std::int64_t x)
  {
    return log ((double)x) / 0.6931471805599453094;
  }

  /**
   * Returns the absolute value of \f$x\f$.
   */
  template <typename Scal>
    inline Scal abs (Scal x)
    {
      if (x < 0) return -x;  else return x;
    }

  /**
   * Returns the sign of `x`. The sign is 1 if `x>0`, 0 if `x==0` and -1 if 
   * `x<0`.
   */
  template <typename T>
    inline std::int64_t sign (const T & x)
    {
      if (x > 0) return 1; else if (x < 0) return -1; else return 0;
    }

  /**
   * Returns the value of x rounded to the NEAREST integer value. (This does not
   * truncate the integer value as is usual in computer arithmetic.)
   */
  template <typename Real>
    inline Real Round (Real x)
    {
      return floor (x + 0.5);
    }

  /**
   * Calculates \f$t!\f$, the factorial of \f$t\f$ and returns it as an
   * std::int64_t.
   */
  std::int64_t Factorial (int t);

  /**
   * @}
   */

  /**
   * \name Division and modular arithmetic
   *
   * \remark **Richard:** Pour certaines fonctions, les résultats sont mis dans
   * les premiers arguments de la fonction pour être compatible avec NTL; pour
   * d’autres, ils sont mis dans les derniers arguments pour être compatible
   * avec notre ancienne version de LatMRG en Modula-2. Plutôt détestable. Je
   * crois qu’il faudra un jour réarranger les arguments des fonctions pour
   * qu’elles suivent toutes la même convention que NTL.
   *
   * @{
   * This module offers function to perform division and find remainders in a
   * standard way. These functions are usefull in the case where one wants to do
   * divisions or find remainders of operations with negative operands. The
   * reason is that NTL and primitive types do not use the same logic when doing
   * calculations on negative numbers.
   *
   * Basically, NTL will always floor a division and C++ will always truncate a
   * division (which effectively means the floor function is replaced by a roof
   * function if the answer is a negative number). When calculating the
   * remainder of x/y, both apply the same logic but get a different result
   * because they do not do the same division. In both representations, we have
   * that
   * \f[
   *    y\cdot(x/y) + x%y = x.
   * \f]
   * It turns out that, with negative values, NTL will return an integer with
   * the same sign as y where C++ will return an integer of opposite sign
   * (but both will return the same number modulo y).
   */

  /**
   * Computes a/b, truncates the fractionnal part and puts the result in q.
   * This function is overloaded to work as specified on NTL::ZZ integers.
   * Example:
   *
   * <center>
   *
   * <table class="LatSoft-table LatSoft-has-hlines">
   * <tr class="bt">
   *   <td class="r bl br">\f$a\f$</td>
   *   <td class="r bl br">\f$b\f$</td>
   *   <td class="r bl br">\f$q\f$</td>
   * </tr><tr class="bt">
   *   <td class="r bl br">\f$5\f$</td>
   *   <td class="r bl br">3</td>
   *   <td class="r bl br">1</td>
   * </tr><tr>
   *   <td class="r bl br">\f$-5\f$</td>
   *   <td class="r bl br">\f$3\f$</td>
   *   <td class="r bl br">\f$-1\f$</td>
   * </tr><tr>
   *   <td class="r bl br">\f$5\f$</td>
   *   <td class="r bl br">\f$-3\f$</td>
   *   <td class="r bl br">\f$-1\f$</td>
   * </tr><tr>
   *   <td class="r bl br">\f$-5\f$</td>
   *   <td class="r bl br">\f$-3\f$</td>
   *   <td class="r bl br">\f$1\f$</td>
   * </tr>
   * </table>
   *
   * </center>
   */
  template <typename Int>
    inline void Quotient (const Int & a, const Int & b, Int & q)
    {
      q = a/b;
    }

  /**
   * \cond
   * This is the NTL overload of the Quotient function.
   */
  inline void Quotient (const NTL::ZZ & a, const NTL::ZZ & b, NTL::ZZ & q)
  {
    NTL::ZZ r;
    DivRem (q, r, a, b);
    if (q < 0 && r != 0)
      ++q;
  }
  /// \endcond

  /**
   * Computes the remainder of a/b and stores its positive equivalent mod b in
   * r. This works with std::int64_t, NTL::ZZ and real valued numbers.
   * */
  template <typename Real>
    inline void Modulo (const Real & a, const Real & b, Real & r)
    {
      //std::cout << "Modulo Real non testé" << std::endl;
      //exit(1);
      r = fmod (a, b);
      if (r < 0) {
        if (b < 0)
          r -= b;
        else
          r += b;
      }
    }

  /**
   * \cond
   * This is the overload of the Modulo function to work with std::int64_t ints.
   */
  inline void Modulo (const std::int64_t & a, const std::int64_t & b, std::int64_t & r)
  {
    r = a % b;
    if (r < 0) {
      if (b < 0)
        r -= b;
      else
        r += b;
    }
  }

  /**
   * This is the overload of the Modulo function to work with NTL:ZZ.
   */
  inline void Modulo (const NTL::ZZ & a, const NTL::ZZ & b, NTL::ZZ & r)
  {
    r = a % b;
    // if r < 0, we know that b < 0 because of the way NTL works
    if (r < 0)
      r -= b;
  }
  /// \endcond

  /**
   * Computes the quotient \f$q = a/b\f$ and remainder \f$r = a
   * \bmod b\f$. Truncates \f$q\f$ to the nearest integer towards 0. One always
   * has \f$a = qb + r\f$ and \f$|r| < |b|\f$. This works with std::int64_t,
   * NTL::ZZ and real numbers.
   *
   * <center>
   *
   * <table class="LatSoft-table LatSoft-has-hlines">
   * <tr class="bt">
   *   <td class="r bl br">\f$a\f$</td>
   *   <td class="r bl br">\f$b\f$</td>
   *   <td class="r bl br">\f$q\f$</td>
   *   <td class="r bl br">\f$r\f$</td>
   * </tr><tr class="bt">
   *   <td class="r bl br">\f$5\f$</td>
   *   <td class="r bl br">3</td>
   *   <td class="r bl br">1</td>
   *   <td class="r bl br">\f$2\f$</td>
   * </tr><tr>
   *   <td class="r bl br">\f$-5\f$</td>
   *   <td class="r bl br">\f$3\f$</td>
   *   <td class="r bl br">\f$-1\f$</td>
   *   <td class="r bl br">\f$-2\f$</td>
   * </tr><tr>
   *   <td class="r bl br">\f$5\f$</td>
   *   <td class="r bl br">\f$-3\f$</td>
   *   <td class="r bl br">\f$-1\f$</td>
   *   <td class="r bl br">\f$2\f$</td>
   * </tr><tr>
   *   <td class="r bl br">\f$-5\f$</td>
   *   <td class="r bl br">\f$-3\f$</td>
   *   <td class="r bl br">\f$1\f$</td>
   *   <td class="r bl br">\f$-2\f$</td>
   * </tr>
   * </table>
   *
   * </center>
   */
  template <typename Real>
    inline void Divide (Real & q, Real & r, const Real & a, const Real & b)
    {
      q = a / b;
      conv (q, trunc(q));
      r = a - q * b;
    }

  /**
   * \cond
   * This is the overload of the Divide function for std::int64_t.
   */
  inline void Divide (std::int64_t & q, std::int64_t & r,
      const std::int64_t & a, const std::int64_t & b)
  {
    ldiv_t z = ldiv (a, b);
    q = z.quot;      // q = a / b;
    r = z.rem;       // r = a % b;
  }

  /**
   * This is the overload of the Divide function for NTL::ZZ.
   */
  inline void Divide (NTL::ZZ & q, NTL::ZZ & r, const NTL::ZZ & a,
      const NTL::ZZ & b)
  {
    DivRem (q, r, a, b);
    if (q < 0 && r != 0) {
      ++q;
      r -= b;
    }
  }
  /// \endcond

  /**
   * Computes \f$a/b\f$, rounds the result to the nearest integer and returns
   * the result in \f$q\f$. This works with std::int64_t, NTL::ZZ and real numbers.
   */
  template <typename Real>
    inline void DivideRound (const Real & a, const Real & b, Real & q)
    {
      q = a/b;
      q = floor(q + 0.5);
    }

  /**
   * \cond
   * This is the overload of the DivideRound function for std::int64_t.
   */
  inline void DivideRound (const std::int64_t & a, const std::int64_t & b, std::int64_t & q)
  {
    bool neg;
    if ((a > 0 && b < 0) || (a < 0 && b > 0))
      neg = true;
    else
      neg = false;
    std::int64_t x = abs(a);
    std::int64_t y = abs(b);
    ldiv_t z = ldiv (x, y);
    q = z.quot;
    std::int64_t r = z.rem;
    r <<= 1;
    if (r > y)
      ++q;
    if (neg)
      q = -q;
  }

  /**
   * This is the overload of the DivideRound function for NTL::ZZ.
   */
  inline void DivideRound (const NTL::ZZ & a, const NTL::ZZ & b, NTL::ZZ & q)
  {
    bool s = false;
    if ((a > 0 && b < 0) || (a < 0 && b > 0))
      s = true;
    NTL::ZZ r, x, y;
    x = abs (a);
    y = abs (b);
    //****** ATTENTION: bug de NTL: DivRem change le signe de a quand a < 0.
    DivRem (q, r, x, y);
    LeftShift (r, r, 1);
    if (r > y)
      ++q;
    if (s)
      q = -q;
  }
  /// \endcond

  /**
   * Returns the value of the greatest common divisor of \f$a\f$ and \f$b\f$ by
   * using Stein's binary GCD algorithm.
   */
  std::int64_t gcd (std::int64_t a, std::int64_t b);

  /**
   * This method computes the greater common divisor of `A` and `B` with Euclid's
   * algorithm. This will store this gcd in `G` and also the linear combination
   * that permits to get `G` from `A` and `B`. This function should work with
   * std::int64_t and NTL::ZZ.
   *
   * For \f$A\f$ and \f$B\f$ this will assign to \f$C\f$, \f$D\f$, \f$E\f$,
   * \f$F\f$ and \f$G\f$ values such that:
   * \f{align*}{
   *    C a + D b & = G = \mbox{GCD } (a,b)\\
   *   E a + F b & = 0.
   * \f}
   */
  template<typename Int>
    void Euclide (const Int & A, const Int & B, Int & C, Int & D,
        Int & E, Int & F, Int & G)
    {

      //Int oldA = A;
      //Int oldB = B;

      Int X, Y, Z;
      G = A;
      Z = B;
      set9 (C);
      clear (D);
      clear (E);
      set9 (F);

      //cout << "   Euclide :" << endl;
      //cout << "      inputs: A = " << oldA << ", B = " << oldB << endl;
      //cout << "      values used: A=" << A << ", B=" << B << ", C=" << C << ", D=" << D << ", E=" << E << ", F=" << F << ", G=" << G << ", (X=" << X << ", Y=" << Y << ", Z=" << Z << ")" << endl;

      if (IsZero(A)) {
        swap9<Int>(G, Z);
        swap9<Int>(C, E);
        swap9<Int>(D, F);
        return;
      }

      while (!IsZero(Z)) {
        swap9<Int>(G, Z);
        swap9<Int>(C, E);
        swap9<Int>(D, F);
        Quotient (Z, G, X);
        X = -X;
        Y = X * G;
        Z += Y;
        Y = X * C;
        E += Y;
        Y = X * D;
        F += Y;
      }

      //cout << "             " << C << " x A + " << D << " x B = " << G << endl;
      //cout << "             " << E << " x A + " << F << " x B = 0" << endl;
      //cout << "      outputs: A=" << A << ", B=" << B << ", C=" << C << ", D=" << D << ", E=" << E << ", F=" << F << ", G=" << G << endl;
      //cout << endl;
    }

  /**
   * @}
   */

  /**
   * \name Vectors
   *
   * @{
   * These are utilities to manipulate vectors ranging from instanciation to 
   * scalar product.
   *
   * \todo Functions of this section take the size of the input vectors as an
   * argument with name d. Some of them assume d is the number of components of
   * the vectors while somme of them assume that d+1 is the number of components
   * of the vector. We need to uniformize that.
   */

  /**
   * Allocates memory to `A` as an array of `Real` of dimension `d` and
   * initializes its elements to 0. `Real has to be a numeric type.
   */
  template <typename Real>
    inline void CreateVect (Real* & A, int d)
    {
      A = new Real[d];
      for (int i = 0; i < d; i++)
        A[i] = 0;
    }

  /**
   * Creates the vector `A` of dimensions `d+1` and initializes its
   * elements to 0. The type `Vect` has to have a `resize(integer_type)` method
   * that sets the size of the instance to the value of the argument.
   */
  template <typename Vect>
    inline void CreateVect (Vect & A, int d)
    {
      A.resize (d+1);
      for (int i = 0; i < (d+1); i++)
        A[i] = 0;
    }

  /**
   * Frees the memory used by the vector `A`. This calls `delete[]` on `A`
   * so trying to access `A after using this is unsafe.
   */
  template <typename Real>
    inline void DeleteVect (Real* & A)
    {
      delete[] A;
      // A = 0;
    }

  /**
   * Frees the memory used by the vector `A`, destroying all the elements it
   * contains. `Vect` type has to have a `clear()` method that deallocates all
   * the elements in the vector.
   */
  template <typename Vect>
    inline void DeleteVect (Vect & A)
    {
      A.clear ();
    }

  /**
   * Sets the first `d` of `A` to 0.
   */
  template <typename Real>
    inline void SetZero (Real* A, int d)
    {
      for (int i = 0; i < d; i++)
        A[i] = 0;
    }

  /**
   * Sets the first `d` components of `A` to 0.
   */
  template <typename Vect>
    inline void SetZero (Vect & A, int d)
    {
      for (int i = 0; i < d; i++)
        A[i] = 0;
    }

  /**
   * Sets the first `d` components of `A` to the value `x`.
   */
  template <typename Real>
    inline void SetValue (Real* A, int d, const Real & x)
    {
      for (int i = 0; i < d; i++)
        A[i] = x;
    }

  /**
   * Returns a string containing `A[c]` to `A[d-1]` formated as 
   * `[A[c]sep...sepA[d-1]]`. In this string, components are separated by string
   * `sep`. By default, `sep` is just a whitespace character.
   */
  template <typename Vect>
    std::string toString (const Vect & A, int c, int d, const char* sep = " ")
    {
      std::ostringstream out;
      out << "[";
      for (int i = c; i < d-1; i++)
        out << A[i] << sep;
      out << A[d-1] << "]";
      return out.str ();
    }

  /**
   * Returns a string containing the first `d` components of the vector `A` as a
   * string. Calls toString(const Vect&, int, int, const char*).
   */
  template <typename Vect>
    std::string toString (const Vect & A, int d)
    {
      return toString<Vect>(A, 0, d);
      /*
       *  std::ostringstream ostr;
       *  ostr << "[";
       *  for (int i = 0; i < d; i++) {
       *    ostr << std::setprecision(2) << std::setw(3) << A[i] <<
       *      std::setw(2) << "  ";
       *  }
       *  ostr << "]";
       *  return ostr.str();
       * */
    }

  /**
   * Computes the scalar product of vectors `A` and `B` truncated to their `n`
   * first components, then puts the result in `D`. There is a lot to consider
   * when passing types to this function. The best for Vect1 to be the same type
   * as Vect2 and for Scal to be the same as Int, and that those types are the
   * ones stored in Vect1 and Vect2.
   * \todo This uses so many types without check about them and also assumes all
   * those types can be converted to each other without problem... The types
   * here should be uniformized.
   * THIS DOES NOT SEEM TO BE A VERY SAFE IMPLEMENTATION
   */
  template <typename Int, typename Vect1, typename Vect2, typename Scal>
    inline void ProdScal (const Vect1 & A, const Vect2 & B, int n, Scal & D)
    {
      // Le produit A[i] * B[i] peut déborder, d'où conv.
      Int C;   C = 0;
      for (int i = 0; i < n; i++)
        C += A[i] * B[i];
      NTL::conv(D, C);
    }



  /**
   * Takes an input vector `A` of dimension `n+1` and fill the vector `B` with
   * the values `[-A[n] -A[n-1] ... -A[1] 1]`. `B` is assumed to be of dimension
   * at least `n+1`.
   *
   * **OLD DOC**:
   * Transforms the polynomial \f$ A_0 + A_1x^1 + \cdots + A_nx^n\f$ into
   * \f$x^n - A_1x^{n-1} - \cdots - A_n\f$. The result is put in `B`.
   * \todo This function does not do what the old doc says. Check if it causes a
   * problem where this function is called.
   */
  template<typename IntVec>
    inline void Invert (const IntVec & A, IntVec & B, int n)
    {
      conv(B[n], 1);
      for(int i = 0; i < n; i++){
        B[i] = -A[n - i - 1];
      }
    }

  /**
   * Computes the `norm` norm of vector `V` trunctated to its `n` first
   * components, and puts the result in `S`. `Scal` has to be a floating point
   * type.
   */
  template <typename Vect, typename Scal>
    inline void CalcNorm (const Vect & V, int n, Scal & S, NormType norm)
    {
      Scal y;
      S = 0;
      switch (norm) {
        case L1NORM:
          for (int i = 0; i < n; i++) {
            NTL::conv (y, V[i]);
            S += abs(y);
          }
          break;

        case L2NORM:
          for (int i = 0; i < n; i++) {
            NTL::conv (y, V[i]);
            S += y*y;
          }
          //S = sqrt(S);
          break;

        case SUPNORM:
          for (int i = 0; i < n; i++) {
            NTL::conv (y, abs(V[i]));
            if (y > S)
              S = y;
          }
          break;

        case ZAREMBANORM:
          S = 1.0;
          for (int i = 0; i < n; i++) {
            NTL::conv (y, abs(V[i]));
            if (y > 1.0)
              S *= y;
          }
          break;
        default:
          ;
      }
    }

  /**
   * Copies the first `n` components of vector `B` into vector `A`.
   */
  template <typename Vect>
    inline void CopyVect (Vect & A, const Vect & B, int n)
    {   
      for (int k = 0; k < n; k++)  A[k] = B[k];
    }

  /**
   * Adds the first `n` components of vector `B` multiplied by `x` to first `n`
   * components of vector `A`. This will modify `A`. This does wierd type
   * convertion and might not work well if different types are used.
   */
  template <typename Vect1, typename Vect2, typename Scal>
    inline void ModifVect (Vect1 & A, const Vect2 & B, Scal x, int n)
    {
      typename Vect2::value_type a;
      NTL::conv (a, x);
      for (int i = 0; i < n; i++)
        A[i] = A[i] + B[i]*a;
    }

  /**
   * Changes the sign (multiplies by -1) the first `n` components of vector `A`.
   */
  template <typename Vect>
    inline void ChangeSign (Vect & A, int n)
    {
      for (int i = 0; i < n; i++) A[i] = -A[i];
    }

  /**
   * Computes the greatest common divisor of `V[k],...,V[n-1]`.
   */
  inline std::int64_t GCD2vect (std::vector<std::int64_t> V, int k, int n)
  {
    int i = k + 1;
    std::int64_t r, d, c;
    d = abs (V[k]);
    while (d != 1 && i < n) {
      c = abs (V[i]);
      // This does Euclid's algorithm with the current value of the gcd and the
      // next vector component
      while (c) {
        r = d % c;
        d = c;
        c = r;
      }
      ++i;
    }
    return d;
  }

  /**
   * @}
   */

  /**
   * \name Matrices
   *
   * @{
   */

  /**
   * Allocates memory to a square matrix `A` of dimensions
   * \f$d \times d\f$ and initializes its elements to 0.
   */
  template <typename Real>
    inline void CreateMatr (Real** & A, int d)
    {
      A = new Real*[d];
      for (int i = 0; i < d; i++) {
        A[i] = new Real[d];
        for (int j = 0; j < d; j++)
          A[i][j] = 0; }
    }

  /**
   * Allocates memory for the matrix `A` of dimensions \f$\text{line} \times
   * \text{col}\f$ and initializes its elements to 0.
   */
  template <typename Real>
    inline void CreateMatr (Real** & A, int line, int col)
    {
      A = new Real *[line];
      for (int i = 0; i < line; i++) {
        A[i] = new Real[col];
        for (int j = 0; j < col; j++)
          A[i][j] = 0;
      }
    }

  /**
   * Resizes the matrix `A` to a square matrix of dimensions `d \times d`
   * and re-initializes its elements to 0.
   */
  template<typename IntMat>
    inline void CreateMatr (IntMat& A, int d)
    {
      A.resize (d, d);
      for (int i = 0; i<d; i++) {
        for (int j = 0; j<d; j++) A[i][j] = 0;
      }
      //clear (A);
    }

  /**
   * Resizes the matrix `A` to a matrix of dimensions \f$line \times col\f$
   * and re-initializes its elements to 0.
   */
  template<typename IntMat>
    inline void CreateMatr (IntMat& A, int line, int col)
    {
      A.resize(line, col);
      for (int i = 0; i<line; i++) {
        for (int j = 0; j<col; j++) A[i][j] = 0;
      }
      //clear (A);
    }

  /**
   * Frees the memory used by the \f$d \times d\f$ matrix `A`. This will not
   * free all the memory allocated to `A` if `A` is of greater dimension and it
   * can cause a memory leak.
   */
  template <typename Real>
    inline void DeleteMatr (Real** & A, int d)
    {
      for (int i = d; i >= 0; --i)
        delete[] A[i];
      delete[] A;
      //   A = 0;
    }

  /**
   * Frees the memory used by the matrix `A` of dimension \f$\text{line} \times
   * \text{col}\f$. This will not free all the memory allocated to `A` if `A` is
   * of greater dimension and it can cause a memory leak.
   */
  template<typename Real>
    inline void DeleteMatr (Real** & A, int line, int col)
    {
      for (int i = line; i >= 0; --i)
        delete [] A[i];
      delete [] A;
      //    A = 0;
    }

  /**
   * Calls the `clear()` method on `A`. `A` has to have a `clear()` method that
   * frees the memory allocated to it.
   */
  template<typename IntMat>
    inline void DeleteMatr (IntMat& A)
    {
      A.clear ();
    }

  /**
   * Copies the \f$n \times n\f$ submatrix of the first lines and columns
   * of `B` into matrix `A`. This function does not check for sizes, so `A` and
   * `B` both have to be at leat \f$n \times n\f$.
   */
  template <typename Matr>
    inline void CopyMatr (Matr & A, const Matr & B, int n)
    {
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          A[i][j] = B[i][j];
    }

  /**
   * Copies the \f$\text{line} \times col\f$ submatrix of the first lines and
   * columns of `B` into matrix `A`. This function does not check for sizes, so
   * `A` and `B` both have to be at leat \f$line \times col\f$.
   */
  template <typename Matr>
    inline void CopyMatr (Matr & A, const Matr & B, int line, int col)
    {
      for (int i = 0; i < line; i++)
        for (int j = 0; j < col; j++)
          A[i][j] = B[i][j];
    }

  /**
   * Returns a string that is a representation of `mat`. This string will
   * represent the \f$d1 \times d2\f$ submatrix of the first lines and colums of
   * `mat`.
   */
  template <typename MatT>
    std::string toStr (const MatT & mat, int d1, int d2)
    {
      std::ostringstream ostr;
      for (int i = 0; i < d1; i++) {
        ostr << "[";
        for (int j = 0; j < d2; j++) {
          ostr << std::setprecision (2) << std::setw (6) << mat[i][j] <<
            std::setw (2) << " ";
        }
        ostr << "]\n";
      }
      return ostr.str ();
    }

  /**
   * Checks that the upper \f$\text{dim} \times \text{dim}\f$ submatrix of `A`
   * is triangular modulo `m`. This will return `true` if all the elements under
   * the diagonal are equal to zero modulo `m` and `false` otherwise.
   * If `m` is `0`, this function simply verifies that the matrix is triangular.
   */
  template <typename Matr, typename Int>
    bool CheckTriangular (const Matr & A, long dim, const Int m)
    {
      for (int i = 1; i < dim; i++) {
        for (int j = 0; j < i; j++) {
          if (m != 0) {
            if (A[i][j] % m != 0) {
              std::cout << "******  CheckTriangular failed for element A[" << i <<
                "][" << j << "]" << std::endl;
              return false;
            }
          } else {
            if (A[i][j] != 0) {
              std::cout << "******  CheckTriangular failed for element A[" << i <<
                "][" << j << "]" << std::endl;
              return false;
            }
          }
        }
      }
      return true;
    }

  /**
   * Performs an integer triangularization operation modulo `m` on the
   * matrix `W` to obtain an upper triangular matrix `V`, dual to
   * `W`. However, the matrix `W` will be transformed too in order to
   * preserve duality. Only the first `lin` lines and the first `col` columns
   * of the matrices will be considered.
   * The main idea is to transform a generating family of vectors into a 
   * basis (removing from the family the vectors that are linear combination 
   * of other vectors).
   *
   * Refer to the article: R. Couture and P L'Ecuyer, Orbits and lattices for
   * linear random number generators with composite moduli, Mathematics of
   * Computation, Volume 65, Number 213, bottom of page 199.
   * \todo Check the algorithm used here and what it does.
   */
  template <typename Matr, typename Int>
    void Triangularization (Matr & W, Matr & V, int lin, int col,
        const Int & m)
    {
      Int T1, T2, T3, T4, T5, T6, T7, T8;

      for (int j = 0; j < col; j++) {

        for (int i = 0; i < lin; i++)
          Modulo (W(i,j), m, W(i,j));

        int r = 0;
        while (r < lin-1) {
          while (IsZero (W(r,j)) && r < lin-1)
            ++r;
          if (r < lin-1) {
            int s = r + 1;
            while (IsZero (W(s,j)) && s < lin-1)
              ++s;
            if (!IsZero (W(s,j))) {

              Int temp;
              Euclide (W(r,j), W(s,j), T1, T2, T3, T4, temp);
              W(s,j) = temp;
              /*
                 -- Remark --
                 Previously, the call to function Euclide was made with
                 the following arguments : 
                 (W(r,j), W(s,j), T1, T2, T3, T4, W(s,j))
                 But because this function uses references, the result 
                 was always the same:
                 W(r,j) = W(r,j)
                 W(s,j) = W(r,j)
                 T1 = 0
                 T2 = 1
                 T3 = -1
                 T4 = 1
                 We then added a temporary variable *temp* so that the
                 calculation of coefficients T1, T2, T3 and T4 is
                 performed properly.
                 */

              clear (W(r,j));

              for (int j1 = j + 1; j1 < col; j1++) {
                T5 = T1 * W(r,j1);
                T6 = T2 * W(s,j1);
                T7 = T3 * W(r,j1);
                T8 = T4 * W(s,j1);
                W(s,j1) = T5 + T6;
                Modulo (W(s,j1), m, W(s,j1));
                W(r,j1) = T7 + T8;
                Modulo (W(r,j1), m, W(r,j1));
              }
            } else {
              for (int j1 = j; j1 < col; j1++) {
                std::swap (W(r,j1), W(s,j1));
              }
            }
            r = s;
          }
        }
        if (IsZero (W(lin-1,j))) {

          for (int j1 = 0; j1 < col; j1++) {
            if (j1 != j)
              clear (V(j,j1));
            else
              V(j,j1) = m;
          }
        } else {
          Euclide (W(lin-1,j), m, T1, T2, T3, T4, V(j,j));
          for (int j1 = 0; j1 < j; j1++)
            clear (V(j,j1));
          for (int j1 = j + 1; j1 < col; j1++) {
            T2 = W(lin-1,j1) * T1;
            Modulo (T2, m, V(j,j1));
          }
          Quotient (m, V(j,j), T1);
          for (int j1 = j + 1; j1 < col; j1++) {
            W(lin-1,j1) *= T1;
            Modulo (W(lin-1,j1), m, W(lin-1,j1));
          }
        }
      }

      //CheckTriangular (V, col, m);
    }


  /**
   * Calculates the `m`-dual of the matrix `A`. The result is placed in the
   * matrix `B`. Only the first `d` lines and columns are considered.
   *
   * The vectors of the basis (lines of A) need to verify the properties (i), (ii)
   * (iii), (iv) as described in the article: R. Couture and P L'Ecuyer, Orbits 
   * and lattices for linear random number generators with composite moduli, 
   * Mathematics of Computation, Volume 65, Number 213, bottom of page 199.
   *
   * \todo Check the algorithm used here and explain how it works and what it
   * takes for it to work.
   */

  template <typename Matr, typename Int>
    void CalcDual (const Matr & A, Matr & B, int d, const Int & m)
    {
      for (int i = 0; i < d; i++) {

        for (int j = i + 1; j < d; j++)
          clear (B(i,j));

        // WARNING:
        // Dans l'original, c'est *Quotient* pour Lac et *DivideRound* pour non-Lac ??

        //Quotient(m, A(i,i), B(i,i));
        DivideRound (m, A(i,i), B(i,i));

        for (int j = i - 1; j >= 0; j--) {

          clear (B(i,j));

          for (int k = j + 1; k <= i; k++)
            B(i,j) += A(j,k) * B(i,k);

          if (B(i,j) != 0)
            B(i,j) = -B(i,j);

          // WARNING:
          // Dans l'original, c'est *Quotient* pour Lac et *DivideRound* pour non-Lac ??
          //Quotient(B(i,j), A(j,j), B(i,j));
          DivideRound (B(i,j), A(j,j), B(i,j));
        }

      }
    }

  /**
   * @}
   */

  /**
   * \name Debugging functions
   *
   * @{
   */

  /**
   * Special exit function. `status` is the status code to return to the
   * system, `msg` is the message to print upon exit.
   */
  void MyExit (int status, std::string msg);

  /**
   * @}
   */

  /**
   * \name Printing functions and operators
   * @{
   */

  /**
   * Streaming operator for maps.
   * Formats a map as: \texttt{\{ key1=>val1, ..., keyN=>valN \}}.
   */
  template <class K, class T, class C, class A>
    std::ostream & operator<< (std::ostream & out, const std::map<K,T,C,A> & x)
    {
      out << "{";
      typename std::map<K,T,C,A>::const_iterator it = x.begin();
      if (it != x.end()) {
        out << it->first << "=>" << it->second;
        while (++it != x.end())
          out << ", " << it->first << "=>" << it->second;
      }
      out << "}";
      return out;
    }

  /**
   * Streaming operator for vectors.
   * Formats a pair as: \c (first,second).
   */
  template <class T1, class T2>
    std::ostream & operator<< (std::ostream & out, const std::pair<T1,T2> & x)
    {
      out << "(" << x.first << "," << x.second << ")";
      return out;
    }


  /**
   * Streaming operator for vectors.
   * Formats a vector as: <tt>[ val1, ..., valN ]</tt>.
   */
  template <class T, class A>
    std::ostream & operator<< (std::ostream & out, const std::vector<T,A> & x)
    {
      out << "[";
      typename std::vector<T,A>::const_iterator it = x.begin();
      if (it != x.end()) {
        out << *it;
        while (++it != x.end())
          out << ", " << *it;
      }
      out << "]";
      return out;
    }


  /**
   * Streaming operator for sets.
   * Formats a set as: \texttt{\{ val1, ..., valN \}}.
   */
  template <class K, class C, class A>
    std::ostream & operator<< (std::ostream & out, const std::set<K,C,A> & x)
    {
      out << "{";
      typename std::set<K,C,A>::const_iterator it = x.begin();
      if (it != x.end()) {
        out << *it;
        while (++it != x.end())
          out << ", " << *it;
      }
      out << "}";
      return out;
    }


  /**
   * @}
   */

  /**
   * \name Miscellaneous functions
   * @{
   * These are the other functions implemented in this module, that have various
   * usage in this program.
   * */
  
  /**
   * Tests if `x` is odd. Returns 1 if it is odd, and 0 if it is even.
   * \todo It is useless that this function is implemented. Look if it is removable.
   */
  inline std::int64_t IsOdd (const std::int64_t & x) {
    return x & 1;
  }

  /**
   * Returns the `bool` resulting of the statement `x == 0`.
   * \todo It is useless that this function is implemented. Look if it is removable.
   */
  inline bool IsZero (const std::int64_t & x)
  {
    return x == 0;
  }

  /**
   * Sets `x` to 0.
   */
  inline void clear (double & x) { x = 0; }

  /**
   * Sets `x` to 0.
   */
  inline void clear (std::int64_t & x) { x = 0; }

  /**
   * Sets `x` to 1.
   */
  inline void set9 (std::int64_t & x)
  {
    x = 1;
  }

  /**
   * Sets `x` to 1.
   */
  inline void set9 (NTL::ZZ & x)
  {
    NTL::set(x);
  }

  /**
   * @}
   */

}     // namespace LatticeTester

#endif
