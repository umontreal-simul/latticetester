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
 * \file latticetester/Util.h
 *
 * This module describes various useful functions as well as functions
 * interfacing with NTL.
 *
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
#include "latticetester/Types.h"
#include "latticetester/Const.h"

#ifdef WITH_NTL
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#endif


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
extern const long TWO_EXP[];

/**
 * Swaps the values of `x` and `y`.
 */
template <typename T>
inline void swap9 (T & x, T & y) { T t = x; x = y; y = t; } // WARNING: can we rename this?

/**
 * \name Random numbers
 *
 * @{
 */
/**
 * Returns a random number in \f$[0, 1)\f$. The number has 53 random bits.
 */
double RandU01();

/**
 * Return a random integer in \f$[i, j]\f$. Numbers \f$i\f$ and \f$j\f$ can
 * occur. Restriction: \f$i < j\f$.
 */
int RandInt (int i, int j);

/**
 * Returns random blocks of \f$s\f$ bits (\f$s\f$-bit integers).
 */
unsigned long RandBits (int s);

/**
 * Sets the seed of the generator. If not called, a default seed will be
 * used.
 */
void SetSeed (unsigned long seed);

/**
 * @}
 */

/**
 * \name NTL compatibility utilities
 *
 * @{
 */

#ifndef WITH_NTL
/**
 * Converts \f$y\f$ to \f$x\f$.
 */
template <typename A, typename B>
void conv (A & x, const B & y) { x = y; }
#endif

/**
 * Returns 1 if \f$x\f$ is odd, and 0 otherwise.
 */
inline long IsOdd (const long & x) {
    return x & 1;
}

/**
 * Return `true` if \f$x = 0\f$.
 */
inline bool IsZero (const long & x)
{
    return x == 0;
}

/**
 * Sets \f$x\f$ to 0.
 */
inline void clear (double & x) { x = 0; }

/**
 * Sets \f$x\f$ to 0.
 */
inline void clear (long & x) { x = 0; }

/**
 * Sets \f$x\f$ to 1.
 */
inline void set9 (long & x)
{
    x = 1;
}

#ifdef WITH_NTL
/**
 * Sets \f$x\f$ to 1.
 */
inline void set9 (NTL::ZZ & x)
{
    NTL::set(x);
}
#endif

/**
 * @}
 */

/**
 * \name Mathematical functions
 *
 * @{
 */

#ifdef WITH_NTL
/**
 * Returns \f$p^i\f$.
 */
inline long power (long p, long i)
{
   return NTL::power_long (p, i);
}
#endif

#ifdef WITH_NTL
/**
 * Sets \f$z = 2^i\f$.
 */
inline void power2 (long & z, long i)
{
   z = NTL::power_long (2, i);
}
/**
 * Sets \f$z = 2^i\f$.
 */
inline void power2 (ZZ& z, long i)
{
   z = NTL::power_ZZ (2, i);
}
#else
/**
 * Sets \f$z = 2^i\f$.
 */
inline void power2 (long & z, long i)
{
   z = 1L << i;
}
#endif

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
 * Returns \f$\sqrt{x}\f$.
 * \remark **Richard:** Cette fonction est-elle encore utilisée?
 */
inline double SqrRoot (double x)
{
    return sqrt (x);
}

/**
 * Logarithm of \f$x\f$ in base 2.
 */
template <typename T>
inline double Lg (const T & x)
{
   return log (x) / 0.6931471805599453094;
}

/**
 * Logarithm of \f$x\f$ in base 2.
 */
inline double Lg (long x)
{
   return log ((double)x) / 0.6931471805599453094;
}

/**
 * Returns the absolute value.
 */
template <typename Scal>
inline Scal abs (Scal x)
{
   if (x < 0) return -x;  else return x;
}

/**
 * Returns 1, 0 or \f$-1\f$ depending on whether \f$x> 0\f$, \f$x= 0\f$ or
 * \f$x< 0\f$ respectively.
 */
template <typename T>
inline long sign (const T & x)
{
    if (x > 0) return 1; else if (x < 0) return -1; else return 0;
}

/**
 * Rounds to the nearest integer value.
 */
template <typename Real>
inline Real Round (Real x)
{
    return floor (x + 0.5);
}

/**
 * Calculates \f$t!\f$, the factorial of \f$t\f$.
 */
long Factorial (int t);

/**
 * @}
 */


/**
 * \name Division and remainder
 *
 * For negative operands, the `/` and <tt>%</tt> operators do not give the
 * same results for NTL large integers `ZZ` and for primitive types `int` and
 * `long`. The negative quotient differs by 1 and the remainder also differs.
 * Thus the following small `inline` functions for division and remainder.
 * \remark **Richard:** Pour certaines fonctions, les résultats sont mis dans
 * les premiers arguments de la fonction pour être compatible avec NTL; pour
 * d’autres, ils sont mis dans les derniers arguments pour être compatible
 * avec notre ancienne version de LatMRG en Modula-2. Plutôt détestable. Je
 * crois qu’il faudra un jour réarranger les arguments des fonctions pour
 * qu’elles suivent toutes la même convention que NTL.
 *
 * @{
 */

/**
 * Computes \f$q = a/b\f$ by dropping the fractionnal part, i.e. truncates
 * towards 0. Example:
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

#ifdef WITH_NTL
/**
 * \copydoc Quotient(const Int&, const Int&, Int*)
 */
inline void Quotient (const NTL::ZZ & a, const NTL::ZZ & b, NTL::ZZ & q)
{
    NTL::ZZ r;
    DivRem (q, r, a, b);
    if (q < 0 && r != 0)
       ++q;
}
#endif

template <typename Real>
inline void Modulo (const Real & a, const Real & b, Real & r)
{
    std::cout << "Modulo Real non testé" << std::endl;
    exit(1);
    r = fmod (a, b);
    if (r < 0) {
        if (b < 0)
            r -= b;
        else
            r += b;
    }
}

/**
 * Computes the "positive" remainder \f$r = a \bmod b\f$, i.e. such that \f$0
 * \le r < b\f$. Example:
 *
 * <center>
 *
 * <table class="LatSoft-table LatSoft-has-hlines">
 * <tr class="bt">
 *   <td class="r bl br">\f$a\f$</td>
 *   <td class="r bl br">\f$b\f$</td>
 *   <td class="c bl br">\f$r\f$</td>
 * </tr><tr class="bt">
 *   <td class="r bl br">\f$5\f$</td>
 *   <td class="r bl br">3</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="r bl br">\f$-5\f$</td>
 *   <td class="r bl br">\f$3\f$</td>
 *   <td class="c bl br">\f$1\f$</td>
 * </tr><tr>
 *   <td class="r bl br">\f$5\f$</td>
 *   <td class="r bl br">\f$-3\f$</td>
 *   <td class="c bl br">\f$2\f$</td>
 * </tr><tr>
 *   <td class="r bl br">\f$-5\f$</td>
 *   <td class="r bl br">\f$-3\f$</td>
 *   <td class="c bl br">\f$1\f$</td>
 * </tr>
 * </table>
 *
 * </center>
 */
inline void Modulo (const long & a, const long & b, long & r)
{
    r = a % b;
    if (r < 0) {
        if (b < 0)
            r -= b;
        else
            r += b;
    }
}

#ifdef WITH_NTL
/**
 * \copydoc Modulo(const long&, const long&, long&)
 */
inline void Modulo (const NTL::ZZ & a, const NTL::ZZ & b, NTL::ZZ & r)
{
    r = a % b;
    if (r < 0)
        r -= b;
}
#endif

/**
 * Computes the quotient \f$q = a/b\f$ and remainder \f$r = a
 * \bmod b\f$. Truncates \f$q\f$ to the nearest integer towards 0. One always
 * has \f$a = qb + r\f$ and \f$|r| < |b|\f$.
 */
template <typename Real>
inline void Divide (Real & q, Real & r, const Real & a, const Real & b)
{
    q = a / b;
    conv (q, trunc(q));
    r = a - q * b;
}

/**
 * Computes the quotient \f$q = a/b\f$ and remainder \f$r = a \bmod b\f$ by
 * truncating \f$q\f$ towards 0. The remainder can be negative. One always
 * has \f$a = qb + r\f$ and \f$|r| < |b|\f$. Example:
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
inline void Divide (long & q, long & r, const long & a, const long & b)
{
    ldiv_t z = ldiv (a, b);
    q = z.quot;      // q = a / b;
    r = z.rem;       // r = a % b;
}

#ifdef WITH_NTL
/**
 * \copydoc Divide(long&, long&, const long&, const long&)
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
#endif

/**
 * Integer division: \f$a = b/d\f$.
 */
inline void div (long & a, const long & b, const long & d)
{
    a = b/d;
}

/**
 * Computes \f$a/b\f$, rounds the result to the nearest integer and returns
 * the result in \f$q\f$.
 */
template <typename Real>
inline void DivideRound (const Real & a, const Real & b, Real & q)
{
    q = a/b;
    q = floor(q + 0.5);
}

/**
 * Computes \f$a/b\f$, rounds the result to the nearest integer and returns
 * the result in \f$q\f$.
 */
inline void DivideRound (const long & a, const long & b, long & q)
{
    bool neg;
    if ((a > 0 && b < 0) || (a < 0 && b > 0))
       neg = true;
    else
       neg = false;
#ifdef WITH_NTL
    long x = abs(a);
    long y = abs(b);
#else
    long x = std::abs(a);
    long y = std::abs(b);
#endif
    ldiv_t z = ldiv (x, y);
    q = z.quot;
    long r = z.rem;
    r <<= 1;
    if (r > y)
       ++q;
    if (neg)
       q = -q;
}

#ifdef WITH_NTL
/**
 * \copydoc DivideRound(const long&, const long&, long&)
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
#endif

/**
 * Returns the value of the greatest common divisor of \f$a\f$ and \f$b\f$.
 * \remark **Richard:** Il y a déjà des fonctions GCD dans NTL, pour les
 * `long` et les `ZZ` (voir fichier ZZ.h)
 */
long gcd (long a, long b);

/**
 * For given \f$a\f$ and \f$b\f$, returns the values \f$C\f$, \f$D\f$,
 * \f$E\f$, \f$F\f$ and \f$G\f$ such that:
 * \f{align*}{
 *    C a + D b
 *    &
 *   =
 *    G = \mbox{GCD } (a,b)
 *  \\
 *   E a + F b
 *    &
 *   =
 *    0.
 * \f}
 */
void Euclide (const MScal & a, const MScal & b, MScal & C, MScal & D,
              MScal & E, MScal & F, MScal & G);

/**
 * @}
 */

/**
 * \name Vectors
 *
 * @{
 */

/**
 * Allocates memory for the vector \f$A\f$ of dimensions \f$d\f$ and
 * initializes its elements to 0.
 */
template <typename Real>
inline void CreateVect (Real* & A, int d)
{
    A = new Real[d];
    for (int i = 0; i < d; i++)
        A[i] = 0;
}

/**
 * Frees the memory used by the vector \f$A\f$.
 */
template <typename Real>
inline void DeleteVect (Real* & A)
{
    delete[] A;
//    A = 0;
}

/**
 * Creates the vector \f$A\f$ of dimensions \f$d+1\f$ and initializes its
 * elements to 0.
 */
template <typename Vect>
inline void CreateVect (Vect & A, int d)
{
    A.resize (1 + d);
    // aaa clear (A);
}

/**
 * Frees the memory used by the vector \f$A\f$.
 */
template <typename Vect>
inline void DeleteVect (Vect & A)
{
    A.clear ();
}

/**
 * Sets components \f$[0..d-1]\f$ of \f$A\f$ to 0.
 */
template <typename Real>
inline void SetZero (Real* A, int d)
{
    for (int i = 0; i < d; i++)
        A[i] = 0;
}

/**
 * Sets components \f$[0..d-1]\f$ of \f$A\f$ to 0.
 */
template <typename Vect>
inline void SetZero (Vect & A, int d)
{
    for (int i = 0; i < d; i++)
        A[i] = 0;
}

/**
 * Sets all components \f$[0..d]\f$ of \f$A\f$ to the value \f$x\f$.
 */
template <typename Real>
inline void SetValue (Real* A, int d, const Real & x)
{
    for (int i = 0; i < d; i++)
        A[i] = x;
}

/**
 * Prints components \f$[c..d]\f$ of vector \f$A\f$ as a string. Components
 * are separated by string `sep`.
 */
template <typename Vect>
std::string toString (const Vect & A, int c, int d, const char* sep)
{
    std::ostringstream out;
    out << "[";
    for (int i = c; i < d; i++)
        out << A[i] << sep;
    out << A[d] << "]";
    return out.str ();
}

/**
 * Prints components \f$[c..d]\f$ of vector \f$A\f$ as a string.
 */
template <typename Vect>
inline std::string toString (const Vect & A, int c, int d)
{
    return toString (A, c, d, "  ");
}

/**
 * Prints components \f$[0..d-1]\f$ of vector \f$A\f$ as a string.
 */
template <typename Vect>
std::string toString (const Vect & A, int d)
{
   std::ostringstream ostr;
   ostr << "[";
   for (int i = 0; i < d; i++) {
      ostr << std::setprecision(2) << std::setw(3) << A[i] <<
              std::setw(2) << "  ";
   }
   ostr << "]";
   return ostr.str();
}

/**
 * Computes the scalar product of vectors \f$A\f$ and \f$B\f$, using
 * components \f$[0..n-1]\f$, and puts the result in \f$D\f$.
 */
template <typename Vect1, typename Vect2, typename Scal>
inline void ProdScal (const Vect1 & A, const Vect2 & B, int n, Scal & D)
{
    // Le produit A[i] * B[i] peut déborder, d'oÃ¹ conv.
    MScal C;   C = 0;
    for (int i = 0; i < n; i++)
        C += A[i] * B[i];
    conv (D, C);
}



/**
 * Transforms the polynomial \f$A_0 + A_1x^1 + \cdots+ A_nx^n\f$ into
 * \f$x^n - A_1x^{n-1} - \cdots- A_n\f$. The result is put in \f$B\f$.
 */
inline void Invert (const MVect & A, MVect & B, int n)
{
    conv(B[n], 1);
    for(int i = 0; i < n; i++){
       B[i] = -A[n - i];
    }
}

/**
 * Computes the `norm` norm of vector \f$V\f$, using components \f$[1..n]\f$,
 * and puts the result in \f$S\f$.
 */
template <typename Vect, typename Scal>
inline void CalcNorm (const Vect & V, int n, Scal & S, NormType norm)
{
    Scal y;
    S = 0;
    switch (norm) {
        case L1NORM:
            for (int i = 0; i < n; i++) {
                conv (y, V[i]);
#ifdef WITH_NTL
                S += abs(y);
#else
                S += std::abs(y);
#endif
            }
            break;

        case L2NORM:
            for (int i = 0; i < n; i++) {
                conv (y, V[i]);
                S += y*y;
            }
            //S = sqrt(S);
            break;

        case SUPNORM:
            for (int i = 0; i < n; i++) {
#ifdef WITH_NTL
                conv (y, abs(V[i]));
#else
                conv (y, std::abs(V[i]));
#endif
                if (y > S)
                    S = y;
            }
            break;

        case ZAREMBANORM:
            S = 1.0;
            for (int i = 0; i < n; i++) {
#ifdef WITH_NTL
                conv (y, abs(V[i]));
#else
                conv (y, std::abs(V[i]));
#endif
                if (y > 1.0)
                   S *= y;
            }
            break;

        default:
            ;
    }
}

/**
 * Copies vector \f$A\f$ into vector \f$B\f$ using components \f$[0..n-1]\f$.
 */
template <typename Vect>
inline void CopyVect (const Vect & A, Vect & B, int n)
{
    for (int k = 0; k < n; k++)  B[k] = A[k];
}

/**
 * Adds vector \f$B\f$ multiplied by \f$x\f$ to vector \f$A\f$ using
 * components \f$[0..n-1]\f$, and puts the result in \f$A\f$.
 */
template <typename Vect1, typename Vect2, typename Scal>
inline void ModifVect (Vect1 & A, const Vect2 & B, Scal x, int n)
{
    typename Vect2::value_type a;
    conv (a, x);
    for (int i = 0; i < n; i++)
        A[i] = A[i] + B[i]*a;
}

/**
 * Changes the sign of the components \f$[0..n-1]\f$ of vector \f$A\f$.
 */
template <typename Vect>
inline void ChangeSign (Vect & A, int n)
{
    for (int i = 0; i < n; i++) A[i] = -A[i];
}

/**
 * Computes the greatest common divisor of \f$V[k],…,V[n-1]\f$.
 */
inline long GCD2vect (std::vector<long> V, int k, int n)
{
    int i = k;
    long r, d, c;
    d = labs (V[k]);
    while (d != 1 && i < n) {
        c = labs (V[i]);
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
 * Allocates memory for the square matrix \f$A\f$ of dimensions
 * \f$(d+1)\times(d+1)\f$. Initializes its elements to 0.
 */
template <typename Real>
inline void CreateMatr (Real** & A, int d)
{
    A = new Real *[d];
    for (int i = 0; i < d; i++) {
        A[i] = new Real[d];
        for (int j = 0; j < d; j++)
            A[i][j] = 0; }
}

/**
 * Frees the memory used by the \f$d\times d\f$ matrix \f$A\f$.
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
 * Allocates memory for the matrix \f$A\f$ of dimensions (<tt>line</tt>)
 * \f$\times\f$ (<tt>col</tt>). Initializes its elements to 0.
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
 * Frees the memory used by the matrix \f$A\f$.
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
 * Creates the square matrix \f$A\f$ of dimensions
 * \f$d\times d\f$ and initializes its elements to 0.
 */
inline void CreateMatr (MMat& A, int d)
{
    A.resize (d, d);
    //clear (A);
}

#ifdef WITH_NTL
/**
 * \copydoc CreateMatr(MMat&, int)
 */
inline void CreateMatr (MMatP & A, int d)
{
    A.SetDims (d, d);   clear (A);
}
#endif

/**
 * Creates the matrix \f$A\f$ of dimensions (<tt>line</tt>)
 * \f$\times\f$ (<tt>col</tt>). Initializes its elements to 0.
 */
inline void CreateMatr (MMat& A, int line, int col)
{
    A.resize(line, col);
    //clear (A);
}

#ifdef WITH_NTL
/**
 * \copydoc CreateMatr(MMat&, int, int)
 */
inline void CreateMatr (MMatP & A, int line, int col)
{
    A.SetDims (line, col);   clear (A);
}
#endif

/**
 * Deletes the matrix \f$A\f$.
 */
inline void DeleteMatr (MMat& A)
{
    A.clear ();
}

#ifdef WITH_NTL
/**
 * As above.
 */
inline void DeleteMatr (MMatP & A)
{
    A.kill ();
}
#endif

/**
 * Copies matrix \f$A\f$ into matrix \f$B\f$.
 */
template <typename Matr>
inline void CopyMatr (const Matr & A, Matr & B, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B[i][j] = A[i][j];
}

/**
 * As above.
 */
template <typename Matr>
inline void CopyMatr (const Matr & A, Matr & B, int line, int col)
{
    for (int i = 0; i < line; i++)
        for (int j = 0; j < col; j++)
            B[i][j] = A[i][j];
}

/**
 * Transforms `mat` into a string. Prints the first \f$d1\f$ rows and
 * \f$d2\f$ columns. Indices start at 1. Elements with index 0 are not
 * printed.
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
 * Checks that square matrix \f$A\f$ is upper triangular (modulo \f$m\f$) for
 * dimensions 1 to `dim`.
 */
template <typename Matr>
bool CheckTriangular (const Matr & A, int dim, const MScal & m)
{
   for (int i = 1; i < dim; i++) {
      for (int j = 0; j < i; j++) {
         if (A(i,j) % m != 0) {
            std::cout << "******  CheckTriangular failed for element A[" << i <<
                 "][" << j << "]" << std::endl;
            return false;
         }
      }
    }
    return true;
}

/**
 * Performs an integer triangularization operation modulo \f$m\f$ on the
 * matrix \f$W\f$ to obtain an upper triangular matrix \f$V\f$, dual to
 * \f$W\f$. However, the matrix \f$W\f$ will be transformed too in order to
 * preserve duality. Only the first `lin` lines and the first `col` columns
 * of the matrices will be considered.
 */


template <typename Matr>
void Triangularization (Matr & W, Matr & V, int lin, int col,
                        const MScal & m)
{
   MScal T1, T2, T3, T4, T5, T6, T7, T8;

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
               Euclide (W(r,j), W(s,j), T1, T2, T3, T4, W(s,j));
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
//   CheckTriangular (V, col, m);
}

/**
 * Calculates the \f$m\f$-dual of the matrix `A`. The result is placed in the
 * matrix `B`. Only the first \f$d\f$ lines and columns are considered.
 */

template <typename Matr>
void CalcDual (const Matr & A, Matr & B, int d, const MScal & m)
{
   for (int i = 0; i < d; i++) {
      for (int j = i + 1; j < d; j++)
         clear (B(i,j));
// Dans l'original, c'est Quotient pour Lac et DivideRound pour non-Lac ??
//        Quotient(m, A(i,i), B(i,i));
      DivideRound (m, A(i,i), B(i,i));
      for (int j = i - 1; j >= 0; j--) {
         clear (B(i,j));
         for (int k = j + 1; k <= i; k++) {
            B(i,j) += A(j,k) * B(i,k);
         }
         if (B(i,j) != 0)
            B(i,j) = -B(i,j);
// Dans l'original, c'est Quotient pour Lac et DivideRound pour non-Lac ??
//         Quotient(B(i,j), A(j,j), B(i,j));
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
 * @}
 */

}     // namespace LatticeTester

#if NTL_TYPES_CODE > 1
#else

#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

template<typename ValType>
ValType det_double(const boost::numeric::ublas::matrix<ValType>& matrix)
{
    // create a working copy of the input
    boost::numeric::ublas::matrix<ValType> mLu(matrix);
    boost::numeric::ublas::permutation_matrix<std::size_t> pivots(matrix.size1());

    auto isSingular = boost::numeric::ublas::lu_factorize(mLu, pivots);
    if (isSingular)
        return static_cast<ValType>(0);

    ValType det = static_cast<ValType>(1);
    for (std::size_t i = 0; i < pivots.size(); ++i)
    {
        if (pivots(i) != i)
            det *= static_cast<ValType>(-1);

        det *= mLu(i, i);
    }

    return det;
}

#endif


#ifdef WITH_NTL

namespace NTL {

/**
 * \name Some other compatibility utilities
 *
 * @{
 */

/**
 * Converts the array of characters (string) `c` into `l`.
 */
inline void conv (long & l, const char* c) {
    l = strtol(c, (char **) NULL, 10);
}

/**
 * Converts the array of characters (string) `c` into `r`.
 */
inline void conv (double & r, const char* c) {
    r = strtod(c, (char **) NULL);
}

/**
 * @}
 */

}     // namespace NTL

#endif

#endif

