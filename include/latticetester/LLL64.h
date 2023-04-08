/**
 * This file is taken from NTL, written and maintained by Victor Shoup.
 * It was slightly modified by Pierre L'Ecuyer to make computations with
 * ordinary 64-bit integers.
 */

#ifndef LATTICETESTER_LLL64__H
#define LATTICETESTER_LLL64__H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <type_traits>

#include "NTL/tools.h"
#include <NTL/fileio.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_double.h>
#include <NTL/ZZ.h>

#include <latticetester/NTLWrap.h>

// namespace LatticeTester {
// namespace NTL {

// NTL_OPEN_NNS

// class LLL64;

// template<typename Int>
// class LLL64 {

NTL_START_IMPL

typedef NTL::matrix<int64_t> matrix64;
typedef NTL::vector<int64_t> vector64;

// int64_t LatticeSolve(vec_int64_t& x, const mat_int64_t& A, const vec_int64_t& y, int64_t reduce=0);

/**
 * On input, the rows of B contain the initial generating vectors.
 * On output, the rows of B form an LLL-reduced lattice basis.
 */
//static int64_t LLL64_FP(matrix64 &B, double delta = 0.999999);


//static int64_t BKZ64_FP(matrix64 &BB, double delta=0.999999, int64_t blockSize=10);

// };

// NTL_CLOSE_NNS


/**********************************************************************/

//static inline void MulSubFrom(int64_t & x, int64_t a, int64_t b)
//   { x -= a * b; }

static inline void CheckFinite64(double *p) {
   if (!IsFinite(p)) ResourceError("LLL64_FP: numbers too big...use LLL_XD");
}

// Returns the inner product of two arrays of double of size n.
static double InnerProductD(double *a, double *b, int64_t n) {
   double s = 0;
   int64_t i;
   for (i = 0; i < n; i++)
      s += a[i]*b[i];
   return s;
}

// Inner product of two vectors of integers a and b, returned in prod.
static void InnerProductV(int64_t &prod, const vector64& a, const vector64& b) {
   int64_t x = 0;
   int64_t n = min(a.length(), b.length());
   int64_t i;
   for (i = 0; i < n; i++) {
      x += a[i] * b[i];
   }
   prod = x;
}

#define TR_BND (NTL_FDOUBLE_PRECISION/2.0)
// Just to be safe!!

static double max_abs64(double *v, int64_t n) {
   int64_t i;
   double res, t;
   res = 0;
   for (i = 0; i < n; i++) {
      t = fabs(v[i]);
      if (t > res) res = t;
   }
   return res;
}

// Check if all the |a[i]| are less than TR_BND.
static void RowTransformStart64(double *a, int64_t *in_a, int64_t& in_float, int64_t n) {
   int64_t i;
   int64_t inf = 1;
   for (i = 0; i < n; i++) {
      in_a[i] = (a[i] < TR_BND && a[i] > -TR_BND);
      inf = inf & in_a[i];
   }
   in_float = inf;    // Returns 1 if no |a[i]| is too large, 0 otherwise.
}


static void RowTransformFinish64(vector64& A, double *a, int64_t *in_a) {
   int64_t n = A.length();
   int64_t i;
   for (i = 0; i < n; i++) {
      if (in_a[i])  {
         conv(A[i], a[i]);
      }
      else {
         conv(a[i], A[i]);
         CheckFinite64(&a[i]);
      }
   }
}


// A = A - B*MU1
static void RowTransformSub64 (vector64& A, vector64& B, int64_t& MU1) {
   register int64_t MU = MU1;
   int64_t n = A.length();
   int64_t i;

   if (MU == 1) {
      for (i = 0; i < n; i++)
         // sub(A[i], A[i], B[i]);
    	 A[i] -= B[i];
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++)
         // add(A[i], A[i], B[i]);
    	 A[i] += B[i];
      return;
   }
   if (MU == 0) return;

   for (i = 0; i < n; i++) {
	  A[i] -= MU * B[i];
      // MulSubFrom(A[i], B[i], MU);
   }
}

// A = A + B*MU
static void RowTransformAdd64 (vector64& A, vector64& B, const int64_t& MU1) {
   register int64_t T, MU = MU1;
   // int64_t k;
   int64_t n = A.length();
   int64_t i;

   if (MU == 1) {
      for (i = 0; i < n; i++)
         add(A[i], A[i], B[i]);
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++)
         sub(A[i], A[i], B[i]);
      return;
   }
   if (MU == 0) return;

   for (i = 0; i < n; i++) {
       // mul(T, B[i], MU);
       // add(A[i], A[i], T);
       T = MU * B[i];  A[i] += T;
   }
}


// A = A - B * MU1
static void RowTransform64(vector64& A, vector64& B, const int64_t& MU1,
                         double *a, double *b, long *in_a,
                         double& max_a, double max_b, int64_t& in_float) {
   register int64_t T, MU;
   // int64_t k;
   double mu;
   conv(mu, MU1);
   CheckFinite64(&mu);
   int64_t n = A.length();
   int64_t i;

   if (in_float) {
      double mu_abs = fabs(mu);
      if (mu_abs > 0 && max_b > 0 && (mu_abs >= TR_BND || max_b >= TR_BND)) {
         in_float = 0;
      }
      else {
         max_a += mu_abs*max_b;
         if (max_a >= TR_BND)
            in_float = 0;
      }
   }
   if (in_float) {
      if (mu == 1) {
         for (i = 0; i < n; i++)
            a[i] -= b[i];
         return;
      }
      if (mu == -1) {
         for (i = 0; i < n; i++)
            a[i] += b[i];
         return;
      }
      if (mu == 0) return;

      for (i = 0; i < n; i++)
         a[i] -= mu * b[i];
      return;
   }

   std::cout << "RowTransform64, not in_float! \n";

   MU = MU1;
   if (MU == 1) {
      for (i = 0; i < n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
               b[i] < TR_BND && b[i] > -TR_BND) {
            a[i] -= b[i];
         }
         else {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }
            sub(A(i), A(i), B(i));
         }
      }
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
                b[i] < TR_BND && b[i] > -TR_BND) {
            a[i] += b[i];
         }
         else {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }
            add(A(i), A(i), B(i));
         }
      }
      return;
   }
   if (MU == 0) return;

   double b_bnd = fabs(TR_BND/mu) - 1;
   if (b_bnd < 0) b_bnd = 0;

/*
   if (NTL::NumTwos(MU) >= NTL_ZZ_NBITS)
      k = MakeOdd(MU);
   else
      k = 0;
   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      if (k > 0) {
         for (i = 1; i <= n; i++) {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }

            mul(T, B(i), mu1);
            LeftShift(T, T, k);
            sub(A(i), A(i), T);
         }
      }
      else {
         for (i = 1; i <= n; i++) {
            if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
                b[i] < b_bnd && b[i] > -b_bnd) {

               a[i] -= b[i]*mu;
            }
            else {
               if (in_a[i]) {
                  conv(A(i), a[i]);
                  in_a[i] = 0;
               }
               MulSubFrom(A(i), B(i), mu1);
            }
         }
      }
   }
   else {
*/

    for (i = 0; i < n; i++) {
         if (in_a[i]) {
            conv(A(i), a[i]);
            in_a[i] = 0;
         }
         mul(T, B(i), MU);
//         if (k > 0) LeftShift(T, T, k);
         sub(A(i), A(i), T);
//      }
   }
}


static void ComputeGS(matrix64& B, double **B1, double **mu, double *b,
               double *c, int64_t k, double bound, int64_t st, double *buf) {
	// k and st are both reduced by 1 compared with NTL version
   int64_t n = B.NumCols();
   int64_t i, j;
   double s, t1, y, t;
   int64_t T1;
   int64_t test;
   double *mu_k = mu[k];

   if (st < k) {
      for (i = 0; i < st; i++)
         buf[i] = mu_k[i]*c[i];
   }
   // std::cout << "ComputeGS, st = " << st << " k = " << k << "\n";

   for (j = st; j < k; j++) {
      s = InnerProductD(B1[k], B1[j], n);
      // std::cout << "ComputeGS, j = " << j << " Inner product s = " << s << "\n";

      // test = b[k]*b[j] >= NTL_FDOUBLE_PRECISION^2
      test = (b[k]/NTL_FDOUBLE_PRECISION >= NTL_FDOUBLE_PRECISION/b[j]);

      // test = test && s^2 <= b[k]*b[j]/bound,
      // but we compute it in a strange way to avoid overflow
      if (test && (y = fabs(s)) != 0) {
         t = y/b[j];
         t1 = b[k]/y;
         // std::cout << "ComputeGS, t1 = " << t1 << "\n";
         if (t <= 1)
            test = (t*bound <= t1);
         else if (t1 >= 1)
            test = (t <= t1/bound);
         else
            test = 0;
      }

      if (test) {
         InnerProductV (T1, B[k], B[j]);
         conv(s, T1);
         // std::cout << "ComputeGS, T1 = s = " << s << "\n";
      }
      double *mu_j = mu[j];

      t1 = 0;
      for (i = 0; i < j; i++) {
         t1 += mu_j[i] * buf[i];
      }
      mu_k[j] = (buf[j] = (s - t1))/c[j];
      // std::cout << "ComputeGS, mu_k[j] = " << mu_k[j] << "\n";
   }

   // Kahan summation
   double c1;
   s = c1 = 0;
   for (j = 0; j < k; j++) {
      y = mu_k[j]*buf[j] - c1;
      t = s+y;
      c1 = t-s;
      c1 = c1-y;
      s = t;
   }
   c[k] = b[k] - s;
}


NTL_CHEAP_THREAD_LOCAL double LLLStatusInterval64 = 900.0;
NTL_CHEAP_THREAD_LOCAL char *LLLDumpFile64 = 0;

static NTL_CHEAP_THREAD_LOCAL double red_fudge64 = 0;
static NTL_CHEAP_THREAD_LOCAL int64_t log_red64 = 0;
// static NTL_CHEAP_THREAD_LOCAL int64_t verbose = 0;

static NTL_CHEAP_THREAD_LOCAL uint64_t NumSwaps64 = 0;
//static NTL_CHEAP_THREAD_LOCAL double RR_GS_time = 0;
//static NTL_CHEAP_THREAD_LOCAL double StartTime = 0;
//static NTL_CHEAP_THREAD_LOCAL double LastTime = 0;

static void init_red_fudge64() {
   int64_t i;
   log_red64 = int64_t(0.50*NTL_DOUBLE_PRECISION);
   red_fudge64 = 1;
   for (i = log_red64; i > 0; i--)
      red_fudge64 = red_fudge64*0.5;
}

static void inc_red_fudge64() {
   red_fudge64 = red_fudge64 * 2;
   log_red64--;
   cerr << "LLL64: warning--relaxing reduction (" << log_red64 << ")\n";
   if (log_red64 < 4)
      ResourceError("LLL64: too much loss of precision...stop!");
}


static int64_t ll_LLL_FP64(matrix64& B, double delta,
           double **B1, double **mu,
           double *b, double *c,
           int64_t m, int64_t init_k, int64_t &quit) {
	// m is the number of rows (generating vectors), same value as in NTL
	// init_k and k are one less compared with NTL.

   int64_t n = B.NumCols();
   int64_t i, j, k, Fc1;
   int64_t MU;
   double mu1;
   double t1;
   double *tp;
      // we tolerate a 15% loss of precision in computing
      // inner products in ComputeGS.
   static double bound = 1;
   for (i = 2*int64_t(0.15*NTL_DOUBLE_PRECISION); i > 0; i--)
         bound = bound * 2;
   double half_plus_fudge = 0.5 + red_fudge64;
   quit = 0;
   k = init_k;

   vector64 st_mem;
   st_mem.SetLength(m+2);
   int64_t *st = st_mem.elts();    // An array of integers.

   for (i = 0; i < k; i++)
      st[i] = i;
   for (i = k; i <= m; i++)
      st[i] = 0;    // With k=0, we do only this.

   UniqueArray<double> buf_store;
   buf_store.SetLength(m);
   double *buf = buf_store.get();

   vector64 in_vec_mem;
   in_vec_mem.SetLength(n);
   int64_t *in_vec = in_vec_mem.elts();

   UniqueArray<double> max_b_store;
   max_b_store.SetLength(m);
   double *max_b = max_b_store.get();

   // cerr << "LLL64: after creating UniqueArray's \n";

   for (i = 0; i < m; i++)
      max_b[i] = max_abs64(B1[i], n);

   int64_t in_float;
   int64_t rst;
   int64_t counter;
   int64_t start_over;
   int64_t trigger_index;
   int64_t small_trigger;
   int64_t cnt;
   // int64_t m_orig = m;
   int64_t rr_st = 0;   // One less than in NTL.
   int64_t max_k = 0;
   int64_t swap_cnt = 0;

   // cerr << "LLL64: before first `while` on k. \n";

   while (k < m) {
      if (k > max_k) {
         max_k = k;
         swap_cnt = 0;
      }
      if (k < rr_st) rr_st = k;  // Both are 1 less than in NTL.

      if (st[k] == k)
         rst = 0;   // 1 less than in NTL.
      else
         rst = k;
      if (st[k] < st[k+1]) st[k+1] = st[k];

      // std::cout << "LLL64: before ComputeGS, B1[k][1] = " << B1[k][1] << "\n";
      ComputeGS(B, B1, mu, b, c, k, bound, st[k], buf);
      CheckFinite64(&c[k]);
      st[k] = k;
      // std::cout << "LLL64: after ComputeGS, mu[k][0] = " << mu[k][0] << "\n";

      if (swap_cnt > 200000) {
         cerr << "LLL64: swap loop?\n";
         // In NTL, there is more stuff here...
      }

      counter = 0;
      trigger_index = k;
      small_trigger = 0;
      cnt = 0;
      int64_t thresh = 10;
      int64_t sz=0, new_sz;

      do {
         // size reduction
         counter++;
     	 // std::cout << "do loop: (k  counter) =  " <<  k << "  " << counter << "  \n";
     	 // std::cout <<  B << "  \n";
         // cerr << " vector b = " << b[0] << "  " << b[1] << "  " << b[2] << "\n";

     	 if ((counter & 127) == 0) {   // Should be 127
            new_sz = 0;
            for (j = 0; j < n; j++)
               new_sz += NumBits(B[k][j]);
            if ((counter >> 7) == 1 || new_sz < sz) {  // Should be 7
               sz = new_sz;
            }
            else {
               cerr << "LLL64: sz not smaller; warning--infinite loop? \n";
               abort();
            }
         }
         Fc1 = 0;
         start_over = 0;

         // std::cout << "rst = " << rst << "\n";
         for (j = rst-1; j >= 0; j--) { // both j and rst are 1 less than in NTL.
            t1 = fabs(mu[k][j]);
        	 // std::cout << "entered for loop: j =  " <<  j << "  \n";
        	 // std::cout << "mu[k,j] =  " <<  mu[k][j] << "  t1 =  " <<  t1 << "  \n";
            if (t1 > half_plus_fudge) {
            	 // std::cout << "we have t1 > half_plus_fudge, j =  " <<  j << "  \n";
               if (!Fc1) {
                  if (j > trigger_index ||
                        (j == trigger_index && small_trigger)) {
                     cnt++;
                     if (cnt > thresh) {
                        if (log_red64 <= 15) {
                           while (log_red64 > 10)
                              inc_red_fudge64();
                           half_plus_fudge = 0.5 + red_fudge64;
                        }
                        else {
                           inc_red_fudge64();
                           half_plus_fudge = 0.5 + red_fudge64;
                           cnt = 0;
                        }
                     }
                  }
                  trigger_index = j;
                  small_trigger = (t1 < 4);
                  Fc1 = 1;
                  if (k < rr_st) rr_st = k;
                  RowTransformStart64(B1[k], in_vec, in_float, n);
               }
               mu1 = mu[k][j];
               if (mu1 >= 0)
                  mu1 = ceil(mu1-0.5);
               else
                  mu1 = floor(mu1+0.5);

               double *mu_k = mu[k];
               double *mu_j = mu[j];

               if (mu1 == 1) {
                  for (i = 0; i < j-1; i++)
                     mu_k[i] -= mu_j[i];
               }
               else if (mu1 == -1) {
                  for (i = 0; i < j-1; i++)
                     mu_k[i] += mu_j[i];
               }
               else {
                  for (i = 0; i < j-1; i++)
                     mu_k[i] -= mu1 * mu_j[i];
               }
               mu_k[j] -= mu1;
               conv(MU, mu1);
               // std::cout << "Before row transform, mu1 = " << mu1 << " \n";

               register int64_t T, MU2 = MU;
               for (i = 0; i < n; i++) {
            	  T = MU2 * B1[j][i];
            	  // mul(T, MU2, B1[j][i]);
            	  B1[k][i] -= T;
               }
               // RowTransformSub (B[k], B[j], MU);
               //RowTransform64 (B[k], B[j], MU, B1[k], B1[j], in_vec,
               //             max_b[k], max_b[j], in_float);
               // std::cout << "After row transform, MU = " << MU << " \n";
           	   // std::cout <<  B << "  \n";
               }
         }

         if (Fc1) {
            vector64 temp = B[k];
            RowTransformFinish64(temp, B1[k], in_vec);
            B[k] = temp;
            max_b[k] = max_abs64(B1[k], n);
            b[k] = InnerProductD (B1[k], B1[k], n);
            CheckFinite64(&b[k]);
            ComputeGS (B, B1, mu, b, c, k, bound, 0, buf);
            CheckFinite64(&c[k]);
            rst = k;
         }
     	 // std::cout << "End of loop, B = " <<  B << "  \n";
      } while (Fc1 || start_over);  // End of `do` loop.

      if (b[k] == 0) {
         for (i = k; i < m-1; i++) {
            // swap i, i+1
            // swap(B[i], B[i+1]);
            B[i].swap(B[i+1]);
            tp = B1[i]; B1[i] = B1[i+1]; B1[i+1] = tp;
            t1 = b[i];  b[i] = b[i+1];   b[i+1] = t1;
            t1 = max_b[i]; max_b[i] = max_b[i+1]; max_b[i+1] = t1;
         }
         for (i = k; i <= m; i++) st[i] = 0;
         if (k < rr_st) rr_st = k;

         m--;
         if (quit) break;
         continue;
      }
      if (quit) break;

      // test LLL reduction condition
      // std::cout << "Test reduction condition  \n";
      if (k > 0 && delta*c[k-1] > c[k] + mu[k][k-1]*mu[k][k-1]*c[k-1]) {
         // swap rows k, k-1
         swap(B[k], B[k-1]);
         tp = B1[k]; B1[k] = B1[k-1]; B1[k-1] = tp;
         tp = mu[k]; mu[k] = mu[k-1]; mu[k-1] = tp;
         t1 = b[k]; b[k] = b[k-1]; b[k-1] = t1;
         t1 = max_b[k]; max_b[k] = max_b[k-1]; max_b[k-1] = t1;
         k--;
         NumSwaps64++;
         swap_cnt++;
         // cout << "-\n";
      }
      else {
         k++;
         // cout << "+\n";
      }
   }
   return m;
}

/* ************************************************ */

static int64_t LLL64_FP (matrix64& B, double delta) {
   int64_t m = B.NumRows();  // m is the number of generating vectors.
   int64_t n = B.NumCols();  // n is the dimension of the vectors.
   int64_t i, j;
   int64_t new_m, dep, quit;
   // int64_t MU;
   // int64_t T1;
   // RR_GS_time = 0;
   NumSwaps64 = 0;
   if (delta < 0.50 || delta >= 1) LogicError("LLL64: bad delta");

   init_red_fudge64();

   //  cerr << "LLL64: starting. \n";

   Unique2DArray<double> B1_store;
   B1_store.SetDims(m, n);
   double **B1 = B1_store.get();  // approximates B

   Unique2DArray<double> mu_store;
   mu_store.SetDims(m, m);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m);
   double *c = c_store.get(); // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m);
   double *b = b_store.get(); // squared lengths of basis vectors

   // cerr << "LLL64: UniqueArray's created. \n";

   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) {
         conv(B1[i][j], B[i][j]);
         CheckFinite64(&B1[i][j]);
      }
   for (i = 0; i < m; i++) {
      b[i] = InnerProductD (B1[i], B1[i], n);
      CheckFinite64(&b[i]);
   }

   // Indices in b and B1 start at 0, which is 1 less than in NTL.
   new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, m, 0, quit);
   dep = m - new_m;   // Number of rows that are dependent on others.
   m = new_m;         // Number of independent rows.

   if (dep > 0) {
      // for consistency, we move all of the zero rows to the front.
      for (i = 1; i <= m; i++) {
         swap(B[m+dep-i], B[m-i]);
      }
   }
   return m;
}


/******************************************************/

static int64_t BKZ64_FP (matrix64& BB, double delta, int64_t blockSize) {
   int64_t m = BB.NumRows();
   int64_t n = BB.NumCols();
   int64_t m_orig = m;
   int64_t i, j, ii;
   int64_t MU;

   double t1;
   // int64_t T1;
   double *tp;

   // RR_GS_time = 0;
   NumSwaps64 = 0;
   if (delta < 0.50 || delta >= 1) LogicError("BKZ_FP: bad delta");
   if (blockSize < 2) LogicError("BKZ_FP: bad block size");

   init_red_fudge64();

   matrix64 B;
   B = BB;

   B.SetDims(m+1, n);

   Unique2DArray<double> B1_store;
   B1_store.SetDimsFrom1(m+1, n);
   double **B1 = B1_store.get();  // approximates B


   Unique2DArray<double> mu_store;
   mu_store.SetDimsFrom1(m+1, m);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m+1);
   double *c = c_store.get(); // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m+1);
   double *b = b_store.get(); // squared lengths of basis vectors

   double cbar;


   UniqueArray<double> ctilda_store;
   ctilda_store.SetLength(m+1);
   double *ctilda = ctilda_store.get();


   UniqueArray<double> vvec_store;
   vvec_store.SetLength(m+1);
   double *vvec = vvec_store.get();

   UniqueArray<double> yvec_store;
   yvec_store.SetLength(m+1);
   double *yvec = yvec_store.get();

   UniqueArray<double> uvec_store;
   uvec_store.SetLength(m+1);
   double *uvec = uvec_store.get();

   UniqueArray<double> utildavec_store;
   utildavec_store.SetLength(m+1);
   double *utildavec = utildavec_store.get();

   UniqueArray<int64_t> Deltavec_store;
   Deltavec_store.SetLength(m+1);
   int64_t *Deltavec = Deltavec_store.get();

   UniqueArray<int64_t> deltavec_store;
   deltavec_store.SetLength(m+1);
   int64_t *deltavec = deltavec_store.get();;

   int64_t quit;
   int64_t new_m;
   int64_t z, jj, kk;
   int64_t s, t;
   int64_t h;
   double eta;

   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) {
         conv(B1[i][j], B[i][j]);
         CheckFinite64(&B1[i][j]);
      }

   for (i = 0; i < m; i++) {
      b[i] = InnerProductD(B1[i], B1[i], n);
      CheckFinite64(&b[i]);
   }
   m = ll_LLL_FP64(B, delta, B1, mu, b, c, m, 0, quit);

   // double tt;
   // double enum_time = 0;
   uint64_t NumIterations = 0;
   uint64_t NumTrivial = 0;
   uint64_t NumNonTrivial = 0;
   uint64_t NumNoOps = 0;
   int64_t clean = 1;

   if (m < m_orig) {
      for (i = m_orig; i >= m+1; i--) {
         // swap i, i-1
         swap(B[i], B[i-1]);
      }
   }

   if (!quit && m > 1) {
      if (blockSize > m) blockSize = m;
      z = 0;
      jj = 0;         // jj is the same as for NTL ?

      while (z < m-1) {
         jj++;
         kk = min(jj+blockSize-2, m-1);      // kk is 1 less than for NTL

         if (jj == m) {
            jj = 1;
            kk = blockSize-1;
            clean = 1;
         }

         cbar = c[jj];
         utildavec[jj] = uvec[jj] = 1;

         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;

         s = t = jj;
         deltavec[jj] = 1;

         for (i = jj; i <= kk+1; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }

         while (t <= kk) {
            ctilda[t] = ctilda[t+1] +
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

            ForceToMem(&ctilda[t]);  // prevents an infinite loop

            eta = 0;
            if (ctilda[t] < cbar - eta) {
               if (t > jj) {
                  t--;
                  t1 = 0;
                  for (i = t+1; i <= s; i++)
                     t1 += utildavec[i]*mu[i][t];
                  yvec[t] = t1;
                  t1 = -t1;
                  if (t1 >= 0)
                     t1 = ceil(t1-0.5);
                  else
                     t1 = floor(t1+0.5);
                  utildavec[t] = vvec[t] = t1;
                  Deltavec[t] = 0;
                  if (utildavec[t] > -yvec[t])
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
               }
               else {
                  cbar = ctilda[jj];
                  for (i = jj; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec[t] = -Deltavec[t];
               if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         }
         NumIterations++;

         h = min(kk, m-1);   // h is 1 less than for NTL

         if ((delta - 8*red_fudge64)*c[jj] > cbar) {

            clean = 0;
            // we treat the case that the new vector is b_s (jj < s <= kk)
            // as a special case that appears to occur most of the time.
            s = 0;
            for (i = jj+1; i <= kk; i++) {
               if (uvec[i] != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
            if (s == 0) LogicError("BKZ_FP: internal error");

            if (s > 0) {
               // special case
               NumTrivial++;
               for (i = s-1; i >= jj; i--) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
               // cerr << "special case\n";
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h, jj, quit);
               if (new_m != h) LogicError("BKZ_FP: internal error");
               if (quit) break;
            }
            else {
               // the general case
               register int64_t T;
               NumNonTrivial++;
               for (i = 0; i < n; i++) conv(B[m][i], 0);
               for (i = jj; i <= kk; i++) {
                  if (uvec[i] == 0) continue;
                  conv(MU, uvec[i]);
                  // RowTransformAdd (B[m], B[i], MU);
                  for (ii = 0; ii < n; ii++) {
                      T = MU * B[i][ii];  B[m][ii] += T;
                  }
               }

               for (i = m; i >= jj+1; i--) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }

               for (i = 0; i < n; i++) {
                  conv(B1[jj][i], B[jj-1][i]);
                  CheckFinite64(&B1[jj][i]);
               }

               b[jj] = InnerProductD(B1[jj], B1[jj], n);
               CheckFinite64(&b[jj]);

               if (b[jj] == 0) LogicError("BKZ_FP: internal error, b[jj]==0");

               // remove linear dependencies

               // cerr << "general case\n";
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, kk+1, jj, quit);

               if (new_m != kk) LogicError("BKZ_FP: internal error, new_m != kk");

               // remove zero vector

               for (i = kk+1; i <= m; i++) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
               if (h > kk) {
                  // extend reduced basis
                  new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h+1, h, quit);

                  if (new_m != h+1) LogicError("BKZ_FP: internal error, new_m != h");
                  if (quit) break;
               }
            }
            z = 0;
         }
         else {
            // LLL_FP
            // cerr << "progress\n";
            NumNoOps++;
            if (!clean) {
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h+1, h, quit);
               if (new_m != h+1) LogicError("BKZ_FP: internal error, new_m != h");
               if (quit) break;
            }
            z++;
         }
      }
   }

   // clean up
   if (m_orig > m) {
      // for consistency, we move zero vectors to the front
      for (i = m; i < m_orig; i++) {
         swap(B[i], B[i+1]);
      }
      for (i = 1; i <= m; i++) {
         swap(B[m_orig-i], B[m-i]);
      }
   }
   B.SetDims(m_orig, n);
   BB = B;
   return m;
}

NTL_END_IMPL

//} // end namespace

#endif
