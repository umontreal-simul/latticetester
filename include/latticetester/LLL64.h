/**
 * This file is taken from NTL, written and maintained by Victor Shoup.
 * It was slightly modified by Pierre L'Ecuyer to make computations with
 * ordinary 64-bit integers.
 */

#ifndef NTL_LLL64__H
#define NTL_LLL64__H

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

//#include "NTL/tools.h"
//#include "NTL/ZZ.h"
//#include "NTL/RR.h"
//#include "NTL/vec_ZZ.h"
//#include "NTL/mat_ZZ.h"
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include "latticetester/NTLWrap.h"

namespace LatticeTester {
// namespace NTL {

typedef NTL::matrix<int64_t> mat_long;
typedef NTL::vector<int64_t> vec_long;
// typedef Mat<std::int64_t> mat_long;
// typedef Vec<std::int64_t> vec_long;

// NTL_OPEN_NNS

class LLL64;

// template<typename Int>
class LLL64 {

public:

// long LatticeSolve(vec_long& x, const mat_long& A, const vec_long& y, long reduce=0);

long LLL64_FP(NTL::matrix<int64_t>& B, double delta = 0.999999);

long BKZ64_FP(NTL::matrix<int64_t>& BB, double delta=0.999999, long BlockSize=10, long prune=0);

};

// NTL_CLOSE_NNS

} // end namespace

#endif
