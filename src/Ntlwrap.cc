#include "latticetester/ntlwrap.h"

namespace NTL {

  void ident(Mat_64& mat, long dim) {
    mat.SetDims(dim, dim);
    for (long i = 0; i < dim; i++) {
      for (long j = 0; j < dim; j++) {
        if (i == j) mat[i][j] = 1;
        else mat[i][j] = 0;
      }
    }
  }

  Vec_64 operator*(const Vec_64& vec, std::int64_t a) {
    Vec_64 ret_vec = vec;
    long size = vec.length();
    for (long i = 0; i < size; i++) {
      ret_vec[i] *= a;
    }
    return ret_vec;
  }

  Vec_64 operator*(std::int64_t a, const Vec_64& vec) {
    Vec_64 ret_vec = vec;
    long size = vec.length();
    for (long i = 0; i < size; i++) {
      ret_vec[i] *= a;
    }
    return ret_vec;
  }

  Vec_64& operator+=(Vec_64& vec1, const Vec_64& vec2){
    long size = vec1.length();
    for (long i = 0; i < size; i++) {
      vec1[i] += vec2[i];
    }
    return vec1;
  }

  Vec_64& operator-=(Vec_64& vec1, const Vec_64& vec2){
    long size = vec1.length();
    for (long i = 0; i < size; i++) {
      vec1[i] -= vec2[i];
    }
    return vec1;
  }

  Vec_64& operator*=(Vec_64& vec, std::int64_t a){
    long size = vec.length();
    for (long i = 0; i < size; i++) {
      vec[i] *= a;
    }
    return vec;
  }

  Mat_64& operator*=(Mat_64& mat, std::int64_t a){
    long row = mat.NumRows();
    long col = mat.NumCols();
    for (long i = 0; i < row; i++) {
      for (long j = 0; j < col; j++) {
        mat[i][j] *= a;
      }
    }
    return mat;
  }

}
