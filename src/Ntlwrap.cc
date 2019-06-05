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

  std::int64_t operator*(const Vec_64& vec1, const Vec_64& vec2) {
    if (vec1.length() != vec2.length()) {
      std::cerr << "Dimensions do not match in Vec_64 operator*(const Vec_64& v"
        "ec1, const Vec_64& vec2)\n";
      exit(1);
    }
    std::int64_t ret = 0;
    for (long i = 0; i < vec1.length(); i++) {
      ret += vec1[i] * vec2[i];
    }
    return ret;
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

  Mat_64 operator*(const Mat_64& mat1, const Mat_64& mat2){
    long row = mat2.NumRows();
    long col = mat1.NumCols();
    if (row != col) {
      std::cerr << "Dimensions do not match in Mat_64 operator*(const Mat_64& m"
        "at1, const Mat_64& mat2)\n";
      exit(1);
    }
    Mat_64 mat;
    mat.SetDims(mat1.NumRows(), mat2.NumCols());
    Mat_64 temp(NTL::transpose(mat2));
    for (long i = 0; i < mat1.NumRows(); i++) {
      for (long j = 0; j < mat2.NumCols(); j++) {
        mat[i][j] = mat1[i]*mat2[j];
      }
    }
    return mat;
  }

  double determinant(const NTL::matrix<std::int64_t>& mat) {
    NTL::matrix<NTL::ZZ> temp_mat;
    temp_mat.SetDims(mat.NumCols(), mat.NumCols());
    for (int i = 0; i < mat.NumCols(); i++) {
      for (int j = 0; j < mat.NumCols(); j++) {
        temp_mat[i][j] = mat[i][j];
      }
    }
    double temp;
    NTL::conv(temp, NTL::determinant(temp_mat));
    return temp;
  }


}
