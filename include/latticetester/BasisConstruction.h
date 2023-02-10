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

#ifndef LATTICETESTER_BASISCONSTRUCTION_H
#define LATTICETESTER_BASISCONSTRUCTION_H

#include "NTL/LLL.h"
#include <NTL/mat_GF2.h>
#include "latticetester/IntLatticeBase.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/Const.h"
#include "latticetester/NTLWrap.h"
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





namespace LatticeTester {

template<typename IntMat>
struct LLLConstr {
	void LLLConstruction(IntMat &matrix);
  void LLLConstruction(IntMat &matrix, double delta);
};

/**
 * This class offers methods to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 * The implementation relies on NTL and uses NTL matrices.
 *
 * NTL already offers a very efficient method to construct an LLL-reduced basis from a set
 * of generating vectors.  This is the most effective way of constructing a basis
 * and it is encapsulated in the `LLLConstruction` method given below.
 * We also offer an alternative that constructs a triangular basis, in `GCDTriangularBasis`.
 * To compute the $m$-dual of a given basis, we have a general method implemented in
 * `mDualComputation`, and a faster method in `mDualTriangular` that works only when
 * the basis is upper-triangular.
 * The methods `Util::Triangularization` and `Util::calcDual` do essentially the same
 * things; however, the methods given here perform more verifications.
 *  ***  We should compare the speeds
 *
 * A few tips about the usage of this class:
 * - Prefer the usage of NTL types when using this module. The methods here do not
 *   have any kind of overflow detection.
 * - Reduce the basis before doing a triangularization. Reducing a basis with
 *   LLL is much faster than the GCDTriangularBasis and seems to make this operation
 *   easier to perform.    ***  To be tested again.
 * - Use specialized methods. With a more in depth knowledge of your problem, it
 *   is possible that there are much more efficient ways to build a basis and its
 *   dual (and/or those matrices may already be triangular).
 */

template<typename Int> class BasisConstruction {

private:
	typedef NTL::vector<Int>  IntVec;
	typedef NTL::matrix<Int>  IntMat;
 	struct LLLConstr<IntMat> spec;   // Why this?

public:

	/**
	 * This functions takes a set of generating vectors of a vector space and
	 * finds a basis for this space by applying LLL reduction, using the NTL implementation.
	 * This is much faster than applying GCDConstruction, but it does not provide a triangular basis.
	 * It is a good idea to use it before computing a triangular basis.
	 */
	void LLLConstruction(IntMat &matrix);


 
	/**
	 * This functions takes a set of generating vectors of a vector space and
	 * finds a basis for this space by applying LLL reduction, using the NTL implementation.
	 * This is much faster than applying GCDConstruction, but it does not provide a triangular basis.
	 * It is a good idea to use it before computing a triangular basis.
	 */
	void LLLConstruction(IntMat &matrix, double delta);




	/**
	 * This function does essentially the same thing as `Util::Triangularization`.
	 * It uses a form of Gaussian elimination to obtain an upper triangular basis
	 * for the smallest lattice that contains the generating vectors which are the
	 * rows of the given matrix `matrix`.
	 * In each column, it applies Euclid's algorithm to the elements under the
	 * diagonal and change the corresponding rows to set all these elements to
	 * zero except the one in the diagonal. The allowed row operations are only
	 *   - Multiply a row by \f$-1\f$,
	 *   - Add an integer multiple of row \f$i\f$ to row \f$j\f$ for \f$i \neq j\f$,
	 *   - Swap row \f$i\f$ with row \f$j\f$.
	 * All these operations can be performed modulo the scaling factor `m`.
	 * After constructing this basis, the algorithm eliminates negative
	 * coefficients in the matrix.
	 * Warning: In this implementation, the numbers below the diagonal can grow
	 * very large, so the method may require a lot of memory.
	 *  ***  But everything should be done modulo m ???
	 */
  void GCDTriangularBasis(IntMat &matrix);




  	/**
	 * This function does essentially the same thing as `Util::Triangularization`.
	 * It uses a form of Gaussian elimination to obtain an upper triangular basis
	 * for the smallest lattice that contains the generating vectors which are the
	 * rows of the given matrix `matrix`.
	 * In each column, it applies Euclid's algorithm to the elements under the
	 * diagonal and change the corresponding rows to set all these elements to
	 * zero except the one in the diagonal. The allowed row operations are only
	 *   - Multiply a row by \f$-1\f$,
	 *   - Add an integer multiple of row \f$i\f$ to row \f$j\f$ for \f$i \neq j\f$,
	 *   - Swap row \f$i\f$ with row \f$j\f$.
	 * All these operations can be performed modulo the scaling factor `m`.
	 * After constructing this basis, the algorithm eliminates negative
	 * coefficients in the matrix.
   * All computational operation is done modulo 'm'. 
	 */
  	void GCDTriangularBasis(IntMat &matrix, Int m);



	/**
	 * This does the same thing as mDualTriangular(), but is much slower. This
	 * is here simply for the sake of comparison and should not be used in practice.
	 * Suppose the basis matrix contains basis vectors on its lines and is
	 * \f$p\times q\f$ where \f$q \geq p\f$. We can compute the \f$m\f$-dual
	 * as follows.
	 * Let's note the basis matrix `V` and the dual matrix `W` and have lines
	 * of `W` also contain dual basis vectors in its lines. We know that
	 * \f$VW^t = mI_{p\times p}\f$ where \f$I\f$ is the identity matrix. Now, in
	 * the case of a \f$m\f$-dual computation, we can assume that all the
	 * arithmetic is done modulo \f$m\f$.

	 void DualSlow(IntMat& matrix, IntMat& dualMatrix, Int& m);
	 */

	/**
	 * This function assumes that `matrix` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular.
	 * It computes and returns the `m`-dual basis in `dualMatrix`.
	 * ***  TO BE IMPLEMENTED **
	 */
	//void mDualComputation(IntMat &matrix, IntMat &dualMatrix, Int m);

	/**
	 * This method constructs a basis for the projection `proj` of the basis `in`.
	 * using the `LLLConstruction` method, and puts it in `out`. This will
	 * overwrite the lattice basis in `out`, changing also the dimension.
	 * The returned basis is not triangular in general.
	 */
	//template<typename Int, typename Real, typename RealRed>
	//template<typename Int>
  template<typename Real, typename RealRed>
	void ProjectionConstruction(IntLatticeBase<Int, Real, RealRed> &in,
			IntLatticeBase<Int, Real, RealRed> &out, const Coordinates &proj);




    /**
   * Takes a set of generating vectors in the matrix `mat` and iteratively
   * transforms it into a lower triangular lattice basis into the matrix `mat2`.
   * `mat` and `mat2` have to have the same number of rows and the same number of columns.
   *  All the computations will be done modulo `mod`, which means that you
   * must know the rescaling factor for the vector system to call this function.
   * After the execution, `mat` will be a matrix containing irrelevant information
   * and `mat2` will contain an upper triangular basis.
   *
   * For more details please look at \cite latTesterGide. This algorithm basically
   * implements what is written in this guide. The matrix
   * `mat` contains the set of vectors that is used and modified at each step to
   * get a new vector from the basis.
   */
  // template<typename IntMat,typename IntVec> 
   void lowerTriangular(IntMat &mat, IntMat &mat2, Int &mod);	



    /**
   * Takes a set of generating vectors in the matrix `mat` and iteratively
   * transforms it into an upper triangular lattice basis into the matrix `mat2`.
   * `mat` and `mat2` have to have the same number of rows and the same number of columns.
   *  All the computations will be done modulo `mod`, which means that you
   * must know the rescaling factor for the vector system to call this function.
   * After the execution, `mat` will be a matrix containing irrelevant information
   * and `mat2` will contain an upper triangular basis.
   *
   * For more details please look at \cite latTesterGide. This algorithm basically
   * implements what is written in this guide. The matrix
   * `mat` contains the set of vectors that is used and modified at each step to
   * get a new vector from the basis.
   */

  //template<typename IntMat,typename IntVec> 
   void upperTriangular(IntMat &mat, IntMat &mat2, Int &mod);		

  /**
   * Takes an upper triangular basis `A` and computes an m-dual lattice basis
   * to this matrix. For this algorithm to work, `A` has to be upper
   * triangular and all the coefficients on the diagonal have to divide `m`.
   *
   * For `B` to be `m`-dual to `A`, we have to have that \f$AB^t = mI\f$. It
   * is quite easy to show that, knowing `A` is upper triangular, `B` will be a
   * lower triangular matrix with `A(i,i)*B(i,i) = m` for all `i` and
   * \f$ A_i \cdot B_j = 0\f$ for \f$i\neq j\f$. To get the second condition,
   * we simply have to recursively take for each line
   * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
   */
   //template <typename IntMat>
   void calcDualUpperTriangular (const IntMat & A, IntMat & B, int d, const Int & m);
     

  	/**
	 * This function assumes that `matrix` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular.
   * Takes a basis `A` and computes an m-dual lattice basis and put it to 'B'.
   * The matrix B is the m-dual basis of A.
   */
  // template <typename IntMat>
    void calcDual (const IntMat & A, IntMat & B, const Int & m);
        
    void calcDual(const NTL::Mat<NTL::ZZ>  & A, NTL::Mat<NTL::ZZ>  & B, const NTL::ZZ & m);   
   /**
	 * This function does essentially the same thing as `CalcDualUpperTriangularBasis`.
	 * It assumes that `matrix` contains a triangular basis of the primal lattice
	 * scaled by the factor `m`.  It computes and returns the `m`-dual basis
	 * in `dualMatrix`.  This function uses the method described in \cite rCOU96a.
	 * Since `A=matrix` is upper triangular, `B=mdualMatrix` will be lower triangular
	 * with `A(i,i)*B(i,i) = m` for all `i` and \f$ A_i \cdot B_j = 0\f$ for
	 * \f$i\neq j\f$. To get the second condition, we simply have to
	 * recursively take for each line
	 * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
	 * This is much faster than a traditional solving of a linear system.
	 */
	void mDualTriangular(IntMat &matrix, IntMat &dualMatrix, Int m);

};

//============================================================================
// Implementation

template<>
void LLLConstr<NTL::matrix<std::int64_t>>::LLLConstruction(
		NTL::matrix<std::int64_t> &matrix);

template<>
void LLLConstr<NTL::matrix<std::int64_t>>::LLLConstruction(
		NTL::matrix<std::int64_t> &matrix, double delta);   

template<>
void LLLConstr<NTL::matrix<NTL::ZZ>>::LLLConstruction(
		NTL::matrix<NTL::ZZ> &matrix);

template<>
void LLLConstr<NTL::matrix<NTL::ZZ>>::LLLConstruction(
		NTL::matrix<NTL::ZZ> &matrix, double delta);    

template<typename Int>
void BasisConstruction<Int>::LLLConstruction(IntMat &matrix) {
	spec.LLLConstruction(matrix);
}


template<typename Int>
void BasisConstruction<Int>::LLLConstruction(IntMat &matrix, double delta) {
	spec.LLLConstruction(matrix,delta);

}



template<typename Int>
void BasisConstruction<Int>::GCDTriangularBasis(IntMat &matrix, Int mod) {
	// It is important to note that the lines of matrix are the basis vectors
	long rows = matrix.NumRows();
	long cols = matrix.NumCols();
	long max_rank = rows < cols ? rows : cols;
	long rank = 0;
	// The basis will have at most max_rank vectors.
	Int q;
	//Int r;
	for (long i = 0; i < max_rank; i++) {
		// We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
		// algorithm and applying transformations to the matrix
		for (long j = i + 1; j < rows; j++) {
			while (matrix[j][i] != 0) {
				matrix[i].swap(matrix[j]);
				q = matrix[j][i] / matrix[i][i];
				matrix[j] -= q * matrix[i];
				for(int k=0;k<max_rank;k++)
				  Modulo(matrix[j][k], mod, matrix[j][k]);
			}
		}
		// We make sure that the coefficients are positive.
		// This is because the algorithms we use work for positive vectors.
		// if (matrix[i][i] < 0) matrix[i] *= -1;
		// for (long j = i-1; j >= 0; j--) {
		//   if (matrix[j][i] < 0) {
		//     if (-matrix[j][i] > matrix[i][i]){
		//       q = -matrix[j][i]/matrix[i][i] + 1;
		//       matrix[j] += q * matrix[i];
		//     } else {
		//       matrix[j] += matrix[i];
		//     }
		//   }
		// }
		if (matrix[i][i] != 0) {
			rank++;
			if (matrix[i][i] < 0)
				matrix[i] *= Int(-1);
		}
	}
	// We remove zero vectors from the basis.
	matrix.SetDims(rank, cols);
}


template<typename Int>
void BasisConstruction<Int>::GCDTriangularBasis(IntMat &matrix) {
	// It is important to note that the lines of matrix are the basis vectors
	long rows = matrix.NumRows();
	long cols = matrix.NumCols();
	long max_rank = rows < cols ? rows : cols;
	long rank = 0;
	// The basis will have at most max_rank vectors.
	Int q;
	//Int r;
	for (long i = 0; i < max_rank; i++) {
		// We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
		// algorithm and applying transformations to the matrix
		for (long j = i + 1; j < rows; j++) {
			while (matrix[j][i] != 0) {
				matrix[i].swap(matrix[j]);
				q = matrix[j][i] / matrix[i][i];
				matrix[j] -= q * matrix[i];
			 
			}
		}
		// We make sure that the coefficients are positive.
		// This is because the algorithms we use work for positive vectors.
		// if (matrix[i][i] < 0) matrix[i] *= -1;
		// for (long j = i-1; j >= 0; j--) {
		//   if (matrix[j][i] < 0) {
		//     if (-matrix[j][i] > matrix[i][i]){
		//       q = -matrix[j][i]/matrix[i][i] + 1;
		//       matrix[j] += q * matrix[i];
		//     } else {
		//       matrix[j] += matrix[i];
		//     }
		//   }
		// }
		if (matrix[i][i] != 0) {
			rank++;
			if (matrix[i][i] < 0)
				matrix[i] *= Int(-1);
		}
	}
	// We remove zero vectors from the basis.
	matrix.SetDims(rank, cols);
}

//============================================================================

/**
 * This algorithm calculates the dual as well as the `m` used for rescaling.
 * It checks if this `m` divides the modulo given to the algorithm.
 * Right now, this assumes the basis is triangular, might need to change it.

 template<typename Int>
 void BasisConstruction<Int>::DualSlow(IntMat& matrix,
 IntMat& dualMatrix, Int& modulo)
 {
 // We need to have a triangular basis matrix
 if (! CheckTriangular(matrix, matrix.NumRows(), modulo))
 GCDConstruction(matrix);
 long dim = matrix.NumRows();
 if (dim != matrix.NumCols()) {
 std::cout << "matrix has to be square, but dimensions do not fit.\n";
 return;
 }
 Int m(1);
 NTL::ident(dualMatrix, dim);
 Int gcd;
 for (long i = dim-1; i>=0; i--) {
 for (long j = i+1; j < dim; j++) {
 dualMatrix[i] -= matrix[i][j] * dualMatrix[j];
 matrix[i] -= matrix[i][j] * matrix[j];
 }
 gcd = matrix[i][i];
 for (long j = i; j < dim; j++) {
 gcd = NTL::GCD(gcd, dualMatrix[i][j]);
 }
 gcd = matrix[i][i] / gcd;
 if (gcd != 1) {
 dualMatrix *= gcd;
 m *= gcd;
 }
 for (long j = i; j < dim; j++) {
 dualMatrix[i][j] /= matrix[i][i];
 }
 matrix[i][i] = 1;
 }
 }
 */

template<typename Int>
void BasisConstruction<Int>::mDualTriangular(IntMat &matrix, IntMat &dualMatrix,
		Int m) {
	// We need to have a triangular basis matrix
	if (!CheckTriangular(matrix, matrix.NumRows(), Int(0)))
		GCDTriangularBasis(matrix,m);
	long dim = matrix.NumRows();
	if (dim != matrix.NumCols()) {
		std::cout << "matrix has to be square, but dimensions do not fit.\n";
		return;
	}
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);
		return;
	}
	dualMatrix.SetDims(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = i + 1; j < dim; j++)
			NTL::clear(dualMatrix(i, j));
		if (!NTL::IsZero(matrix(i, i))) {
			Int gcd = NTL::GCD(m, matrix(i, i));
			m *= matrix(i, i) / gcd;
			dualMatrix *= matrix(i, i) / gcd;
		}

		DivideRound(m, matrix(i, i), dualMatrix(i, i));
		for (int j = i - 1; j >= 0; j--) {
			NTL::clear(dualMatrix(i, j));
			for (int k = j + 1; k <= i; k++)
				dualMatrix(i, j) += matrix(j, k) * dualMatrix(i, k);
			if (dualMatrix(i, j) != 0)
				dualMatrix(i, j) = -dualMatrix(i, j);
			if (!NTL::IsZero(dualMatrix(i, j) % matrix(j, j))) {
				Int gcd = NTL::GCD(dualMatrix(i, j), matrix(j, j));
				m *= matrix(j, j) / gcd;
				dualMatrix *= matrix(j, j) / gcd;
			}
			DivideRound(dualMatrix(i, j), matrix(j, j), dualMatrix(i, j));
		}
	}
}

//============================================================================

/*template<typename Int>
void BasisConstruction<Int>::mDualComputation(IntMat &matrix,
		IntMat &dualMatrix, Int m) {
	// **  TO BE IMPLEMENTED  **
}*/

//============================================================================

//template<typename Int>s
//template<typename Int, typename Real, typename RealRed>
template<typename Int>
template<typename Real, typename RealRed>
void BasisConstruction<Int>::ProjectionConstruction(
		IntLatticeBase<Int, Real, RealRed>& in,
		IntLatticeBase<Int, Real, RealRed> &out, const Coordinates& proj) {
	std::size_t dim = proj.size();
	unsigned int lat_dim = in.getDim();
	if (dim > lat_dim)
		MyExit(1, "Coordinates do not match the dimensions of `in`.");
	IntMat new_basis, tmp(NTL::transpose(in.getBasis()));
	new_basis.SetDims(dim, tmp.NumRows());
	tmp = NTL::transpose(tmp);
	auto it = proj.cbegin();
	for (std::size_t i = 0; i < dim; i++) {
		if (*it <= lat_dim)
			new_basis[i] = tmp[*it];
		else
			MyExit(1, "Coordinates do not match the dimensions of `in`.");
		it++;
	}
	new_basis = NTL::transpose(new_basis);
	LLLConstruction(new_basis);
	out = IntLatticeBase<Int, Real, RealRed>(new_basis, dim, in.getNormType());
}



 /**
   * Takes a set of generating vectors in the matrix `mat` and iteratively
   * transforms it into an upper triangular lattice basis into the matrix `mat2`.
   * `mat` and `mat2` have to have the same number of rows and the same number of columns.
   *  All the computations will be done modulo `mod`, which means that you
   * must know the rescaling factor for the vector system to call this function.
   * After the execution, `mat` will be a matrix containing irrelevant information
   * and `mat2` will contain an upper triangular basis.
   *
   * For more details please look at \cite latTesterGide. This algorithm basically
   * implements what is written in this guide. The matrix
   * `mat` contains the set of vectors that is used and modified at each step to
   * get a new vector from the basis.
   */

   //template<typename IntMat,typename IntVec, typename Int> 
   template<typename Int> 
   void BasisConstruction<Int>::upperTriangular(IntMat &mat, IntMat &mat2, Int &mod){
     IntVec coeff, vl,v2; 
     Int C, D, val, gcd;  
     int pc, pl, k;
     int dim1=mat.NumRows();
     int dim2=mat.NumCols();

     pl=0;
     pc=0;
     while(pl<dim1 && pc<dim2){
           for(int i=0;i<dim1;i++)
             Modulo (mat(i,pc), mod, mat(i,pc));
                
            coeff.SetLength(dim2);
            k=0;     
            while( k<dim1 && mat(k,pc)==0)
             { coeff[k]=0; 
               k++;
             }
                
           if(k<dim1)
            { gcd=mat(k,pc);
              coeff[k]=1;
              val=gcd;
             for(int i=k+1;i<dim1; i++){
               if(mat(i,pc)==0)
               { coeff[i]= 0;
                 continue;
                }           
              Euclide (val, mat(i,pc), C, D , gcd);
              coeff[i]= D;
              for(int j=0;j<i;j++) 
                  coeff[j]*=C;
              val=gcd; 
              }       
              int coeffN[dim2];
              int nb=0;
              for(int a=0;a<dim1;a++) 
              { if(coeff[a]!=0)
                 { coeffN[nb]=a;
                   nb++;
                  }
               }    
            vl.SetLength(dim2);
            int ind=0;
            for(int j=0;j<dim2;j++) {
              for(int i=0;i<nb;i++)
              { ind=coeffN[i];
                 vl[j]=vl[j]+coeff[ind]*mat(ind,j);      
              } 
              Modulo (vl[j], mod, vl[j]);  
             }
             for(int i=0;i<dim1;i++)
             {  if(mat(i,pc)!=0){
                v2= (mat(i,pc)/gcd)*vl;
                for(int j=pc;j<dim2;j++)
                    Modulo (v2[j], mod, v2[j]);
                for(int j=pc;j<dim2;j++)
                 {   
                   mat(i,j)=mat(i,j)-v2[j];  
                   Modulo (mat(i,j), mod, mat(i,j));
                 } 
                }    
             }
             mat2[pl]=vl; 
          }
          else
          {  for (int j1 = 0; j1 < dim2; j1++) {
             if (j1 != pl)
               NTL::clear (mat2(pl,j1));
             else
               mat2(pl,j1) = mod;
             }   
          }
          coeff.clear();
          vl.clear();
          pl++; 
          pc++;
       }
    }


    /**
   * Takes a set of generating vectors in the matrix `mat` and iteratively
   * transforms it into a lower triangular lattice basis into the matrix `mat2`.
   * `mat` and `mat2` have to have the same number of rows and the same number of columns.
   *  All the computations will be done modulo `mod`, which means that you
   * must know the rescaling factor for the vector system to call this function.
   * After the execution, `mat` will be a matrix containing irrelevant information
   * and `mat2` will contain an upper triangular basis.
   *
   * For more details please look at \cite latTesterGide. This algorithm basically
   * implements what is written in this guide. The matrix
   * `mat` contains the set of vectors that is used and modified at each step to
   * get a new vector from the basis.
   */
  // template<typename IntMat,typename IntVec,typename Int > 
   template<typename Int> 
   void  BasisConstruction<Int>::lowerTriangular(IntMat &mat, IntMat &mat2, Int &mod){
     IntVec coeff, vl,v2; 
     Int C, D, val, gcd;  
     int pc, pl, k;
     int dim1=mat.NumRows();
     int dim2=mat.NumCols();   
     pl=dim1-1;
     pc=dim2-1;
     while(pl>=0 && pc>=0){
           for(int i=0;i<dim1;i++)
              Modulo (mat(i,pc), mod, mat(i,pc));  
            coeff.SetLength(dim2);
            k=0;     
            while( k<dim1 && mat(k,pc)==0)
             { coeff[k]=0; 
               k++;
             }     
           if(k<dim1)
            { gcd=mat(k,pc);
              coeff[k]=1;
              val=gcd;      
             for(int i=k+1;i<dim1; i++){
               if(mat(i,pc)==0)
               { coeff[i]= 0;
                 continue;
                }   
              Euclide (val, mat(i,pc), C, D , gcd);
              coeff[i]= D;
              for(int j=0;j<i;j++) 
                  coeff[j]*=C;
              val=gcd; 
              }      
              int coeffN[dim2];
              int nb=0;
              for(int a=0;a<dim1;a++) 
              { if(coeff[a]!=0)
                 { coeffN[nb]=a;
                   nb++;
                  }
               }        
            vl.SetLength(dim2);
            int ind=0;
            for(int j=0;j<dim2;j++) {
              for(int i=0;i<nb;i++)
              { ind=coeffN[i];
                 vl[j]=vl[j]+coeff[ind]*mat(ind,j);       
              } 
              Modulo (vl[j], mod, vl[j]);  
             }
             for(int i=0;i<dim1;i++)
             {  if(mat(i,pc)!=0){
                v2= (mat(i,pc)/gcd)*vl;
                for(int j=0;j<dim2;j++)
                    Modulo (v2[j], mod, v2[j]);
                for(int j=0;j<dim2;j++)
                 {   
                   mat(i,j)=mat(i,j)-v2[j];  
                   Modulo (mat(i,j), mod, mat(i,j));
                 } 
                }    
             }
             mat2[pl]=vl; 
          }
          else
          {  for (int j1 = 0; j1 < dim2; j1++) {
             if (j1 != pl)
               NTL::clear (mat2(pl,j1));
             else
               mat2(pl,j1) = mod;
             }   
          }
          coeff.clear();
          vl.clear();
          pl--; 
          pc--;
       }
   }


    /**
   * Takes a basis `A` and computes an m-dual lattice basis B.
   * The matrix B is the m-dual basis of A.
   */
   
   // template <typename Int>
  //  template <typename Int>
    template <typename Matr,typename Int >
    void calcDual (Matr & A, Matr & B,  Int & m) {
      Int  d;
      Matr C;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      C.SetDims(dim1, dim2);
      inv(d,B,A);
      transpose(C,B);
      for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++)
           B(i,j)= (m*C(i,j))/d;     
        }
     }
           
   // template <typename Int>
  /**
   template<typename Int> 
   void BasisConstruction<Int>::calcDual(const NTL::Mat<NTL::GF2>  & A, NTL::Mat<NTL::GF2>  & B, const NTL::GF2 & m) {
      NTL::GF2 d;
    //  Int d;// mult;
    //  IntMat C;
      NTL::Mat<NTL::GF2> C;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      C.SetDims(dim1, dim2);
      inv(d,B,A);
      transpose(C,B);
      for (int i = 1; i < dim1; i++) {
        for (int j = 1; j < dim2; j++)
           B(i,j)= (m*C(i,j))/d;
          
        }
     }
   */

   template<typename Int> 
   void BasisConstruction<Int>::calcDual(const NTL::Mat<NTL::ZZ>  & A, NTL::Mat<NTL::ZZ>  & B, const NTL::ZZ & m) {
      NTL::ZZ d;
    //  Int d;// mult;
    //  IntMat C;
      NTL::Mat<NTL::ZZ> C;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      C.SetDims(dim1, dim2);
      inv(d,B,A);
      transpose(C,B);
      for (int i = 1; i < dim1; i++) {
        for (int j = 1; j < dim2; j++)
           B(i,j)= (m*C(i,j))/d;
          
        }
     }
   


   //  void BasisConstruction<std::int64_t>::
 /**   void calcDual (const NTL::matrix<std::int64_t>  & A, NTL::matrix<std::int64_t>  & B, const std::int64_t & m) {
      std::int64_t d;
     // Int d;// mult;
     // IntMat C;
      NTL::matrix<std::int64_t>  C;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      C.SetDims(dim1, dim2);
      inv(d,B,A);
      transpose(C,B);
      for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++)
           B(i,j)= (m*C(i,j))/d;
          
        }
     }
  **/



  /**
   * Takes an upper triangular basis `A` and computes an m-dual lattice basis
   * to this matrix. For this algorithm to work, `A` has to be upper
   * triangular and all the coefficients on the diagonal have to divide `m`.
   *
   * For `B` to be `m`-dual to `A`, we have to have that \f$AB^t = mI\f$. It
   * is quite easy to show that, knowing `A` is upper triangular, `B` will be a
   * lower triangular matrix with `A(i,i)*B(i,i) = m` for all `i` and
   * \f$ A_i \cdot B_j = 0\f$ for \f$i\neq j\f$. To get the second condition,
   * we simply have to recursively take for each line
   * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
   */
  // template <typename Matr, typename Int>

   template<typename Int> 
   void BasisConstruction<Int>::calcDualUpperTriangular (const IntMat & A, IntMat & B, int d, const Int & m) {
      for (int i = 0; i < d; i++) {
        for (int j = i + 1; j < d; j++)
          NTL::clear (B(i,j));
        DivideRound (m, A(i,i), B(i,i));
        for (int j = i - 1; j >= 0; j--) {
          NTL::clear (B(i,j));
          for (int k = j + 1; k <= i; k++)
            B(i,j) += A(j,k) * B(i,k);
          if (B(i,j) != 0)
            B(i,j) = -B(i,j);
          DivideRound (B(i,j), A(j,j), B(i,j));
        }
      }
    }


extern template class BasisConstruction<std::int64_t> ;
extern template class BasisConstruction<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
