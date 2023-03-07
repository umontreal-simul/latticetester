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
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/IntLattice.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"

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
	 * Same as `LLLConstruction(IntMat &matrix, 0.999999)` (the default value of `delta`).  ?????  *****
	 */
	void LLLConstruction(IntMat &matrix);

	/**
	 * This functions takes a set of generating vectors of a vector space and
	 * finds a basis for this space by applying LLL reduction with the given value of `delta`,
	 * using the NTL implementation. It is implemented only for the \texttt{NTL::ZZ} type.
	 * This function *does not* assume that all vectors m e_i belong to the lattice, and
	 * it may return a basis matrix that has fewer rows than columns!   TRUE?    ***********
	 * If we want to make sure that these vectors belong to the lattice, we can add them
	 * explicitly to the set of generating vectors.
	 */
	void LLLConstruction(IntMat &matrix, double delta);

  	/**
	 * This is an old implementationSame that uses a form of Gaussian elimination to
	 * obtain an upper triangular basis for the smallest lattice that contains the generating
	 * vectors which are the rows of the given `matrix`. It returns the basis in the same `matrix`.
	 * This function *does not* assume that all vectors `m e_i` belong to the lattice, and
	 * it may return a basis matrix that has fewer rows than columns!
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
	 */
	//void mDualComputation(IntMat &matrix, IntMat &dualMatrix, Int m);

    /**
     * THE NAMES OF PARAMETERS SHOULD BE SIGNIFICANT AND UNIFORM ACROSS ALL METHODS and ALL CLASSES,
     * as much as possible!   *************
     * Takes a set of generating vectors in the matrix `gen` and iteratively
     * transforms it into a lower triangular lattice basis into the matrix `basis`.
     * `gen` and `basis` must have the same number of rows and the same number of columns.
     * All the computations are done modulo the scaling factor `m`.
     * After the execution, `gen` will contain irrelevant information (garbage)
     * and `basis` will contain an upper triangular basis.
     * Perhaps with zero rows at the end, in general, unless we assume implicitly that
     * all vectors of the form m e_i are in the generating set.  ???  *****************
     * The algorithm is explained in the \lattester{} guide.
     * Why pass &m instead of just m as elsewhere?  We cannot pass sometimes &m and sometimes m,
     * this is confusing and dangerous.                     **************
     */
    void lowerTriangularBasis(IntMat &gen, IntMat &basis, Int &m);

    /**
     * Similar to `lowerTriangularBasis`, except that the returned basis is upper triangular.
     */
    void upperTriangularBasis(IntMat &gen, IntMat &basis, Int &m);

    /**
     * Takes an upper triangular basis matrix `matrix` and computes the m-dual basis `dualMatrix`.
     * The method assumes that each coefficient on the diagonal of `matrix` is nonzero and divides `m`.
     * The algorithm is described in the Lattice Tester guide \cite iLEC22l.
     * What is `d`?  Note that the dual exists only if `matrix` is invertible.
     * So we must assume that there are no zeros on the diagonal.         *************
     */
    void mDualUpperTriangular (const IntMat &matrix, IntMat &matrixDual, int d, const Int &m);

    /**
	 * This function does essentially the same thing as `mDualUpperTriangular`, but it is
	 * slightly slower. It uses the method described in \cite rCOU96a.
	 */
	void mDualUpperTriangular96(IntMat &matrix, IntMat &matrixDual, Int m);

	/**
	 * This function assumes that `matrix` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular, and it returns in `matrixDual`
	 * the m-dual of matrix.
     */
    void mDualBasis (const IntMat &matrix, IntMat &dualMatrix, const Int & m);
        
    // Why this function here ???   This is only a special case.   **********************
    void mDualBasis(const NTL::Mat<NTL::ZZ>  & A, NTL::Mat<NTL::ZZ>  & B, const NTL::ZZ & m);

	/**
	 * Constructs a basis for the projection `proj` of the lattice `in`,
	 * using LLLConstruction, and puts it in `out`. The basis is not triangular.
	 * This will overwrite the lattice basis in `out` and change the dimension.
	 * It does not update the dual.
	 */
   template<typename Real>
    void projectionConstructionLLL(IntLattice<Int, Real> &in,
			IntLattice<Int, Real> &out, const Coordinates &proj);

};

//============================================================================
// Implementation

// Do we need all of these?  Why?   ***************
//
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
	spec.LLLConstruction(matrix, delta);
}

template<typename Int>
void BasisConstruction<Int>::GCDTriangularBasis(IntMat &matrix, Int m) {
	// On exit, the rows of matrix are the basis vectors.
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
				  Modulo(matrix[j][k], m, matrix[j][k]);
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

//===================================================================

   template<typename Int> 
   void BasisConstruction<Int>::upperTriangularBasis(IntMat &gen, IntMat &basis, Int &m) {
     IntVec coeff, vl,v2; 
     Int C, D, val, gcd;  
     int pc, pl, k;
     int dim1=gen.NumRows();
     int dim2=gen.NumCols();

     pl=0;
     pc=0;
     while(pl<dim1 && pc<dim2){
           for(int i=0;i<dim1;i++)
             Modulo (gen(i,pc), m, gen(i,pc));
                
            coeff.SetLength(dim2);
            k=0;     
            while( k<dim1 && gen(k,pc)==0)
             { coeff[k]=0; 
               k++;
             }
                
           if(k<dim1)
            { gcd=gen(k,pc);
              coeff[k]=1;
              val=gcd;
             for(int i=k+1;i<dim1; i++){
               if(gen(i,pc)==0)
               { coeff[i]= 0;
                 continue;
                }           
              Euclide (val, gen(i,pc), C, D , gcd);
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
                 vl[j]=vl[j]+coeff[ind]*gen(ind,j);
              } 
              Modulo (vl[j], m, vl[j]);
             }
             for(int i=0;i<dim1;i++)
             {  if(gen(i,pc)!=0){
                v2= (gen(i,pc)/gcd)*vl;
                for(int j=pc;j<dim2;j++)
                    Modulo (v2[j], m, v2[j]);
                for(int j=pc;j<dim2;j++)
                 {   
                   gen(i,j)=gen(i,j)-v2[j];
                   Modulo (gen(i,j), m, gen(i,j));
                 } 
                }    
             }
             basis[pl]=vl;
          }
          else
          {  for (int j1 = 0; j1 < dim2; j1++) {
             if (j1 != pl)
               NTL::clear (basis(pl,j1));
             else
               basis(pl,j1) = m;
             }   
          }
          coeff.clear();
          vl.clear();
          pl++; 
          pc++;
       }
    }

//==============================================================================

   template<typename Int> 
   void  BasisConstruction<Int>::lowerTriangularBasis(IntMat &gen, IntMat &basis, Int &m){
     IntVec coeff, vl, v2;
     Int C, D, val, gcd;  
     int pc, pl, k;
     int dim1=gen.NumRows();
     int dim2=gen.NumCols();
     pl=dim1-1;
     pc=dim2-1;
     while(pl>=0 && pc>=0) {
           for(int i=0;i<dim1;i++)
              Modulo (gen(i,pc), m, gen(i,pc));
            coeff.SetLength(dim2);
            k=0;     
            while( k<dim1 && gen(k,pc)==0)
             { coeff[k]=0; 
               k++;
             }     
           if(k<dim1)
            { gcd=gen(k,pc);
              coeff[k]=1;
              val=gcd;      
             for(int i=k+1;i<dim1; i++){
               if(gen(i,pc)==0)
               { coeff[i]= 0;
                 continue;
                }   
              Euclide (val, gen(i,pc), C, D , gcd);
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
                 vl[j]=vl[j]+coeff[ind]*gen(ind,j);
              } 
              Modulo (vl[j], m, vl[j]);
             }
             for(int i=0;i<dim1;i++)
             {  if(gen(i,pc)!=0){
                v2= (gen(i,pc)/gcd)*vl;
                for(int j=0;j<dim2;j++)
                    Modulo (v2[j], m, v2[j]);
                for(int j=0;j<dim2;j++)
                 {   
                   gen(i,j)=gen(i,j)-v2[j];
                   Modulo (gen(i,j), m, gen(i,j));
                 } 
                }    
             }
             basis[pl]=vl;
          }
          else
          {  for (int j1 = 0; j1 < dim2; j1++) {
             if (j1 != pl)
               NTL::clear (basis(pl,j1));
             else
               basis(pl,j1) = m;
             }   
          }
          coeff.clear();
          vl.clear();
          pl--; 
          pc--;
       }
   }

   //======================================================

   // This is the old version from Couture and L'Ecuyer (1996).
   template<typename Int>
   void BasisConstruction<Int>::mDualUpperTriangular96(IntMat &matrix, IntMat &dualMatrix,
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

   //===================================================

   /**
    * This is the version that we recommend.    ****************
    *
   * For `B` to be `m`-dual to `A`, we have to have that \f$AB^t = mI\f$. It
   * is quite easy to show that, knowing `A` is upper triangular, `B` will be a
   * lower triangular matrix with `A(i,i)*B(i,i) = m` for all `i` and
   * \f$ A_i \cdot B_j = 0\f$ for \f$i\neq j\f$. To get the second condition,
   * we simply have to recursively take for each line
   * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
   */
   template<typename Int>
   void BasisConstruction<Int>::mDualUpperTriangular
          (const IntMat & A, IntMat & B, int d, const Int & m) {
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

   
    //===================================================

    // This one seems to be too specific!    ***************
    /*
    template<typename Int>    // There is no Int in this function!
    void BasisConstruction<Int>::mDualBasis(const NTL::Mat<NTL::ZZ>  & A, NTL::Mat<NTL::ZZ>  & B, const NTL::ZZ & m) {
      NTL::ZZ d;
      //  Int d;// mult;
      //  IntMat C;
      NTL::Mat<NTL::ZZ> C;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      C.SetDims(dim1, dim2);
      inv(d, B, A);
      transpose(C, B);
      for (int i = 1; i < dim1; i++) {
        for (int j = 1; j < dim2; j++)
           B(i,j)= (m*C(i,j))/d;
          
        }
     }

    */
   //============================================================================

   template<typename Int> 
   void BasisConstruction<Int>::mDualBasis(const IntMat & A, IntMat & B, const Int & m) {
  // switch (typeof(Int)) { 

  // case NTL::ZZ :    
      NTL::ZZ d;
      int dim1=A.NumRows();
      int dim2=A.NumCols();
      NTL::Mat<NTL::ZZ> A1,C1,B1,B2;
     // NTL::Mat<NTL::ZZ> B2;
      A1.SetDims(dim1, dim2);
      B1.SetDims(dim1, dim2);
      B2.SetDims(dim1, dim2);
      C1.SetDims(dim1, dim2);
      for(int i=0;i<dim1;i++){
       for(int j=0;j<dim2;j++)
        A1[i][j]=NTL::conv<NTL::ZZ>(A[i][j]);
      }  
      inv(d,B1,A1);
      transpose(C1,B1);
      NTL::ZZ mm=NTL::conv<NTL::ZZ>(m);
      for (int i = 0; i < dim1; i++) {
       for (int j = 0; j < dim2; j++){
        B[i][j]= NTL::conv<Int>(mm*C1[i][j]/d);}   
      }
    //  break;
  //  case std::int64_t :   
   //      std::cout << " Error  stg::int64_t not supported \n";
    //     break;

   //  }  

     }       
   //=================================================================================

   
   template<typename Int>
   template<typename Real>
   void BasisConstruction<Int>::projectionConstructionLLL(
   		IntLattice<Int, Real>& in,
   		IntLattice<Int, Real>& out, const Coordinates& proj) {
	  std::size_t size = proj.size();
   	  unsigned int lat_dim = in.getDim();
   	  if (size > lat_dim)
   		 MyExit(1, "More projection coordinates than the dimension of `in`.");
   	  IntMat new_basis, tmp(NTL::transpose(in.getBasis()));
   	  new_basis.SetDims(size, tmp.NumRows());
   	  tmp = NTL::transpose(tmp);
   	  auto it = proj.cbegin();
   	  for (std::size_t i = 0; i < size; i++) {
   		if (*it <= lat_dim)
   			new_basis[i] = tmp[*it];
   		else
   			MyExit(1, "Coordinates do not match the dimensions of `in`.");
   		it++;
    	}
   	  new_basis = NTL::transpose(new_basis);
   	  LLLConstruction(new_basis);
   	  out = IntLattice<Int, Real>(new_basis, size, in.getNormType());
      }


extern template class BasisConstruction<std::int64_t> ;
extern template class BasisConstruction<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
