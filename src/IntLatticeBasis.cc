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

#include "latticetester/Util.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaPalpha.h"

#include <fstream>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <typeinfo>

#ifdef WITH_NTL
#include "NTL/quad_float.h"
#include "NTL/RR.h"
using namespace NTL;
#else
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
using namespace boost::numeric::ublas;
#endif

using namespace std;


namespace LatticeTester
{

IntLatticeBasis::IntLatticeBasis (const int dim, NormType norm):
   m_dim (dim),
   m_norm (norm),
   m_modulo(0),
   m_withDual(false),
   m_xx(0)

{
   m_basis.resize(dim,dim);
   /*
#ifdef WITH_NTL
   ident(m_basis, dim);
#else
   m_basis = identity_matrix<long>(dim);
#endif
 */
   m_vecNorm.resize (dim);
   initVecNorm();
}

/*=========================================================================*/

IntLatticeBasis::IntLatticeBasis (const BMat basis, const int dim, NormType norm):
   m_basis (basis),
   m_dim (dim),
   m_norm (norm),
   m_modulo(0),
   m_withDual(false),
   m_xx(0)
{
   m_vecNorm.resize (dim);
   initVecNorm();
}

/*=========================================================================*/

IntLatticeBasis::IntLatticeBasis (
   const BMat primalbasis,
   const BMat dualbasis,
   const MScal modulo,
   const int dim,
   NormType norm):
   IntLatticeBasis(primalbasis, dim, norm)
{
   m_dualbasis = dualbasis;
   m_withDual = true;
   m_dualvecNorm.resize (dim);
   setDualNegativeNorm();
   m_modulo = modulo;
}

/*=========================================================================*/

IntLatticeBasis::IntLatticeBasis (const IntLatticeBasis & lat):
   m_dim (lat.getDim ()),
   m_norm (lat.getNorm ()),
   m_xx(0)
{
   copyBasis (lat);
}

/*=========================================================================*/


IntLatticeBasis::~IntLatticeBasis ()
{
   kill();
   m_basis.BMat::clear ();
   m_dualbasis.BMat::clear ();
   m_vecNorm.clear ();
   m_dualvecNorm.clear ();
}

/*=========================================================================*/


void IntLatticeBasis::kill ()
{
   delete [] m_xx;
   m_xx = 0;
}

/*=========================================================================*/

void IntLatticeBasis::copyBasis (const IntLatticeBasis & lat)
{
   if(m_dim == lat.m_dim)
      m_basis = lat.m_basis;
      m_dualbasis = lat.m_dualbasis;
      m_vecNorm = lat.m_vecNorm;
      m_dualvecNorm = lat.m_dualvecNorm;
      m_withDual = lat.m_withDual;
      m_modulo = lat.m_modulo;
      m_xx = new bool[m_dim];
      for (int i = 0; i < m_dim; i++)
         m_xx[i] = lat.getXX(i);
}

/*=========================================================================*/

void IntLatticeBasis::copyBasis (const IntLatticeBasis & lat, int n)
{
   if(m_dim == n) {
      CopyMatr(m_basis, lat.m_basis, n);
      CopyVect(m_vecNorm, lat.m_vecNorm, n);
      m_withDual = lat.m_withDual;
      if(m_withDual){
         m_dualbasis.resize(m_basis.size1(), m_basis.size1());
         m_dualvecNorm.resize(m_basis.size1());
         CopyMatr(m_dualbasis, lat.m_dualbasis, n);
         CopyVect(m_dualvecNorm, lat.m_dualvecNorm, n);
      }
      m_modulo = lat.m_modulo;
      m_xx = new bool[n];
      for (int i = 0; i < n; i++)
         m_xx[i] = lat.getXX(i);
      }
}

/*=========================================================================*/

void IntLatticeBasis::initVecNorm ()
{
   m_xx = new bool[m_dim];
   for(int i = 0; i < m_dim; i++){
      m_vecNorm[i] = -1;
      m_xx[i] = true;
   }
}


/*=========================================================================*/

void IntLatticeBasis::setNegativeNorm ()
{
   for (int i = 0; i<m_dim; i++){
      m_vecNorm[i] = -1;
   }
}

/*=========================================================================*/

void IntLatticeBasis::setDualNegativeNorm ()
{
   for(int i = 0; i < m_dim; i++){
      m_dualvecNorm[i] = -1;
   }
}

/*=========================================================================*/

void IntLatticeBasis::updateVecNorm ()
{
   updateVecNorm (0);
}

/*=========================================================================*/

void IntLatticeBasis::updateVecNorm (const int & d)
{
    assert (d >= 0);

   for (int i = d; i < m_dim; i++) {
      if (m_vecNorm[i] < 0) {
         matrix_row<BMat> row(m_basis, i);
         if (m_norm == L2NORM) {
            ProdScal (row, row, m_dim, m_vecNorm[i]);
         } else {
            CalcNorm <BVect, NScal> (row, m_dim, m_vecNorm[i], m_norm);
         }
      }
   }
}

/*=========================================================================*/

void IntLatticeBasis::updateDualVecNorm ()
{
   updateDualVecNorm (0);
}

/*=========================================================================*/

void IntLatticeBasis::updateDualVecNorm (const int & d)
{
    assert (d >= 0);

   for (int i = d; i < m_dim; i++) {
      if (m_dualvecNorm[i] < 0) {
         matrix_row<BMat> row(m_dualbasis, i);
         if (m_norm == L2NORM) {
            ProdScal (row, row, m_dim, m_dualvecNorm[i]);
         } else {
            CalcNorm <BVect, NScal> (row, m_dim, m_dualvecNorm[i], m_norm);
         }
      }
   }
}

/*=========================================================================*/

void IntLatticeBasis::updateScalL2Norm (const int i)
{
   matrix_row<BMat> row(m_basis, i);
   ProdScal (row, row, m_dim, m_vecNorm[i]);
}

/*=========================================================================*/

void IntLatticeBasis::updateScalL2Norm (const int k1, const int k2)
{
   for (int i = k1; i < k2; i++) {
      updateScalL2Norm(i);
   }
}

/*=========================================================================*/

void IntLatticeBasis::updateDualScalL2Norm (const int i)
{
   matrix_row<BMat> row(m_dualbasis, i);
   ProdScal (row, row, m_dim, m_dualvecNorm[i]);
}

/*=========================================================================*/

void IntLatticeBasis::updateDualScalL2Norm (const int k1, const int k2)
{
   for (int i = k1; i < k2; i++) {
      updateDualScalL2Norm(i);
   }
}


/*=========================================================================*/

void IntLatticeBasis::permute (int i, int j)
{
   if (i == j)
      return ;
   for (int k = 0; k < m_dim; k++){
      swap9 (m_basis(j,k), m_basis(i,k));
      if(m_withDual){
          swap9(m_dualbasis(j,k), m_dualbasis(i,k));
      }
   }
   swap9 (m_vecNorm[i], m_vecNorm[j]);
   if(m_withDual){
      swap9 (m_dualvecNorm[i], m_dualvecNorm[j]);
   }
   bool b = m_xx[j];
   m_xx[j] = m_xx[i];
   m_xx[i] = b;
}


/*=========================================================================*/

bool IntLatticeBasis::checkDuality ()
{
   if(!m_withDual) {
      cout << "DO NOT USE IntLatticeBasis::checkDuality without dual" << endl;
      return false;
   }
   BScal S;
   int dim = getDim ();

   for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
         matrix_row<const BMat> row1(m_basis, i);
         matrix_row<const BMat> row2(m_dualbasis, j);
         ProdScal (row1, row2, dim, S);
         if (j != i) {
            if (S != 0) {
               cout << "******  checkDuality failed for V[" << i <<
                    "] and W[" << j << "]" << endl;
               return false;
            }
         } else if (S != m_modulo) {
            cout << "******  checkDuality failed for i, j = " << i << " , " <<
                     j << endl;
            return false;
         }
      }
   }
   return true;


}

/*=========================================================================*/

void IntLatticeBasis::sort (int d)
/*
 * We assume that the square lengths are already updated.
 * This gives flexibility to the user to put something else than
 * the square Euclidean length in V.vecNorm, W.vecNorm, etc.
 */
{
   int dim = getDim ();
   for (int i = 0; i < dim; i++){
      if (getVecNorm(i) < 0) {
         cout << "\n***** ERROR: sort   Negative norm for i = " << i <<
            ",  dim = " << dim << endl;
      }
   }

   for (int i = d; i < dim; i++){
      int k = i;
      for (int j = i + 1; j < dim; j++) {
         if (getVecNorm (j) < getVecNorm (k))
            k = j;
      }
      if (i != k)
         permute (i, k);
   }
}

/*=========================================================================*/


void IntLatticeBasis::write () const
{
   cout << "Dim = " << m_dim << " \n \n";
   for (int i = 0; i < m_dim; i++) {
      cout << "  |     ";
      for (int j = 0; j < m_dim; j++) {
         cout << setprecision (15) << m_basis(i,j) << "\t";
      }
      cout << "    |";
      if(m_withDual){
         cout <<"  |     ";
         for (int j = 0; j < m_dim; j++) {
            cout << setprecision (15) <<
            m_dualbasis(i,j) << "\t";
         }
         cout << "    |";
      }
      cout << "\n";

   }
   cout << "\n";
   cout << "Norm used : " << toStringNorm(m_norm) << "\n" << endl;
   cout << "Norm of each Basis vector : \n";
   cout << " Primal     ";
   if(m_withDual)
      cout << "\t Dual \n";
   cout << "\n";

   for (int i = 0; i < m_dim; i++) {
      cout << "   ";
      if (m_vecNorm[i] < 0) {
         cout << "NaN OR Not computed";
      } else {
         cout << m_vecNorm[i];
      }
      if(m_withDual){
         cout << "\t \t \t ";
         if (m_dualvecNorm[i] < 0)
            cout << "NaN OR Not computed";
         else
            cout << m_dualvecNorm[i];
      }
      cout << "\n";
   }
}

/*=========================================================================*/

string IntLatticeBasis::toStringBasis () const
{
   ostringstream os;
   os << "Primal Basis:\n";
   os << "  Dim = " << m_dim << " \n";
   for (int i = 0; i < m_dim; i++) {
      os << "    [";
      for (int j = 0; j < m_dim; j++)
         os << " " <<  setprecision (15) << m_basis(i,j);
      os << " ]\n";
   }

   os << "  Norm:\n";
   for (int i = 0; i < m_dim; i++) {
      os << "    ";
      if (m_vecNorm[i] < 0) {
         os << "-1" << endl;
      } else {
         os << m_vecNorm[i] << endl;
      }
   }
   os << endl;
   return os.str ();
}

/*=========================================================================*/

string IntLatticeBasis::toStringDualBasis () const
{
   ostringstream os;
   os << "Dual Basis:\n";
   os << "  Dim = " << m_dim << " \n";
   for (int i = 0; i < m_dim; i++) {
      os << "    [";
      for (int j = 0; j < m_dim; j++)
         os << " " <<  setprecision (15) << m_dualbasis(i,j);
      os << " ]\n";
   }

   os << "  Norm:\n";
   for (int i = 0; i < m_dim; i++) {
      os << "    ";
      if (m_dualvecNorm[i] < 0) {
         os << "-1" << endl;
      } else {
         os << m_dualvecNorm[i] << endl;
      }
   }
   os << endl;
   return os.str ();
}


} //namespace LatticeTester

























