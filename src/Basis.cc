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

#include <iomanip>
#include <cassert>
#include <ostream>

#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latticetester/Basis.h"
#include "latticetester/Util.h"


#ifdef WITH_NTL
using namespace NTL;
#else
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;
#endif

using namespace std;


namespace LatticeTester
{

//===========================================================================

Basis::Basis (int d, int maxDim, NormType norm)
{
   assert (d > 0);
   assert (maxDim > 0);

   if (d > maxDim)
      d = maxDim;

   m_dim = d;
   m_maxDim = maxDim;
   m_norm = norm;

   BMat::resize (1 + maxDim, 1 + maxDim);
   m_vecNorm.resize (1 + maxDim);
   m_negFlag = new bool[1 + maxDim];
   for (int i = 0; i < maxDim; i++) {
      m_negFlag[i] = true;
      for (int j = 0; j < maxDim; j++) {
         (*this)(i,j) = 0;
      }
   }
}

//===========================================================================
// Erwan

Basis::Basis (int dim, NormType norm)
{
   assert (dim > 0);

   m_dim = dim;
   m_maxDim = dim;
   m_norm = norm;

   BMat::resize (dim, dim);
   m_vecNorm.resize (dim);
   m_negFlag = new bool[1 + dim];
   for (int i = 0; i < dim; i++) {
      m_negFlag[i] = true;
      for (int j = 0; j < dim; j++) {
         (*this)(i,j) = 0;
      }
   }
}


/*=========================================================================*/

Basis::Basis (const Basis & base): BMat (base)
{
   m_dim = base.m_dim;
   m_maxDim = base.m_maxDim;
   BMat::operator=(base);
   m_vecNorm = base.m_vecNorm;

   int d = m_maxDim + 1;
   m_negFlag = new bool[d];
   for (int i = 0; i < d; i++)
      m_negFlag[i] = base.m_negFlag[i];
   m_norm = base.m_norm;
}


/*=========================================================================*/

Basis::~Basis ()
{
   kill ();
}


/*=========================================================================*/

void Basis::kill ()
{
   if (m_negFlag != 0) {
      delete[] m_negFlag;
      m_negFlag = 0;
   }
   m_vecNorm.clear ();
   BMat::clear ();
}

/*=========================================================================*/

Basis & Basis::operator= (const Basis & base)
{
   if (this == &base)
      return * this;
   kill();
   m_norm = base.m_norm;
   m_dim = base.m_dim;
   m_maxDim = base.m_maxDim;
   BMat::operator=(base);
   m_vecNorm = base.m_vecNorm;
   m_negFlag = new bool[m_maxDim + 1];
   for (int i = 0; i < m_dim; i++)
      m_negFlag[i] = base.m_negFlag[i];

   return *this;
}


/*=========================================================================*/

void Basis::swap (Basis & b)
{
   LatticeTester::swap9 (m_dim, b.m_dim);
   LatticeTester::swap9 (m_maxDim, b.m_maxDim);
   LatticeTester::swap9 (m_vecNorm, b.m_vecNorm);
   LatticeTester::swap9 (*this, b);
   LatticeTester::swap9 (m_negFlag, b.m_negFlag);
   LatticeTester::swap9 (m_norm, b.m_norm);
}


/*=========================================================================*/

void Basis::permute (int i, int j)
{
   if (i == j)
      return ;
   for (int k = 0; k < m_dim; k++) {
      LatticeTester::swap9 ((*this)(j,k), (*this)(i,k));
   }
   LatticeTester::swap9 (m_vecNorm[i], m_vecNorm[j]);
   LatticeTester::swap9 (m_negFlag[i], m_negFlag[j]);
}


/*=========================================================================*/

void Basis::setDim (int d)
{
   assert (d > 0);

   if (d > m_maxDim) {
      m_dim = m_maxDim;
   } else {
      m_dim = d;
   }
}


/*=========================================================================*/

void Basis::setNorm (NormType norm)
{
   if (m_norm == norm)
      return ;
   m_norm = norm;
   setNegativeNorm (true);
}


/*=========================================================================*/

void Basis::setVecNorm (NScal & value, int i)
{
      assert (value >= 0);
      assert (i >= 0 && i < m_maxDim);
   m_vecNorm[i] = value;
   m_negFlag[i] = false;
}


/*=========================================================================*/

void Basis::updateVecNorm ()
{
   updateVecNorm (0);
}


/*=========================================================================*/

void Basis::updateVecNorm (int d)
{
   assert (d >= 0);

   for (int i = d; i < m_dim; i++) {
      if (m_negFlag[i]) {
         matrix_row<BMat> row(*this, i);
         if (m_norm == L2NORM) {
            ProdScal (row, row, m_dim, m_vecNorm[i]);
         } else {
            CalcNorm <BVect, NScal> (row, m_dim, m_vecNorm[i], m_norm);
         }
        m_negFlag[i] = false;
      }

      if (m_vecNorm[i] < 0) {
         cout << "\n***** ERROR:   Negative norm for i = " << i <<
         ",  dim = " << m_dim << endl;
         cout << "Basis vector:   ";
         matrix_row<BMat> row2(*this, i);
         cout << LatticeTester::toString (row2, 0, m_dim, "  ") << endl;
         exit(1);
      }
   }
}


/*=========================================================================*/

void Basis::updateScalL2Norm (int i)
{
   if (m_negFlag[i]) {
      matrix_row<BMat> row(*this, i);
      ProdScal (row, row, m_dim, m_vecNorm[i]);
   }
   m_negFlag[i] = false;
   if (m_vecNorm[i] < 0) {
      cout << "\n***** ERROR:   Negative norm L2 for i = " << i <<
      ",  dim = " << m_dim << endl;
      cout << "Basis vector:   ";
      matrix_row<BMat> row2(*this, i);
      cout << LatticeTester::toString (row2, 0, m_dim, "  ") << endl;
      exit(1);
   }
}


/*=========================================================================*/

void Basis::updateScalL2Norm (int k1, int k2)
{
   for (int i = k1; i < k2; i++) {
      if (m_negFlag[i]) {
         matrix_row<BMat> row(*this, i);
         ProdScal (row, row, m_dim, m_vecNorm[i]);

      }
      m_negFlag[i] = false;
      if (m_vecNorm[i] < 0) {
         cout << "\n***** ERROR:   Negative norm L2 for i = " << i <<
         ",  dim = " << m_dim << endl;
         cout << "Basis vector:   ";

         matrix_row<BMat> row2(*this, i);
         cout << LatticeTester::toString (row2, 0, m_dim, "  ") << endl;
         exit(1);
      }
   }
}


/*=========================================================================*/

void Basis::setNegativeNorm (bool flag)
{
   for (int i = 0; i < m_dim; i++) {
      m_negFlag[i] = flag;
   }
}


/*=========================================================================*/

string Basis::toString() const
{
   ostringstream os;
   os << "Basis:\n";
   os << "  Dim = " << m_dim << " \n";
   for (int i = 0; i < m_dim; i++) {
      os << "    [";
      for (int j = 0; j < m_dim; j++)
         os << " " <<  setprecision (15) << (*this)(i,j);
      os << " ]\n";
   }

   os << "  Norm:\n";
   for (int i = 0; i < m_dim; i++) {
      os << "    ";
      if (m_negFlag[i]) {
         os << "-1" << endl;
      } else if (m_vecNorm[i] < 0) {
         os << "Erreur:  norm < 0 and negFlag not set" << endl;
      } else {
         os << m_vecNorm[i] << endl;
      }
   }
   os << endl;
   return os.str ();
}


/*=========================================================================*/

string Basis::toString(int i) const
{
   ostringstream os;
   os << "Dim = " << m_dim << " \n";
   os << "   Basis[" << i << "] =  [";
   for (int j = 0; j < m_dim; j++)
      os << "  " << (*this)(i,j);
   os << " ]\n";

   if (L2NORM == m_norm)
      os << "   Length^2 =   ";
   else
      os << "   Length =   ";
   if (m_negFlag[i]) {
      os << "-1" << endl;
   } else if (m_vecNorm[i] < 0) {
      os << "Erreur:  norm < 0 and negFlag not set" << endl;
   } else {
      os << setprecision (15) << m_vecNorm[i] << endl;
      NScal x = 1.0 / m_vecNorm[i];
      if (L2NORM == m_norm)
         x = sqrt(x);
      os << "   1/Length =   " << x << endl;
   }
   os << endl;
   return os.str ();
}


/*=========================================================================*/

void Basis::write () const
{
   cout << "Dim = " << m_dim << " \n";
   for (int i = 0; i < m_dim; i++) {
      cout << "   | ";
      for (int j = 0; j < m_dim; j++) {
         cout << setprecision (15) << (*this)(i,j) << "\t";
      }
      cout << " |" << endl;
   }
   // return;
   cout << "Norm:\n";

   for (int i = 0; i < m_dim; i++) {
      cout << "   ";
      if (m_negFlag[i]) {
         cout << "-1" << endl;
      } else if (m_vecNorm[i] < 0) {
         cout << "Erreur:  norm < 0 and negFlag not set" << endl;
      } else {
         cout << m_vecNorm[i] << endl;
      }
   }
   cout << endl;
}


/*=========================================================================*/

void Basis::write (int i) const
{
   cout << "Dim = " << m_dim << " \n";
   cout << "   Basis[" << i << "] =  [ ";
   for (int j = 0; j < m_dim; j++) {
      cout << setprecision (15) << (*this)(i,j) << "  ";
   }
   cout << "]" << endl;

   if (L2NORM == m_norm)
      cout << "   Length^2 =   ";
   else
      cout << "   Length =   ";
   if (m_negFlag[i]) {
      cout << "-1" << endl;
   } else if (m_vecNorm[i] < 0) {
      cout << "Erreur:  norm < 0 and negFlag not set" << endl;
   } else {
      cout << m_vecNorm[i] << endl;
      NScal x = 1.0 / m_vecNorm[i];
      if (L2NORM == m_norm)
         x = sqrt(x);
      cout << "   1/Length =   " << x << endl;
   }
   cout << endl;
}

}
