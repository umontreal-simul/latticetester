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

#include <limits>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <typeinfo>
#include <stdexcept>

#include "latticetester/Util.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Reducer.h"

#ifdef WITH_NTL
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include "NTL/vec_ZZ.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#else
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;
#endif


using namespace std;


namespace LatticeTester
{

// Initialization of non-const static members
bool Reducer::PreRedDieterSV = false;
bool Reducer::PreRedLLLSV = false;
bool Reducer::PreRedLLLRM = false;
bool Reducer::PreRedBKZ = true;
long Reducer::maxNodesBB = 10000000;



//=========================================================================

Reducer::Reducer (IntLatticeBasis & lat)
{
   // Squared length of vectors must not overflow max(double)
   m_lat = &lat;
   const int dim1 = m_lat->getDim ();
   int dim2 = dim1;
   if (dim2 <= 2)
      dim2++;

   m_c0.resize (dim1, dim1);
   m_c2.resize (dim1, dim1);
   m_cho2.resize (dim2, dim2);
   m_gramVD.resize (dim2, dim2);
   // Indices in Cholesky go as high as m_lat->getMaxDim() + 2
   m_IC = new int[2 + dim1];

   m_nv.resize (dim1);
   m_bv.resize (dim1);
   m_bw.resize (dim1);
   m_n2.resize (dim1);
   m_zLR.resize (dim1);
   m_zLI.resize (dim1);
   m_zShort.resize (dim1);
   m_dc2.resize (dim1);
   m_BoundL2.resize (dim1);

   m_lMin = std::numeric_limits <double>::max ();
   m_lMin2 = m_lMin;
   for (int i = 0; i < dim1; i++) {
      m_zLI[i] = -1;
      m_zShort[i] = -1;
      m_BoundL2[i] = -1;
      m_IC[i] = -1;
   }
   m_countNodes = 0;
}


//=========================================================================

Reducer::Reducer (const Reducer & red)
{
   copy (red);
}


//=========================================================================

Reducer & Reducer::operator= (const Reducer & red)
{
   if (this != &red)
      copy (red);
   return *this;
}


//=========================================================================

void Reducer::copy (const Reducer & red)
{
   m_lat = red.m_lat;
   m_c0 = red.m_c0;
   m_c2 = red.m_c2;
   m_dc2 = red.m_dc2;
   m_nv = red.m_nv;
   m_bv = red.m_bv;
   m_bw = red.m_bw;
   m_n2 = red.m_n2;
   m_zLR = red.m_zLR;
   m_zLI = red.m_zLI;
   m_zShort = red.m_zShort;
   m_cho2 = red.m_cho2;
   m_gramVD = red.m_gramVD;
   m_lMin = red.m_lMin;
   m_lMin2 = red.m_lMin2;
   m_BoundL2 = red.m_BoundL2;
   m_IC = new int[3 + m_lat->getDim ()];
   for (int i = 0; i < 2 + m_lat->getDim (); i++)
      m_IC[i] = red.m_IC[i];
}


//=========================================================================

Reducer::~Reducer ()
{
   m_c0.clear();
   m_c2.clear();
   m_cho2.clear();
   m_gramVD.clear();
   m_nv.clear();
   m_bv.clear();
   m_bw.clear();
   m_n2.clear();
   m_zLR.clear();
   m_zLI.clear();
   m_zShort.clear();
   m_dc2.clear();
   m_BoundL2.clear();
   delete[] m_IC;
}


//=========================================================================

void Reducer::setBoundL2 (const NVect & V, int dim1, int dim2)
{
   m_BoundL2.resize(dim2);
   for (int i = dim1; i < dim2; i++)
      m_BoundL2[i] = V[i];
}


//=========================================================================

void Reducer::permuteGramVD (int i, int j, int n)
/*
 * Permute la i-ieme et la j-ieme ligne, et la i-ieme et j-ieme colonne
 * de la matrice des produits scalaires.   Pour redLLL.
 */
{
   int k;
   for (k = 0; k < n; k++)
   {
      swap(m_gramVD(i,k), m_gramVD(j,k));
   }
   for (k = 0; k < n; k++)
   {
      swap(m_gramVD(k,i), m_gramVD(k,j));
   }
}


//=========================================================================

void Reducer::calculCholeski2LLL (int n, int j)
/*
 * Recalcule les n premieres entrees de la colonne j de la matrice de
 * Choleski d'ordre 2.   Pour redLLL.
 */
{
   m_cho2(0,j) = m_gramVD(0,j);
   for (int i = 1; i <= n; i++) {
      m_cho2(i,j) = m_gramVD(i,j);
      for (int k = 0; k < i; k++) {
         m_cho2(i,j) -= (m_cho2(k,j) / m_cho2(k,k)) * m_cho2(k,i);
      }
   }
}


//=========================================================================

inline void Reducer::calculCholeski2Ele (int i, int j)
// Recalcule l'entree (i,j) de la matrice de Choleski d'ordre 2
{
   m_cho2(i,j) = m_gramVD(i,j);
   for (int k = 0; k < i; k++) {
      m_cho2(i,j) -= m_cho2(k,i) * (m_cho2(k,j) / m_cho2(k,k));
   }
}


//=========================================================================

void negativeCholeski ()
{
   cout << "\n***** Negative diagonal element in Choleski Decomposition\n"
        << endl;
}


//=========================================================================

inline void Reducer::calculGramVD ()
/*
  Retourne dans m_gramVD la matrice des produits scalaires m_lat->V[i]*m_lat->V[j].
  Rem.: m_lat->V.vecNorm ne contient que les m_lat->V[i]*m_lat->V[i].
  Reviens à faire m_lat->V * transpose(m_lat->V).
  Utilise pour redLLL.
 */
{
   const int dim = m_lat->getDim ();

   for (int i = 0; i < dim; i++) {
      for (int j = i; j < dim; j++) {
         matrix_row<BMat> row1(m_lat->getBasis(), i);
         matrix_row<BMat> row2(m_lat->getBasis(), j);
         ProdScal (row1, row2, dim, m_gramVD(i,j));
         m_gramVD(j,i) = m_gramVD(i,j);
      }
   }
}


//=========================================================================

inline void Reducer::miseAJourGramVD (int j)
/*
 * Recalcule la j-ieme ligne et la j-ieme colonne de la matrice des
 * produits scalaires.  Pour redLLL.
 */
{
   const int dim = m_lat->getDim ();
   for (int i = 0; i < dim; i++) {
      matrix_row<BMat> row1(m_lat->getBasis(), i);
      matrix_row<BMat> row2(m_lat->getBasis(), j);
      ProdScal (row1, row2, dim, m_gramVD(i,j));
      m_gramVD(j,i) = m_gramVD(i,j);
   }
}


//=========================================================================

bool Reducer::calculCholeski (RVect & DC2, RMat & C0)
/*
 * Returns in C0 the elements of the upper triangular matrix of the
 * Choleski decomposition that are above the diagonal. Returns in DC2 the
 * squared elements of the diagonal. These elements are rescaled by EchVV
 * when SISquares= true.
 */
{
   const int dim = m_lat->getDim ();

   int k, j, i;
   RScal m2;
   // C2(i,j) = C0(i,j) * C2(i,i) if i != j.
   // C2(i,i) = DC2[i].
   conv (m2, m_lat->getModulo ());
   m2 = m2*m2;
   int d = dim;
   if(m_lat->withDual())
      d = dim / 2; // If we use the Dual, we compute Choleski
                   // with the Dual

   // Compute the d first lines of C0 with the primal Basis.
   for (i = 0; i < d; i++) {
      m_lat->updateScalL2Norm (i);
      for (j = i; j < dim; j++) {
         if (j == i)
            conv (m_c2(i,i), m_lat->getVecNorm (i));
         else {
            matrix_row<BMat> row1(m_lat->getBasis(), i);
            matrix_row<BMat> row2(m_lat->getBasis(), j);
            ProdScal (row1, row2, dim, m_c2(i,j));
         }
         for (k = 0; k < i; k++)
            m_c2(i,j) -= C0(k,i) * m_c2(k,j);
         if (i == j) {
            DC2[i] = m_c2(i,i);
            if (DC2[i] < 0.0) {
               negativeCholeski ();
               return false;
            }
         } else
            C0(i,j) = m_c2(i,j) / DC2[i];
      }
   }

   // Compute the d last lines of C0 with the dual Basis.
   /* This operation with the dual is needed in case of high dimension
    * and large number (30 bits). The choleski decomposition
    * convert numbers in double which is not sufficient in that case.
    * You need to use RR of the NTL library for this calcul.
   */
   if(m_lat->withDual()){
      for (i = dim-1; i >= d; i--)
      {
         m_lat->updateDualScalL2Norm (i);
         for (j = i; j >= 0; j--) {
            if (j == i)
               conv (m_c2(i,i), m_lat->getDualVecNorm (i));
            else {
               matrix_row<BMat> row1(m_lat->getDualBasis(), i);
               matrix_row<BMat> row2(m_lat->getDualBasis(), j);
               ProdScal (row1, row2, dim, m_c2(i,j));
            }
            for (k = i + 1; k < dim; k++){
               m_c2(i,j) -= C0(k,i) * m_c2(k,j);
            }
            if (i != j)
               C0(i,j) = m_c2(i,j) / m_c2(i,i);
         }
         DC2[i] = m2 / m_c2(i,i);
         if (DC2[i] < 0.0) {
            negativeCholeski();
            return false;
         }
         for (j = i + 1; j < dim; j++) {
            C0(i,j) = -C0(j,i);
            for (k = i + 1; k < j; k++) {
               C0(i,j) -= C0(k,i) * C0(k,j);
            }
         }
      }
   }
   return true;
}


//=========================================================================

void Reducer::trace (char *mess)
{
   cout << endl << "================================= " << mess << endl;
   //cout << "dim = " << m_lat->getDim () << endl;
   m_lat->setNegativeNorm();
   m_lat->setDualNegativeNorm();
   m_lat->updateVecNorm ();
   m_lat->updateDualVecNorm();
   m_lat->sort(0);
   m_lat->write();
}


//===========================================================================

void Reducer::pairwiseRedPrimal (int i, int d)
{
   const int dim = m_lat->getDim ();
   ++m_countDieter;
   m_lat->updateScalL2Norm (i);
   bool modifFlag;

   for (int j = d; j < dim; j++) {
      if (i == j)
         continue;
      modifFlag = false;
      {
         matrix_row<BMat> row1(m_lat->getBasis(), i);
         matrix_row<BMat> row2(m_lat->getBasis(), j);
         ProdScal (row1, row2, dim, m_ns);
      }
      DivideRound <NScal> (m_ns, m_lat->getVecNorm (i),
                           m_ns);
      if (m_ns == 0)
         continue;
      conv (m_bs, m_ns);
      if (m_ns < 1000 && m_ns > -1000) {
         m_lat->updateScalL2Norm (j);
         {
            matrix_row<BMat> row1(m_lat->getBasis(), j);
            matrix_row<BMat> row2(m_lat->getBasis(), i);
            ModifVect (row1, row2, -m_bs, dim);
         }

         // Verify that m_lat->getBasis()[j] is really shorter
         {
         matrix_row<BMat> row1(m_lat->getBasis(), j);
         ProdScal (row1, row1, dim, m_ns);
         }
         if (m_ns >= m_lat->getVecNorm (j)) {
            matrix_row<BMat> row1(m_lat->getBasis(), j);
            matrix_row<BMat> row2(m_lat->getBasis(), i);
            ModifVect (row1, row2, m_bs, dim);
         } else {
            modifFlag = true;
            m_lat->setVecNorm (m_ns, j);
         }
      } else {
         matrix_row<BMat> row1(m_lat->getBasis(), j);
         matrix_row<BMat> row2(m_lat->getBasis(), i);
         ModifVect (row1, row2, -m_bs, dim);
         m_lat->setNegativeNorm (j);
         modifFlag = true;
      }

      if (modifFlag) {
         m_countDieter = 0;
         ++m_cpt;
         if(m_lat->withDual()){
            matrix_row<BMat> row1(m_lat->getDualBasis(), i);
            matrix_row<BMat> row2(m_lat->getDualBasis(), j);
            ModifVect (row1, row2, m_bs, dim);
            m_lat->setDualNegativeNorm(i);

         }
         m_lat->setXX (false, i);
         m_lat->setXX (false, j);
      }
   }
}

//=========================================================================

void Reducer::pairwiseRedDual (int i)
{
   int j;
   const int dim = m_lat->getDim ();

   ++m_countDieter;
   m_lat->updateDualScalL2Norm (i);
   matrix_row<BMat> row9(m_lat->getBasis(), i);
   m_bv = row9;
   for (j = 0; j < dim; j++) {
      if (i != j) {
         matrix_row<BMat> row1(m_lat->getDualBasis(), i);
         matrix_row<BMat> row2(m_lat->getDualBasis(), j);
         ProdScal (row1, row2, dim, m_ns);
         // ProdScal (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
         //           dim, m_ns);
         DivideRound <NScal> (m_ns, m_lat->getDualVecNorm (i),
                              m_nv[j]);
         if (m_nv[j] != 0) {
            conv (m_bs, m_nv[j]);
            matrix_row<BMat> row7(m_lat->getBasis(), j);
            ModifVect (m_bv, row7, m_bs, dim);
         }
      }
   }

   m_lat->updateScalL2Norm (i);
   ProdScal (m_bv, m_bv, dim, m_ns);
   if (m_ns < m_lat->getVecNorm (i)) {
      ++m_cpt;
      m_countDieter = 0;
      matrix_row<BMat> row6(m_lat->getBasis(), i);
      for (j = 0; j < dim; j++)
         row6(j) = m_bv[j];
      m_lat->setNegativeNorm (i);
      m_lat->setXX (false, i);
      m_lat->setVecNorm (m_ns, i);
      for (j = 0; j < dim; j++) {
         if (i != j && m_nv[j] != 0) {
            conv (m_bs, -m_nv[j]);
            matrix_row<BMat> row1(m_lat->getDualBasis(), j);
            matrix_row<BMat> row2(m_lat->getDualBasis(), i);
            ModifVect (row1, row2, m_bs, dim);
            //  ModifVect (m_lat->getDualBasis ()[j], m_lat->getDualBasis ()[i],
            //            m_bs, dim);
            m_lat->setDualNegativeNorm (j);
            m_lat->setXX (false, j);
         }
      }
   }
}

//=========================================================================

void Reducer::preRedDieter(int d)
{
   long BoundCount;
   const int dim = m_lat->getDim ();
   bool withDual = m_lat->withDual();

   m_lat->updateScalL2Norm (d, dim);
   m_lat->sort (d);
   int i = dim-1;
   m_cpt = 0;
   m_countDieter = 0;
   BoundCount = 2 * dim - d;
   do {
      pairwiseRedPrimal (i, d);
      if (i > d && withDual)
         pairwiseRedDual (i);
      if (i < 1)
         i = dim-1;
      else
         --i;
   } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
}

//=========================================================================

void Reducer::preRedDieterRandomized (int d, int seed)
{
   long BoundCount;
   const int dim = m_lat->getDim ();
   bool withDual = m_lat->withDual();

   m_lat->updateScalL2Norm (d, dim);
   //m_lat->getDualBasis ().updateScalL2Norm (d, dim);
   m_lat->sort (d);
   int i = dim-1;
   m_cpt = 0;
   m_countDieter = 0;
   BoundCount = 2 * dim - d;
   srand(seed);
   do {
      pairwiseRedPrimal (rand() % dim, d);
      if (i > d && withDual)
         pairwiseRedDual (i);
      if (i < 1)
         i = dim-1;
      else
         --i;
   } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
}



//=========================================================================

/**
 * Reduce the Choleski matrix with adding a multiple of the i-th vector
 * to the j-th vector. It updates the Gram Schmidt matrix
 */
void Reducer::reductionFaible (int i, int j)
/*
 * Reduit la matrice de Choleski (d'ordre 2) en ajoutant un multiple du
 * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
 * consequence et remet a jour la matrice des produits scalaires.
 * Utilise par redLLL.
 */
{
   RScal cte;
   long cteLI;
   bool withDual = m_lat->withDual();
   cte = m_cho2(i,j) / m_cho2(i,i);

   const int dim = m_lat->getDim ();

   if (fabs (cte) < std::numeric_limits<double>::max()) {
      // On peut representer cte en LONGINT.
      if (fabs (cte) > 0.5) {
         conv (cteLI, Round (cte));
         matrix_row<BMat> row1(m_lat->getBasis(), j);
         matrix_row<BMat> row2(m_lat->getBasis(), i);
         ModifVect (row1, row2, -cteLI, dim);
         if(withDual){
            matrix_row<BMat> row3(m_lat->getDualBasis(), i);
            matrix_row<BMat> row4(m_lat->getDualBasis(), j);
            ModifVect (row3, row4, cteLI, dim);
         }

      } else
         return;

   } else {
      // On represente cte en double.
      if (fabs (cte) < std::numeric_limits<long double>::max())
         cte = Round (cte);
      matrix_row<BMat> row1(m_lat->getBasis(), j);
      matrix_row<BMat> row2(m_lat->getBasis(), i);
      ModifVect (row1, row2, -cte, dim);
      if(withDual){
         matrix_row<BMat> row3(m_lat->getDualBasis(), i);
         matrix_row<BMat> row4(m_lat->getDualBasis(), j);
         ModifVect (row3, row4, cte, dim);
      }

   }
   m_lat->setNegativeNorm (j);
   m_lat->updateVecNorm(j);
   if(withDual){
      m_lat->setDualNegativeNorm (i);
      m_lat->updateDualVecNorm(i);
   }


   miseAJourGramVD (j);
   calculCholeski2LLL (i, j);
}




//=========================================================================

#ifdef WITH_NTL

void Reducer::redLLLNTLExact(double fact){

#if NTL_TYPES_CODE > 1

   bool withDual = m_lat->withDual();
   if (withDual) {
      mat_ZZ U;
      U.SetDims(m_lat->getBasis().NumRows(), m_lat->getBasis().NumCols());
      ZZ det(0);
      LLL(det, m_lat->getBasis(), U, 99999, 100000);
      m_lat->getDualBasis() = transpose(inv(U)) * m_lat->getDualBasis();
   } else{
      ZZ det(0);
      LLL(det, m_lat->getBasis(), 99999, 100000);
   }

#else

   cout << "**** WARNING redLLLNTLExact cannot be use with Double Precision (LLDD)\n";
   cout << "**** in NTL nor without NTL. LLL reduction is used with our algorithm \n";
   cout << "which can be slower.";
   redLLL(fact, 1000000, m_lat->getDim ());

#endif
}



//=========================================================================

void Reducer::redBKZ(double fact, long blocksize, PrecisionType precision, int dim) {

   IntLatticeBasis *lattmp = 0;
   if(dim >0){
      lattmp = new IntLatticeBasis(dim, m_lat->getNorm());
      lattmp->copyBasis(*m_lat, dim);
   }
   else
      lattmp = m_lat;

#if NTL_TYPES_CODE > 1

   bool withDual = m_lat->withDual();
   mat_ZZ U;
   U.SetDims(lattmp->getBasis().size1(), lattmp->getBasis().size2());

   switch(precision){
      case DOUBLE:
         BKZ_FP(lattmp->getBasis(), U, fact, blocksize);
         break;
      case QUADRUPLE:
         BKZ_QP(lattmp->getBasis(), U, fact, blocksize);
         break;
      case EXPONENT:
         BKZ_XD(lattmp->getBasis(), U, fact, blocksize);
         break;
      case ARBITRARY:
         BKZ_RR(lattmp->getBasis(), U, fact, blocksize);
         break;
      default:
         MyExit(1, "BKZ PrecisionType:   NO SUCH CASE");
   }

   if (withDual)
      lattmp->getDualBasis() = transpose(inv(U)) * lattmp->getDualBasis();

   m_lat->copyBasis(*lattmp, dim);

#else

   cout << "\n**** WARNING redBKZ cannot be use with long type for integers\n";
   cout << "**** (LLDD) and requires ZZ type. Instead, LLL reduction is\n";
   cout << "**** performed with our algorithm which can be much slower.\n";
   cout << endl;
   redLLL(fact, 1000000, m_lat->getDim ());

#endif

   if (dim>0)
      delete lattmp;

}


//=========================================================================

void Reducer::redLLLNTL(double fact, PrecisionType precision, int dim) {

   IntLatticeBasis *lattmp = 0;
   if(dim >0){
      lattmp = new IntLatticeBasis(dim, m_lat->getNorm());
      lattmp->copyBasis(*m_lat, dim);
   }
   else
      lattmp = m_lat;

#if NTL_TYPES_CODE > 1

   if(precision == EXACT)
      redLLLNTLExact(fact);
   else{
      bool withDual = m_lat->withDual();
      mat_ZZ U;
      U.SetDims(lattmp->getBasis().NumRows(), lattmp->getBasis().NumCols());

      switch(precision){
         case DOUBLE:
            LLL_FP(lattmp->getBasis(), U, fact, 0, 0);
            break;
         case QUADRUPLE:
            LLL_QP(lattmp->getBasis(), U, fact, 0, 0);
            break;
         case EXPONENT:
            LLL_XD(lattmp->getBasis(), U, fact, 0, 0);
            break;
         case ARBITRARY:
            LLL_RR(lattmp->getBasis(), U, fact, 0, 0);
            break;
         case EXACT:
            break;
         default:
            MyExit(1, "LLL PrecisionType:   NO SUCH CASE");
      }

      if (withDual)
         lattmp->getDualBasis() = transpose(inv(U)) * lattmp->getDualBasis();

      m_lat->copyBasis(*lattmp, dim);

   }

#else

   cout << "\n**** WARNING redLLLNTL cannot be use with long type for integers\n";
   cout << "**** (LLDD) and requires ZZ type. Instead, LLL reduction is performed\n";
   cout << "**** with our algorithm which can be much slower.\n";
   cout << endl;
   redLLL(fact, 1000000, m_lat->getDim ());

#endif

   if (dim>0)
      delete lattmp;
}

//=========================================================================

#endif

//=========================================================================

void Reducer::redLLL (double fact, long maxcpt, int Max)
/*
 * Effectue la pre-reduction de B au sens de Lenstra-Lenstra-Lovasz. N'utilise
 * pas les vecteurs m_lat->getBasis().vecNorm, Wm_lat->getDualBasis().
 */
{
   bool withDual = m_lat->withDual();
   const int REDBAS_e = 40;
   int i, j, k, h;
   RScal Cho0ij;
   RScal limite;
   long cpt;

   const int dim = m_lat->getDim ();

   cpt = 0;
   calculGramVD ();
   limite = 1.0;
   for (k = 1; k <= REDBAS_e; k++)
      limite *= 2.0;
   limite *= dim;
   m_cho2(0,0) = m_gramVD(0,0);
   m_cho2(0,1) = m_gramVD(0,1);
   m_IC[0] = 1;
   m_cho2(1,1) = m_gramVD(1,1) - m_cho2(0,1) * (m_cho2(0,1) / m_cho2(0,0));
   m_IC[1] = 1;
   for (i = 2; i < dim; i++)
      m_IC[i] = -1;
   h = 0;

   while (h < Max-1 && cpt < maxcpt) {
      if (m_gramVD(h + 1,h + 1) > limite) {
         for (i = h; i >= 0; i--)
            reductionFaible (i, h + 1);
      } else
         reductionFaible (h, h + 1);

      calculCholeski2Ele (h + 1, h + 1);
      if (m_IC[h + 1] == -1)
         m_IC[h + 1] = h + 1;

      if (m_cho2(h + 1,h + 1)/m_cho2(h,h) + (m_cho2(h,h + 1)/m_cho2(h,h))
            * (m_cho2(h,h + 1) / m_cho2(h,h)) < fact) {
         ++cpt;

         m_lat->permute (h, h + 1);
         permuteGramVD (h, h + 1, dim);
         m_cho2(h,h) = m_gramVD(h,h);
         for (i = 0; i < h; i++) {
            swap(m_cho2(i,h), m_cho2(i,h + 1));
            m_cho2(h,h) -= m_cho2(i,h) * (m_cho2(i,h) / m_cho2(i,i));
         }
         if (h == 0) {
            Cho0ij = m_cho2(0,1) / m_cho2(0,0);
            if (fabs (Cho0ij) > 0.5) {
               m_IC[0] = 1;
               m_IC[1] = -1;
               h = 0;
            } else {
               m_cho2(1,1) = m_gramVD(1,1) -
                  m_cho2(0,1) * m_cho2(0,1) / m_cho2(0,0);
               calculCholeski2LLL (2, 2);
               m_IC[0] = 2;
               m_IC[1] = 2;
               m_IC[2] = 2;
               h = 1;
            }
         } else {
            m_IC[h] = h + 1;
            m_IC[h + 1] = -1;
            --h;
         }

      } else {
         for (i = 0; i <= h + 2; i++) {
            if (h + 2 > m_IC[i]) {
               if (h + 2 < dim)
                  calculCholeski2Ele (i, h + 2);
               m_IC[i] = h + 2;
            }
         }
         ++h;
      }
   }

   if(cpt == maxcpt){
      cout << "***** in redLLL cpt > maxcpt = " << maxcpt << endl;
   }

   for (j = 2; j < Max; j++) {
      for (i = j - 2; i >= 0; i--)
         reductionFaible (i, j);
   }
   m_lat->setNegativeNorm ();
   if(withDual){
      m_lat->setDualNegativeNorm ();
   }
}


//=========================================================================

void Reducer::transformStage3 (std::vector<long> & z, int & k)
{
   int j, i;
   long q;
   const int dim = m_lat->getDim ();
   bool withDual = m_lat->withDual();

   j = dim-1;
   while (z[j] == 0)
      --j;
   while (labs (z[j]) > 1) {
      i = j - 1;
      while (z[i] == 0)
         --i;
      // On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
      while (z[j]) {
         // Troncature du quotient vers 0
         q = z[i] / z[j];
         if (q) {
            // On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
            z[i] -= q * z[j];
            matrix_row<BMat> row2(m_lat->getBasis(), i);
            matrix_row<BMat> row1(m_lat->getBasis(), j);
            //    ModifVect (m_lat->getBasis ()[j], m_lat->getBasis ()[i],
            //            q, dim);
            ModifVect (row1, row2, q, dim);
            m_lat->setNegativeNorm (j);

            if(withDual){
               matrix_row<BMat> row3(m_lat->getDualBasis(), i);
               matrix_row<BMat> row4(m_lat->getDualBasis(), j);
               //    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
               //             -q, dim);
               ModifVect (row3, row4, -q, dim);
               m_lat->setDualNegativeNorm (i);
            }
         }
         // Permutation.
         swap <long>(z[i], z[j]);
         m_lat->permute (i, j);
      }
      j = i;
   }
   k = j;
}


//=========================================================================

bool Reducer::tryZ (int j, int i, int Stage, bool & smaller, const BMat & WTemp)
// Si m_countNodes > MaxNodes retourne false, sinon retourne true.
{
   long max0, min0;
   RScal x, dc;
   RScal center;
   long zhigh, zlow, h;
   bool high;
   int k;
   RScal S1, S2, S3, S4, mR;

   const int dim = m_lat->getDim ();
   conv (mR, m_lat->getModulo());


   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      cout << "-------- m_countNodes > maxNodesBB = " << maxNodesBB << endl;
      return false;
   }


   // Calcul d'un intervalle contenant les valeurs admissibles de zj.
   center = 0.0;
   if (j < dim-1) {
      // Calcul du centre de l'intervalle.
      for (k = j + 1; k < dim; k++)
         center = center - m_c0(j,k) * m_zLR[k];

      // Distance du centre aux extremites de l'intervalle.
      dc = sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j]);

      /* Calcul de deux entiers ayant la propriete qu'un entier */
      /* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
      /* n'appartient pas a l'intervalle.  */
      x = center - dc;
      conv (min0, trunc (x));
      if (x > 0.0)
         ++min0;

      x = center + dc;
      conv (max0, trunc (x));
      if (x < 0.0)
         --max0;

      // En vue du choix initial de zj. On determine zlow et zhigh.
      if (min0 > max0){
         return true;
      }
      if (min0 == max0) {
         zlow = min0;
         zhigh = max0 + 1;
         high = false;
      } else if (center >= 0.0) {
         conv (zlow, trunc (center));
         zhigh = zlow + 1;
         conv (h, trunc (2.0 * center));
         high = h & 1;
      } else {
         conv (zhigh, trunc (center));
         zlow = zhigh - 1;
         conv (h, -trunc (2.0 * center));
         high = (h & 1) == 0;
      }

   } else {  // j = dim-1

      zlow = 0;
      high = true;
      if (Stage == 2) {
         min0 = 1;
         max0 = 1;
         zhigh = 1;
      } else {
         min0 = 2;
         zhigh = 2;
         conv (max0, trunc (sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j])));
      }
   }


   NScal temp;
   /* On essaie maintenant chacun des z[j] dans l'intervalle, en */
   /* commencant par le centre puis en alternant d'un cote a l'autre. */
   while (zlow >= min0 || zhigh <= max0) {

      if (high) {
         m_zLI[j] = zhigh;
      } else {
         m_zLI[j] = zlow;
      }
      m_zLR[j] = m_zLI[j];

      // Calcul de m_n2[j-1].
      x = m_zLR[j] - center;

      if (j == 0) {
         RScal tmps_n2 = m_n2[0] + x * x * m_dc2[0];
         if (tmps_n2 < m_lMin2) {
            // On verifie si on a vraiment trouve un vecteur plus court
            matrix_row<const BMat> row1(m_lat->getBasis(), dim-1);
            m_bv = row1;
            for (k = 0; k < dim-1; k++) {
               if (m_zLI[k] != 0) {
                  matrix_row<const BMat> row1(m_lat->getBasis(), k);
                  ModifVect (m_bv, row1, m_zLI[k], dim);
               }
            }
            if (Stage == 3) {
               matrix_row<const BMat> row1(m_lat->getBasis(), dim-1);
               ModifVect (m_bv, row1, m_zLR[dim-1] - 1.0, dim);
            }

            ProdScal (m_bv, m_bv, dim, S1);
            conv (S4, m_lat->getVecNorm (dim-1));
            if (S1 < S4) {
               if (Stage == 2) {
                  smaller = true;
                  if (!PreRedLLLRM)
                     m_zShort = m_zLI;
                  else {
                     for (k = 1; k < dim; k++) {
                        matrix_row<const BMat> row1(WTemp, k);
                        ProdScal (m_bv, row1, dim, S2);
                        Quotient (S2, mR, S3);
                        conv (m_zShort[k], S3);
                     }
                     m_zShort[dim-1] = 1;
                  }
               } else if (Stage == 3 && !PreRedLLLRM) {
                  if (GCD2vect (m_zLI, i, dim) == 1) {
                     m_zShort = m_zLI;
                     smaller = true;
                  }
               } else {
                  for (k = 0; k < dim; k++) {
                     matrix_row<const BMat> row1(WTemp, k);
                     ProdScal (m_bv, row1, dim, S2);
                     Quotient (S2, mR, S3);
                     conv (m_zShort[k], S3);
                  }
                  if (GCD2vect (m_zShort, i, dim) == 1) {
                     smaller = true;
                  }
               }
               if (smaller) {
                  conv (temp, S1);
                  m_lat->setVecNorm (temp, dim-1);
                  return true;
               }
            }
         }
      } else { // j > 0
         m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
         if (m_lMin2 >= m_n2[j - 1]) {
            if (!tryZ (j - 1, i, Stage, smaller, WTemp))
               return false;
         // Des qu'on a trouve quelque chose, on sort de la recursion
         // et on retourne dans reductMinkowski.
            if (smaller)
               return true;
         }
     }
      if (high) {
         ++zhigh;
         if (zlow >= min0)
            high = false;
      } else {
         --zlow;
         if (zhigh <= max0)
            high = true;
      }
   }
   return true;
}


//=========================================================================

bool Reducer::redBB (int i, int d, int Stage, bool & smaller)
/*
 * Tries to shorten m_lat->getBasis()[i] using branch-and-bound.
 * Used in Minkowski Reduction.
 * Stage is 2 or 3.
 * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
 * Stops and returns false if not finished after examining MaxNodes
 * nodes in the branch-and-bound tree.  When succeeds, returns true.
 * Assumes that the norm is Euclidean.
 */
{
   bool withDual = m_lat->withDual();
   const int dim = m_lat->getDim ();
   BMat VTemp (dim, dim), WTemp (dim, dim);
   bool XXTemp[dim];
   NScal tmp;
// trace( "AVANT redBB");
   smaller = false;

   // Approximation du carre de la longueur de Vi.
   if (m_lat->getVecNorm(i) < 0) {
      matrix_row<BMat> row1(m_lat->getBasis(), i);
      ProdScal (row1, row1, dim, tmp);
      //  ProdScal (m_lat->getBasis()[i], m_lat->getBasis()[i],
      //            dim, tmp);
      m_lat->setVecNorm (tmp, i);
   }
   conv (m_lMin2, m_lat->getVecNorm (i));

   if (Stage == 3 && withDual)
   {
      if (m_lat->getDualVecNorm (i) < 0) {
         matrix_row<BMat> row1(m_lat->getDualBasis(), i);
         ProdScal (row1, row1, dim, tmp);
         //   ProdScal (m_lat->getDualBasis()[i], m_lat->getDualBasis()[i],
         //            dim, tmp);
         m_lat->setDualVecNorm (tmp, i);
      }
      NScal m2;
      conv(m2, m_lat->getModulo());
      m2 = m2*m2;
      if (m_lMin2 * m_lat->getDualVecNorm (i) < 4*m2)
         return true; // if the angle between the basis vector i and the dual
                      // basis vector i is between -Pi/3 and Pi/3
   }
   if(withDual){
      m_lat->updateDualVecNorm ();
   }
   m_lat->updateVecNorm ();
   m_lat->permute (i, dim-1);

   int k, h;

   if (PreRedLLLRM)
   {
      // On memorise la base courante.
      VTemp = m_lat->getBasis ();
      if(withDual){
         WTemp = m_lat->getDualBasis ();
      }
      for (h = 0; h < dim; h++)
         XXTemp[h] = true;
      redLLL (1.0, 1000000, dim - 1);
      m_lat->updateVecNorm ();
      if(withDual){
         m_lat->updateDualVecNorm ();
      }
   }
   if (!calculCholeski (m_dc2, m_c0))
      return false;
   m_countNodes = 0;
   m_n2[dim-1] = 0.0;
   if (!tryZ (dim-1, i, Stage, smaller, WTemp))
      return false;

   if (PreRedLLLRM)
   {
      /* On remet l'anciennne base, celle d'avant LLL, avant de considerer
         la prochaine m_lat->dimension.  */
      m_lat->getBasis () = VTemp;
      m_lat->updateVecNorm ();
      if(withDual){
         m_lat->getDualBasis () = WTemp;
         m_lat->updateDualVecNorm ();
      }
      for (h = 0; h < dim; h++)
         m_lat->setXX (XXTemp[h], h);
   }
   if (smaller)
   {
      /* On a trouve un plus court vecteur.  On ameliore
         m_lat->getBasis()[k].  */
      if (Stage == 2)
         k = dim-1;
      else
         transformStage3 (m_zShort, k);
      matrix_row<BMat> row1(m_lat->getBasis(), k);
      for (h = 0; h < dim; h++)
         row1(h) = m_bv[h];
      //  m_lat->getBasis ()[k] = m_bv;
      m_lat->setNegativeNorm (k);
      if (m_zShort[k] < 0 && withDual) {
         matrix_row<BMat> row2(m_lat->getDualBasis(), k);
         ChangeSign (row2, dim);
      }
      /* Mise a jour des vecteurs de la base duale selon le nouveau
         m_lat->getBasis()[k] */
      for (h = 0; h < dim; h++) {
         if ((m_zShort[h] != 0) && (h != k)) {
            if(withDual){
               matrix_row<BMat> row1(m_lat->getDualBasis(), h);
               matrix_row<BMat> row2(m_lat->getDualBasis(), k);
               ModifVect (row1, row2, -m_zShort[h], dim);
               m_lat->setDualNegativeNorm (h);
            }
            if (Stage == 2) {
               if (h >= d)
                  m_lat->setXX (false, h);
            }
         }
      }
   } else if (Stage == 2)
       m_lat->setXX(true, dim-1);

   m_lat->permute (i, dim-1);
// trace( "APRES redBB");
   return true;
}


//=========================================================================

bool Reducer::tryZ0 (int j, bool & smaller)
/*
 * Procedure recursive implantant le BB pour shortestVector (redBB0).
 * Si m_countNodes > MaxNodes, retourne false.
 * Sinon, retourne true: on a trouve le plus court vecteur.
 * Le carre de la longueur du plus court vecteur a date est dans m_lMin2.
 * Si j=1 et on trouve un plus court vecteur, on retrouve les zj
 * correspondants dans m_zShort.
 */
{
   /* Note: Pour une implantation non recursive, ces variables devraient
      etre des tableaux indices par j. */
   RScal dc, x, center;
   long min0, max0;            // Bornes de l'intervalle pour les z_j.
   long zlow, zhigh;           // Valeur courante examinee a gauche et a
                               // droite du centre.
   bool high;                  // Indique si on est a d. ou g. du centre.
   int k;
   long temp;

   const int dim = m_lat->getDim ();

   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      cout << "*****   m_countNodes > maxNodesBB = " << maxNodesBB << endl;
      return false;
   }
   /* Calcul d'un intervalle contenant les valeurs admissibles de zj. */
   /* Calcul du centre de l'intervalle.  */
   center = 0.0;
   for (k = j + 1; k < dim; ++k)
      center -= m_c0(j,k) * m_zLR[k];

   // Distance du centre aux extremites de l'intervalle.

   dc = sqrt ((m_lMin2 - m_n2[j]) / m_dc2[j]);

   /* Calcul de deux entiers ayant la propriete qu'un entier */
   /* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
   /* n'appartient pas a l'intervalle.  */
   /* Si NOT m_foundZero on pose min egal a zero par symetrie.  */
   if (!m_foundZero)
      min0 = 0;
   else {
      x = center - dc;
      conv (min0, trunc (x));
      if (x > 0.0)
         ++min0;
   }
   x = center + dc;
   conv (max0, trunc (x));
   if (x < 0.0)
      --max0;

   // En vue du choix initial de zj. On determine zlow et zhigh.
   if (min0 > max0)
      return true;
   if (min0 == max0) {
      zlow = min0;
      zhigh = max0 + 1;
      high = false;
   } else if (center >= 0.0) {
      conv (zlow, trunc (center));
      zhigh = zlow + 1;
      conv (temp, trunc (2.0 * center));
      high = (temp & 1);
   } else {
      conv (zhigh, trunc (center));
      zlow = zhigh - 1;
      conv (temp, -trunc (2.0 * center));
      high = (temp & 1) == 0;
   }

   // On essaie maintenant chacun des z[j] dans l'intervalle, en
   // commencant par le centre puis en alternant d'un cote a l'autre.
   while (zlow >= min0 || zhigh <= max0) {
      if (high)
         m_zLI[j] = zhigh;
      else
         m_zLI[j] = zlow;
      m_zLR[j] = m_zLI[j];

      // Calcul de m_n2[j-1], qui est normalise dans le cas SISquares.
      x = m_zLR[j] - center;


      if (j == 0) {
         // Tous les zj sont choisis: on a un vecteur a tester.
         if (m_lMin2 > m_n2[0] + x * x * m_dc2[0]) {
            /* On verifie si on a vraiment trouve un vecteur */
            /* non-nul et plus court.  */
            if (!m_foundZero) {
               // Le premier vecteur trouve sera zero.
               m_foundZero = true;
            } else {
               SetZero (m_bv, dim);
               for (k = 0; k < dim; k++) {
                  if (m_zLI[k] != 0) {
                     matrix_row<BMat> row1(m_lat->getBasis(), k);
                     ModifVect (m_bv, row1, m_zLI[k], dim);
                  }
               }
               // Le nouveau vecteur trouve est maintenant dans m_bv.
               if (m_lat->getNorm () == L2NORM) {
                  ProdScal (m_bv, m_bv, dim, x);
               } else {
                  CalcNorm <BVect, RScal> (m_bv, dim, x, m_lat->getNorm ());
                  x = x * x;
               }
               if (x < m_lMin2) {
    /* La condition suivante ralentit le programme; je l'ai mise donc en
       commentaire. Il est très rare qu'elle soit effective et permette de sortir
       prématurément, et seulement pour de grandes dimensions dans le cas L2NORM.
       Mais on doit la tester à chaque passage ici, et la longueur des candidats
       plus courts diminue lentement.
       Il se pourrait que ce soit parfois plus rapide dans le cas L1NORM,
       dépendant de la dimension. Mais L2NORM est primordial. */
                  /*
                  if (x <= m_BoundL2[dim]) {
                     return false;
                  }
                  */

                  // Il est plus court!
                  smaller = true;
                  conv (m_lMin2, x);
                  m_zShort = m_zLI;
                  m_bw = m_bv;
               }
            }
         }
      } else if (m_lMin2 > m_n2[j] + x * x * m_dc2[j]) {
         // Encore de l'espoir: recursion.
         m_n2[j-1] = m_n2[j] + x * x * m_dc2[j];
         if (!tryZ0 (j-1, smaller))
            return false;
      } else
         m_n2[j-1] = m_n2[j] + x * x * m_dc2[j];

      if (high) {
         ++zhigh;
         if (zlow >= min0)
            high = false;
      } else {
         --zlow;
         if (zhigh <= max0)
            high = true;
      }
   }
   return true;
}


//=========================================================================

bool Reducer::redBB0 (NormType norm)
/*
 * Finds shortest non-zero vector, using branch-and-bound. Stops and returns
 * false if not finished after examining MaxNodes nodes in the
 * branch-and-bound tree.  When succeeds, returns true, and the shortest
 * vector square length will be in m_lMin2.
 */
{
      
   bool withDual = m_lat->withDual();
   bool smaller = false;
   int k, h;
   const int dim = m_lat->getDim ();

   RScal x;
   /* We sort here to get same results as in xds98 version. Otherwise,
      Cholesky will fail more rapidly due to floating-point errors. We do
      not sort after redLLL because this greatly slows down the program in
      the case of the L2 Norm. */

   m_lat->updateScalL2Norm (0, dim);
   if(withDual){
      m_lat->updateDualScalL2Norm (0, dim);
   }
   if (m_countNodes < SHORT_LLL)
      m_lat->sort (0);

   /* Approximation de la norme du plus court vecteur. */
   if (norm == L2NORM) {
      conv (m_lMin2, m_lat->getVecNorm (0));
   } else {
      // Looking for the min of of the vecteur according to the considered norm
      matrix_row<BMat> row1(m_lat->getBasis(), 0);
      CalcNorm <BVect, NScal> (row1, dim, m_lMin, norm);
      for (k = 1; k < dim; k++) {
         matrix_row<BMat> row2(m_lat->getBasis(), k);
         CalcNorm <BVect, NScal> (row2, dim, x, norm);
         if (x < m_lMin)
            m_lMin = x;
      }
      m_lMin2 = m_lMin * m_lMin;
   }

   if (m_lMin2 <= m_BoundL2[dim-1]) {

      /* If there is in this lattice a vector that are shorter than the shortest
         of the current lattice, we don't need to find the shortest vector in this
         lattice. Use only in seek function

         S'il existe dans ce réseau un vecteur de longueur plus courte que celle
         du meilleur réseau trouvé à date, il n'est pas nécessaire de trouver le
         plus court vecteur de ce réseau: on peut l'éliminer immédiatement. */

      return false;
   }

   if (!calculCholeski (m_dc2, m_c0))
      return false;

   /* On effectue le branch and bound.  */
   /* m_n2[j] sera la somme des termes |z*k|^2 ||v*k||^2 pour k > j.  */
   m_n2[dim-1] = 0.0;
   m_countNodes = 0;
   smaller = false;
   m_foundZero = false;
   if (!tryZ0 (dim-1, smaller))
      return false;

   if (smaller) {
      // We have found a shorter vector. Its square length is in m_lMin2.
      transformStage3 (m_zShort, k);
      matrix_row<BMat> row1(m_lat->getBasis(), k);
      for (h = 0; h < dim; h++)
         row1(h) = m_bw[h];
      m_lat->setNegativeNorm (k);

      /* The new shorter vector is now in m_lat->getBasis()[k].  */
      /*  We update the vectors of the dual basis. */
      if (withDual){
         if (m_zShort[k] < 0) {
            matrix_row<BMat> row2(m_lat->getDualBasis(), k);
            ChangeSign (row2, dim);
         }
         for (h = 0; h < dim; h++) {
            if (m_zShort[h] && h != k) {
               matrix_row<BMat> row3(m_lat->getDualBasis(), h);
               matrix_row<BMat> row4(m_lat->getDualBasis(), k);
               ModifVect (row3, row4, -m_zShort[h], dim);
               m_lat->setDualNegativeNorm (h);
            }
         }
      }

      /* The new candidate for a shortest vector will be in
         m_lat->getBasis()(0). */
      /* In the case of L1NORM or others, check that it is really smaller.  */
      if (norm == L2NORM)
         m_lat->permute (k, 0);
      else {
         matrix_row<BMat> row5(m_lat->getBasis(), k);
         CalcNorm (row5, dim, x, norm);
         if (x < m_lMin) {
            m_lMin = x;
            m_lMin2 = m_lMin * m_lMin;
            m_lat->permute (k, 0);
         }
      }
   }
   return true;
}




//=========================================================================

bool Reducer::reductMinkowski (int d)
{
   bool withDual = m_lat->withDual();
   const int dim = m_lat->getDim ();
   int i;
   long totalNodes = 0;
   bool found;
   bool smaller;               // A smaller vector has been found

   //redBKZ(0.999999, 10);
   do {
      // The first d vectors should not be modified.
      for (i = 0; i < d; i++)
         m_lat->setXX (true, i);
      for (i = d; i < dim; i++)
         m_lat->setXX (false, i);


      found = false; 

      do {
         preRedDieter (d);

         m_lat->setNegativeNorm(d);
         m_lat->updateVecNorm (d);
         if(withDual){
            m_lat->setDualNegativeNorm(d);
            m_lat->updateDualVecNorm (d);
         }
         m_lat->sort (d);
         found = false;

         for (i = 0; i < dim; i++) {
            if (!m_lat->getXX (i)) {
               // On essaie de reduire le i-ieme vecteur.
               if (!redBB (i, d, 2, smaller))
                  return false;
               totalNodes += m_countNodes;
               if (smaller){
                  found = true;
               }
            }
         }

      } while (found);
      // Stage 3
      if (dim > 7) {
         for (i = d; i < dim; i++) {
            if (!redBB (i, d, 3, smaller))
               return false;
            totalNodes += m_countNodes;
            if (smaller)
               found = true;
         }
      }
   } while (found);

   /*
   if (totalNodes > MINK_LLL) {
      PreRedLLLRM = true;
   }
   */
   
   m_lat->setNegativeNorm();
   m_lat->updateScalL2Norm (0, dim);
   if(withDual){
      m_lat->setDualNegativeNorm();
      m_lat->updateDualScalL2Norm (0, dim);
   }
   m_lMin2 = m_lat->getVecNorm(0);
   return true;
}


//=========================================================================

bool Reducer::shortestVector (NormType norm)
// Square length of shortest vector can be recovered in m_lMin2
{

   if (norm != L2NORM) {
      m_lat->setNegativeNorm ();
      if(m_lat->withDual())
         m_lat->setDualNegativeNorm ();
   }

   /* Find the shortest vector for the selected norm.  */
   /* The L2 norm is used for the Choleski decomposition and BB bounds. */
   bool ok;
   if (norm == L1NORM || norm == L2NORM || norm == ZAREMBANORM) {
      ok = redBB0 (norm);
   } else {
      ok = false;
      cerr << "RedLattice::shortestVector:   wrong norm";
      exit (3);
   }

   m_lat->updateVecNorm();
   m_lat->sort(0);
   if(m_lat->withDual())
      m_lat->updateDualVecNorm();

   return ok;
}

//=========================================================================

bool Reducer::redDieter (NormType norm)
/*
 * Finds shortest non-zero vector with specified norm.
 * Does not modify the basis.
 * Stops and returns false if not finished after examining possib.
 * When succeeds, returns true, and returns length in m_lMin.
 * Uses the algorithm of Dieter (1975), given in Knuth (1981).
 */
{
   if(!m_lat->withDual()){
      cout << "***************************************************\n";
      cout << "Dieter Reduction needs the Dual. Evaluation Aborted\n";
      cout << "***************************************************\n" << endl;
      return false;
   }
   const int dim = m_lat->getDim ();
   int k;
   RScal x;
   NVect supz, z;

   supz.resize (dim);
   z.resize (dim);

   /* Compute the bounds in eq.(7) of Dieter (1975). */
   /* For this, VV and WW temporarily hold the vector norms. */
   m_lat->setNegativeNorm();
   m_lat->setDualNegativeNorm();
   m_lat->updateVecNorm ();
   m_lat->updateDualVecNorm ();
   m_lat->sort (0);

   // Basis is now sorted according to the desired norm.
   conv (m_lMin, m_lat->getVecNorm (0));
   m_lMin2 = m_lMin * m_lMin;
   NScal temp1, temp2;
   for (k = 0; k < dim; k++)
   {
      conv (temp1, m_lMin);
      conv (temp2, m_lat->getModulo ());
      supz[k] =
         1.0E-9 + m_lat->getDualVecNorm (k) * temp1 / temp2;
      z[k] = 0;
   }

   // Algorithm given in Knuth (1981).  Memoire de F. Blouin, p. 25.
   m_countNodes = 0;
   for (k = 0; k < dim; k++)
      m_bv[k] = 0;

   k = dim-1;
   while (k >= 0)
   {
      if (z[k] < supz[k]) {
         ++z[k];
         matrix_row<BMat> row1(m_lat->getBasis(), k);
         ModifVect (m_bv, row1, 1, dim);
         while (k < dim-1) {
            ++k;
            z[k] = -supz[k];
            matrix_row<BMat> row2(m_lat->getBasis(), k);
            ModifVect (m_bv, row2, 2 * z[k], dim);
         }
         CalcNorm <BVect, RScal> (m_bv, dim, x, norm);
         if (x < m_lMin) {
            conv (m_lMin, x);
         }
      } else
         --k;
   }
   if (norm == L2NORM)
      m_lMin2 = m_lMin * m_lMin;

   m_lat->setNegativeNorm ();
   m_lat->setDualNegativeNorm ();

   supz.clear ();
   z.clear ();

   return true;
}

//=========================================================================

bool Reducer::shortestVectorDieter (NormType norm)
// Length of shortest vector can be recovered in m_lMin.
{
   // Perform pre-reductions using L2 norm.
   if (PreRedDieterSV)
      preRedDieter (0);
   if (PreRedLLLSV)
      redLLL (0.999999, 1000000, m_lat->getDim ());

   // Find the shortest vector for the selected norm.
   bool ok = redDieter (norm);
   if (m_countNodes > SHORT_DIET)
      PreRedDieterSV = true;
   if (m_countNodes > SHORT_LLL)
      PreRedLLLSV = true;
   // UpdateVVWW (0);
   // m_lat->getBasis().Sort (0);
   return ok;
}


//=========================================================================

}                                 // end namespace
