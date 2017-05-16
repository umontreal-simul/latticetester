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


#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"


#ifdef WITH_NTL
#else
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
using namespace boost::numeric::ublas;
#endif

using namespace std;

using namespace LatticeTester;


/**
 * Maximum number of transformations in the method `PreRedDieter`.
 * After <tt>MAX_PRE_RED</tt> successful transformations have been
 * performed, the prereduction is stopped.
 */
static const long MAX_PRE_RED = 1000000;


/*=========================================================================*/

/**
 * Returns the matrix of the scalar product (the Gram-Schmidt matrix) :
 * mat_gram(i,j) = basis[i] * T(basis [j])
 * Used in redLLL
 */
static RMat calculGram ( IntLatticeBasis & lat )
{
   const int dim = lat.getDim();
   RMat mat_gram(dim, dim);

   for (int i = 0; i < dim; i++) {
      for (int j = i; j < dim; j++) {
         matrix_row<BMat> row1(lat.getBasis(), i);
         matrix_row<BMat> row2(lat.getBasis(), j);
         ProdScal (row1, row2, dim, mat_gram(i,j));
         mat_gram(j,i) = mat_gram(i,j);
      }
   }
   return mat_gram;
}

/*=========================================================================*/

/**
 * Recompute the j-th line and the j-th of the Gram-Schmidt matrix
 * Used in redLLL
 */
static void majGram ( RMat & mat_gram, IntLatticeBasis & lat, const int & j)
{
   const int dim = lat.getDim();
   for (int i = 0; i < dim; i++) {
      matrix_row<BMat> row1(lat.getBasis(), i);
      matrix_row<BMat> row2(lat.getBasis(), j);
      ProdScal (row1, row2, dim, mat_gram(i,j));
      mat_gram(j,i) = mat_gram(i,j);
   }
}


/*=========================================================================*/

/**
 * update the (i,j) element of the Choleski matrix.
 * It needs the Gram-Schmidt decomposition
 */
static void calculCholeskiele (
   RMat & mat_cho,
   const RMat & mat_gram,
   const int & i,
   const int & j )
{
   mat_cho(j,i) = mat_gram(i,j);  //mat_cho(i,i);
   for (int k = 0; k < i; k++) {
      mat_cho(j,i) -= mat_cho(i,k) * mat_cho(j,k);
   }
   if(i==j){
      mat_cho(j,i) = sqrt(mat_cho(j,i));
   }
   else{
      mat_cho(j,i) /= mat_cho(i,i);
   }
}

//=========================================================================

/**
 * update the j last elements of the Choleski matrix of the j-th column.
 * It needs the Gram-Schmidt decomposition
 */
static void calculCholeskiuntiln (
   RMat & mat_cho,
   const RMat & mat_gram,
   const int & dim,
   const int & i )
{
   //mat_cho(0,j) = sqrt(mat_gram(0,j));
   calculCholeskiele(mat_cho, mat_gram, i, i);
   for (int j = 0; j < i; j++) {
      calculCholeskiele(mat_cho, mat_gram, i, j);
   }
}


//=========================================================================

/**
 * Compute the decomposition of Choleski
 */
static RMat calculCholeski (IntLatticeBasis & lat, RMat & mat_gram)
{
   const int dim = lat.getDim();
   RMat mat_cho(dim, dim);
   mat_cho(0,0) = sqrt(mat_gram(0,0));
   for(int j = 1; j<dim; j++){
      mat_cho(j,0) = mat_gram(j,0)/mat_cho(0,0);
   }
   for(int i = 1; i < dim; i++){
      mat_cho(i,i) = mat_gram(i,i);
      for(int k = 0; k<i; k++){
         mat_cho(i,i) -= mat_cho(i,k)*mat_cho(i,k);
      }
      mat_cho(i,i) = sqrt(mat_cho(i,i));
      for(int j = i+1; j<dim; j++){
         mat_cho(j,i) = mat_gram(i,j);
         for(int k = 0; k<i; k++){
            mat_cho(j,i) -= mat_cho(i,k)*mat_cho(j,k);
         }
         mat_cho(j,i) /= mat_cho(i,i);
      }
      for(int j = 0; j<i; j++){
         mat_cho(j,i) = 0;
      }
   }
   return mat_cho;
}


/*=========================================================================*/

/**
 * Performs pairwise reductions. This method tries to reduce each basis
 * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
 * adding to it a multiple of the \f$i\f$-th vector. Always uses the
 * Euclidean norm.
 */

static void pairwiseRed (
   IntLatticeBasis & lat,
   const int & i,
   const int & d,
   int & countModif,
   int & countDieter)
{
 // trace( "AVANT pairwiseRedPrimal");
   const int dim = lat.getDim ();
   ++countDieter;
   lat.setNorm(L2NORM);
   lat.updateVecNorm();

   bool modifFlag;

   NScal ns;
   BScal bs;

   for (int j = d; j < dim; j++) {
      if (i == j)
         continue;
      modifFlag = false;
      {
         matrix_row<BMat> row1(lat.getBasis(), i);
         matrix_row<BMat> row2(lat.getBasis(), j);
         //cout << "la ligne : " << row1[2] << endl;
         ProdScal (row1, row2, dim, ns);
      }
      DivideRound <NScal> (ns, lat.getVecNorm (i), ns); // donne le int le plus proche
      if (ns == 0)
         continue;
      conv (bs, ns);
      if (ns < 1000 && ns > -1000) {
         lat.updateScalL2Norm (j);
         {
            matrix_row<BMat> row1(lat.getBasis(), j);
            matrix_row<BMat> row2(lat.getBasis(), i);
            ModifVect (row1, row2, -bs, dim);
         }

         // Verify that m_lat->getPrimalBasis()[j] is really shorter
         {
         matrix_row<BMat> row1(lat.getBasis(), j);
         ProdScal (row1, row1, dim, ns);
         }
         if (ns >= lat.getVecNorm (j)) {
            matrix_row<BMat> row1(lat.getBasis(), j);
            matrix_row<BMat> row2(lat.getBasis(), i);
            ModifVect (row1, row2, bs, dim);
         } else {
            modifFlag = true;
            lat.setVecNorm (ns, j);
         }
      } else {
         matrix_row<BMat> row1(lat.getBasis(), j);
         matrix_row<BMat> row2(lat.getBasis(), i);
         ModifVect (row1, row2, -bs, dim);
         //   ModifVect (m_lat->getPrimalBasis ()[j],
         //             m_lat->getPrimalBasis ()[i], -m_bs, dim);
         //m_lat->getPrimalBasis ().setNegativeNorm (true, j);
         modifFlag = true;
      }
      // modification du dual dans Reducer.h
      if (modifFlag) {
         countDieter = 0;
         ++countModif;
      /*
         matrix_row<BMat> row1(lat.getBasis(), i);
         matrix_row<BMat> row2(lat.getBasis(), j);
         ModifVect (row1, row2, bs, dim);

         //    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
         //              m_bs, dim);
         m_lat->getDualBasis ().setNegativeNorm (true, i);
         m_lat->setXX (false, i);
         m_lat->setXX (false, j);
         */

      }

   }
 // trace( "APRES pairwiseRedPrimal");
}

//=========================================================================

/*
 * Performs pairwise reductions. This method tries to reduce each basis
 * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
 * adding to it a multiple of the \f$i\f$-th vector. Always uses the
 * Euclidean norm.
 */
static IntLatticeBasis preRedDieter (const IntLatticeBasis & lattice, int d)
{
    // trace( "AVANT preRedDieter");
   IntLatticeBasis lat(lattice);
   long BoundCount;
   long dim = lat.getDim();

   lat.updateScalL2Norm (d, dim);
   //m_lat->getDualBasis ().updateScalL2Norm (d+1, dim);
   lat.sort (d);
   int i = dim-1;
   int countModif = 0; //count the nb of modifications
   int countDieter = 0; //count the nb of the call of the fonction
   // pairwiseRed with no modification
   BoundCount = 2 * dim - d;
   do {
      pairwiseRed (lat, i, d, countModif, countDieter);
      //if (i > d)
      //   pairwiseRedDual (i);
      if (i < 1)
         i = dim-1;
      else
         --i;
   } while (!(countDieter >= BoundCount || countModif > MAX_PRE_RED));
   cout << "Nombre de modif : " << countModif << endl;
   cout << "Nombre sans modif : " << countDieter << endl;
   return lat;
}


//=========================================================================
/*
static void reductionFaible (
   IntLatticeBasis & lat,
   RMat & cho2,
   RMat & mat_gram,
   int i,
   int j)
/*
 * Reduit la matrice de Choleski (d'ordre 2) en ajoutant un multiple du
 * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
 * consequence et remet a jour la matrice des produits scalaires.
 * Utilise par redLLL.
 {
   RScal cte;
   long cteLI;
   cte = cho2(i,j) / cho2(i,i);
 // trace( "AVANT reductionFaible");

   const int dim = lat.getDim ();

   if (fabs (cte) < std::numeric_limits<double>::max()) {
      // On peut representer cte en LONGINT.
      if (fabs (cte) > 0.5) {
         conv (cteLI, Round (cte));
         matrix_row<BMat> row1(lat.getBasis(), j);
         matrix_row<BMat> row2(lat.getBasis(), i);
         ModifVect (row1, row2, -cteLI, dim);
         //  ModifVect (m_lat->getPrimalBasis ()[j], m_lat->getPrimalBasis ()[i],
         //            -cteLI, dim);

         // DUAL
         //matrix_row<Basis> row3(m_lat->getDualBasis(), i);
         //matrix_row<Basis> row4(m_lat->getDualBasis(), j);
         //ModifVect (row3, row4, cteLI, dim);
         //  ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
         //            cteLI, dim);
      } else
         return;

   } else {
      // On represente cte en double.
      if (fabs (cte) < std::numeric_limits<long double>::max())
         cte = Round (cte);
      matrix_row<BMat> row1(lat.getBasis(), j);
      matrix_row<BMat> row2(lat.getBasis(), i);
      ModifVect (row1, row2, -cte, dim);
         //      ModifVect (m_lat->getPrimalBasis ()[j], m_lat->getPrimalBasis ()[i],
         //          -cte, dim);

      // DUAL
      //matrix_row<Basis> row3(m_lat->getDualBasis(), i);
      //matrix_row<Basis> row4(m_lat->getDualBasis(), j);
      //ModifVect (row3, row4, cte, dim);
      //  ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j], cte, dim);
   }
   lat.updateVecNorm(j);
   //m_lat->getPrimalBasis ().setNegativeNorm (true, j);
   //m_lat->getDualBasis ().setNegativeNorm (true, i);

   miseAJourGramVD (j);
   calculCholeski2LLL (i, j);
 // trace( "APRES reductionFaible");
}
*/


























