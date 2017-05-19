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

/**
 * Maximum number of Nodes in the Branch and Bound algorithm.
 */
static const long maxNodesBB = 10000000;

/**
 * Choose the type of prereduction
 */
static const bool PreRedLLLRM = false;



static void permutecol ( RMat & mat, int i, int j, int n ){
   int k = 0;
   for (k=0; k<n; k++){
      swap(mat(k,i), mat(k,j));
   }
}

static void permutelin ( RMat & mat, int i, int j, int n ){
   int k = 0;
   for (k=0; k<n; k++){
      swap(mat(i,k), mat(j,k));
   }
}


/*=========================================================================*/

/**
 * Returns the matrix of the scalar product :
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
 * Recompute the j-th line and the j-th of the matrix of the scalar product.
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
 * It needs the matrix of the scalar product.
 */
static void calculCholeskiele (
   RMat & mat_cho,
   const RMat & mat_gram,
   const int & i,
   const int & j )
{
   mat_cho(i,j) = mat_gram(i,j);
   for (int k = 0; k < i; k++) {
      mat_cho(i,j) -= mat_cho(k,i) * (mat_cho(k,j) / mat_cho(k,k));
   }
}

//=========================================================================

/**
 * update the n first elements of the Choleski matrix of the j-th column.
 * It needs the Gram-Schmidt decomposition
 */
static void calculCholeskiuntiln (
   RMat & mat_cho,
   const RMat & mat_gram,
   const int & n,
   const int & j )
{
   mat_cho(0,j) = mat_gram(0,j);
   for (int i = 1; i <= n; i++) {
      mat_cho(i,j) = mat_gram(i,j);
      for (int k = 0; k < i; k++) {
         mat_cho(i,j) -= (mat_cho(k,j) / mat_cho(k,k)) * mat_cho(k,i);
      }
   }
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
   return lat;
}


static IntLatticeBasis preRedDieter2 (const IntLatticeBasis & lattice, int d)
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
      pairwiseRed (lat, rand()%dim, d, countModif, countDieter);
      //if (rand()%dim == 0)
      if (i < 1)
         i = dim-1;
      else
         --i;
   } while (!(countDieter >= BoundCount || countModif > MAX_PRE_RED));
   return lat;
}


//=========================================================================


/**
 * Reduce the Choleski matrix with adding a multiple of the i-th vector
 * to the j-th vector. It updates the Gram Schmidt matrix
 */
static void reductionFaible (
   IntLatticeBasis & lat,
   RMat & mat_cho,
   RMat & mat_gram,
   const int & i,
   const int & j)
/*
 * Reduit la matrice de Choleski (d'ordre 2) en ajoutant un multiple du
 * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
 * consequence et remet a jour la matrice des produits scalaires.
 * Utilise par redLLL.
 */
 {
   RScal cte;
   long cteLI;
   cte = mat_cho(i,j) / mat_cho(i,i);
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
   lat.setNegativeNorm(j);
   lat.updateVecNorm(j);
   //m_lat->getPrimalBasis ().setNegativeNorm (true, j);
   //m_lat->getDualBasis ().setNegativeNorm (true, i);

   majGram (mat_gram, lat, j);
   calculCholeskiuntiln (mat_cho, mat_gram, i, j);
 // trace( "APRES reductionFaible");
}



static void redLLL (
   IntLatticeBasis & lat,
   double fact, long maxcpt, int Max)
/*
 * Effectue la pre-reduction de B au sens de Lenstra-Lenstra-Lovasz. N'utilise
 * pas les vecteurs m_lat->getPrimalBasis().vecNorm, Wm_lat->getDualBasis().
 */
{
   const int REDBAS_e = 40;
   int i, j, k, h;
   RScal Cho0ij;
   RScal limite;
   long cpt;


   const int dim = lat.getDim ();

   int *IC = new int[Max];
   IC[0] = 1;
   IC[1] = 1;
   for (i = 2; i < dim; i++)
      IC[i] = -1;

   RMat mat_cho(dim, dim);

   cpt = 0;
   RMat mat_gram = calculGram (lat);
   limite = 1.0;
   for (k = 1; k <= REDBAS_e; k++)
      limite *= 2.0;
   limite *= dim;
   mat_cho(0,0) = mat_gram(0,0);
   mat_cho(0,1) = mat_gram(0,1);


   //cout << "" << endl;
   mat_cho(1,1) = mat_gram(1,1) - mat_cho(0,1) * (mat_cho(0,1) / mat_cho(0,0));

   h = 0;
   while (h < Max-1 && cpt < maxcpt) {
      if (mat_gram(h + 1,h + 1) > limite) { // Si la norme du vecteur v(h+1,h+1) est trop élevée
         for (i = h; i >= 0; i--) // On essaye de le réduire avec une combinaison linéaire (et entière) des vecteurs précédents
            reductionFaible (lat, mat_cho, mat_gram, i, h + 1);
      } else
         // Sinon on effectue la pairwise reduction de LLL
         reductionFaible (lat, mat_cho, mat_gram, h, h + 1);

      calculCholeskiele (mat_cho, mat_gram, h + 1, h + 1);
      if (IC[h + 1] == -1)
         IC[h + 1] = h + 1;

      if (mat_cho(h + 1,h + 1)/mat_cho(h,h) + (mat_cho(h,h + 1)/mat_cho(h,h))
            * (mat_cho(h,h + 1) / mat_cho(h,h)) < fact) {
         // If the Lovasz-condition is not true, we swap the element h and h+1
         // in the basis, the matrix of scalar product and in the cholesky matrix
         ++cpt;
         lat.permute (h, h + 1);
         permutelin ( mat_gram, h, h + 1, dim);
         permutecol ( mat_gram, h, h + 1, dim);
         mat_cho(h,h) = mat_gram(h,h);

         for (i = 0; i < h; i++) {
            swap(mat_cho(i,h), mat_cho(i,h + 1));
            //mat_cho(i,0) = mat_cho(i,h);
            //mat_cho(i,h) = mat_cho(i,h + 1);
            //mat_cho(i,h + 1) = mat_cho(i,0);
            mat_cho(h,h) -= mat_cho(i,h) * (mat_cho(i,h) / mat_cho(i,i));
         }
         if (h == 0) {
            Cho0ij = mat_cho(0,1) / mat_cho(0,0);
            if (fabs (Cho0ij) > 0.5) {
               IC[0] = 1;
               IC[1] = -1;
               h = 0;
            } else {
               mat_cho(1,1) = mat_gram(1,1) -
                  mat_cho(0,1) * mat_cho(0,1) / mat_cho(0,0);
               calculCholeskiuntiln (mat_cho, mat_gram, 2, 2);
               IC[0] = 2;
               IC[1] = 2;
               IC[2] = 2;
               h = 1;
            }
         } else {
            IC[h] = h + 1;
            IC[h + 1] = -1;
            --h;
         }

      } else {
         for (i = 0; i <= h + 2; i++) {
            if (h + 2 > IC[i]) {
               if (h + 2 < dim)
                  calculCholeskiele (mat_cho, mat_gram, i, h + 2);
               IC[i] = h + 2;
            }
         }
         ++h;
      }
   }

   for (j = 2; j < Max; j++) {
      for (i = j - 2; i >= 0; i--)
         reductionFaible (lat, mat_cho, mat_gram, i, j);
   }
   lat.setNegativeNorm();
   lat.updateVecNorm();
   //m_lat->getPrimalBasis ().setNegativeNorm (true);
   //m_lat->getDualBasis ().setNegativeNorm (true);
}


//=========================================================================

/**
 * Returns in C0 the elements of the upper triangular matrix of the
 * Choleski decomposition. Returns in DC2 the squared elements of the
 * diagonal.
 */
static bool calculCholeski (IntLatticeBasis & lat, RMat & mat_gram, RVect & diag_cho, RMat & mat_cho)
{
   const int dim = lat.getDim ();

   int d, k, j, i;
   RScal m2r; // TO COMMENT
   RMat tmp_cho; // TO COMMENT
   // C2(i,j) = C0(i,j) * C2(i,i) if i != j.
   // C2(i,i) = DC2[i].
   //conv (m2r, m_lat->getM2 ());
   d = dim / 2;

   for (i = 0; i < dim; i++) {
      // Compute the i-th line of the cholesky matrix c2
      lat.updateScalL2Norm (i);
      for (j = i; j < dim; j++) {
         if (j == i)
            conv (tmp_cho(i,i), lat.getVecNorm (i));
         else {
            matrix_row<BMat> row1(lat.getBasis(), i);
            matrix_row<BMat> row2(lat.getBasis(), j);
            ProdScal (row1, row2, dim, tmp_cho(i,j));
         }
         for (k = 0; k < i; k++)
            tmp_cho(i,j) -= mat_cho(k,i) * tmp_cho(k,j);
         if (i == j) {
            diag_cho[i] = tmp_cho(i,i);
            if (diag_cho[i] < 0.0) {
               cout << "\n***** Negative diagonal element in Choleski Decomposition\n"
                    << endl;
               return false;
            }
         } else
            mat_cho(i,j) = tmp_cho(i,j) / diag_cho[i];
      }
   }
   /*
    * Previously, in LatCommun, we use the duale basis to compute the second half
    * of the Cholesky matrix. We now complete all the cholesky basis with the primal
    * basis.
    *
   // Calcul des m_lat->dim-d dernieres lignes de C0 via la base duale.
   for (i = dim-1; i > d; i--)
   {
      m_lat->getDualBasis ().updateScalL2Norm (i);
      for (j = i; j >= 0; j--) {
         if (j == i)
            conv (m_c2(i,i), m_lat->getDualBasis ().getVecNorm (i));
         else {
            matrix_row<Basis> row1(m_lat->getDualBasis(), i);
            matrix_row<Basis> row2(m_lat->getDualBasis(), j);
            ProdScal (row1, row2, dim, m_c2(i,j));
         }
         for (k = i + 1; k < dim; k++)
            m_c2(i,j) -= C0(k,i) * m_c2(k,j);
         if (i != j)
            C0(i,j) = m_c2(i,j) / m_c2(i,i);
      }

      DC2[i] = m2r / m_c2(i,i);
      if (DC2[i] < 0.0) {
            cout << "\n***** Negative diagonal element in Choleski Decomposition\n"
                 << endl;
         return false;
      }
      for (j = i + 1; j < dim; j++) {
         C0(i,j) = -C0(j,i);
         for (k = i + 1; k < j; k++) {
            C0(i,j) -= C0(k,i) * C0(k,j);
         }
      }
   }
   */
   return true;
}


//=========================================================================




/**
    * Method used in `reductMinkowski` to perform a transformation of
    * stage 3 described in \cite rAFF85a&thinsp;. Also used in
    * `ShortestVector`. Assumes that \f$\sum_{i=1}^t z_i V_i\f$ is a
    * short vector that will enter the basis. Tries to reduce some vectors
    * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
    * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
    * \f$V_j\f$ when this happens. Returns in \f$k\f$ the last index
    * \f$j\f$ such that \f$|z_j|=1\f$.
    */
static void transformStage3 (IntLatticeBasis & lat, std::vector<long> & z, int & k)
{
   int j, i;
   long q;
   const int dim = lat.getDim();

   j = dim;
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
            // On ajoute q * v[i] au vecteur m_lat->getPrimalBasis()[j]
            z[i] -= q * z[j];
            matrix_row<BMat> row2(lat.getBasis(), i);
            matrix_row<BMat> row1(lat.getBasis(), j);
            //    ModifVect (m_lat->getPrimalBasis ()[j], m_lat->getPrimalBasis ()[i],
            //            q, dim);
            ModifVect (row1, row2, q, dim);

            //matrix_row<Basis> row3(m_lat->getDualBasis(), i);
            //matrix_row<Basis> row4(m_lat->getDualBasis(), j);
            //    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
            //             -q, dim);
            //ModifVect (row3, row4, -q, dim);
            //m_lat->getPrimalBasis ().setNegativeNorm (true, j);
            lat.setNegativeNorm (i);
         }
         // Permutation.
         swap <long>(z[i], z[j]);
         lat.permute (i, j);
      }
      j = i;
   }
   k = j;
}



#if 0

//=========================================================================

static bool redBB (IntLatticeBasis & lat, int i, int d, int Stage, bool & smaller)
/*
 * Tries to shorten m_lat->getPrimalBasis()[i] using branch-and-bound.
 * Stage is 2 or 3.
 * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
 * Stops and returns false if not finished after examining MaxNodes
 * nodes in the branch-and-bound tree.  When succeeds, returns true.
 * Assumes that the norm is Euclidean.
 */
{
   const int dim = lat.getDim ();
   //const int maxDim = m_lat->getMaxDim ();
   BMat VTemp; //, WTemp (dim, dim)
   bool XXTemp[1 + dim];
   NScal tmp, lMin2;
   RMat mat_cho;
   RVect diag_cho;
// trace( "AVANT redBB");
   smaller = false;

   // Approximation du carre de la longueur de Vi.
   if (lat.getVecNorm()[i] < 0) {
      matrix_row<BMat> row1(lat.getBasis(), i);
      ProdScal (row1, row1, dim, tmp);
      //  ProdScal (m_lat->getPrimalBasis()[i], m_lat->getPrimalBasis()[i],
      //            dim, tmp);
      lat.setVecNorm (tmp, i);
   }
   conv (lMin2, lat.getVecNorm (i));

   /*
   if (Stage == 3)
   {
      if (m_lat->getDualBasis ().isNegativeNorm (i)) {
         matrix_row<Basis> row1(m_lat->getDualBasis(), i);
         ProdScal (row1, row1, dim, tmp);
         //   ProdScal (m_lat->getDualBasis()[i], m_lat->getDualBasis()[i],
         //            dim, tmp);
         m_lat->getDualBasis ().setVecNorm (tmp, i);
      }
      if (m_lMin2 * m_lat->getDualBasis ().getVecNorm (i) < 4 * m_lat->getM2 ())
         return true;
   }
   */

   lat.updateVecNorm ();
   //m_lat->getDualBasis ().updateVecNorm ();
   lat.permute (i, dim);

   int k, h;

   if (PreRedLLLRM)
   {
      // On memorise la base courante.
       VTemp.BMat::operator=(lat.getBasis());
      //WTemp = m_lat->getDualBasis ();
      //for (h = 1; h <= dim; h++)
      //   XXTemp[h] = m_lat->getXX (h);
      redLLL (lat, 1.0, 1000000, dim);
      lat.updateVecNorm ();
      //m_lat->getDualBasis ().updateVecNorm ();
   }
   if (!calculCholeski (lat, diag_cho, mat_cho))
      return false;
   m_countNodes = 0;
   m_n2[dim] = 0.0;
   if (!tryZ (dim, i, Stage, smaller, WTemp))
      return false;

   if (PreRedLLLRM)
   {
      /* On remet l'anciennne base, celle d'avant LLL, avant de considerer
         la prochaine dimension.  */
      m_lat->getPrimalBasis () = VTemp;
      m_lat->getDualBasis () = WTemp;
      m_lat->getPrimalBasis ().updateVecNorm ();
      m_lat->getDualBasis ().updateVecNorm ();
      for (h = 1; h <= dim; h++)
         m_lat->setXX (XXTemp[h], h);
   }

   if (smaller)
   {
      /* On a trouve un plus court vecteur.  On ameliore
         m_lat->getPrimalBasis()[k].  */
      if (Stage == 2)
         k = dim;
      else
         transformStage3 (m_zShort, k);
      matrix_row<Basis> row1(m_lat->getPrimalBasis(), k);
      for (h = 1; h <= dim; h++)
         row1(h) = m_bv[h];
      //  m_lat->getPrimalBasis ()[k] = m_bv;
      m_lat->getPrimalBasis ().setNegativeNorm (true, k);
      if (m_zShort[k] < 0) {
         matrix_row<Basis> row2(m_lat->getDualBasis(), k);
         ChangeSign (row2, dim);
      }
      /* Mise a jour des vecteurs de la base duale selon le nouveau
         m_lat->getPrimalBasis()[k] */
      for (h = 1; h <= dim; h++) {
         if ((m_zShort[h] != 0) && (h != k)) {
            matrix_row<Basis> row1(m_lat->getDualBasis(), h);
            matrix_row<Basis> row2(m_lat->getDualBasis(), k);
            ModifVect (row1, row2, -m_zShort[h], dim);
            m_lat->getDualBasis ().setNegativeNorm (true, h);
            if (Stage == 2) {
               if (h > d)
                  m_lat->setXX (false, h);
            }
         }
      }
   } else if (Stage == 2)
      m_lat->setXX (true, dim);

   m_lat->permute (i, dim);
// trace( "APRES redBB");
   return true;
}

#endif









































#if 0


//=========================================================================



/**
 * Tries to find shorter vectors in `reductMinkowski`.
 */
static bool tryZ (
   IntLatticeBasis & lat,
   RMat & mat_cho,
   int j,
   int i,
   int Stage,
   bool & smaller,
   const BMat & WTemp,
   long & countNodes)
// Si m_countNodes > MaxNodes retourne false, sinon retourne true.
{
   long max0, min0;
   RScal x, dc;
   RScal center;
   long zhigh, zlow, h;
   bool high;
   int k;
   RScal S1, S2, S3, S4, mR;
   NScal lMin2;
   conv (lMin2, lat.getVecNorm (i));
// trace( "AVANT tryZ");
//  cout << j << "  " << i << "  " << Stage << "  " << smaller << endl;

   const int dim = lat.getDim ();
   //conv (mR, m_lat->getM ());

   ++countNodes;
   if (countNodes > maxNodesBB) {
      cout << "*****  m_countNodes > maxNodesBB = " << maxNodesBB << endl;
      return false;
   }
   // Calcul d'un intervalle contenant les valeurs admissibles de zj.
   center = 0.0;
   if (j < dim) {
      // Calcul du centre de l'intervalle.
      for (k = j + 1; k <= dim; k++)
         center = center - mat_cho(j,k) * m_zLR[k];

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
      if (min0 > max0)
         return true;
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

   } else {                    // j = dim
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
      m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];

      if (j == 1) {
         if (m_n2[0] < m_lMin2) {
            // On verifie si on a vraiment trouve un vecteur plus court
            matrix_row<const Base> row1(m_lat->getPrimalBasis(), dim);
            m_bv = row1;
            //    m_bv = m_lat->getPrimalBasis ()[dim];
            for (k = 1; k < dim; k++) {
               if (m_zLI[k] != 0) {
                  matrix_row<const Base> row1(m_lat->getPrimalBasis(), k);
                  ModifVect (m_bv, row1, m_zLI[k], dim);
               }
            }
            if (Stage == 3) {
               matrix_row<const Base> row1(m_lat->getPrimalBasis(), dim);
               ModifVect (m_bv, row1, m_zLR[dim] - 1.0, dim);
            }

            ProdScal (m_bv, m_bv, dim, S1);
            conv (S4, m_lat->getPrimalBasis ().getVecNorm (dim));
            if (S1 < S4) {
               if (Stage == 2) {
                  smaller = true;
                  if (!PreRedLLLRM)
                     m_zShort = m_zLI;
                  else {
                     for (k = 1; k < dim; k++) {
                        matrix_row<const Base> row1(WTemp, k);
                        ProdScal (m_bv, row1, dim, S2);
                        Quotient (S2, mR, S3);
                        conv (m_zShort[k], S3);
                     }
                     m_zShort[dim] = 1;
                  }
               } else if (Stage == 3 && !PreRedLLLRM) {
                  if (GCD2vect (m_zLI, i, dim) == 1) {
                     m_zShort = m_zLI;
                     smaller = true;
                  }
               } else {
                  for (k = 1; k <= dim; k++) {
                     matrix_row<const Base> row1(WTemp, k);
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
                  m_lat->getPrimalBasis ().setVecNorm (temp, dim);
                  return true;
               }
            }
         }
      } else { // j > 1
         if (m_lMin2 >= m_n2[j - 1]) {
            if (!tryZ (j - 1, i, Stage, smaller, WTemp))
               return false;
            // Des qu'on a trouve quelque chose, on sort de la recursion */
            // et on retourne dans reductMinkowski.  */
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

 // trace( "APRES tryZ");
   return true;
}






#endif
















