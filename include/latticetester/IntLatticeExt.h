// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
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


#ifndef LATTICETESTER_INTLATTICE_H
#define LATTICETESTER_INTLATTICE_H

#include "latticetester/IntLattice.h"
/*
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaMinkL2.h"
#include "latticetester/Normalizer.h"
*/
#include "latticetester/Coordinates.h"
#include "latticetester/Lacunary.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"
#include <cassert>


namespace LatticeTester {

  /**
   * This abstract class extends `IntLattice` and is a skeleton for the
   * specialized subclasses that define specific types of lattices.
   * It is not intended to be used directly, but only via subclasses.
   * An `IntLatticeExt` object is an `IntLattice` with additional (virtual) methods
   * that must be implemented in subclasses.
   *
   * There are virtual methods to construct a basis and a dual basis,
   * to extend a basis by one coordinate, to construct the lattice defined as the
   * projection of the full lattice on a subset  \f$I\f$  of coordinates indices
   * defined by a `Coordinates` object, and to
   * recompute a basis for different numbers of dimensions and subsets of coordinates.
   * An `IntLatticeExt` object keeps a current projection of the whole lattice on a subset
   * of coordinates, as well as a basis (and perhaps its dual basis) for this current projection.
   * When computing figures of merit, this projection is frequently changed and the basis must
   * be updated.  The method `buildProjection` takes care of that.
   *
   * **REMOVE ?**
   * (I do not think we need to have this here, but only where we compute the FOMs.)
   * The lattices considered here are assumed to have a special structure, which is used
   * for the computation of the lattice density and the normalization constants in the
   * figures of merit.  It is assumed that the lattice has rank \f$k\f$ and that the
   * rescaling was done by multiplying the primal basis vectors by \f$m\f$.
   * All the lattices considered in the LatMRG and LatNet Builder software tools 
   * have this property. 
   * A lattice of rank \f$k\f$ with integer vectors modulo \f$m\f$ contains
   * \f$m^k\f$ distinct vectors (modulo $m$). If we divide the basis vectors by \f$m\f$,
   * this gives \f$m^k\f$ vectors per unit of volume, so \f$m^k\f$ is the density of the 
   * original (non-scaled) lattice. This number is used to obtain bounds on the shortest vector length,
   * which are used to normalize the shortest vector length in the spectral test.
   * This class offers methods to compute and store the constants
   * \f$ \log_2(m^{2i}) \f$ for \f$ 1 \leq i \leq k \f$,
   * which are used elsewhere to compute the normalization constants.
   */

  template<typename Int, typename Real>
  class IntLatticeExt : public IntLattice<Int, Real> {
     private:
		typedef NTL::vector<Int> IntVec;
		typedef NTL::matrix<Int> IntMat;
		typedef NTL::vector<Real> RealVec;
		typedef NTL::matrix<Real> RealMat;

     public:

          /**
           * A constructor that initializes the primal and dual bases with the
           * identity matrix. The dimension of the lattice is set to `maxDim` 
           * and the norm type is set to `norm`.
           * @param m The scaling factor `m` for the integer coordinates
           * @param maxDim The maximal dimension for which this lattice can be
           * expanded/tested
           * @param withDual Specifies whether this object contains a dual or not
           * @param norm  The norm type to measure the vector lengths.
           */
          IntLatticeExt (Int m, int maxDim, bool withDual,
              NormType norm = L2NORM);

          /**
           * Copy constructor that makes a copy of `lat`. The maximal dimension 
           * of the created basis is set equal to the current dimension in `lat`.
           */
          IntLatticeExt (const IntLatticeExt<Int, Real> & lat);

          /**
           * Makes a copy of `lattice` into this object. It copies the
           * internal vectors and matrices using the NTL assign operator =
           * (see https://libntl.org/doc/matrix.cpp.html).
           */
          void copy (const IntLatticeExt<Int, Real> & lattice);

          /**
           * Destructor. Depends on the specific subclass.
           */
          virtual ~IntLatticeExt ();

          /**
           * This returns the rank (order) of the lattice.
           */
          // int getOrder() const { return m_order; }

          /**
           * Increments the dimension of the basis and dual basis vectors by one.
           * This implementation initializes the added components to `0` and does not
           * compute the value taken by the added components and vector. It also
           * resets the norm vectors to -1. The implementation in this
           * class is meant to be overriden by subclasses.
           */
          virtual void incDim ();

          /**
           * Computes and stores the logarithm of the normalization factors
           * (<tt>m_lgVolDual2</tt>) in all dimensions up to `MaxDim`, for this
           * lattice. Here, `lgm2` must be \f$\log_g m^2\f$ and the computed values are
           * those returned by `getLgVolDual2` below.
           * These normalization contants are for the Euclidean norm only. 
           */
          // void calcLgVolDual2 (double lgm2);

          /**
           * Returns \f$\log m^{2i}\f = i \log m^2$  for \f$1\le i \le k\f$,
           * and \f$\log m^{2k}\f$ otherwise, where \f$k\f$ is the lattice rank (or order).
           */
          // double getLgVolDual2 (int i) const { return m_lgVolDual2[i]; }

          /**
           * REMOVE: This method precomputes the log of the lattice density (or of its dual),
           * as a function of the dimension.  These values are part of the normalization constants used to get
           * the normalized merit from the shortest vectors in the lattice. If
           * `dualF` is `true`, the values are computed for the m-dual
           * lattice, otherwise they are computed for the primal lattice.
           *
           * Currently, this only computes the log of m^(k/dim) or its inverse. 
           * ** Done only once in a search? **  
           */
          // void computeNormalConstants(bool dualF);
          //    void fixLatticeNormalization (bool dualF);

          /**
           * Builds the basis (and perhaps m-dual basis) for the projection `proj` for this
           * lattice. The result is placed in the `lattice` object. The LLL algorithm is 
           * applied to recover a proper basis for the primal, and the riangular method for the dual.
           * WHY ?   ***************
           */
          virtual void buildProjection (IntLatticeExt<Int, Real>* lattice,
              const Coordinates & proj);

          /**
           * This virtual method builds a basis for the lattice in `dim` dimensions.
           */
          virtual void buildBasis (int dim);

          /**
           * REMOVE: This depends on the lattice only via the density.
           * This method is used nowhere else inside lattice tester.
           * In LatMRG, it is used only once in progs/Test.h
           * In Latnet builder, it seems to be used nowhere.
           *
           * Creates and returns a Normalizer corresponding to the normalization
           * type `norma` and the density of the current lattice. 
           * The argument `alpha` = \f$\alpha\f$ is used only for the 
           * \f$P_{\alpha}\f$ measure. For all other cases, it is unused.
           * The returned Normalizer is returned but not stored in this object.
           * It contains the complete normalization constants for the number of dimensions
           * of this lattice. 
           */
          // LatticeTester::Normalizer * getNormalizer (NormaType norma,
          //    int alpha, bool dualF);

          /**
           * Selects and stores a vector of indices with lacunary values.
           */
          virtual void setLac (const Lacunary<Int> &) {};

          /**
           * Returns a string describing the lattice. 
           */
          virtual std::string toString() const;

        protected:

          /**
           * \copydoc LatticeTester::IntLattice::kill()
           *  ** USEFUL ? **
           */
          virtual void kill ();

          /**
           * Allocates space to the vectors m_vSI and m_wSI used internally to store
           * the bases for the current projection, and updates the norms.
           * This should not be called directly by the user.
           * It is called by the constructors and the copy method.
           */
          void init ();

          /**
           * REMOVE?  The order (rank) of the basis. Usually defined in subclasses.
           */
          // int m_order;

          /**
           * The maximum Dimension for the basis (for tests)
           */
          int m_maxDim;

          /**
           * A vector of normalization constants.  See `calcLgVolDual2`.
           */
          // double *m_lgVolDual2;

          /**
           * \f$\log_2 (m^2)\f$.
           */
          // double m_lgm2;

          /**
           * The m-dual basis of the current projection.
           */
          IntMat m_wSI;

          /**
           * The primal basis of the current projection.
           */
          IntMat m_vSI;

          /**
           * Working Variables used in MRGLattice.h
           * **  WHAT ARE THEY DOING HERE? MOVE THEM TO WHERE THEY BELONG.  **
           */
          // Int m_t1, m_t2, m_t3;

      }; // Class IntLatticeExt

  //===========================================================================

  template<typename Int, typename Real>
       IntLatticeExt<Int, Real>::IntLatticeExt (Int m,
          int maxDim, bool withDual, NormType norm): 
       IntLattice<Int, Real>(maxDim, norm)  {
    this->m_dim = maxDim;
    this->m_withDual = withDual;
    this->m_modulo = m;
    // m_order = k;
    init ();
    this->m_basis.resize(this->m_dim, this->m_dim);
    this->m_vecNorm.resize(this->m_dim);
    this->setNegativeNorm();
    if (withDual) {
      this->m_dualbasis.resize(this->m_dim, this->m_dim);
      this->m_dualvecNorm.resize(this->m_dim);
      this->setDualNegativeNorm();
    }
  }

  //===========================================================================

  template<typename Int, typename Real>
  IntLatticeExt<Int, Real>::IntLatticeExt (const IntLatticeExt<Int, Real> & lat):
      IntLattice<Int, Real>(lat)  {
    this->m_withDual = lat.withDual();
    // m_order = lat.m_order;
    init ();
    m_vSI = lat.m_vSI;
    if (this->m_withDual){
      this->setDualNegativeNorm();
      m_wSI = lat.m_wSI;
    }
  }

  //===========================================================================

  template<typename Int, typename Real>
  void IntLatticeExt<Int, Real>::init ()  {
      int dim = this->getDim ();
      IntLattice<Int, Real>::initVecNorm();
      // double temp;   // Used only for m_lgVolDual2.
      // NTL::conv (temp, this->m_modulo);
      m_vSI.resize(dim, dim);
      if (this->m_withDual) {
        // m_lgVolDual2 = new double[dim+1];
        // m_lgm2 = 2.0 * Lg (temp);
        // m_lgVolDual2[1] = m_lgm2;
        m_wSI.resize(dim, dim);
      }
    }

  //===========================================================================

  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::kill ()
    {
     // IntLattice<Int, Real>::kill();
      // m_vSI.clear();
    }


  //===========================================================================

  template<typename Int, typename Real>
      IntLatticeExt<Int, Real>::~IntLatticeExt ()
    {
      //sIntLatticeBase<Int, Real>::kill();
    }

  //===========================================================================

  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::copy (
          const IntLatticeExt<Int, Real> & lat) {
	  // Uses the NTL assignment operator = to make a copy of the bases.
      // m_order = lat.getOrder();
      this->m_modulo = lat.m_modulo;
      //m_m2 = lat.m_m2;
      this->m_basis = lat.m_basis;
      if (lat.withDual())
         this->m_dualbasis = lat.m_dualbasis;
      init ();  // Resizes the bases to the dimensions of the current object.
    }

  //===========================================================================

  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::incDim () {
      IntLatticeExt<Int, Real> lattmp (*this);
      int dim = this->getDim();
      this->m_basis.resize(dim+1, dim+1);
      this->m_vecNorm.resize(dim+1);
      if (this->m_withDual) {
        this->m_dualbasis.resize(dim+1, dim+1);
        this->m_dualvecNorm.resize(dim+1);
      }
      for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
          this->m_basis(i,j) = lattmp.m_basis(i,j);
          if (this->m_withDual)
            this->m_dualbasis(i,j) = lattmp.m_dualbasis(i,j);
        }
        this->m_vecNorm(i) = lattmp.m_vecNorm(i);
        if (this->m_withDual)
          this->m_dualvecNorm(i) = lattmp.m_dualvecNorm(i);
      }
      this->setNegativeNorm(dim);
      if (this->m_withDual)
        this->setDualNegativeNorm(dim);
      this->setDim(dim+1);
      return;
    }

  //===========================================================================
  /**
  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::calcLgVolDual2 (double lgm2)
    {
      if(!(this->m_withDual)) return;
      int dim = this->getDim();
      int rmax = std::min(m_order, dim);

      m_lgVolDual2[1] = lgm2;
      for (int r = 2; r <= rmax; r++)
        m_lgVolDual2[r] = m_lgVolDual2[r - 1] + lgm2;
      // WARNING [David]: one version had `m_order` instead of `rmax`.
      for (int r = rmax + 1; r <= dim; r++)
        m_lgVolDual2[r] = m_lgVolDual2[r - 1];
    }
  */

  //===========================================================================

  /**
  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::computeNormalConstants(
          bool dualF)
    {
      // Normalization factor: dual to primal : m^(k/dim) -> 1/m^(k/dim)
	  // This is the part of the normalization that depends on the lattice density.
      if (( dualF && m_lgVolDual2[1] < 0.0) ||
          (!dualF && m_lgVolDual2[1] > 0.0)) {
        for (int i = 0; i < this->getDim(); i++)
          m_lgVolDual2[i] = -m_lgVolDual2[i];
      }
      //   for (int i = 1; i <= getMaxDim(); i++)
      //      std::cout << " fix  " << m_lgVolDual2[i] << endl;
    }ss
  */

  //===========================================================================

  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::buildProjection (
          IntLatticeExt<Int, Real>* lattice, const Coordinates  & proj)
    {
      const int dim = this->getDim ();
      //  std::cout << "      ESPION_2\n";  getPrimalBasis ().write();
      int i = 0;
      IntMat temp;
      temp.SetDims(dim, dim);
      for (auto iter = proj.begin(); iter != proj.end(); ++iter) {
        for (int j = 0; j < dim; j++){
          temp(j, i) = this->m_basis(j, (*iter));
        }
        ++i;
      }
      // Constructs a basis for the projection `lattice` using the LLL construction.
      lattice->setDim (static_cast<int>(proj.size()));
      // lattice->m_order = m_order;
      BasisConstruction<Int> constr;
      constr.LLLConstruction(temp);
      temp.SetDims(lattice->getDim(), lattice->getDim());
      lattice->setNegativeNorm ();
      lattice->m_basis = temp;

      lattice->m_withDual = this->m_withDual;
      if (this->m_withDual) {
    	 // For the dual, we use the triangular method?  Why?   **************
        constr.mDualTriangular(lattice->m_basis, lattice->m_dualbasis, this->m_modulo);
        lattice->setDualNegativeNorm ();
      }

      //Triangularization<IntMat> (lattice->m_dualbasis, lattice->m_basis, dim,
      //    static_cast<int>(proj.size()), this->m_modulo);
      // lattice->trace("\nESPION_4");
      /* std::cout << "  ***** build 2\n";
         lattice->getPrimalBasis ().setNegativeNorm (true);
         lattice->getPrimalBasis ().updateScalL2Norm (1,proj.size());
         lattice->getPrimalBasis ().write();*/
      // calcDual<IntMat> (lattice->m_basis, lattice->m_dualbasis,
      //     static_cast<int>(proj.size()), this->m_modulo);
      /*
         std::cout << "  ***** build 3\n";
         lattice->getDualBasis ().setNegativeNorm (true);
         lattice->getDualBasis ().updateScalL2Norm (1,proj.size());
         lattice->getDualBasis ().write();
         */

      //lattice->updateDualScalL2Norm (0, proj.size());
      //lattice->updateScalL2Norm (0,proj.size());
      //lattice->setNegativeNorm ();
    }

  //===========================================================================

  template<typename Int, typename Real>
      void IntLatticeExt<Int, Real>::buildBasis (int d) {
	  // To be re-implemented in subclasses.
      MyExit(1, " buildBasis(d) does nothing, must be inmplemented in subclass");
      d++;  // eliminates compiler warning
    }

  //===========================================================================

  /**
  template<typename Int, typename Real>
      Normalizer<Real> * IntLatticeExt<Int, Real>::getNormalizer(
          NormaType norma, int alpha, bool dualF)
    {
      int dim = this->getDim();
      Normalizer<Real> *normal;

      Real logDensity;

      // The primal lattice density is assumed to be m^k, and m^{-k} for the dual.
      if (dualF) // dual basis 
        logDensity = - m_order * NTL::log(this->m_modulo);
      else // primal basis
        logDensity = m_order * NTL::log(this->m_modulo);

      // We create a normalizer normal for the given density, and return it.
      // This normalizer is not stored in this object.
      switch (norma) {
        case BESTLAT:
          normal = new NormaBestLat<Real> (logDensity, dim);
          break;
        case BESTBOUND:
          normal = new NormaBestBound<Real> (logDensity, dim);
          break;
        case LAMINATED:
          normal = new NormaLaminated<Real> (logDensity, dim);
          break;
        case ROGERS:
          normal = new NormaRogers<Real> (logDensity, dim);
          break;
        case MINKL1:
          normal = new NormaMinkL1<Real> (logDensity, dim);
          break;
        case MINK:
          normal = new NormaMinkL2<Real> (logDensity, dim);
          break;
        case NONE:
          normal = new Normalizer<Real> (logDensity, dim, "Norma_generic");
          break;
        default:
          std::cout << "normalizer:   no such case";
          exit (2);
      }
      return normal;
    }
    */

  //===========================================================================

  template<typename Int, typename Real>
      std::string IntLatticeExt<Int, Real>::toString() const {
      // To be re-implemented in subclasses.
	  assert (0);
      return std::string();
    }

  //===========================================================================

  extern template class IntLatticeExt<std::int64_t, double>;
  extern template class IntLatticeExt<NTL::ZZ, double>;
  extern template class IntLatticeExt<NTL::ZZ, NTL::RR>;

} // End namespace LatticeTester

#endif
