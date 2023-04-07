// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
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

#include "latticetester/IntLatticeBasis.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/Normalizer.h"
#include "latticetester/Coordinates.h"
#include "latticetester/Lacunary.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"

#include <cassert>

namespace LatticeTester {

  /**
   * This class is a skeleton for the implementation of different types of 
   * lattices. This class is not really intended to be used directly, hence the
   * lack of constructor allowing the specification of a basis.
   * 
   * This class can store a lattice with or without dual and contains a few
   * virtual methods to perform common computations on lattices.
   * This class contains a method to compute lattices of projections of 
   * \f$\{x_i\}_{ 0 \leq i}\f$, a method to exchange the basis and the dual 
   * basis, and a virtual method that can be implemented in subclasses to 
   * recompute the basis for different dimensions.
   *
   * A lattice of rank \f$k\f$ with integer vectors modulo \f$m\f$ contains
   * \f$m^k\f$ distinct vectors. This number, the density, can then be used to 
   * compute bounds on the spectral test. This class implements methods to 
   * compute \f$ \log_2(m^{2i}) \f$ for \f$ 1 \leq i \leq k \f$ to help with the 
   * computation of such bounds. 
   */
  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      class IntLattice : public IntLatticeBasis<Int, BasInt, Dbl, RedDbl> {
        private:
          typedef NTL::vector<BasInt> BasIntVec;
          typedef NTL::matrix<BasInt> BasIntMat;
          typedef NTL::vector<Dbl> DblVec;
        public:

          /**
           * Constructor initializing the primal and the dual basis with the 
           * identity matrix. The dimension of the lattice is set to `maxDim` 
           * and the norm used for reduction to `norm`.
           * @param modulo The modulo of the integer coordinates
           * @param k The rank of the lattice to be constructed
           * @param maxDim The maximal dimension for which this lattice can be
           * expanded/tested (?)
           * @param withDual Specifies wether this object contains a dual or not
           * @param norm The norm to use in for reduction
           */
          IntLattice (Int modulo, int k, int maxDim, bool withDual,
              NormType norm = L2NORM);

          /**
           * Copy constructor that makes a copy of `Lat`. The maximal dimension 
           * of the created basis is set equal to `Lat`’s current dimension.
           */
          IntLattice (const IntLattice<Int, BasInt, Dbl, RedDbl> & Lat);

          /**
           * Copies `lattice` into this object. This should be equivalent to
           * the creation of a new IntLattice using the copy constructor with
           * `lattice` as an argument.
           */
          void copy (const IntLattice<Int, BasInt, Dbl, RedDbl> & lattice);

          /**
           * Destructor.
           */
          virtual ~IntLattice ();

          /**
           * Allocates space to vectors used internally. This should probably be
           * private or protected because it should not be needed to call it 
           * directly (the constructors and copy already call it).
           */
          void init ();

          /**
           * This returns the rank of the lattice.
           */
          int getOrder() const { return m_order; }

          /**
           * Increments the dimension of the basis and dual basis vectors by 
           * one. This initializes the added components to `0` and does not 
           * compute the value taken by the added components and vector. It also
           * resets vectors containing the norms. The implementation in this
           * class is meant to be overriden by subclasses and probably should 
           * not be used.
           */
          virtual void incDim ();

          /**
           * Computes the logarithm of the normalization factor
           * (<tt>m_lgVolDual2</tt>) in all dimensions \f$\leq\f$ `MaxDim` for
           * the lattice. `lgm2` is the logarithm in base 2 of \f$m^2\f$.
           */
          void calcLgVolDual2 (double lgm2);

          /**
           * Gives the log of m^(2*i) if i < order, else gives the log of m^(2*i)
           */
          double getLgVolDual2 (int i) const { return m_lgVolDual2[i]; }

          /**
           * Exchange the primal basis and the dual basis.
           */
          void dualize ();

          /**
           * This function is called to fix the normalization constants to get
           * the normalized merit from the shortest distance in the lattice. If
           * `dualF` is `true`, the normalization constant is reset for the dual
           * lattice, otherwise it is reset for the primal lattice.
           */
          void fixLatticeNormalization (bool dualF);

          /**
           * Builds the basis (and dual basis) of the projection `proj` for this
           * lattice. The result is placed in the `lattice` lattice. The basis is
           * triangularized to form a proper basis.
           */
          virtual void buildProjection (IntLattice<Int, BasInt, Dbl, RedDbl>* lattice,
              const Coordinates & proj);

          /**
           * Builds the basis for the lattice in dimension `d`. This function is
           * not implemented for this class. The general basis construction for
           * a lattice such as this one is located in BasisConstruction.
           */
          virtual void buildBasis (int d);

          /**
           * Creates and returns the normalizer corresponding to criterion
           * `norma`. In the case of the \f$P_{\alpha}\f$ test, the argument
           * `alpha` = \f$\alpha\f$. In all other cases, it is unused.
           */
          LatticeTester::Normalizer<RedDbl> * getNormalizer (NormaType norma,
              int alpha, bool dualF);

          /**
           * A utility method to store a vector of indices with lacunary values
           * in subclasses of this one. This method has no implementation in
           * this base class.
           */
          virtual void setLac (const Lacunary<BasInt> &) {};

          /**
           * Returns a string describing the lattice. 
           */
          virtual std::string toString() const;

        protected:

          /**
           * \copydoc LatticeTester::IntLatticeBasis::kill()
           */
          virtual void kill ();

          /**
           * The order of the basis.
           */
          int m_order;

          /*
           * The maximum Dimension for the test
           */
          int m_maxDim;

          /**
           * Represente sur dual along the diagonal?? ERWAN
           */
          double *m_lgVolDual2;

          /**
           * The logarithm \f$\log (m^2)\f$.
           */
          double m_lgm2;

          /**
           * The dual basis of the current projection.
           */
          BasIntMat m_wSI;

          /**
           * The primal basis of the current projection.
           */
          BasIntMat m_vSI;

          /**
           * Working Variables used in MRGLattice.h
           */
          Int m_t1, m_t2, m_t3;

      }; // Class IntLattice

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      IntLattice<Int, BasInt, Dbl, RedDbl>::IntLattice ( Int modulo, int k,
          int maxDim, bool withDual, NormType norm): 
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>(maxDim, norm)
  {
    this->m_dim = maxDim;
    this->m_withDual = withDual;
    this->m_modulo = modulo;
    m_order = k;
    init ();
    this->m_basis.resize(this->m_dim,this->m_dim);
    this->m_vecNorm.resize(this->m_dim);
    this->setNegativeNorm();
    if (withDual) {
      this->m_dualbasis.resize(this->m_dim,this->m_dim);
      this->m_dualvecNorm.resize(this->m_dim);
      this->setDualNegativeNorm();
    }
  }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      IntLattice<Int, BasInt, Dbl, RedDbl>::IntLattice (
          const IntLattice<Int, BasInt, Dbl, RedDbl> & Lat):
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>(Lat)
  {
    this->m_withDual = Lat.withDual();
    m_order = Lat.m_order;
    init ();
    m_vSI = Lat.m_vSI;
    if (this->m_withDual){
      this->setDualNegativeNorm();
      m_wSI = Lat.m_wSI;
    }
  }

  //===========================================================================


  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::init ()
    {
      int dim = this->getDim ();
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::initVecNorm();
      double temp;
      NTL::conv (temp, this->m_modulo);
      m_vSI.resize(dim, dim);

      if (this->m_withDual) {
        m_lgVolDual2 = new double[dim+1];
        m_lgm2 = 2.0 * Lg (temp);
        m_lgVolDual2[1] = m_lgm2;
        m_wSI.resize(dim, dim);
      }

    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::kill ()
    {
      IntLatticeBasis<Int, BasInt, Dbl, RedDbl>::kill();

      if (this->m_withDual){
        if (m_lgVolDual2 == 0)
          return;
        delete [] m_lgVolDual2;
        m_lgVolDual2 = 0;
      }
      // m_vSI.clear();

    }


  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      IntLattice<Int, BasInt, Dbl, RedDbl>::~IntLattice ()
    {
      kill ();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::incDim ()
    {
      IntLattice<Int, BasInt, Dbl, RedDbl> lattmp (*this);
      int dim = this->getDim();

      // std::int64_t sizemat = m_basis.size1();
      // declared as an "unused variable" by the compiler

      this->m_basis.resize(dim+1, dim+1);
      this->m_vecNorm.resize(dim+1);

      if (this->m_withDual) {
        if(this->m_lgVolDual2 != 0)
          delete[] this->m_lgVolDual2;
        this->m_lgVolDual2 = new double[dim+2]();
        this->calcLgVolDual2 (m_lgm2);
        this->m_dualbasis.resize(dim+1, dim+1);
        this->m_dualvecNorm.resize(dim+1);
      }

      for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
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

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::calcLgVolDual2 (double lgm2)
    {
      if(!(this->m_withDual)) return;
      int dim = this->getDim();
      int rmax = std::min(m_order, dim);

      m_lgVolDual2[1] = lgm2;
      for (int r = 2; r <= rmax; r++)
        m_lgVolDual2[r] = m_lgVolDual2[r - 1] + lgm2;
      // WARNING [David]: one version had `m_order` instead of `rmax`.
      // I am not sure which is the fix and which is the bug.
      for (int r = rmax + 1; r <= dim; r++)
        m_lgVolDual2[r] = m_lgVolDual2[r - 1];
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::dualize ()
    {
      if(!(this->m_withDual)) return;
      std::swap(this->m_basis, this->m_dualbasis);
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::fixLatticeNormalization(
          bool dualF)
    {
      // Normalization factor: dual to primal : M^(k/dim) -> 1/M^(k/dim)
      if (( dualF && m_lgVolDual2[1] < 0.0) ||
          (!dualF && m_lgVolDual2[1] > 0.0)) {
        for (int i = 0; i < this->getDim(); i++)
          m_lgVolDual2[i] = -m_lgVolDual2[i];
      }
      //   for (int i = 1; i <= getMaxDim(); i++)
      //      std::cout << " fix  " << m_lgVolDual2[i] << endl;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::buildProjection (
          IntLattice<Int, BasInt, Dbl, RedDbl>* lattice, const Coordinates & proj)
    {
      const int dim = this->getDim ();
      //  std::cout << "      ESPION_2\n";  getPrimalBasis ().write();
      int i = 0;
      BasIntMat temp;
      temp.SetDims(dim, dim);
      for (auto iter = proj.begin(); iter != proj.end(); ++iter) {
        for (int j = 0; j < dim; j++){
          temp(j, i) = this->m_basis(j, (*iter));
        }
        ++i;
      }

      lattice->setDim (static_cast<int>(proj.size()));
      lattice->m_order = m_order;
      BasisConstruction<BasInt> constr;
      constr.LLLConstruction(temp);
      temp.SetDims(lattice->getDim(), lattice->getDim());
      lattice->setNegativeNorm ();
      lattice->m_basis = temp;

      lattice->m_withDual = this->m_withDual;
      if (this->m_withDual) {
        constr.DualConstruction(lattice->m_basis, lattice->m_dualbasis, this->m_modulo);
        lattice->setDualNegativeNorm ();
      }

      //Triangularization<BasIntMat> (lattice->m_dualbasis, lattice->m_basis, dim,
      //    static_cast<int>(proj.size()), this->m_modulo);
      // lattice->trace("\nESPION_4");
      /* std::cout << "  ***** build 2\n";
         lattice->getPrimalBasis ().setNegativeNorm (true);
         lattice->getPrimalBasis ().updateScalL2Norm (1,proj.size());
         lattice->getPrimalBasis ().write();*/
      // CalcDual<BasIntMat> (lattice->m_basis, lattice->m_dualbasis,
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

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::buildBasis (int d)
    {
      MyExit(1, " buildBasis does nothing");
      d++;  // eliminates compiler warning
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      void IntLattice<Int, BasInt, Dbl, RedDbl>::copy (
          const IntLattice<Int, BasInt, Dbl, RedDbl> & lat)
    {
      m_order = lat.getOrder();
      this->m_modulo = lat.m_modulo;
      //m_m2 = lat.m_m2;
      this->m_basis = lat.m_basis;
      if(lat.withDual())
        this->m_dualbasis = lat.m_dualbasis;
      init ();
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      Normalizer<RedDbl> * IntLattice<Int, BasInt, Dbl, RedDbl>::getNormalizer(
          NormaType norma, int alpha, bool dualF)
    {
      int dim = this->getDim();
      Normalizer<RedDbl> *normal;

      RedDbl logDensity;

      if (dualF) // dual basis 
        logDensity = - m_order * NTL::log(this->m_modulo);
      else // primal basis
        logDensity = m_order * NTL::log(this->m_modulo);

      switch (norma) {
        case BESTLAT:
          normal = new NormaBestLat<RedDbl> (logDensity, dim);
          break;
        case BESTBOUND:
          normal = new NormaBestBound<RedDbl> (logDensity, dim);
          break;
        case LAMINATED:
          normal = new NormaLaminated<RedDbl> (logDensity, dim);
          break;
        case ROGERS:
          normal = new NormaRogers<RedDbl> (logDensity, dim);
          break;
        case MINKL1:
          normal = new NormaMinkL1<RedDbl> (logDensity, dim);
          break;
        case MINK:
          normal = new NormaMinkowski<RedDbl> (logDensity, dim);
          break;
        case NONE:
          normal = new Normalizer<RedDbl> (logDensity, dim, "Norma_generic");
          break;
        default:
          std::cout << "normalizer:   no such case";
          exit (2);
      }
      return normal;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename Dbl, typename RedDbl>
      std::string IntLattice<Int, BasInt, Dbl, RedDbl>::toString() const
    {
      assert (0);
      return std::string();
    }
  extern template class IntLattice<std::int64_t, std::int64_t, double, double>;
  extern template class IntLattice<NTL::ZZ, NTL::ZZ, double, double>;
  extern template class IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>;

} // End namespace LatticeTester

#endif
