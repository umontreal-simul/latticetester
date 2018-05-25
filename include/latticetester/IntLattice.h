#ifndef LATTICETESTER__INTLATTICE_H
#define LATTICETESTER__INTLATTICE_H

#include <cassert>

#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/Lacunary.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
// #include "latticetester/MRGComponent.h"

namespace LatticeTester {

  /**
   * \copydoc LatticeTester::IntLatticeBasis
   * Beside, it contains fonction use in LatticeTester
   */
  class IntLattice : public IntLatticeBasis {
    public:

      /**
       * \copydoc LatticeTester::IntLatticeBasis::IntLatticeBasis (
       const BMat primalbasis,
       const BMat dualbasis,
       const MScal modulo,
       const int dim,
       NormType norm = L2NORM);
       */
      IntLattice (MScal modulo, int k, int maxDim, LatticeTester::NormType norm = LatticeTester::L2NORM);

      /**
       * \copydoc LatticeTester::IntLatticeBasis::IntLatticeBasis (const IntLatticeBasis & Lat)
       */
      IntLattice (const IntLattice & Lat);

      /**
       * Same as assignment operator above.
       */
      void copy (const IntLattice & lattice);

      /**
       * Destructor.
       */
      virtual ~IntLattice ();

      /**
       * Init the matrix
       */
      void init ();

      /**
       * Return the order of the matrix
       */
      int getOrder() const { return m_order; }

      /**
       * Increment the dimension of the element of the lattice
       */
      virtual void incDim ();

      /**
       * Computes the logarithm of the normalization factor
       * (<tt>m_lgVolDual2</tt>) in all dimensions \f$\le\f$ `MaxDim` for
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
      void buildProjection (IntLattice* lattice, const Coordinates & proj);

      /**
       * Builds the basis for the lattice in dimension `d`.
       * This fonction is implemented in subclasses
       */
      virtual void buildBasis (int d);

      /**
       * Creates and returns the normalizer corresponding to criterion
       * `norma`. In the case of the \f$P_{\alpha}\f$ test, the argument
       * `alpha` = \f$\alpha\f$. In all other cases, it is unused.
       */
      LatticeTester::Normalizer * getNormalizer (LatticeTester::NormaType norma, 
          int alpha, bool dualF);

      /**
       * Does nothing is this base class.
       */
      virtual void setLac (const Lacunary &) {
        LatticeTester::MyExit(1, "IntLattice.setLac is empty"); }

      /**
       * Returns the vector of multipliers (or coefficients) \f$A\f$ as a
       * string.
       */
      virtual std::string toStringCoef() const;


      // #ifdef WITH_NTL
      //    /**
      //     * The components of the lattice when it is built out of more than one
      //     * component. When there is only one component, it is unused as the
      //     * parameters are the same as above.
      //     */
      // //    std::vector<MRGComponent *> comp;
      // #endif


    protected:

      /**
       * \copydoc LatticeTester::IntLatticeBasis::kill()
       */
      virtual void kill ();

      /**
       * The order of the basis.
       */
      int m_order;

      /**
       * The maximum Dimension for the test
       */
      //int m_maxDim;

      /*
       * Represente sur dual along the diagonal?? ERWAN
       */
      double *m_lgVolDual2;

      /**
       * The logarithm \f$\log_2 (m^2)\f$.
       */
      double m_lgm2;

      /**
       * Use to save the dual basis and the basis in some works.
       */
      BMat m_wSI, m_vSI;

      /**
       * Working Variables use in MRGLattice.h
       */
      MScal m_t1, m_t2, m_t3;

  };

}
#endif
