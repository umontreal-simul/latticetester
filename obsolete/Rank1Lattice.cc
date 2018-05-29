#include "latticetester/Rank1Lattice.h"
#include "latticetester/Util.h"
#include <cassert>

using namespace std;
using namespace LatticeTester;


namespace LatticeTester
{

  Rank1Lattice::Rank1Lattice (const MScal & n, const MVect & a, int maxDim,
      NormType norm):
    IntLattice<MScal> (n, 1, maxDim, norm)
  {
    this->m_a = a;
    init();
  }


  //=========================================================================

  Rank1Lattice::~Rank1Lattice()
  {
    this->m_a.clear ();
  }


  //=========================================================================

  void Rank1Lattice::init()
  {
    IntLattice::init();
    for (int r = 1; r < getDim(); r++)
      this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
  }


  //=========================================================================

  Rank1Lattice & Rank1Lattice::operator= (const Rank1Lattice & lat)
  {
    if (this == &lat)
      return * this;
    copy (lat);
    init ();
    this->m_a = lat.m_a;
    return *this;
  }


  //=========================================================================

  Rank1Lattice::Rank1Lattice (const Rank1Lattice & lat):
    IntLattice<MScal> (lat.m_modulo, lat.getOrder (),
        lat.getDim (), lat.getNorm ())
  {
    // MyExit (1, "Rank1Lattice:: constructeur n'est pas terminé " );
    init ();
    this->m_a = lat.m_a;
  }


  //===========================================================================

  std::string Rank1Lattice::toStringCoef ()const
  {
    return toString (this->m_a, 0, getDim ());
  }


  //=========================================================================

  void Rank1Lattice::incDim ()
  {
    // kill();
    buildBasis (1 + getDim ());
    setNegativeNorm ();
    setDualNegativeNorm ();
  }


  //===========================================================================

  void Rank1Lattice::buildBasis (int d)
  {
    // assert(d <= getMaxDim());
    setDim (d);

    // conv(m_v[1][1], 1);

    for (int j = 0; j < d; j++) {
      this->m_basis (0, j) = this->m_a[j];
    }

    for (int i = 1; i < d; i++) {
      for (int j = 0; j < d; j++) {
        if (i == j) {
          this->m_basis (i, j) = this->m_modulo;
        } else {
          this->m_basis (i, j) = 0;
        }
      }
    }

    // if a[0] != 1, the basis must be triangularized
    if (this->m_basis (0, 0) != 1) {
      Triangularization < BMat > (this->m_basis, this->m_dualbasis, d, d, this->m_modulo);
      dualize ();
    }
    CalcDual < BMat > (this->m_basis, this->m_dualbasis, d, this->m_modulo);
    setNegativeNorm ();
    setDualNegativeNorm ();
  }

  //===========================================================================

  void Rank1Lattice::dualize ()
  {
    BMat tmps(this->m_basis);   this->m_basis = this->m_dualbasis;   this->m_dualbasis = tmps;
  }


} //namespace
