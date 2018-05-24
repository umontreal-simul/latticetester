#include "latticetester/Rank1Lattice.h"
#include "latticetester/Util.h"
#include <cassert>

#ifdef WITH_NTL
#else
using namespace boost::numeric::ublas;
#endif

using namespace std;
using namespace LatticeTester;


namespace LatticeTester
{

Rank1Lattice::Rank1Lattice (const MScal & n, const MVect & a, int maxDim,
                           NormType norm):
                  IntLattice::IntLattice (n, maxDim, norm)
{
   m_a.resize(maxDim);
   m_a = a;
   init();
}


//=========================================================================

Rank1Lattice::~Rank1Lattice()
{
   m_a.clear ();
}


//=========================================================================

void Rank1Lattice::init()
{
   IntLattice::init();
   for (int r = 1; r < getDim(); r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}


//=========================================================================

Rank1Lattice & Rank1Lattice::operator= (const Rank1Lattice & lat)
{
   if (this == &lat)
      return * this;
   copy (lat);
   init ();
   int dim = lat.getDim ();
   m_a.resize (dim);
   m_a = lat.m_a;
   return *this;
}


//=========================================================================

Rank1Lattice::Rank1Lattice (const Rank1Lattice & lat):
      IntLattice::IntLattice (lat.m_modulo, lat.getOrder (),
                              lat.getDim (), lat.getNorm ())
{
   // MyExit (1, "Rank1Lattice:: constructeur n'est pas terminé " );
   init ();
   int maxDim = lat.getDim ();
   m_a.resize (maxDim);
   m_a = lat.m_a;
}


//===========================================================================

std::string Rank1Lattice::toStringCoef ()const
{
   return toString (m_a, 1, getDim ());
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
      m_basis (0, j) = m_a[j];
   }

   for (int i = 1; i < d; i++) {
      for (int j = 0; j < d; j++) {
         if (i == j) {
            m_basis (i, j) = m_modulo;
         } else {
            m_basis (i, j) = 0;
         }
      }
   }

   // if a[0] != 1, the basis must be triangularized
   if (m_basis (0, 0) != 1) {
      Triangularization < BMat > (m_basis, m_dualbasis, d, d, m_modulo);
      dualize ();
   }
   CalcDual < BMat > (m_basis, m_dualbasis, d, m_modulo);
   setNegativeNorm ();
   setDualNegativeNorm (true);
}

//===========================================================================

void Rank1Lattice::dualize ()
{
   BMat tmps(m_basis);   m_basis = m_dualbasis;   m_dualbasis = tmps;
}


} //namespace
