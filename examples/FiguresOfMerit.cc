/**
 * This example should showcase the usage of the Normalizer and Weights classes.
 * */
/** 
 * Suppose that we have a `m`-lattice and its `m`-dual. We can get the figure of
 * merit called the normalized spectral test by dividing the length of the
 * shortest vector in the dual lattice with an upper bound. LatticeTester
 * provides several modules which are different approximations of this
 * upper bound. This example introduces them and presents the syntax to use
 * them.
 * */

#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/ParamReader.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"

#include "latticetester/CoordinateSets.h"
#include "latticetester/UniformWeights.h"
#include "latticetester/BasisConstruction.h"

#include "latticetester/Types.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"

using namespace LatticeTester;

int main() {
  //! Reading a matrix to use the normalizer class
  int min_dim = 0, max_dim = 10;
  ParamReader<MScal, BScal, RScal> reader("./44matrixEx.dat");
  reader.getLines();
  BMat matrix(max_dim,max_dim);
  unsigned int ln = 0;
  reader.readBMat(matrix, ln, 0, max_dim);
  IntLatticeBasis<MScal, BScal, NScal, RScal> lat_basis(matrix, max_dim);
  double merit1 = 1.0, merit2 = 1.0;

  UniformWeights weights(1.0);
  BasisConstruction<BScal> constructor;
  IntLatticeBasis<MScal, BScal, NScal, RScal> proj_basis(max_dim);
  CoordinateSets::FromRanges coord(min_dim+1, max_dim, min_dim, max_dim-1);
  for(auto it = coord.begin(); it != coord.end(); it++){
    constructor.ProjectionConstruction(lat_basis, proj_basis, *it);

    //! Computing the shortest vector in the lattice spanned by matrix
    proj_basis.updateVecNorm();
    proj_basis.sort(0);
    Reducer<MScal, BScal, NScal, RScal> red(proj_basis);
    red.redBKZ();
    red.shortestVector(L2NORM);
    double shortest = NTL::conv<double>(red.getMinLength());
    // We now instanciate the different Normalization modules and print the
    // "different" results (the results won't actually be different because the
    // bound is exact in dimension 4 and all the modules have the same).
    RScal log_density=-log(abs(NTL::determinant(proj_basis.getBasis())));
    // As described in Normalizer.h, we need to create a pointer to a Normalizer
    // object and dynamically allocate memory/objects to it.
    Normalizer<RScal>* norma = new NormaBestLat<RScal>(log_density, max_dim);
    double merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    if (merit < merit1) {
      merit1 = merit;
      std::cout << "New worst projection is: " << *it << std::endl;
    }
    delete norma;
    norma = new NormaBestBound<RScal>(log_density, max_dim);
    merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    if (merit < merit2) {
      merit2 = merit;
      std::cout << "New worst projection is: " << *it << std::endl;
    }
    delete norma;
  }

  std::cout << "Figure of merit with BestLat: " << merit1 << "\n";
  std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;

  //! Printing the results of the figure of merit calculation.
  return 0;
}
