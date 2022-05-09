/**
 * This example presents how to use LatticeTester to make a figure of merit that
 * is a weighted sum, maximum, minimum or average of a mesure on projections on
 * a lattice. This examples shows how to use the weight classes, how to build
 * projections of a basis, but also how to normalize a computation. Typically,
 * when building figures of merit, mesures need to be rescalled to the same
 * interval to be compared with one another, this is what is called
 * normalization.
 * 
 * This example computes a simple spectral test on all projections of a lattice
 * in dimension 10, normalizes it between 0 and 1 and then takes the minimal
 * value observed as a figure of merit for that lattice.
 *
 * This example outputs only the figure of merit for the lattice for two
 * different normalizers. This is not really interesting in itself, hence it is
 * not included here. To get interesting informations on a figure of merit like
 * this one, it would be possible to store an few of the worst projections and
 * print them after the test. 
 * */

#define NTL_TYPES_CODE 2

#include <iostream>

#include "latticetester/Types.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Reducer.h"
#include "latticetester/ParamReader.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/UniformWeights.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;

int main() {
  //! Reading a matrix to use the normalizer class
  int min_dim = 0, max_dim = 10;
  ParamReader<MScal, BScal, RScal> reader("./44matrixEx.dat");
  reader.getLines();
  BMat matrix(max_dim,max_dim);
  unsigned int ln = 0;
  reader.readBMat(matrix, ln, 0, max_dim);
  IntLatticeBase<MScal, BScal, NScal, RScal> lat_basis(matrix, max_dim);
  double merit1 = 1.0, merit2 = 1.0;

  // The variables specific to the construction of a figure of merit
  UniformWeights weights(1.0); // This just puts a weight of 1 to everything
  BasisConstruction<BScal> constructor; // Computes projections basis
  IntLatticeBase<MScal, BScal, NScal, RScal> proj_basis(max_dim); // To store projections
  // CoordinateSets namespace contains classes to create iterators on sets of coordinates
  CoordinateSets::FromRanges coord(min_dim+1, max_dim, min_dim, max_dim-1);

  // For loop on the iterator built previously
  for(auto it = coord.begin(); it != coord.end(); it++){
    // Computing the projection
    constructor.ProjectionConstruction(lat_basis, proj_basis, *it);

    //! Computing the shortest vector in the lattice spanned by matrix
    proj_basis.updateVecNorm();
    proj_basis.sort(0);
    Reducer<MScal, BScal, NScal, RScal> red(proj_basis);
    red.redBKZ();
    red.shortestVector(L2NORM);
    double shortest = NTL::conv<double>(red.getMinLength());

    // Instanciating the normalizers
    // The prefered way of doing this is descibed in Normalizer documentation
    RScal log_density=-log(abs(NTL::determinant(proj_basis.getBasis())));
    Normalizer<RScal>* norma = new NormaBestLat<RScal>(log_density, max_dim);

    // Computing the figure of merit for this projection
    double merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    // Testing if it is the minimum as of now
    if (merit < merit1) merit1 = merit;
    delete norma;
    norma = new NormaBestBound<RScal>(log_density, max_dim);
    merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    if (merit < merit2) merit2 = merit;
    delete norma;
  }

  //! Printing the results in three simple lines
  std::cout << "Figure of merit with BestLat: " << merit1 << "\n";
  std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
    " weights and projections choices\n";

  return 0;
}
