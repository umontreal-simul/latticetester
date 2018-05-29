


#ifndef LATTICETESTERROUTINES_H
#define LATTICETESTERROUTINES_H

#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/LatticeAnalysis.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaRogers.h"

#include "latticetester/Reducer.h"

namespace LatticeTester {

  /**
   * This function allows computation of the shortest non-zero vector in a lattice,
   * according to the selected norm. Many parameters can bet set by the user, otherwise
   * the function work with default values.
   * Returns -1.0 if there was an error in Branch-and-Bound procedure. Return the length
   * of the shortest non-zero vector otherwise.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double ShortestVector(BasIntMat matrix, NormType norm, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      int dimension;
      if (matrix.size1() != matrix.size2()) {
        MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
        exit(1);
        // C'est pas un peu nul ça ?
      } else 
        dimension = (int) matrix.size1();

      // creating objects needed to perform the test
      IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl> basis (matrix, dimension, norm);
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat> red (basis);

      // performing pre-reduction
      switch (preRed) {
        case BKZ:
          red.redBKZ(fact, blocksize, doublePrecision);
          break;
        case LenstraLL:
          red.redLLLNTL(fact, doublePrecision);
          break;
        case PreRedDieter:
          red.preRedDieter(0);
          break;
        case NOPRERED:
          cout << "WARNING: no pre-reduction is performed before applying Branch-and-Bound";
          cout << " procedure. Running time could increase dramaticaly with the matrix dimension.";
          cout << endl;
          break;
        default:
          MyExit(1, "LatticeTesterRoutines::ShortestVector:   NO SUCH CASE FOR PreReductionType");
          exit(1);
      }

      // performing the Branch-and-Bound procedure to find the shortest non-zero vector.
      // red.shortestVector(norm) is a bool set to *true* if the Branch-and-Bound algorithm
      // terminates without any error.
      if (red.shortestVector(norm))
        return conv<double>(red.getMinLength());
      else 
        return -1.0;
    }

  /**
   * Same thing as before but with the possibility to set a different value for
   * the variable maxNodesBB.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double ShortestVector(BasIntMat matrix, NormType norm, long maxNodesBB,
        PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE,
        double fact = 0.999, int blocksize = 20)
    {
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat>::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
      return ShortestVector<Int, BasInt>(matrix, norm, preRed, doublePrecision, fact, blocksize);
    }

  /**
   * This function compute the Figure of Merit to a given matrix, according to a
   * normalization criteria. It first computes the shortest non-zero vector using the
   * above functions. It then normalizes this value.
   * Returns -1.0 if there was an error in Branch-and-Bound procedure while calculating
   * the length of shortest non-zero vector. Return the figure of merit otherwise.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double FigureOfMerit(BasIntMat matrix, NormaType normalizerType, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      double merit;

      int dimension;
      if (matrix.size1() != matrix.size2()) {
        MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
        exit(1);
        // C'est pas un peu nul ça ?
      } else 
        dimension = (int) matrix.size1();

      // creation of the norm and normalizer objects
      NormType norm;
      Normalizer<RedDbl>* normalizer;

      // calculation of the log-density of the matrix used to initialize the normalizer
      RedDbl logDensity;
#if NTL_TYPES_CODE > 1
      logDensity = - log( abs(NTL::determinant(matrix)) );
#else 
      // NTL library does not support matrix with double: we then use the boost library
      boost::numeric::ublas::matrix<long>  mat_tmps;
      mat_tmps.resize(dimension, dimension);
      for (unsigned int i = 0; i < dimension; i++) {
        for (unsigned int j = 0; j < dimension; j++) {
          mat_tmps(i,j) = matrix(i,j);
        }
      }
      logDensity = -log( abs(det_double(mat_tmps)) );
#endif

      // we initialize the norm and normalizer objects according to the normalizerType used
      switch (normalizerType) {
        case BESTLAT:
          norm = L2NORM;
          normalizer = new NormaBestLat<RedDbl> (logDensity, dimension);
          break;
        case LAMINATED:
          norm = L2NORM;
          normalizer = new NormaLaminated<RedDbl> (logDensity, dimension);
          break;
        case ROGERS:
          norm = L2NORM;
          normalizer = new NormaRogers<RedDbl> (logDensity, dimension);
          break;
        case MINKOWSKI:
          norm = L2NORM;
          normalizer = new NormaMinkowski<RedDbl> (logDensity, dimension);
          break;
        case MINKL1:
          norm = L1NORM;
          normalizer = new NormaMinkL1<RedDbl> (logDensity, dimension);
          break;
        default: //PALPHA_N, NORMA_GENERIC, L1, L2
          MyExit(1, "LatticeTesterRoutines::FigureOfMerit:   NO SUCH CASE FOR *normalization type*");
          exit(1);
      }

      // compute the shortest non-zero vector
      merit = ShortestVector<Int, BasInt>(matrix, norm, preRed, doublePrecision, fact, blocksize);

      if (merit == -1.0)
        return -1.0; // the BB procedure didn't terminated well

      // normalization step
      merit /= conv<double>(normalizer->getPreComputedBound(dimension));

      delete normalizer;
      return merit; 
    }

  /**
   * Same thing as before but with the possibility to set a different value for
   * the variable maxNodesBB.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double FigureOfMerit(BasIntMat matrix, NormaType normalizerType, long maxNodesBB,
        PreReductionType preRed = BKZ, PrecisionType doublePrecision = DOUBLE,
        double fact = 0.999, int blocksize = 20)
    {
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat>::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
      return FigureOfMerit<Int, BasInt>(matrix, normalizerType, preRed, doublePrecision, fact, blocksize);
    }

  /**
   * This function reduces a basis to a Minkowski-reduced basis. Such basis has strong
   * properties regarding the length of its vectors but will require a huge running time,
   * especially when the dimension of the basis increases. Such Minkowski-reduced basis
   * is usefull, for example, to calculate a Beyer quotient (as implemented in
   * FigureOfMerit()).
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    bool MinkowskiReduction(BasIntMat & matrix, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      int dimension;
      if (matrix.size1() != matrix.size2()) {
        MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
        exit(1);
        // C'est pas un peu nul ça ?
      } else 
        dimension = (int) matrix.size1();

      // creating objects needed to perform the test
      IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl> basis (matrix, dimension, L2NORM);
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat> red (basis);

      // performing pre-reduction
      red.preRedDieter(0);

      // performing the Minkowski reduction. Returns *false* if the algorithm didn't terminated well, 
      // returns *true* if it a success.
      bool reductionSuccess = red.reductMinkowski (0);
      matrix = red.getIntLatticeBasis().getBasis();
      return reductionSuccess;
    }

  /**
   * Same thing as before but with the possibility to set a different value for
   * the variable maxNodesBB.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    bool MinkowskiReduction(BasIntMat & matrix, long maxNodesBB, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat>::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
      return MinkowskiReduction<Int, BasInt>(matrix, preRed, doublePrecision, fact, blocksize);
    }

  /**
   * This function compute the Figure of Merit to a given matrix, according to the
   * Beyer criteria. It first computes the Minkowski-reduced basis of the lattice 
   * and then makes the quotient of shortest over longest vector.
   * Returns -1.0 if there was an error in Branch-and-Bound procedure while calculating
   * the Minkowski-reduced basis. Return the figure of merit otherwise.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double FigureOfMeritBeyer(BasIntMat matrix, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      int dimension;
      if (matrix.size1() != matrix.size2()) {
        MyExit(1, "LatticeTesterRoutines::ShortestVector:   NEED A SQUARE MATRIX");
        exit(1);
        // C'est pas un peu nul ça ?
      } else 
        dimension = (int) matrix.size1();

      bool reductionSuccess;
      //m_lat->dualize ();
      reductionSuccess = MinkowskiReduction<Int, BasInt>(matrix, preRed, doublePrecision, fact, blocksize);
      //m_lat->dualize ();

      //if (m_dualF)
      //m_lat->dualize ();

      if (reductionSuccess) {
        IntLatticeBasis<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl> basis (matrix, dimension, L2NORM);
        basis.updateScalL2Norm (0);
        basis.updateScalL2Norm (dimension-1);
        double x1, x2; // maybe using something else than double (xdouble, RR?)
        conv (x1, basis.getVecNorm (0));
        conv (x2, basis.getVecNorm (dimension-1));
        return sqrt(x1 / x2);
      } else 
        return -1.0;
    }

  /**
   * Same thing as before but with the possibility to set a different value for
   * the variable maxNodesBB.
   */
  template<typename Int, typename BasInt, typename BasIntVec, typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl, typename RedDblVec, typename RedDblMat>
    double FigureOfMeritBeyer(BasIntMat matrix, long maxNodesBB, PreReductionType preRed = BKZ,
        PrecisionType doublePrecision = DOUBLE, double fact = 0.999, int blocksize = 20)
    {
      Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RedDblVec, RedDblMat>::maxNodesBB = maxNodesBB; // setting the number of nodes visited in the BB procedure
      return FigureOfMeritBeyer<Int, BasInt>(matrix, preRed, doublePrecision, fact, blocksize);
    }

} // end namespace LatticeTester

#endif
