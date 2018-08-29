#include "latticetester/ntlwrap.h"
#include "latticetester/Util.h"

namespace LatticeTester {
  /**
   * **WARNING** THIS MODULE IS UNDER CONSTRUCTION.
   *
   * This class implements general methods to perform a lattice basis
   * construction from a set of vectors, as well as general methods to obtain
   * the dual lattice basis depending on the current lattice basis.
   *
   * This module still has to be implemented.
   *
   * This module should work only with NTL matrices because those are objects
   * aware of their shape. I think that this is important since we use object
   * oriented stuff.
   * */
  template<typename BasInt> class BasisConstruction{

    private:
      typedef NTL::vector<BasInt> BasIntVec;
      typedef NTL::matrix<BasInt> BasIntMat;

    public:
      /**
       * */
      BasisConstruction(){};

      /**
       * Destructor.
       * */
      ~BasisConstruction(){};

      /**
       * This functions takes a set of generating vectors of a vector space and
       * finds a basis for this space whilst applying LLL reduction.
       * */
      void LLLConstruction();

      /**
       * This function does some kind of Gaussian elimination on
       * \f$\text{span}(\text{matrix}^t)\f$ as a vector space on \f$\mathbb{Z}\f$ to
       * obtain a basis of this vector space (lattice).
       * This takes a matrix `matrix` where each row `v[0], ..., v[n]` is
       * considered as a vector. This basically performs, for each column,
       * Euclid's algorithm on the rows under the diagonal to set all
       * coefficients to zero except the one in the diagonal. Finding the GCD
       * of all the coefficients from the diagonal and under allows to perform
       * Gaussian elimination using only the 3 allowed operations that do not
       * change the span of the vectors:
       *   - Multiply a line by \f$-1\f$,
       *   - Add an integer multiple of line \f$i\f$ to line \f$j\f$ for \f$i \neq j\f$,
       *   - Swap line \f$i\f$ and line \f$j\f$.
       * After constructing this basis, the algorithm eliminates negative
       * coefficients in the matrix.
       *
       * \todo Give a running time. **Marc-Antoine**: This runs in
       * O(min(row, col)*lg(n)*row + row^2) time I think, but it might be a bit
       * more because we do operations on vectors so we might need to multiply
       * this by col.
       * */
      void GCDConstruction(BasIntMat& matrix);

      /**
       * Suppose the basis matrix contains basis vectors on its lines and is
       * \f$p\times q\f$ where \f$q \geq p\f$. We can compute the \f$m\f$-dual
       * as follows.
       * Let's note the basis matrix `V` and the dual matrix `W` and have lines
       * of `W` also contain dual basis vectors in its lines. We know that
       * \f$VW^t = mI_{p\times p}\f$ where \f$I\f$ is the identity matrix. Now, in
       * the case of a \f$m\f$-dual computation, we can assume that all the
       * arithmetic is done modulo \f$m\f$ and that \f$m\f$ is prime??
       * */
      void DualConstruction(BasIntMat& matrix, BasInt modulo = BasInt(0));

    private:

      /**
       * This is the matrix this object is working on. This matrix is modified
       * as the program advances.
       * */
      BasIntMat m_mat;

      /**
       * This stores the matrix the object was initialized with. This is usefull
       * for implementing methods that verify this class does not change the
       * lattice. This is also usefull to compare different reduction methods.
       * Finally, this can be kept for the sake of comparison.
       * */
      BasIntMat m_old;

  };

  template<typename BasInt>
    void BasisConstruction<BasInt>::GCDConstruction(BasIntMat& matrix)
  {
    // It is important to note that the lines of matrix are the basis vectors
    std::cout << "First matrix:\n" << matrix;
    long rows = matrix.NumRows();
    long cols = matrix.NumCols();
    long max_rank = rows < cols ? rows : cols;
    long rank = 0;
    // The basis will have at most max_rank vectors.
    BasInt q;
    for (long i = 0; i < max_rank; i++) {
      // We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
      // algorithm and applying transformations to the matrix
      for (long j = i+1; j < rows; j++) {
        while (matrix[j][i] != 0) {
          matrix[i].swap(matrix[j]);
          q = matrix[j][i]/matrix[i][i];
          matrix[j] -= q * matrix[i];
          //std::cout << "Matrix:\n" << matrix;
        }
      }
      // We make sure that the coefficients that can be not zero are positive.
      // This is because the algorithms we use work for positive vectors.
      if (matrix[i][i] < 0) matrix[i] *= -1;
      for (long j = i-1; j >= 0; j--) {
        if (matrix[j][i] < 0) {
          if (-matrix[j][i] > matrix[i][i]){
            q = -matrix[j][i]/matrix[i][i] + 1;
            matrix[j] += q * matrix[i];
          } else {
            matrix[j] += matrix[i];
          }
        }
      }
      if (matrix[i][i] != 0) rank++;
    }
    //std::cout << "Matrix before truncating\n" << matrix;
    // We remove zero vectors from the basis.
    matrix.SetDims(rank, cols);
    std::cout << "Final matrix\n" << matrix;
  }

  //============================================================================

  /**
   * This algorithm calculates the dual as well as the `m` used for rescalling.
   * It checks if this `m` divides the modulo given to the algorithm.
   * Right now, this assumes the basis is triangular, might need to change it.
   * */
  template<typename BasInt>
    void BasisConstruction<BasInt>::DualConstruction(BasIntMat& matrix,
         BasInt modulo)
  {
    // We first make the matrix triangular. I will change this step.
    if (! CheckTriangular(matrix, matrix.NumRows(), 0)) GCDConstruction(matrix);
    long dim = matrix.NumRows();
    if (dim != matrix.NumCols()) {
      std::cout << "matrix has to be square, but dimensions do not fit.\n";
      return;
    }
    BasInt m;
    m = 1;
    BasIntMat result;
    NTL::ident(result, dim);
    BasInt gcd;
    for (long i = dim-1; i>=0; i--) {
      for (long j = i+1; j < dim; j++) {
        result[i] -= matrix[i][j] * result[j];
        matrix[i] -= matrix[i][j] * matrix[j];
      }
      gcd = matrix[i][i];
      for (long j = i; j < dim; j++) {
        gcd = NTL::GCD(gcd, result[i][j]);
      }
      gcd = matrix[i][i] / gcd;
      if (gcd != 1) {
        result *= gcd;
        m *= gcd;
      }
      for (long j = i; j < dim; j++) {
        result[i][j] /= matrix[i][i];
      }
      matrix[i][i] = 1;
    }
    std::cout << matrix << std::endl;
    std::cout << "Dual to triangular basis:\n " << result << std::endl;
    std::cout << "m:" << m << std::endl;
  }

} // end namespace LatticeTester
