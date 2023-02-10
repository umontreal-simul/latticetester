// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER_REDUCER_H
#define LATTICETESTER_REDUCER_H

#include "NTL/LLL.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/Const.h"
#include "latticetester/Util.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/Types.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <ctime>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <type_traits>





namespace LatticeTester {

/**
 * This `Reducer` class provides facilities to reduce the basis of a lattice
 * (an `IntLatticeBase` object) in different ways (pairwise, LLL, BKZ, Minkowski
 * \cite rDIE75a, \cite mLEN82a, \cite mSCH91a),
 * and find a shortest nonzero vector in the lattice using a BB algorithm \cite rFIN85a.
 * Most of the methods do not use or change the m-dual lattice.
 * To reduce the m-dual basis or find a shortest nonzero vector in it,
 * one should first dualize the lattice; the method `IntLattice::dualize` does that.
 * Some of the lattice reduction methods are NTL wraps.
 * For those, the `Int` type can only be `ZZ`, because NTL offers no other option.
 * For LLL, we have both our simple implementation and implementations from
 * NTL, which are usually more efficient; see redLLLNTL and redLLLNTLExact.
 * The method redBKZ is also a wrapper for the NTL algorithm for BKZ reduction.
 *
 * The `shortestVector` method does not apply any pre-reduction by itself.
 * Before calling it, one should always reduce the basis separately beforehand
 * with an LLL or BKZ reduction, because it drastically reduces the size of
 * the BB search.
 *
 * To use these facilities, one should create an instance of `Reducer` by passing a
 * `IntLatticeBase` object to the constructor. This reducer will always work on this
 * `IntLatticeBase` object, which is accessible internally via a pointer.
 * This pointer is also accessible via `getIntLatticeBase`.
 * The methods of the `Reducer` object will modify this `IntLatticeBase` object.
 */

template<typename Int, typename Real, typename RealRed>
class Reducer;

/// \cond specReducerDec
/*
 * This struct specializes some of the functions in a `Reducer`. This is a
 * workaround needed, implementation-wise, since we cannot specialize member
 * functions of a class without specializing the whole class.
 *
 * What this does is that this structure contains specific functions of the
 * original class Reducer that will act differently depending on the types of
 * the arguments. Instead of specializing all the methods of Reducer for all
 * our use cases (since some of them depend on the type) we can specialize
 * only this structure.
 *
 */
template<typename Int, typename Real, typename RealRed>
struct specReducer {
	void redLLLNTL(Reducer<Int, Real, RealRed> &red, double delta,
			PrecisionType precision, int dim);
	void redLLLNTLExact(Reducer<Int, Real, RealRed> &red, double delta);
	void redBKZ(Reducer<Int, Real, RealRed> &red, double delta,
			std::int64_t blocksize, PrecisionType precision, int dim);
};
/// \endcond

template<typename Int, typename Real, typename RealRed>
class Reducer {

private:
	// Local typedefs for matrix and vector types needed in the class.
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;
	typedef NTL::vector<RealRed> RealRedVec;
	typedef NTL::matrix<RealRed> RealRedMat;

public:

	/**
	 * Constructor that initializes the reducer to work on the lattice `lat`.
	 */
	Reducer(IntLatticeBase<Int, Real, RealRed> &lat);

	/**
	 * Copy constructor.
	 */
	Reducer(const Reducer<Int, Real, RealRed> &red);

	/**
	 * Destructor.
	 */
	~Reducer();

	/**
	 * Assignment operator that makes a deep copy of `red`
	 * into the current object, using `copy`.
	 */
	Reducer<Int, Real, RealRed>& operator=(
			const Reducer<Int, Real, RealRed> &red);

	/**
	 * Copies `red` into the current object.
	 */
	void copy(const Reducer<Int, Real, RealRed> &red);

	/**
	 * Computes a shortest non-zero vector in the `IntLatticeBase` lattice stored
	 * in this object, with respect to norm type `norm`, using the BB
	 * algorithm described in \cite rLEC97c.
	 * The admissible norm types here are `L1NORM` and `L2NORM`.
	 * The `NormType` attribute of this `IntLatticeBase`
	 * object will be changed to `norm`. If `MaxNodesBB` is exceeded
	 * during the branch-and-bound, the method aborts and returns
	 * `false`. Otherwise, it returns `true`. If the reduction was
	 * successful, the new reduced basis can be accessed via
	 * `getIntLatticeBase()`.
	 *
	 * It is strongly recommended to use `redBKZ` or `redLLLNTL` to pre-reduce
	 * the basis before invoking this method; this is not done automatically.
	 * The param 'decomp' acepted are 'cholesky' or 'triangular'. It specify 
	 * the decompistion used in the Branch-and-Bound algorithm.
	 */
	bool shortestVector(NormType norm,std::string decomp);

	/**
	 * This method performs pairwise reduction sequentially on all vectors
	 * of the basis whose indices are greater of equal to `dim >=0`.
	 * The boolean vector `xx[]` is used internally in Minkowski reduction.
	 */
	void redDieter(int dim, bool xx[] = NULL);

	/**
	 * Same as `redDieter(dim)` but the choice of
	 * vectors on which to perform pairwise reduction is randomized,
	 * using a simple RNG from the standard C library, with the given integer seed.
	 */
	void redDieterRandomized(int dim, int seed);

	/**
	 * Performs a LLL basis reduction with factor `delta` \cite iLEC22l.
	 * The reduction is applied to the first `dim` basis vectors and coordinates
	 * when `dim > 0`, and to the entire basis (all vectors) when `dim=0`.
	 * It always uses the Euclidean norm.
	 * The factor `delta` must be between 1/2 and 1. The closer it is to 1,
	 * the more the basis is reduced, in the sense that the LLL
	 * algorithm will enforce tighter conditions on the basis,
	 * but this will require more work. The reduction algorithm is
	 * applied until `maxcpt` successful transformations have been done,
	 * or until the basis is correctly LLL reduced.
	 *
	 * This is our simple implementation of the LLL reduction.
	 * For NTL types, it is considerably slower than what is
	 * available through NTL, but it performs well on standard C++ types.
	 */
	void redLLL(double delta = 0.999999, std::int64_t maxcpt = 1000000000,
			int dim = 0);

	/**
	 * This is the NTL implementation of the floating point version of the
	 * LLL reduction algorithm presented in \cite mSCH91a.
	 * The factor `delta` has the same meaning as in `redLLL`.
	 * The parameter `precision` specifies the precision of the floating point numbers
	 * that the algorithm will use. `Const.h` provides a list of the possible values,
	 * and their description is done in the module `LLL` of NTL.
	 * The reduction is applied to the first `dim` basis vectors when `dim > 0`,
	 * and to the entire basis (all vectors and coordinates) when `dim=0`.
	 */
	void redLLLNTL(double delta = 0.999999, PrecisionType precision = QUADRUPLE,
			int dim = 0);

	/**
	 * This implements an exact algorithm to perform the original LLL reduction.
	 * This is slower than `redLLLNTL`, but more accurate.
	 */
	void redLLLNTLExact(double delta);

	/**
	 * This calls the NTL implementation of the floating point version of the
	 * BKZ reduction algorithm presented in \cite mSCH91a,
	 * with reduction factor `delta` and block size `blocksize`; see \cite iLEC22l.
	 * The factor `delta` has a similar meaning as in `redLLL`.
	 * The `precision` and `dim` parameters have the same meaning as in `redLLLNTL`.
	 * The parameter `blocksize` gives the size of the blocks in the BKZ
	 * reduction. Roughly, larger blocks means a stronger condition.
	 * A `blocksize` of 2 is equivalent to LLL reduction.
	 */
	void redBKZ(double delta = 0.999999, int blocksize = 10,
			PrecisionType precision = QUADRUPLE, int dim = 0);

	/**
	 * Reduces the current basis to a Minkowski-reduced basis with respect
	 * to the Euclidean norm, assuming that the first \f$d\f$ vectors are
	 * already reduced and sorted. If `MaxNodesBB` is exceeded during one
	 * of the branch-and-bound step, the method aborts and returns `false`.
	 * Otherwise it returns `true`, the basis is reduced and sorted by
	 * increasing vector lengths.
	 */
	bool reductMinkowski(int d);

	/**
	 * Returns the length of the current shortest basis vector in the lattice.
	 */
	Real getMinLength() {
		if (m_lat->getNormType() == L2NORM)
			return sqrt(m_lMin2);
		else
			return m_lMin;
	}

	/**
	 * Returns the length of the current longest basis vector in the lattice.
	 */
	Real getMaxLength() {
		if (m_lat->getNormType() == L2NORM)
			return sqrt(m_lat->getVecNorm(m_lat->getDim() - 1));
		else
			return m_lat->getVecNorm(m_lat->getDim() - 1);
	}

	/**
	 * Sets the vector of bounds on the square of the acceptable shortest
	 * vector lengths (for the Euclidean norm),
	 * in dimensions from `dim1+1` to `dim2`. `thresholds[i]` must
	 * contain a lower bound on the square of the length of the shortest
	 * vector in the lattice in dimension `i+1`. Such a bound, if it is
	 * set, will be used during during the Branch-and-Bound step when
	 * searching the shortest vector of a lattice. If the Branch-and-Bound
	 * finds that the shortest vector of the lattice is smaller than the
	 * bound, the algorithm will stop. This is useful in the case where
	 * the user is searching for a good lattice with the spectral test
	 * since this can cut off some work.
	 *
	 * ***  This is used only in `RedBB0`.  Perhaps the value should be passed
	 * as a parameter to `shortestVector`?  Or even the value should be checked
	 * just before and after calling `shortestVector`, for simplicity.  **
	 */
	void setBoundL2(const RealVec &thresholds, int dim1, int dim2);

	/**
	 * Returns the lattice that this object is working on.
	 * */
	IntLatticeBase<Int, Real, RealRed>* getIntLatticeBase() {
		return m_lat;
	}

	/**
	 * The maximum number of nodes in the branch-and-bound tree when
	 * calling `shortestVector` or `reductMinkowski`. When this number is
	 * exceeded, the method aborts and returns `false`.
	 *
	 * *** Perhaps this should be a parameter of these functions?
	 */
	static std::int64_t maxNodesBB;

private:

	/**
	 * The lattice that this object is working on.
	 */
	IntLatticeBase<Int, Real, RealRed> *m_lat;
//	IntLatticeBase<Int, Real, RealRed> *m_lat2;

	/**
	 * Contains specialized implementations of member methods depending on
	 * the types given to the `Reducer` template. This is private
	 * because it is possible to call specialized method at an higher level.
	 * We do this because NTL offers reductions only for ZZ type.
	 *  ** P.L.: I really do not like this.  Perhaps NTL should be able to promote
	 *     the int64_t to ZZ to apply the methods.   **
	 */
	struct specReducer<Int, Real, RealRed> spec;

	/*
	 * Method used in `shortestVector` to perform a transformation of
	 * stage 3 described in \cite rAFF85a.
	 * We call this when we have found a shorter vector which is a linear combination of the
	 * previous basis vectors with coefficients given in z. This procedure updates the basis
	 * so that the new shortest vector is in first and the other vectors are adjusted accordingly.
	 * With the L2 norm, an alternative might be to add the new shortest vector to
	 * the basis to get dim+1 generating vectors, and apply LLL to recover a new basis.
	 *  ***  We should test this and compare the speeds.
	 * Here, we assume that \f$\sum_{i=1}^t z_i V_i\f$ is a
	 * short vector that will enter the basis. Then we try to reduce some vectors
	 * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
	 * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
	 * \f$V_j\f$ when this happens. We return in \f$k\f$ the last index
	 * \f$j\f$ such that \f$|z_j|=1\f$.
	 * This method does not change the m-dual.
	 */
	void transformStage3ShortVec(std::vector<std::int64_t> &z, int &k);

	/**
	 * Method used in `reductMinkowski` to perform a transformation of
	 * stage 3 described in \cite rAFF85a. Assumes that
	 * \f$\sum_{i=1}^t z_i V_i\f$ is a
	 * short vector that will enter the basis. Tries to reduce some vectors
	 * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
	 * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
	 * \f$V_j\f$ when this happens. Returns in \f$k\f$ the last index
	 * \f$j\f$ such that \f$|z_j|=1\f$.    **This method updates the m-dual.**
	 */
	void transformStage3Mink(std::vector<std::int64_t> &z, int &k);

	/**
	 * Performs pairwise reductions. This method tries to reduce each basis
	 * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
	 * adding to it a multiple of the \f$i\f$-th vector. Always uses the
	 * Euclidean norm.
	 */
	void pairwiseRedPrimal(int i, int d, bool xx[] = NULL);

	/**
	 * Performs pairwise reductions, trying to reduce every other vector of
	 * the *dual* basis by adding multiples of the \f$i\f$-th vector. That
	 * may change the \f$i\f$-th vector in the primal basis. Each such dual
	 * reduction is actually performed only if that does not increase the
	 * length of vector \f$i\f$ in the primal basis. Always uses the
	 * Euclidean norm.
	 */
	void pairwiseRedDual(int i, bool xx[] = NULL);

	/**
	 * Tries to shorten the vectors of the primal basis using
	 * branch-and-bound, in `reductMinkowski`.
	 */
	bool redBBMink(int i, int d, int Stage, bool &smaller, bool xx[] = NULL);

	/**
	 * Tries to shorten the smallest vector of the primal basis using
	 * branch-and-bound.  Used in `shortestVector`.
	 * The param 'decomp' accepted are only 'cholesky' and 'triangular' decompostion. 
	 * The branch-and-bound is more faster with 'cholesky' than  'triangular'.
	 */
	bool redBBShortVec(NormType norm,std::string decomp);

	/**
	 * Tries to find shorter vectors; recursive procedure used in `reductMinkowski`.
	 */
	bool tryZMink(int j, int i, int Stage, bool &smaller, const IntMat &WTemp);

	/**
	 * Tries to find a shorter vector; recursive procedure used in `shortestVector`.
	 * The param 'decomp' accepted are only 'cholesky' and 'triangular'. The better
	 * The 'cholesky' decomposition is more faster than the triangular decomposition .  
	 */
	bool tryZShortVec(int j, bool &smaller, NormType norm, std::string decomp);

	/**
	 * Computes a Cholesky decomposition of the basis. Returns in `C0` the
	 * elements of the upper triangular matrix of the Cholesky
	 * decomposition that are above the diagonal. Returns in `DC2` the
	 * squared elements of the diagonal.
	 */
	bool calculCholesky(RealRedVec &DC2, RealRedMat &C0);

	/**
	 * Permutes the \f$i^{th}\f$ and the \f$j^{th}\f$ line, and the
	 * \f$i^{th}\f$ and the \f$j^{th}\f$ column of the scalar products matrix.
	 */
	void permuteGramVD(int i, int j, int n);

	/**
	 * Recalculates the first \f$n\f$ entries of the \f$j^{th}\f$ column of
	 * the Cholesky matrix of order 2.
	 */
	void calculCholesky2LLL(int n, int j);

	/**
	 * Recalculates the entry (\f$i\f$, \f$j\f$) of the Cholesky matrix of order 2.
	 */
	void calculCholesky2Ele(int i, int j);

	/*
	 * Recomputes the element in row `j` and column `j` of the matrix of scalar products.
	 */
	void miseAJourGramVD(int j);

	/*
	 * Computes and stores in `m_gramVD` the matrix of scalar products `m_lat->V[i]*m_lat->V[j]`.
	 * Equivalent to computing `m_lat->V * transpose(m_lat->V)`.  Used in redLLL.
	 */
	void calculGramVD();

	/**
	 * Reduces the Cholesky matrix by adding a multiple of the i-th vector
	 * to the j-th vector. It updates the Gram-Schmidt matrix.
	 */
	void reductionFaible(int i, int j);

	/**
	 * Debug function that sorts and prints the primal and dual bases
	 * to standard output, using the `write` function.
	 */
	void tracePrintBases(char *mess);

	/**
	 * A vector that contains a lower bound on the acceptable (squared) length
	 * of the shortest vector, for each number of dimensions.
	 * If any vector of the lattice is shorter than this bound,
	 * we stop the reduction immediately and reject this lattice since
	 * we already know that its shortest vector is too small.
	 * This is used in `RedBBShortVec`, and useful for the seek procedure in LatMRG.
	 * We could also just check this after the BKZ reduction and after the BB,
	 * in the seek procedures.
	 */
	RealVec m_BoundL2;

	/*
	 * NOT USED.
	 *
	 * Whenever the number of nodes in the branch-and-bound tree exceeds
	 * <tt>SHORT_DIET</tt> in the method `ShortestVector`, `PreRedDieterSV` is
	 * automatically set to `true` for the next call; otherwise it is set to
	 * `false`.
	 *
	 static const std::int64_t SHORT_DIET = 1000;

	 *
	 * Whenever the number of nodes in the branch-and-bound tree exceeds
	 * <tt>SHORT_LLL</tt> in the method `ShortestVector`, `PreRedLLLSV` is
	 * automatically set to `true` for the next call; otherwise it is set
	 * to `false`.
	 *
	 static const std::int64_t SHORT_LLL = 1000;

	 *
	 * Maximum number of transformations in the method `PreRedDieter`.
	 * After <tt>MAX_PRE_RED</tt> successful transformations have been
	 * performed, the prereduction is stopped.
	 * 
	 * */
	 static const std::int64_t MAX_PRE_RED = 1000000;
	

	/**
	 * Whenever the number of nodes in the branch-and-bound tree exceeds
	 * <tt>MINK_LLL</tt> in the method <tt>reductMinkowski</tt>,
	 * `PreRedLLLMink` is automatically set to `true` for the next call;
	 * otherwise it is set to `false`.
	 */
	static const std::int64_t MINK_LLL = 500000;

	/**
	 * Pre-reduction flag for `reductMinkowski`.
	 * When true, LLL is performed automatically at certain steps of the reduction.
	 */
	static bool PreRedLLLMink;

	/*
	 * Some local working variables for this class.
	 * They are used inside the basis reduction and short vector methods, and
	 * are declared here to avoid passing them as parameters across the methods.
	 */
	Int m_bs;
	IntVec m_bv;       // Saves current shorter vector in primal basis
	IntVec m_bw;       // Saves current shorter vector in dual basis
	Real m_lMin;       // The norm of the shortest vector in the primal basis
					   // according to the norm considered
	Real m_lMin2;  // Squared L2-norm of the shortest vector in the primal basis
	Real m_ns;
	RealVec m_nv;

	// RealRed m_rs;
	RealRedVec m_zLR, m_n2, m_dc2;
	RealRedMat m_c0, m_c2, m_cho2, m_gramVD;
	int *m_IC;             // Indices in Cholesky

	std::vector<std::int64_t> m_zLI;
	std::vector<std::int64_t> m_zShort;
	std::int64_t m_countNodes;  // Counts number of nodes in the BB tree
	std::int64_t m_countDieter; // Number of attempts since last successful
								// Dieter transformation
	std::int64_t m_cpt;  // Number of successes in pre-reduction transformations
	bool m_foundZero;    // = true -> the zero vector has been found

};
// End class Reducer



/// \cond specReducerSpec
// Structure specialization

// Specialization for the case where Int is int64_t
//
template<typename Real, typename RealRed>
struct specReducer<std::int64_t, Real, RealRed> {

	void redLLLNTLExact(Reducer<std::int64_t, Real, RealRed> &red,
			double delta) {
		std::cout
				<< "** WARNING: redLLLNTLExact cannot be used with std::int64_t.\n";
		std::cout << "** We are now using redLLL instead.\n";
		red.redLLL(delta, 1000000, red.getIntLatticeBase()->getDim());
	}

	void redBKZ(Reducer<std::int64_t, Real, RealRed> &red, double delta,
			std::int64_t blocksize, PrecisionType precision, int dim) {
		IntLatticeBase<std::int64_t, Real, RealRed> *lattmp = 0;
		if (dim > 0) {
			lattmp = new IntLatticeBase<std::int64_t, Real,
					RealRed>(dim, red.getIntLatticeBase()->getNormType());
			lattmp->copyLattice(*red.getIntLatticeBase(), dim);
		} else
			lattmp = red.getIntLatticeBase();
		std::cout
				<< "\n** WARNING: redBKZ cannot be used with std::int64_t integers;\n";
		std::cout
				<< "** it requires the ZZ type. We are now using redLLL instead.\n";
		std::cout << "** It does not do the same thing.\n";
		std::cout << std::endl;
		red.redLLL(delta, 1000000, red.getIntLatticeBase()->getDim());
		if (dim > 0)
			delete lattmp;
	}

	/* This is a test to see if making the program promote to NTL types and
	 * performs NTL LLL reduction is faster than our own LLL.
	 * If it is, we should re-implement our LLL to match what is done in NTL.
	 */
	void redLLLNTL(Reducer<std::int64_t, Real, RealRed> &red, double delta,
			PrecisionType precision, int dim) {
		// Here we create a new temporary lattice, new NTL matrices, and then destroy them.
		// That's a lot of object creations, probably very ineffective.
		IntLatticeBase<NTL::ZZ, Real, RealRed> *lattmp = 0;
		NTL::matrix<std::int64_t> basis = red.getIntLatticeBase()->getBasis();
		// We copy the basis in B, for all the dimensions!
		NTL::mat_ZZ B;
		B.SetDims(basis.NumRows(), basis.NumCols());
		for (int i = 0; i < basis.NumRows(); i++) {
			for (int j = 0; j < basis.NumCols(); j++) {
				B[i][j] = basis[i][j];
			}
		}
		lattmp = new IntLatticeBase<NTL::ZZ, Real, RealRed>(B, dim,
				red.getIntLatticeBase()->getNormType());
		B.kill();

		switch (precision) {
		case DOUBLE:
			LLL_FP(lattmp->getBasis(), delta);
			break;
		case QUADRUPLE:
			LLL_QP(lattmp->getBasis(), delta);
			break;
		case XDOUBLE:
			LLL_XD(lattmp->getBasis(), delta);
			break;
		case RR:
			LLL_RR(lattmp->getBasis(), delta);
			break;
		//case EXACT:
		//	break;
		default:
			MyExit(1, "Undefined PrecisionType for LLL");
		}
		// Here we replace the old basis for all the dimensions!  This looks incorrect.
		for (int i = 0; i < basis.NumRows(); i++) {
			for (int j = 0; j < basis.NumCols(); j++) {
				red.getIntLatticeBase()->getBasis()[i][j] = NTL::trunc_long(
						lattmp->getBasis()[i][j], 63);
				red.getIntLatticeBase()->getBasis()[i][j] *= NTL::sign(
						lattmp->getBasis()[i][j]);
			}
		}
		red.getIntLatticeBase()->updateVecNorm();
		delete lattmp;
	}

	/*
	 * void redLLLNTL(Reducer<std::int64_t, Real, RealRed>& red,
	 *     double delta, PrecisionType precision, int dim)
	 * {
	 *   IntLatticeBase<std::int64_t, Real, RealRed> *lattmp = 0;
	 *   if(dim > 0){
	 *     lattmp = new IntLatticeBase<std::int64_t, Real, RealRed>(
	 *                dim, red.getIntLatticeBase()->getNorm(Type));
	 *     lattmp->copyLattice(*red.getIntLatticeBase(), dim);
	 *   }
	 *   else
	 *     lattmp = red.getIntLatticeBase();
	 *   std::cout << "\n**** WARNING redLLLNTL cannot be use with std::int64_t integers\n";
	 *   std::cout << "** It requires the ZZ type. Instead, LLL reduction is performed\n";
	 *   std::cout << "** with our algorithm, which can be slower.\n";
	 *   std::cout << std::endl;
	 *   red.redLLL(delta, 1000000, red.getIntLatticeBase()->getDim ());
	 *   if (dim>0) delete lattmp;
	 * }
	 * */
};

// Specialization for the case where Int is NTL::ZZ
//
template<typename Real, typename RealRed>
struct specReducer<NTL::ZZ, Real, RealRed> {

	void redLLLNTLExact(Reducer<NTL::ZZ, Real, RealRed> &red, double delta) {
		NTL::ZZ det(0);
		int64_t denum;
		denum = round(1.0 / (1.0 - delta)); // We want (denum-1)/denum \approx delta.
		NTL::LLL(det, red.getIntLatticeBase()->getBasis(), denum - 1, denum);
	}

	void redBKZ(Reducer<NTL::ZZ, Real, RealRed> &red, double delta,
			std::int64_t blocksize, PrecisionType precision, int dim) {
		IntLatticeBase<NTL::ZZ, Real, RealRed> *lattmp = 0;
		if (dim > 0) {
			lattmp = new IntLatticeBase<NTL::ZZ, Real, RealRed>(dim,
					red.getIntLatticeBase()->getNormType());
			lattmp->copyLattice(*red.getIntLatticeBase(), dim);
		} else
			lattmp = red.getIntLatticeBase(); // dim=0 means we use the full dimension.

		switch (precision) {
		case DOUBLE:
			NTL::BKZ_FP(lattmp->getBasis(), delta, blocksize);
			break;
		case QUADRUPLE:
			NTL::BKZ_QP(lattmp->getBasis(), delta, blocksize);
			break;
		case XDOUBLE:
			NTL::BKZ_XD(lattmp->getBasis(), delta, blocksize);
			break;
		case RR:
			NTL::BKZ_RR(lattmp->getBasis(), delta, blocksize);
			break;
		default:
			MyExit(1, "Undefined precision type for BKZ");
		}
		red.getIntLatticeBase()->copyLattice(*lattmp, dim);
		if (dim > 0)
			delete lattmp;
	}

	void redLLLNTL(Reducer<NTL::ZZ, Real, RealRed> &red, double delta,
			PrecisionType precision, int dim) {
		IntLatticeBase<NTL::ZZ, Real, RealRed> *lattmp = 0;
		if (dim > 0) {
			// We should copy only the basis matrix, not the whole IntLatticeBase object !   ******
			lattmp = new IntLatticeBase<NTL::ZZ, Real, RealRed>(dim,
					red.getIntLatticeBase()->getNormType());
			lattmp->copyLattice(*red.getIntLatticeBase(), dim);
		} else
			lattmp = red.getIntLatticeBase();
		//if (precision == EXACT)
		//	red.redLLLNTLExact(delta);
		//else {
			switch (precision) {
			case DOUBLE:
				LLL_FP(lattmp->getBasis(), delta, 0, 0);
				break;
			case QUADRUPLE:
				LLL_QP(lattmp->getBasis(), delta, 0, 0);
				break;
			case XDOUBLE:
				LLL_XD(lattmp->getBasis(), delta, 0, 0);
				break;
			case RR:
				LLL_RR(lattmp->getBasis(), delta, 0, 0);
				break;
			//case EXACT:
			//	break;
			default:
				MyExit(1, "LLL PrecisionType:   NO SUCH CASE");
			}
			red.getIntLatticeBase()->copyLattice(*lattmp, dim);
		//}
		if (dim > 0)
			delete lattmp;

	}
};
/// \endcond

//===========================================================================

// Initialization of non-const static members
//
template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::PreRedLLLMink = false;
template<typename Int, typename Real, typename RealRed>
std::int64_t Reducer<Int, Real, RealRed>::maxNodesBB = 1000000000;

//=========================================================================

template<typename Int, typename Real, typename RealRed>
Reducer<Int, Real, RealRed>::Reducer(IntLatticeBase<Int, Real, RealRed> &lat) {
	// Squared length of vectors must not overflow max(double)
	m_lat = &lat;
	const int dim1 = m_lat->getDim();
	int dim2 = dim1;
	if (dim2 <= 2)
		dim2++;

	m_c0.resize(dim1, dim1);
	m_c2.resize(dim1, dim1);
	m_cho2.resize(dim2, dim2);
	m_gramVD.resize(dim2, dim2);
	// Indices in Cholesky go as high as m_lat->getMaxDim() + 2
	m_IC = new int[2 + dim1];

	m_nv.resize(dim1);
	m_bv.resize(dim1);
	m_bw.resize(dim1);
	m_n2.resize(dim1);
	m_zLR.resize(dim1);
	m_zLI.resize(dim1);
	m_zShort.resize(dim1);
	m_dc2.resize(dim1);
	m_BoundL2.resize(dim1);

	m_lMin = std::numeric_limits<double>::max();
	m_lMin2 = m_lMin;
	for (int i = 0; i < dim1; i++) {
		m_zLI[i] = -1;
		m_zShort[i] = -1;
		m_BoundL2[i] = -1;
		m_IC[i] = -1;
	}
	m_countNodes = 0;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
Reducer<Int, Real, RealRed>::Reducer(const Reducer<Int, Real, RealRed> &red) {
	copy(red);
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
Reducer<Int, Real, RealRed>&
Reducer<Int, Real, RealRed>::operator=(const Reducer<Int, Real, RealRed> &red) {
	if (this != &red)
		copy(red);
	return *this;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::copy(const Reducer<Int, Real, RealRed> &red) {
	m_lat = red.m_lat;
	m_c0 = red.m_c0;
	m_c2 = red.m_c2;
	m_dc2 = red.m_dc2;
	m_nv = red.m_nv;
	m_bv = red.m_bv;
	m_bw = red.m_bw;
	m_n2 = red.m_n2;
	m_zLR = red.m_zLR;
	m_zLI = red.m_zLI;
	m_zShort = red.m_zShort;
	m_cho2 = red.m_cho2;
	m_gramVD = red.m_gramVD;
	m_lMin = red.m_lMin;
	m_lMin2 = red.m_lMin2;
	m_BoundL2 = red.m_BoundL2;
	if (m_IC != 0)
		delete[] m_IC;
	m_IC = new int[3 + m_lat->getDim()];
	for (int i = 0; i < 2 + m_lat->getDim(); i++)
		m_IC[i] = red.m_IC[i];
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
Reducer<Int, Real, RealRed>::~Reducer() {
	m_c0.clear();
	m_c2.clear();
	m_cho2.clear();
	m_gramVD.clear();
	m_nv.clear();
	m_bv.clear();
	m_bw.clear();
	m_n2.clear();
	m_zLR.clear();
	m_zLI.clear();
	m_zShort.clear();
	m_dc2.clear();
	m_BoundL2.clear();
	delete[] m_IC;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::setBoundL2(const RealVec &thresholds,
		int dim1, int dim2) {
	m_BoundL2.resize(dim2);
	for (int i = dim1; i < dim2; i++)
		m_BoundL2[i] = thresholds[i];
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::permuteGramVD(int i, int j, int n) {
	int k;
	for (k = 0; k < n; k++) {
		std::swap(m_gramVD(i, k), m_gramVD(j, k));
	}
	for (k = 0; k < n; k++) {
		std::swap(m_gramVD(k, i), m_gramVD(k, j));
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::calculCholesky2LLL(int n, int j) {
	// Recalcule les n premieres entrees de la colonne j de la matrice de
	// Cholesky d'ordre 2.   Pour redLLL.
	m_cho2(0, j) = m_gramVD(0, j);
	for (int i = 1; i <= n; i++) {
		m_cho2(i, j) = m_gramVD(i, j);
		for (int k = 0; k < i; k++) {
			m_cho2(i, j) -= (m_cho2(k, j) / m_cho2(k, k)) * m_cho2(k, i);
		}
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
inline void Reducer<Int, Real, RealRed>::calculCholesky2Ele(int i, int j) {
	// Recalcule l'entree (i,j) de la matrice de Cholesky d'ordre 2.
	m_cho2(i, j) = m_gramVD(i, j);
	for (int k = 0; k < i; k++) {
		m_cho2(i, j) -= m_cho2(k, i) * (m_cho2(k, j) / m_cho2(k, k));
	}
}

//=========================================================================

//  void negativeCholesky();   ???

//=========================================================================

template<typename Int, typename Real, typename RealRed>
inline void Reducer<Int, Real, RealRed>::calculGramVD() {
	// Retourne dans m_gramVD la matrice des produits scalaires m_lat->V[i]*m_lat->V[j].
	// Note: m_lat->V.vecNorm ne contient que les m_lat->V[i]*m_lat->V[i].
	const int dim = m_lat->getDim();
	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
			NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
			ProdScal<Int>(row1, row2, dim, m_gramVD(i, j));
			m_gramVD(j, i) = m_gramVD(i, j);
		}
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
inline void Reducer<Int, Real, RealRed>::miseAJourGramVD(int j) {
	// Recalcule la j-ieme ligne et la j-ieme colonne de la matrice des
	// produits scalaires.  Pour redLLL.
	const int dim = m_lat->getDim();
	for (int i = 0; i < dim; i++) {
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
		NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
		ProdScal<Int>(row1, row2, dim, m_gramVD(i, j));
		m_gramVD(j, i) = m_gramVD(i, j);
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::calculCholesky(RealRedVec &DC2,
		RealRedMat &C0) {
	/*
	 * Returns in C0 the elements of the upper triangular matrix of the
	 * Cholesky decomposition that are above the diagonal. Returns in DC2 the
	 * squared elements of the diagonal. These elements are rescaled by EchVV
	 * when SISquares= true.
	 */
	const int dim = m_lat->getDim();
	int k, j, i;
	RealRed m2;
	// C2(i,j) = C0(i,j) * C2(i,i) if i != j.
	// C2(i,i) = DC2[i].
	NTL::conv(m2, m_lat->getModulo());
	m2 = m2 * m2;
	int d = dim;
	//if(m_lat->withDual())
	//  d = dim / 2; // If we use the Dual, we compute Cholesky
	// with the Dual


	// Compute the d first lines of C0 with the primal Basis.
	for (i = 0; i < d; i++) {
		m_lat->updateScalL2Norm(i);
		for (j = i; j < dim; j++) {
			if (j == i)
				NTL::conv(m_c2(i, i), m_lat->getVecNorm(i));
			else {
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
				NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
				ProdScal<Int>(row1, row2, dim, m_c2(i, j));
			}
			
			for (k = 0; k < i; k++)
				m_c2(i, j) -= C0(k, i) * m_c2(k, j);
			if (i == j) {
				DC2[i] = m_c2(i, i);
				if (DC2[i] < 0.0) {
					//negativeCholesky();
					return false;
				}
			} else
				C0(i, j) = m_c2(i, j) / DC2[i];
		//add for test		
		/**  if(i!=j && i<j)		
	       std::cout<< "C0("<<i<<","<<j<<")="<<C0(i,j)<<" ";	
		  else if (i==j) 	
		    std::cout<< "C0("<<i<<","<<j<<")="<<DC2[i]<<" ";	
		  else
		     std::cout<< "C0("<<i<<","<<j<<")="<<"0"<<" ";	*/		
		}
	  // std::cout<<""<<std::endl;
	}

	// Compute the d last lines of C0 with the dual Basis.
	/* This operation with the dual is needed in case of high dimension
	 * and large number (30 bits). The choleski decomposition
	 * convert numbers in double which is not sufficient in that case.
	 * You need to use RR of the NTL library for this calculation.
	 */
	//if(m_lat->withDual()){
	//  for (i = dim-1; i >= d; i--)
	//  {
	//    m_lat->updateDualScalL2Norm (i);
	//    for (j = i; j >= 0; j--) {
	//      if (j == i)
	//        NTL::conv (m_c2(i,i), m_lat->getDualVecNorm (i));
	//      else {
	//        NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
	//        NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
	//        ProdScal<Int> (row1, row2, dim, m_c2(i,j));
	//      }
	//      for (k = i + 1; k < dim; k++){
	//        m_c2(i,j) -= C0(k,i) * m_c2(k,j);
	//      }
	//      if (i != j)
	//        C0(i,j) = m_c2(i,j) / m_c2(i,i);
	//    }
	//    DC2[i] = m2 / m_c2(i,i);
	//    if (DC2[i] < 0.0) {
	//      negativeCholesky();
	//      return false;
	//    }
	//    for (j = i + 1; j < dim; j++) {
	//      C0(i,j) = -C0(j,i);
	//      for (k = i + 1; k < j; k++) {
	//        C0(i,j) -= C0(k,i) * C0(k,j);
	//      }
	//    }
	//  }
	//}
	return true;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::tracePrintBases(char *message) {
	std::cout << std::endl << "================================= " << message
			<< std::endl;
	//std::cout << "dim = " << m_lat->getDim () << std::endl;
	m_lat->setNegativeNorm();
	//m_lat->setDualNegativeNorm();
	m_lat->updateVecNorm();
	//m_lat->updateDualVecNorm();
	m_lat->sort(0);
	m_lat->write();
}

//===========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::pairwiseRedPrimal(int i, int d, bool xx[]) {
	const int dim = m_lat->getDim();
	++m_countDieter;
	m_lat->updateScalL2Norm(i);
	bool modifFlag;

	for (int j = d; j < dim; j++) {
		if (i == j)
			continue;
		modifFlag = false;
		{
			NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
			NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
			ProdScal<Int>(row1, row2, dim, m_ns);
		}
		DivideRound<Real>(m_ns, m_lat->getVecNorm(i), m_ns);
		if (m_ns == 0)
			continue;
		NTL::conv(m_bs, m_ns);
		if (m_ns < 1000 && m_ns > -1000) {
			m_lat->updateScalL2Norm(j);
			{
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
				NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
				ModifVect(row1, row2, -m_bs, dim);
			}

			// Verify that m_lat->getBasis()[j] is really shorter
			{
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
				ProdScal<Int>(row1, row1, dim, m_ns);
			}
			if (m_ns >= m_lat->getVecNorm(j)) {
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
				NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
				ModifVect(row1, row2, m_bs, dim);
			} else {
				modifFlag = true;
				m_lat->setVecNorm(m_ns, j);
			}
		} else {
			NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
			NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
			ModifVect(row1, row2, -m_bs, dim);
			m_lat->setNegativeNorm(j);
			modifFlag = true;
		}

		if (modifFlag) {
			m_countDieter = 0;
			++m_cpt;
			if (m_lat->withDual()) {
				NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
				NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
				ModifVect(row1, row2, m_bs, dim);
				m_lat->setDualNegativeNorm(i);

			}
			if (xx) {
				xx[i] = false;
				xx[j] = false;
			}
		}
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::pairwiseRedDual(int i, bool xx[]) {
	int j;
	const int dim = m_lat->getDim();

	++m_countDieter;
	m_lat->updateDualScalL2Norm(i);
	NTL::matrix_row<IntMat> row9(m_lat->getBasis(), i);
	m_bv = row9;
	for (j = 0; j < dim; j++) {
		if (i != j) {
			NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
			NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
			ProdScal<Int>(row1, row2, dim, m_ns);
			// ProdScal<Int> (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
			//           dim, m_ns);
			DivideRound<Real>(m_ns, m_lat->getDualVecNorm(i), m_nv[j]);
			if (m_nv[j] != 0) {
				NTL::conv(m_bs, m_nv[j]);
				NTL::matrix_row<IntMat> row7(m_lat->getBasis(), j);
				ModifVect(m_bv, row7, m_bs, dim);
			}
		}
	}

	m_lat->updateScalL2Norm(i);
	ProdScal<Int>(m_bv, m_bv, dim, m_ns);
	if (m_ns < m_lat->getVecNorm(i)) {
		++m_cpt;
		m_countDieter = 0;
		NTL::matrix_row<IntMat> row6(m_lat->getBasis(), i);
		for (j = 0; j < dim; j++)
			row6(j) = m_bv[j];
		m_lat->setNegativeNorm(i);
		if (xx)
			xx[i] = false;
		m_lat->setVecNorm(m_ns, i);
		for (j = 0; j < dim; j++) {
			if (i != j && m_nv[j] != 0) {
				NTL::conv(m_bs, -m_nv[j]);
				NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), j);
				NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), i);
				ModifVect(row1, row2, m_bs, dim);
				//  ModifVect (m_lat->getDualBasis ()[j], m_lat->getDualBasis ()[i],
				//            m_bs, dim);
				m_lat->setDualNegativeNorm(j);
				if (xx)
					xx[j] = false;
			}
		}
	}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redDieter(int d, bool xx[]) {
	std::int64_t BoundCount;
	const int dim = m_lat->getDim();
	bool withDual = m_lat->withDual();

	m_lat->updateScalL2Norm(d, dim);
	m_lat->sort(d);
	int i = dim - 1;
	m_cpt = 0;
	m_countDieter = 0;
	BoundCount = 2 * dim - d;
	do {
		pairwiseRedPrimal(i, d, xx);
		if (i > d && withDual)
			pairwiseRedDual(i, xx);
		if (i < 1)
			i = dim - 1;
		else
			--i;
	} while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redDieterRandomized(int d, int seed) {
	std::int64_t BoundCount;
	const int dim = m_lat->getDim();
	bool withDual = m_lat->withDual();

	m_lat->updateScalL2Norm(d, dim);
	//m_lat->getDualBasis ().updateScalL2Norm (d, dim);
	m_lat->sort(d);
	int i = dim - 1;
	m_cpt = 0;
	m_countDieter = 0;
	BoundCount = 2 * dim - d;
	srand(seed);   // We use a simple RNG from the standard C library.
	do {
		pairwiseRedPrimal(rand() % dim, d);
		if (i > d && withDual)
			pairwiseRedDual(i);
		if (i < 1)
			i = dim - 1;
		else
			--i;
	} while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
}

//=========================================================================

/**
 * Reduce the Cholesky matrix with adding a multiple of the i-th vector
 * to the j-th vector. It updates the Gram Schmidt matrix.
 */
template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::reductionFaible(int i, int j)
/*
 * Reduit la matrice de Cholesky (d'ordre 2) en ajoutant un multiple du
 * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
 * consequence et remet a jour la matrice des produits scalaires.
 * Utilise par redLLL.
 */
{
	RealRed cte;
	std::int64_t cteLI;
	//bool withDual = m_lat->withDual();
	cte = m_cho2(i, j) / m_cho2(i, i);

	const int dim = m_lat->getDim();

	if (abs(cte) < std::numeric_limits<double>::max()) {
		// On peut representer cte en LONGINT.
		if (abs(cte) > 0.5) {
			NTL::conv(cteLI, Round(cte));
			NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
			NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
			ModifVect(row1, row2, -cteLI, dim);
			//if(withDual){
			//  NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
			//  NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
			//  ModifVect (row3, row4, cteLI, dim);
			//}

		} else
			return;

	} else {
		// On represente cte en double.
		if (abs(cte) < std::numeric_limits<long double>::max())
			cte = Round(cte);
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
		NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
		ModifVect(row1, row2, -cte, dim);
		//if(withDual){
		//  NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
		//  NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
		//  ModifVect (row3, row4, cte, dim);
		//}

	}
	m_lat->setNegativeNorm(j);
	m_lat->updateVecNorm(j);
	//if(withDual){
	//  m_lat->setDualNegativeNorm (i);
	//  m_lat->updateDualVecNorm(i);
	//}

	miseAJourGramVD(j);
	calculCholesky2LLL(i, j);
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redLLLNTLExact(double delta) {
	spec.redLLLNTLExact(*this, delta);
}

//=========================================================================
//std::int64_t 
template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redBKZ(double delta, int blocksize,
		PrecisionType precision, int dim) {
	spec.redBKZ(*this, delta, blocksize, precision, dim);
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redLLLNTL(double delta,
		PrecisionType precision, int dim) {
	spec.redLLLNTL(*this, delta, precision, dim);
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::redLLL(double delta, std::int64_t maxcpt,
		int Max) {
	// Effectue la pre-reduction de B au sens de Lenstra-Lenstra-Lovasz. N'utilise
	// pas les vecteurs m_lat->getBasis().vecNorm, Wm_lat->getDualBasis().

	//bool withDual = m_lat->withDual();
	const int REDBAS_e = 40;
	int i, j, k, h;
	RealRed Cho0ij;
	RealRed limite;
	std::int64_t cpt;

	const int dim = m_lat->getDim();
	if (Max == 0)
		Max = dim;

	cpt = 0;
	calculGramVD();
	limite = 1.0;
	for (k = 1; k <= REDBAS_e; k++)
		limite *= 2.0;
	limite *= dim;
	m_cho2(0, 0) = m_gramVD(0, 0);
	m_cho2(0, 1) = m_gramVD(0, 1);
	m_IC[0] = 1;
	m_cho2(1, 1) = m_gramVD(1, 1)
			- m_cho2(0, 1) * (m_cho2(0, 1) / m_cho2(0, 0));
	m_IC[1] = 1;
	for (i = 2; i < dim; i++)
		m_IC[i] = -1;
	h = 0;

	while (h < Max - 1 && cpt < maxcpt) {
		if (m_gramVD(h + 1, h + 1) > limite) {
			for (i = h; i >= 0; i--)
				reductionFaible(i, h + 1);
		} else
			reductionFaible(h, h + 1);

		calculCholesky2Ele(h + 1, h + 1);
		if (m_IC[h + 1] == -1)
			m_IC[h + 1] = h + 1;

		if (m_cho2(h + 1, h + 1) / m_cho2(h, h)
				+ (m_cho2(h, h + 1) / m_cho2(h, h))
						* (m_cho2(h, h + 1) / m_cho2(h, h)) < delta) {
			++cpt;

			m_lat->permuteNoDual(h, h + 1);
			permuteGramVD(h, h + 1, dim);
			m_cho2(h, h) = m_gramVD(h, h);
			for (i = 0; i < h; i++) {
				std::swap(m_cho2(i, h), m_cho2(i, h + 1));
				m_cho2(h, h) -= m_cho2(i, h) * (m_cho2(i, h) / m_cho2(i, i));
			}
			if (h == 0) {
				Cho0ij = m_cho2(0, 1) / m_cho2(0, 0);
				if (abs(Cho0ij) > 0.5) {
					m_IC[0] = 1;
					m_IC[1] = -1;
					h = 0;
				} else {
					m_cho2(1, 1) = m_gramVD(1, 1)
							- m_cho2(0, 1) * m_cho2(0, 1) / m_cho2(0, 0);
					calculCholesky2LLL(2, 2);
					m_IC[0] = 2;
					m_IC[1] = 2;
					m_IC[2] = 2;
					h = 1;
				}
			} else {
				m_IC[h] = h + 1;
				m_IC[h + 1] = -1;
				--h;
			}

		} else {
			for (i = 0; i <= h + 2; i++) {
				if (h + 2 > m_IC[i]) {
					if (h + 2 < dim)
						calculCholesky2Ele(i, h + 2);
					m_IC[i] = h + 2;
				}
			}
			++h;
		}
	}

	if (cpt == maxcpt) {
		std::cout << "***** in redLLL cpt > maxcpt = " << maxcpt << std::endl;
	}

	for (j = 2; j < Max; j++) {
		for (i = j - 2; i >= 0; i--)
			reductionFaible(i, j);
	}
	m_lat->setNegativeNorm();
	//if(withDual){
	//  m_lat->setDualNegativeNorm ();
	//}
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::transformStage3Mink(
		std::vector<std::int64_t> &z, int &k) {
	int j, i;
	std::int64_t q;
	const int dim = m_lat->getDim();
	bool withDual = m_lat->withDual();

	j = dim - 1;
	while (z[j] == 0)
		--j;
	while (abs(z[j]) > 1) {
		i = j - 1;
		while (z[i] == 0)
			--i;
		// On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
		while (z[j]) {
			// Troncature du quotient vers 0
			q = z[i] / z[j];
			if (q) {
				// On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
				z[i] -= q * z[j];
				NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
				//    ModifVect (m_lat->getBasis ()[j], m_lat->getBasis ()[i],
				//            q, dim);
				ModifVect(row1, row2, q, dim);
				m_lat->setNegativeNorm(j);

				if (withDual) {
					NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
					NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
					//    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
					//             -q, dim);
					ModifVect(row3, row4, -q, dim);
					m_lat->setDualNegativeNorm(i);
				}
			}
			// Permutation.
			std::swap < std::int64_t > (z[i], z[j]);
			m_lat->permute(i, j);
		}
		j = i;
	}
	k = j;
}

//=========================================================================

/*
 * We call this when we have found a shorter vector which is a linear combination of the
 * previous basis vectors with coefficients given in z.  This procedure updates the basis
 * so that the new shortest vector in first and the other vectors are adjusted accordingly.
 * Note that an alternative for the L2 norm might be to add the new shortest vector to
 * the basis to get dim+1 generating vectors, and apply LLL to recover a new basis.
 */
template<typename Int, typename Real, typename RealRed>
void Reducer<Int, Real, RealRed>::transformStage3ShortVec(
		std::vector<std::int64_t> &z, int &k) {
	int j, i;
	std::int64_t q;
	const int dim = m_lat->getDim();
	j = dim - 1;
	while (z[j] == 0)
		--j;
	while (abs(z[j]) > 1) {
		i = j - 1;
		while (z[i] == 0)
			--i;
		// On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
		while (z[j]) {
			// Troncature du quotient vers 0
			q = z[i] / z[j];
			if (q) {
				// On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
				z[i] -= q * z[j];
				NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
				NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
				ModifVect(row1, row2, q, dim);
				m_lat->setNegativeNorm(j);
			}
			// Permutation.
			std::swap < std::int64_t > (z[i], z[j]);
			m_lat->permuteNoDual(i, j);
		}
		j = i;
	}
	k = j;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::tryZMink(int j, int i, int Stage,
		bool &smaller, const IntMat &WTemp)
		// If m_countNodes > MaxNodesBB, returns false, otherwise returns true.
		{
	std::int64_t max0, min0;
	RealRed x, dc;
	RealRed center;
	std::int64_t zhigh, zlow, h;
	bool high;
	int k;
	RealRed S1, S2, S3, S4, mR;

	const int dim = m_lat->getDim();
	NTL::conv(mR, m_lat->getModulo());

	++m_countNodes;
	if (m_countNodes > maxNodesBB) {
		std::cerr << "-------- m_countNodes > maxNodesBB = " << maxNodesBB
				<< std::endl;
		return false;
	}

	// Calcul d'un intervalle contenant les valeurs admissibles de zj.
	center = 0.0;
	if (j < dim - 1) {
		// Calcul du centre de l'intervalle.
		for (k = j + 1; k < dim; k++)
			center = center - m_c0(j, k) * m_zLR[k];

		// Distance du centre aux extremites de l'intervalle.
		// We use the L2 norm for this, since other norms are not allowed.
		dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);

		/* Calcul de deux entiers ayant la propriete qu'un entier */
		/* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
		/* n'appartient pas a l'intervalle.  */
		x = center - dc;
		NTL::conv(min0, trunc(x));
		if (x > 0.0)
			++min0;

		x = center + dc;
		NTL::conv(max0, trunc(x));
		if (x < 0.0)
			--max0;

		// En vue du choix initial de zj. On determine zlow et zhigh.
		if (min0 > max0) {
			return true;
		}
		if (min0 == max0) {
			zlow = min0;
			zhigh = max0 + 1;
			high = false;
		} else if (center >= 0.0) {
			NTL::conv(zlow, trunc(center));
			zhigh = zlow + 1;
			NTL::conv(h, trunc(2.0 * center));
			high = h & 1;
		} else {
			NTL::conv(zhigh, trunc(center));
			zlow = zhigh - 1;
			NTL::conv(h, -trunc(2.0 * center));
			high = (h & 1) == 0;
		}

	} else {  // j = dim-1

		zlow = 0;
		high = true;
		if (Stage == 2) {
			min0 = 1;
			max0 = 1;
			zhigh = 1;
		} else {
			min0 = 2;
			zhigh = 2;
			NTL::conv(max0, trunc(sqrt((m_lMin2 - m_n2[j]) / m_dc2[j])));
		}
	}

	Real temp;
	/* On essaie maintenant chacun des z[j] dans l'intervalle, en      */
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

		if (j == 0) {
			RealRed tmps_n2 = m_n2[0] + x * x * m_dc2[0];
			if (tmps_n2 < m_lMin2) {
				// On verifie si on a vraiment trouve un vecteur plus court
				NTL::matrix_row<const IntMat> row1(m_lat->getBasis(), dim - 1);
				m_bv = row1;
				for (k = 0; k < dim - 1; k++) {
					if (m_zLI[k] != 0) {
						NTL::matrix_row<const IntMat> row1(m_lat->getBasis(),
								k);
						ModifVect(m_bv, row1, m_zLI[k], dim);
					}
				}
				if (Stage == 3) {
					NTL::matrix_row<const IntMat> row1(m_lat->getBasis(),
							dim - 1);
					ModifVect(m_bv, row1, m_zLR[dim - 1] - 1.0, dim);
				}

				ProdScal<Int>(m_bv, m_bv, dim, S1);
				NTL::conv(S4, m_lat->getVecNorm(dim - 1));
				if (S1 < S4) {
					if (Stage == 2) {
						smaller = true;
						if (!PreRedLLLMink)
							m_zShort = m_zLI;
						else {
							for (k = 1; k < dim; k++) {
								NTL::matrix_row<const IntMat> row1(WTemp, k);
								ProdScal<Int>(m_bv, row1, dim, S2);
								Quotient(S2, mR, S3);
								NTL::conv(m_zShort[k], S3);
							}
							m_zShort[dim - 1] = 1;
						}
					} else if (Stage == 3 && !PreRedLLLMink) {
						if (GCD2vect(m_zLI, i, dim) == 1) {
							m_zShort = m_zLI;
							smaller = true;
						}
					} else {
						for (k = 0; k < dim; k++) {
							NTL::matrix_row<const IntMat> row1(WTemp, k);
							ProdScal<Int>(m_bv, row1, dim, S2);
							Quotient(S2, mR, S3);
							NTL::conv(m_zShort[k], S3);
						}
						if (GCD2vect(m_zShort, i, dim) == 1) {
							smaller = true;
						}
					}
					if (smaller) {
						NTL::conv(temp, S1);
						m_lat->setVecNorm(temp, dim - 1);
						return true;
					}
				}
			}
		} else { // j > 0
			m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
			if (m_lMin2 >= m_n2[j - 1]) {
				if (!tryZMink(j - 1, i, Stage, smaller, WTemp))
					return false;
				// Des qu'on a trouve quelque chose, on sort de la recursion
				// et on retourne dans reductMinkowski.
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
	return true;
}

//=========================================================================

template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::redBBMink(int i, int d, int Stage,
		bool &smaller, bool xx[])
		/*
		 * Tries to shorten m_lat->getBasis()[i] using branch-and-bound.
		 * Used in Minkowski Reduction.
		 * Stage is 2 or 3.
		 * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
		 * Stops and returns false if not finished after examining MaxNodesBB
		 * nodes in the branch-and-bound tree.  When succeeds, returns true.
		 * Assumes that the norm is Euclidean.
		 */
		{
	bool withDual = m_lat->withDual();
	const int dim = m_lat->getDim();
	IntMat VTemp(dim, dim), WTemp(dim, dim);
	bool XXTemp[dim];
	Real tmp;
	// trace( "AVANT redBBMink");
	smaller = false;

	// Approximation du carre de la longueur de Vi.
	if (m_lat->getVecNorm(i) < 0) {
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
		ProdScal<Int>(row1, row1, dim, tmp);
		//  ProdScal<Int> (m_lat->getBasis()[i], m_lat->getBasis()[i],
		//            dim, tmp);
		m_lat->setVecNorm(tmp, i);
	}
	NTL::conv(m_lMin2, m_lat->getVecNorm(i));

	if (Stage == 3 && withDual) {
		if (m_lat->getDualVecNorm(i) < 0) {
			NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
			ProdScal<Int>(row1, row1, dim, tmp);
			//   ProdScal<Int> (m_lat->getDualBasis()[i], m_lat->getDualBasis()[i],
			//            dim, tmp);
			m_lat->setDualVecNorm(tmp, i);
		}
		Real m2;
		NTL::conv(m2, m_lat->getModulo());
		m2 = m2 * m2;
		if (m_lMin2 * m_lat->getDualVecNorm(i) < 4 * m2)
			return true; // if the angle between the basis vector i and the dual
		// basis vector i is between -Pi/3 and Pi/3
	}
	if (withDual) {
		m_lat->updateDualVecNorm();
	}
	m_lat->updateVecNorm();
	m_lat->permute(i, dim - 1);

	int k, h;

	if (PreRedLLLMink) {
		// On memorise la base courante.
		VTemp = m_lat->getBasis();
		if (withDual) {
			WTemp = m_lat->getDualBasis();
		}
		for (h = 0; h < dim; h++)
			XXTemp[h] = true;
		redLLL(1.0, 1000000, dim - 1);
		m_lat->updateVecNorm();
		if (withDual) {
			m_lat->updateDualVecNorm();
		}
	}
	if (!calculCholesky(m_dc2, m_c0))
		return false;
	m_countNodes = 0;
	m_n2[dim - 1] = 0.0;
	if (!tryZMink(dim - 1, i, Stage, smaller, WTemp))
		return false;

	if (PreRedLLLMink) {
		/* On remet l'anciennne base, celle d'avant LLL, avant de considerer
		 la prochaine m_lat->dimension.  */
		m_lat->getBasis() = VTemp;
		m_lat->updateVecNorm();
		if (withDual) {
			m_lat->getDualBasis() = WTemp;
			m_lat->updateDualVecNorm();
		}
		for (h = 0; h < dim; h++)
			xx[h] = XXTemp[h];
	}
	if (smaller) {
		/* On a trouve un plus court vecteur.  On ameliore
		 m_lat->getBasis()[k].  */
		if (Stage == 2)
			k = dim - 1;
		else
			transformStage3Mink(m_zShort, k);
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
		for (h = 0; h < dim; h++)
			row1(h) = m_bv[h];
		//  m_lat->getBasis ()[k] = m_bv;
		m_lat->setNegativeNorm(k);
		if (m_zShort[k] < 0 && withDual) {
			NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), k);
			ChangeSign(row2, dim);
		}
		/* Mise a jour des vecteurs de la base duale selon le nouveau
		 m_lat->getBasis()[k] */
		for (h = 0; h < dim; h++) {
			if ((m_zShort[h] != 0) && (h != k)) {
				if (withDual) {
					NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), h);
					NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), k);
					ModifVect(row1, row2, -m_zShort[h], dim);
					m_lat->setDualNegativeNorm(h);
				}
				if (Stage == 2) {
					if (h >= d)
						xx[h] = false;
				}
			}
		}
	} else if (Stage == 2)
		xx[dim - 1] = true;

	m_lat->permute(i, dim - 1);
	// trace( "APRES redBBMink");
	return true;
}

//=========================================================================

/**
* @decomp take std::string Choles("cholesky") or std::string Triang("triangular");
* @norm take value 'L1NORM' or 'L2NORM'
*/
template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::tryZShortVec(int j, bool &smaller, NormType norm, std::string decomp) {
	/*
	 * This recursive procedure uses the Cholesky decomposition and works for the
	 * L1NORM and L2NORM.  It is called initially with j=dim-1 by redBBShortVec
	 * to find a shortest vector via the BB procedure.
	 * If `m_countNodes > MaxNodesBB`, it returns `false`, otherwise it returns `true`,
	 * meaning that a shortest vector has been found.
	 * The square length of the current shortest vector is maintained in `m_lMin2`.
	 * When j=1 and we find a shorter vector, we return the corresponding `z_j`'s in `m_zShort`.
	 */

	/* Pour une implantation non recursive, ces variables devraient
	 etre des tableaux indices par j. */
	RealRed dc, x, center;
	std::int64_t min0, max0;     // Bornes de l'intervalle pour les z_j.
	std::int64_t zlow, zhigh; // Valeur courante a gauche et a droite du centre.
	bool high; // Indicates if we are on the right (true) or the left of the center.
	int k;
	std::int64_t temp;
	//
	std::string Choles("cholesky");
	std::string Triang("triangular");


	const int dim = m_lat->getDim();
	++m_countNodes;
	if (m_countNodes > maxNodesBB) {
		std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB
				<< std::endl;
		return false;
	}
	
	
	/* Compute an interval that contains the admissible values of zj. */
	/* This computation is for the L2 norm, but also works for the L1 norm. */
	/* 1. Compute the center of the interval.  */
	if (decomp==Choles){
	   center = 0.0;
	   for (k = j + 1; k < dim; ++k)
		  center -= m_c0(j, k) * m_zLR[k];
		  
	// Distance from the center to the boundaries.
	// m_lMin2 contains the square length of current shortest vector with the selected norm.
	    dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);
	 }
	 if(decomp==Triang && norm==L2NORM){    
        center = 0.0;
	    for (k = j + 1; k < dim; ++k)
	 		center -= m_c0(k,j) * m_zLR[k];
       	  dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);
	  
	  }	
     if(decomp==Triang && norm==L1NORM){
           
       center = 0.0;
	   for (k = j + 1; k < dim; ++k)
		  center -= m_c0(k, j) * m_zLR[k];
     
	    dc = (m_lMin - m_n2[j])  / m_c0(j,j);
	  
	  }	
	// Compute two integers min0 and max0 that are the min and max integers in the interval.
	if (!m_foundZero)
		min0 = 0;     // We are at the beginning, min will be zero.
	else {
		x = center - dc;
		NTL::conv(min0, trunc(x));
		if (x > 0.0)
			++min0;
	}
	x = center + dc;
	NTL::conv(max0, trunc(x));
	if (x < 0.0)
		--max0;


    //   std::cout << "Borne Min="<<min0<<"     Borne Max="<<max0<<"     j="<<j<<std::endl;
	// Compute initial values of zlow and zhigh, the search pointers on each side of the interval.
	if (min0 > max0)
		return true;
	if (min0 == max0) {
		zlow = min0;
		zhigh = max0 + 1;
		high = false;
	} else if (center >= 0.0) {
		NTL::conv(zlow, trunc(center));
		zhigh = zlow + 1;
		NTL::conv(temp, trunc(2.0 * center));
		high = (temp & 1);
	} else {
		NTL::conv(zhigh, trunc(center));
		zlow = zhigh - 1;
		NTL::conv(temp, -trunc(2.0 * center));
		high = (temp & 1) == 0;
	}

	// We try each zj in the interval, starting in the center and alternating between left and right.

	while (zlow >= min0 || zhigh <= max0) {
		if (high)
			m_zLI[j] = zhigh;
		else
			m_zLI[j] = zlow;
		m_zLR[j] = m_zLI[j];

		// Computing m_n2[j-1].
		x = m_zLR[j] - center;
		if (j == 0) {
			// All the zj have been selected: we have a vector to test.
			if (m_lMin2 > m_n2[0] + x * x * m_dc2[0]) {
				/* Check if we have a shorter nonzero vector. */
				if (!m_foundZero) {
					// Le premier vecteur trouve sera zero.
					m_foundZero = true;
				} else {
					SetZero(m_bv, dim);
					for (k = 0; k < dim; k++) {
						if (m_zLI[k] != 0) {
							NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
							ModifVect(m_bv, row1, m_zLI[k], dim);
						}
					}
					// The new vector is now in m_bv.
					if (m_lat->getNormType() == L2NORM) {
						ProdScal<Int>(m_bv, m_bv, dim, x);
					} else {
						// Compute the square length for the L1 norm.
						CalcNorm<IntVec, RealRed>(m_bv, dim, x,
								m_lat->getNormType());
						x = x * x;
					}
					if (x < m_lMin2) {
						// The new vector is shorter.
						smaller = true;
						NTL::conv(m_lMin2, x);
						m_zShort = m_zLI;
						m_bw = m_bv;
					}
				}
			}
		} else if (m_lMin2 > m_n2[j] + x * x * m_dc2[j]) {
			// There is still hope; we continue the recursion.
			m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
			if (!tryZShortVec(j - 1, smaller,norm, decomp))
				return false;
		} else
			m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
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
	
	return true;
}

//=========================================================================


/**
* @decomp take std::string Choles("cholesky") or std::string Triang("triangular");
* @norm take value 'L1NORM' or 'L2NORM'
*/

template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::redBBShortVec(NormType norm, std::string decomp) {
	/*
	 * Finds shortest non-zero vector, using branch-and-bound, with L1 or L2 norm.
	 * Stops and returns false if not finished after examining MaxNodesBB nodes in the
	 * branch-and-bound tree.  When succeeds, returns true, and the squared shortest
	 * vector length will be in m_lMin2.
	 */
     
    std::string cholesky("cholesky");
	std::string triangular("triangular");


	if ((norm != L1NORM) & (norm != L2NORM)) {
		std::cerr << "RedBBShortVec: only L1 and L2 norms are supported";
		return false;
	}
	bool smaller = false;  // Found a smaller vector?
	int k, h;
	const int dim = m_lat->getDim();  // Lattice dimension
	RealRed x;
    // Here we sort the basis, otherwise Cholesky will fail more rapidly
    // due to floating-point errors.
	m_lat->updateScalL2Norm(0, dim);
	m_lat->sortNoDual(0);

    // Approximate the square norm of the current shortest vector.
	if (norm == L2NORM) {
		NTL::conv(m_lMin2, m_lat->getVecNorm(0));
	} else {
		// Looking for the shortest vector in basis according to the considered norm
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), 0);
		CalcNorm<IntVec, Real>(row1, dim, m_lMin, norm);
		for (k = 1; k < dim; k++) {
			NTL::matrix_row<IntMat> row2(m_lat->getBasis(), k);
			CalcNorm<IntVec, Real>(row2, dim, x, norm);
			if (x < m_lMin)
				m_lMin = x;
		}
		m_lMin2 = m_lMin * m_lMin;  // Squared shortest length with L1 norm.
	}

	/* If we already have a shorter vector than the minimum threshold, we stop right away.
	 * This is useful for the seek programs in LatMRG.sss
	 */
	if (m_lMin2 <= m_BoundL2[dim - 1])
		return false;


    if(decomp==cholesky){
    /* Perform the Cholesky decomposition; if it fails we exit. */
         if (!calculCholesky(m_dc2, m_c0))
           return false;    
    }
    else if(decomp==triangular){
		/* Perform a triangular decomposition */
		 BasisConstruction<Int> constr;
		 IntMat m_v, m_v2;
		 m_v.resize(dim, dim);
		 m_v2.resize(dim, dim);
		 Int mod=m_lat->getModulo();
	  	 CopyMatr(m_v,m_lat->getBasis(), dim, dim);
        // TriangularizationLower<IntMat,IntVec,Int>(m_v, m_v2 ,mod);
	 	 constr.lowerTriangular(m_v, m_v2 ,mod);	
		// CopyMatr(m_lat->getBasis(), m_v2,dim, dim);
	     for (int i = 0; i < dim; i++){
		   for (int j = 0; j < dim; j++){
			if(i!=j){
			  m_c0(i,j)=NTL::conv<RealRed>(m_v2(i,j))/NTL::conv<RealRed>(m_v2(i,i));		
		     }
		    else{
			   m_c0(i,j)=NTL::conv<RealRed>(m_v2(i,j));
		     }
		    
		   } 
        }
  
        for (int i = 0; i < dim; i++) {
          m_dc2[i] = m_c0(i, i)*m_c0(i, i);
		  }
  
   }
   else{
        std::cerr << "RedBBShortVec:decomp only 'cholesky' and 'triangular' decomposition are supported";
	    return false;
   }
       
    /* Perform the branch and bound.  */
	/* m_n2[j] will be the sum of terms |z*k|^2 ||v*k||^2 for k > j.  */
	m_n2[dim - 1] = 0.0;
	m_countNodes = 0;
	smaller = false;
	m_foundZero = false;
	if (!tryZShortVec(dim - 1, smaller,norm, decomp))   // Here we search for a shortest vector.
		return false;
	if (smaller) {
		// We found a shorter vector. Its square length is in m_lMin2.
		transformStage3ShortVec(m_zShort, k); // Is this useful and OK for L1 ???
		NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
		for (h = 0; h < dim; h++)
			row1(h) = m_bw[h];
		m_lat->setNegativeNorm(k);

		/* The new current shortest vector will be in
		 m_lat->getBasis()(0). */
		/* In the case of L1NORM, check that it is really smaller.  */
		if (norm == L2NORM)
			m_lat->permuteNoDual(k, 0);
		else {
			NTL::matrix_row<IntMat> row5(m_lat->getBasis(), k);
			CalcNorm(row5, dim, x, norm);
			if (x < m_lMin) {
				m_lMin = x;
				m_lMin2 = m_lMin * m_lMin;
				m_lat->permuteNoDual(k, 0);
			}
		}
	}
	return true;
}

//=========================================================================

// This works only for the L2 norm.
template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::reductMinkowski(int d) {
	bool withDual = m_lat->withDual();
	const int dim = m_lat->getDim();
	int i;
	std::int64_t totalNodes = 0;
	bool found;
	bool smaller;          // A smaller vector has been found
	bool xx[dim];      // xx[i]=true means that vector i should not be modified.

	do {
		// The first d vectors should not be modified.
		for (i = 0; i < d; i++)
			xx[i] = true;
		for (i = d; i < dim; i++)
			xx[i] = false;
		found = false;
		do {
			redDieter(d, xx);  // Not so good... should be changed to BKZ!
			m_lat->setNegativeNorm(d);
			m_lat->updateVecNorm(d);
			if (withDual) {
				m_lat->setDualNegativeNorm(d);
				m_lat->updateDualVecNorm(d);
			}
			m_lat->sort(d);
			found = false;

			for (i = 0; i < dim; i++) {
				if (!xx[i]) {
					// On essaie de reduire le i-eme vecteur.
					if (!redBBMink(i, d, 2, smaller, xx))
						return false;
					totalNodes += m_countNodes;
					if (smaller) {
						found = true;
					}
				}
			}
		} while (found);
		// Stage 3
		if (dim > 7) {
			for (i = d; i < dim; i++) {
				if (!redBBMink(i, d, 3, smaller, xx))
					return false;
				totalNodes += m_countNodes;
				if (smaller)
					found = true;
			}
		}
	} while (found);
	m_lat->setNegativeNorm();
	m_lat->updateScalL2Norm(0, dim);
	if (withDual) {
		m_lat->setDualNegativeNorm();
		m_lat->updateDualScalL2Norm(0, dim);
	}
	m_lMin2 = m_lat->getVecNorm(0);
	return true;
}

//=========================================================================

// On successful exit, the basis vectors are sorted by lengths and the norms are updated.
// The square length of the shortest vector can be recovered in m_lMin2.
// The acepted param for 'decomp' are 'cholesky' or 'triangular'. 
// It specify the decompisition basis matrix using in Brach-and-Bound.
template<typename Int, typename Real, typename RealRed>
bool Reducer<Int, Real, RealRed>::shortestVector(NormType norm,std::string decomp) {

	if (norm != L2NORM) {
		m_lat->setNegativeNorm();
	}
	/* Find the shortest vector for the selected norm, L1 or L2.  */
	/* The L2 norm is used for the Cholesky decomposition and BB bounds. */
	bool ok;
	if (norm == L1NORM || norm == L2NORM) {
		ok = redBBShortVec (norm,decomp);
	} else {
		ok = false;
		std::cerr << "RedLattice::shortestVector: this norm is not supported";
		exit(3);
	}
	m_lat->updateVecNorm();
	m_lat->sortNoDual(0);
	return ok;
}

//============================================================================

extern template class Reducer< std::int64_t, double, double> ;
extern template class Reducer< NTL::ZZ, double, double> ;
extern template class Reducer< NTL::ZZ, NTL::RR, NTL::RR> ;

}     // namespace LatticeTester

#endif 