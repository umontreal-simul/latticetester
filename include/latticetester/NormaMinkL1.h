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

#ifndef LATTICETESTER_NORMAMINKL1_H
#define LATTICETESTER_NORMAMINKL1_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

	/**
	 * This class implements theoretical bounds on the length of the shortest
	 * nonzero vector in a lattice, using the \f${\mathcal{L}}_1\f$ norm.
	 * Here, the length of the shortest nonzero vector is related to the minimal number
	 * of hyperplanes that cover all the points of the dual lattice aasociated.
	 * The following upper bound in this case was provided by Marsaglia
	 * \cite rMAR68a; by applying the general convex body theorem of Minkowski:
	 * \f[
	 *   \ell_t^* = (t!)^{1/t}*(n)^{-1/t}) = \gamma_t^{1/2} \eta^{-1/t},
	 * \f]
	 * for a lattice containing \f$eta\f$ points per unit volume, in dimension \f$t\f$.
	 * The Hermite constants are $\gamma_t = (t!)^{2/t}\f$.
	 */

	class NormaMinkL1 : public Normalizer {
	public:

		/**
		 * Constructs a `NormaMinkL1` for up to `maxDim` dimensions, by assuming that the
		 * log density is `logDensity` in all dimensions.
		 */
		NormaMinkL1 (double logDensity, int maxDim);

		/**
		 * This constructor assumes that the primal lattice has scaling factor \f$m\f$
		 * and order \f$k\f$, so its density is \f$m^k\f$ for \f$t\geq k\f$, and cannot
		 * exceed  \f$m^s\f$ for projections in \f$s < k\f$ dimensions.
		 */
		NormaMinkL1 (double logm, int k, int maxDim);

		/**
		 * Constructs a `NormaMinkL1` for up to `maxDim` dimensions, without computing the bounds.
		 */
		NormaMinkL1 (int maxDim);

		/**
		 * Destructor.
		 */
		~NormaMinkL1();

		/**
		 * Returns the value of the lattice constant \f$\gamma_j\f$ in
		 * dimension \f$j\f$.
		 */
		double getGamma (int j) const;

	private:

		/**
		 * The lattice constants \f$\gamma_j\f$ for the Minkowski-Marsaglia bounds in each
		 * dimension \f$j\f$.
		 */
		double *m_gamma;

		/**
		 * Precomputed lattice constants \f$\gamma_j\f$ for the Minkowski bounds
		 * in each dimension \f$j \le 48\f$.
		 */
		static const double m_gamma0[1 + Normalizer::MAX_DIM];

		/**
		 * Computes the MinkL1 bound in dimension \f$d\f$.
		 */
		double calcGamma (int d);
	}; // End class NormaMinkL1

	//===========================================================================

	const double NormaMinkL1::m_gamma0[] =
	{         //  GamMinkL1[t] = (t!)^{2/t}
		/* GamMinkL1[0] = */0.00000000000000,
		/* GamMinkL1[1] = */1.00000000000000,
		/* GamMinkL1[2] = */2.00000000000000,
		/* GamMinkL1[3] = */3.301927248894626,
		/* GamMinkL1[4] = */4.898979485566356,
		/* GamMinkL1[5] = */6.786916380543178,
		/* GamMinkL1[6] = */8.962809493114328,
		/* GamMinkL1[7] = */11.42450247602496,
		/* GamMinkL1[8] = */14.17033543597957,
		/* GamMinkL1[9] = */17.19898810749518,
		/* GamMinkL1[10] = */20.50938353057181,
		/* GamMinkL1[11] = */24.10062539497529,
		/* GamMinkL1[12] = */27.97195541552144,
		/* GamMinkL1[13] = */32.12272326936208,
		/* GamMinkL1[14] = */36.55236475376784,
		/* GamMinkL1[15] = */41.2603855160832,
		/* GamMinkL1[16] = */46.24634867434057,
		/* GamMinkL1[17] = */51.50986522412917,
		/* GamMinkL1[18] = */57.05058648500341,
		/* GamMinkL1[19] = */62.868198068697,
		/* GamMinkL1[20] = */68.96241500217117,
		/* GamMinkL1[21] = */75.33297774027004,
		/* GamMinkL1[22] = */81.97964887293777,
		/* GamMinkL1[23] = */88.90221038131129,
		/* GamMinkL1[24] = */96.10046133233949,
		/* GamMinkL1[25] = */103.5742159272685,
		/* GamMinkL1[26] = */111.3233018382898,
		/* GamMinkL1[27] = */119.3475587818154,
		/* GamMinkL1[28] = */127.64683728756,
		/* GamMinkL1[29] = */136.2209976308053,
		/* GamMinkL1[30] = */145.069908901562,
		/* GamMinkL1[31] = */154.1934481892748,
		/* GamMinkL1[32] = */163.5914998656133,
		/* GamMinkL1[33] = */173.2639549509683,
		/* GamMinkL1[34] = */183.2107105527401,
		/* GamMinkL1[35] = */193.4316693654913,
		/* GamMinkL1[36] = */203.9267392246444,
		/* GamMinkL1[37] = */214.6958327067107,
		/* GamMinkL1[38] = */225.7388667701241,
		/* GamMinkL1[39] = */237.0557624316304,
		/* GamMinkL1[40] = */248.6464444739227,
		/* GamMinkL1[41] = */260.5108411808306,
		/* GamMinkL1[42] = */272.6488840968785,
		/* GamMinkL1[43] = */285.0605078084658,
		/* GamMinkL1[44] = */297.7456497442783,
		/* GamMinkL1[45] = */310.7042499928673,
		/* GamMinkL1[46] = */323.9362511355723,
		/* GamMinkL1[47] = */337.44159809321,
		/* GamMinkL1[48] = */351.220237985132,
	};

	/*=========================================================================*/

	double NormaMinkL1::calcGamma (int dim) {
		double gamma = 0.0;
		for (int i = 1; i <= dim; i++)
		gamma += log(i);
		gamma *= 2.0 / dim;
		return exp(gamma);
	}

	/*=========================================================================*/

	NormaMinkL1::NormaMinkL1 (double logDensity, int maxDim)
	: Normalizer ( maxDim, "MinkL1", L1NORM) {
		m_gamma = new double[maxDim + 1];
		int t0 = maxDim;
		if (t0 > this->MAX_DIM)
		t0 = this->MAX_DIM;
		for (int i = 0; i <= t0; i++)
		m_gamma[i] = m_gamma0[i];
		for (int i = t0 + 1; i <= maxDim; i++)
		m_gamma[i] = calcGamma(i);
		Normalizer::computeBounds (logDensity);
	}

	/*=========================================================================*/

	NormaMinkL1::NormaMinkL1 (double logm, int k, int maxDim)
	: Normalizer (maxDim, "MinkL1", L1NORM) {
		m_gamma = new double[maxDim + 1];
		int t0 = maxDim;
		if (t0 > this->MAX_DIM)
		t0 = this->MAX_DIM;
		for (int i = 0; i <= t0; i++)
		m_gamma[i] = m_gamma0[i];
		for (int i = t0 + 1; i <= maxDim; i++)
		m_gamma[i] = calcGamma(i);
		Normalizer::computeBounds (logm, k);
	}

	/*=========================================================================*/

	NormaMinkL1::~NormaMinkL1() {
		delete[] m_gamma;
	}

	/*=========================================================================*/

	inline double NormaMinkL1::getGamma (int j) const {
		if (j < 1 || j > this->m_maxDim)
		throw std::out_of_range("NormaMinkL1::getGamma");
		return m_gamma[j];
	}

} // End namespace LatticeTester

#endif

