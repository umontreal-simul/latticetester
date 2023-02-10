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

#ifndef LATTICETESTER_NORMABESTBOUND_H
#define LATTICETESTER_NORMABESTBOUND_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

	/**
	 * In this normalizer, the Hermite constants \f$\gamma_s\f$ are approximated using the
	 * best upper bounds that are available.
	 * For dimensions 0 through 8, the Hermite constant are known exactly.
	 * For dimensions 9 through 36, they are computed using the article
	 * "A Conceptual breakthrough in sphere packing" from Henry Cohn (2016).
	 * see https://doi.org/10.4007/annals.2003.157.689.
	 * For dimensions 37 through 48 we use Rogers's bound.
	 * This class is to be used with the L2NORM (the Euclidean norm) exclusively.
	 */
	class NormaBestBound : public Normalizer {
	public:

		/**
		 * Constructs a `NormaBestBound` for up to `maxDim` dimensions, by assuming that the
		 * log density is `logDensity` in all dimensions.
		 * Restriction: `maxDim`\f$ \le 48\f$.
		 */
		NormaBestBound (double logDensity, int maxDim);

    	/**
    	 * This constructor assumes that the primal lattice has scaling factor \f$m\f$
    	 * and order \f$k\f$, so its density is \f$m^k\f$ for \f$t\geq k\f$, and cannot
    	 * exceed  \f$m^s\f$ for projections in \f$s < k\f$ dimensions.
    	 */
    	NormaBestBound (double logm, int k, int maxDim);

    	/**
		 * Destructor.
		 */
		~NormaBestBound();

		/**
		 * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
		 * in dimension \f$j\f$.
		 */
		double getGamma (int j) const;

	private:

		/**
		 * Precomputed lattice constants \f$\gamma_j\f$ for the best bound
		 * in each dimension \f$j \le 48\f$.
		 */
		static const double m_gamma[1 + Normalizer::MAX_DIM];

	}; // End class NormaBestBound

	//===========================================================================

	/*
	 * Values 1 through 36 are calculated from the article "A Conceptual
	 * Breakthrough in sphere packing" from Henry Cohn (2016). They are computed
	 * by taking gamma[n] = 4 * ( bound / V_n )^(n/2) where V_n is the volume of
	 * an n dimensional sphere of radius 1.
	 */
//<<<<<<< HEAD
	//template<typename RealRed>
//=======
//>>>>>>> b23681ea0112bce9ba98c1463251528b775075a4
	const double NormaBestBound::m_gamma[ ] =
	{
		/* GammaBestBound[0] = */0.0,
		/* GammaBestBound[1] = */1.0,
		/* GammaBestBound[2] = */1.1547005395033834,
		/* GammaBestBound[3] = */1.304077209522991,
		/* GammaBestBound[4] = */1.4491512970689033,
		/* GammaBestBound[5] = */1.5907345457787878,
		/* GammaBestBound[6] = */1.7294425778558589,
		/* GammaBestBound[7] = */1.8657438977611798,
		/* GammaBestBound[8] = */2.0000000001950413,
		/* GammaBestBound[9] = */2.1324942979385826,
		/* GammaBestBound[10] = */2.263452600069123,
		/* GammaBestBound[11] = */2.3930576417841327,
		/* GammaBestBound[12] = */2.5214594724429005,
		/* GammaBestBound[13] = */2.6487829527745377,
		/* GammaBestBound[14] = */2.7751332142995158,
		/* GammaBestBound[15] = */2.9005997190604442,
		/* GammaBestBound[16] = */3.025259312863684,
		/* GammaBestBound[17] = */3.149178557041412,
		/* GammaBestBound[18] = */3.2724155399984154,
		/* GammaBestBound[19] = */3.395021298677218,
		/* GammaBestBound[20] = */3.517040950462465,
		/* GammaBestBound[21] = */3.6385145917472883,
		/* GammaBestBound[22] = */3.7594780693416627,
		/* GammaBestBound[23] = */3.8799635383208293,
		/* GammaBestBound[24] = */4.000000015650443,
		/* GammaBestBound[25] = */4.119613703661868,
		/* GammaBestBound[26] = */4.238828491527929,
		/* GammaBestBound[27] = */4.357666099112316,
		/* GammaBestBound[28] = */4.476146404322391,
		/* GammaBestBound[29] = */4.594287637369152,
		/* GammaBestBound[30] = */4.7121065670603075,
		/* GammaBestBound[31] = */4.8296186570429365,
		/* GammaBestBound[32] = */4.94683820125262,
		/* GammaBestBound[33] = */5.063778445850659,
		/* GammaBestBound[34] = */5.180451691100012,
		/* GammaBestBound[35] = */5.296869382898868,
		/* GammaBestBound[36] = */5.413042186941963,
		/* GammaBestBound[37] = */5.7020718581143,
		/* GammaBestBound[38] = */5.8255656070255,
		/* GammaBestBound[39] = */5.9489276473284,
		/* GammaBestBound[40] = */6.0721635670068,
		/* GammaBestBound[41] = */6.1952785955803,
		/* GammaBestBound[42] = */6.3182776348,
		/* GammaBestBound[43] = */6.4411652860615,
		/* GammaBestBound[44] = */6.5639458749555,
		/* GammaBestBound[45] = */6.6866234733141,
		/* GammaBestBound[46] = */6.8092019190592,
		/* GammaBestBound[47] = */6.9316848341156,
		/* GammaBestBound[48] = */7.0540756406128
	};

	/*=======================================================================*/

	NormaBestBound::NormaBestBound (double logDensity, int maxDim)
	: Normalizer (maxDim, "Best", L2NORM) { //Normalizer (logDensity, maxDim, "Best", L2NORM
		Normalizer::computeBounds (logDensity);
	}

    /*=========================================================================*/

    NormaBestBound::NormaBestBound (double logm, int k, int maxDim)
      : Normalizer (maxDim, "BestLat", L2NORM)
      {
        if (maxDim > this->MAX_DIM)
          throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
        Normalizer::computeBounds (logm, k);
      }

    /*=========================================================================*/

	NormaBestBound::~NormaBestBound() {
	}

	/*=========================================================================*/

	inline double NormaBestBound::getGamma (int j) const
	{
		if (j < 1 || j > this->m_maxDim)
		    throw std::out_of_range("NormaBestBound::getGamma");
		return m_gamma[j];
	}

}

#endif
