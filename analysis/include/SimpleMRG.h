
//
// Class for simple MRG
//


#ifndef SIMPLEMRG_H
#define SIMPLEMRG_H

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>


#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

using namespace std;
using namespace NTL;


class SimpleMRG {

public:

	// Constructor. Create the matrix m_A with the given coefficients (a_i)
	// and initialize the current m_state to the input state.
	SimpleMRG (const ZZ modulus, const int k, const vec_ZZ a, vec_ZZ state)
	{
	m_modulus = modulus;
	ZZ_p::init(modulus);
	m_k = k;
	m_A.SetDims(k, k);
	// update last lign of the companion matrix
	for (int j = 0; j < k; j++)
		m_A[k-1][j] = conv<ZZ_p>( a[k-1-j] );
	// update upper diagonal of the companion matrix
	for (int i = 0; i < k-1; i++)
		m_A[i][i+1] = 1;

	m_stateIndex = 0;
	m_state.SetLength(k);
	m_state = conv<vec_ZZ_p> (state);
	}

	// Destructor
	~SimpleMRG () {}

	// Update m_state from current m_index to new index;
	void goToIndex (int index)
	{
		while (m_stateIndex < index) {
			m_state = m_A * m_state;
			m_stateIndex++;
		}
	}

	// Give the last value calculated
	ZZ_p getLastValue () const { return m_state[m_k-1]; }

	// Update m_state to the next step and give the new value calculated
	ZZ_p getNextValue ()
	{
		// update the state
		m_state = m_A * m_state;
		m_stateIndex++;

		// get the value
		return m_state[m_k-1];
	}

	// access function
	ZZ getModulus () const { return m_modulus; }

	// access function
	int getOrder () const { return m_k; }

	// access function
	mat_ZZ_p getCompanionMatrix () const { return m_A; }

	// access function
	int getStateIndex () const { return m_stateIndex; }

	// access function
	vec_ZZ_p getState () const { return m_state; }




private:

	// modulus of the MRG
	ZZ m_modulus;

	// order of the MRG
	int m_k;

	// representation matrix A with (a_i) coefficients
	mat_ZZ_p m_A;

	// index current state
	int m_stateIndex;

	// Values of current state
	vec_ZZ_p m_state;

};

#endif



















