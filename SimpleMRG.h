
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
	SimpleMRG (const ZZ modulus, const int k, const vec_ZZ a, vec_ZZ state);

	// Destructor
	~SimpleMRG ();

	// Update m_state from current m_index to new index;
	void goToIndex (int index);

	// Give the last value calculated
	ZZ_p getLastValue ();

	// Update m_state to the next step and give the new value calculated
	ZZ_p getNextValue ();

	// access function
	ZZ getModulus ();

	// access function
	int getOrder ();

	// access function
	mat_ZZ_p getCompanionMatrix ();

	// access function
	int getStateIndex ();

	// access function
	vec_ZZ_p getState ();




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



















