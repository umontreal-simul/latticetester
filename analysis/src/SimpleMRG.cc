
//
// Class for simple MRG
//

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

#include "SimpleMRG.h"

using namespace std;
using namespace NTL;


//===========================================================

SimpleMRG::SimpleMRG (const ZZ modulus, const int k, const vec_ZZ a, vec_ZZ state)
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


//===========================================================

void SimpleMRG::goToIndex (int index)
{
	while (m_stateIndex < index) {
		m_state = m_A * m_state;
		m_stateIndex++;
	}
}

//===========================================================

ZZ_p SimpleMRG::getNextValue ()
{
	// update the state
	m_state = m_A * m_state;
	m_stateIndex++;

	// get the value
	return m_state[m_k-1];
}












