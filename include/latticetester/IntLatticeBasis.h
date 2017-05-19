// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
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

#ifndef LATTICETESTER__INTLATTICEBASIS_H
#define LATTICETESTER__INTLATTICEBASIS_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"
//#include "latticetester/Normalizer.h"
//#include "latticetester/CoordinateSets.h"
//#include "latticetester/Lacunary.h"
#include <string>

namespace LatticeTester {

/**
 * This class offers tools to manipulate lattice bases. Each lattice is
 * represented by a basis \f$V\f$ and a Norm.
 */

/* Erwan
* Cette classe représente un réseau. Une lattice est caractérisée par
* une base et une norme. Cette classe est née de la fusion de IntLattice
* et Basis. Elle permet une meilleure lisibilitée d'un réseau.
*/

class IntLatticeBasis {
public:

    /**
     * Constructor. The primal basis is initialize with identity,
     * the dimension of the lattice with dim and the Norm used for
     * reduction with norm.
     */
    IntLatticeBasis (const int dim, NormType norm = L2NORM);

    /**
     * Constructor. The primal basis is initialize with \f$basis\f$,
     * the dimension of the lattice with dim and the Norm used for
     * reduction with norm.
     */
    IntLatticeBasis (const BMat basis, const int dim, NormType norm = L2NORM);

    /**
     * Copy constructor. The maximal dimension of the created basis is set
     * equal to <tt>Lat</tt>’s current dimension.
     */
    IntLatticeBasis (const IntLatticeBasis & Lat);

    /**
     * Destructor
     */
    ~IntLatticeBasis ();

    /**
     * Copy the lattice
     */
    void copyBasis (const IntLatticeBasis & lat);

    /**
     * Return the basis, a BMat-type matrix
     */
    BMat & getBasis () { return m_basis; }

    /**
     * Return the dimension of the lattice, a int type
     */
    int getDim () const { return m_dim; }

    /**
     * Return the Norm type used by the lattice
     */
    NormType getNorm () const { return m_norm; }

    /**
     * Return the norm of the i-th line of the basis
     */
    NScal getVecNorm (const int & i) { return m_vecNorm[i]; }

    /**
     * Return the norm of each line of the basis
     */
    NVect getVecNorm () const { return m_vecNorm; }

    /**
     * Set the dimension of the basis
     */
    void setDim (const int & d) { if(d>0) m_dim = d;}

    /**
     * Set the norm used by the lattice
     */
    void setNorm (const NormType & norm) { m_norm = norm; }

    /**
     * Set the basis of the lattice
     */
    void setBasis (const BMat & basis) { m_basis = basis; }

    /**
     * Set all the norm of the matrix to -1
     */
    void initVecNorm ();

    /**
     * Set all the norm egal to the value of NScal
     */
    void setVecNorm ( const NVect & vecnorm ){ m_vecNorm = vecnorm; }

    /**
     * Set the norm i egal to the value of NScal
     */
    void setVecNorm ( const NScal & value, const int & i){ m_vecNorm[i] = value; }

    /**
     * Recalculates the norm of each vector in the basis of
     * the lattice
     */
    void updateVecNorm ();

    /**
     * Get and det m_xx but I don't now what means m_xx
     */
    bool getXX (int i) const { return m_xx[i]; } //???
    void setXX (bool val, int i) { m_xx[i] = val; } //??

    /**
     * Recalculates the norm of each vector in the basis of
     * the lattice from d to dim
     */
    void updateVecNorm (const int & d);

    /**
     * Updates the norm of vector at dimension `d` using the `L2NORM`.
     */
    void updateScalL2Norm (const int i);

    /**
     * Set the norm of all vectors to -1
     */
    void setNegativeNorm ();

    /**
     * Set the norm of the i-eme vector to -1
     */
    void setNegativeNorm (const int & i){ m_vecNorm[i] = -1; }

    /**
     * Updates the norm of all basis vectors from dimensions `d1` to `d2`
     * (inclusive) using the `L2NORM`.
     */
    void updateScalL2Norm (const int k1, const int k2);

    /**
     * Exchanges vectors \f$i\f$ and \f$j\f$ in the basis.
     */
    void permute (int i, int j);

    /**
     * Sorts the basis vectors with indices from \f$d\f$ to the dimension
     * of the basis by increasing length. The dual vectors are permuted
     * accordingly. Assumes that the lengths of the corresponding basis
     * vectors are up to date.
     */
    void sort (int d);

    /**
     * Writes the lattice and the parameters on standard output
     */
    void write () const;

private:

    /**
    * Represent the basis of the lattice. BMat is defined in Types.h
    * This is a matrix where the type of the number can change.
    */
    BMat m_basis;

    /**
    * The norm used in the reduction and for this lattice
    */
    NormType m_norm;

    /**
    * The dimension of the space.
    */
    int m_dim;

    /**
     * The norm of each vector in the basis.
     */
    NVect m_vecNorm;

    /**
     * Implemented in LatCommun but I've no idea what it does
     */
    bool *m_xx;


};

} // namespace LatticeTester

#endif


