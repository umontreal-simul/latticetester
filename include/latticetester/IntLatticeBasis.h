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
     * Return the dimension of the lattice, a long type
     */
    long getDim () const { return m_dim; }

    /**
     * Return the Norm type used by the lattice
     */
    NormType getNorm () const { return m_norm; }

    /**
     * Return the norm of each line of lattice
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
     * Recalculates the norm of each vector in the basis of
     * the lattice from d to m_dim
     */
    void updateVecNorm (const int & d);

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
    long m_dim;

    /**
     * The norm of each vector in the basis.
     */
    NVect m_vecNorm;


};

} // namespace LatticeTester

#endif


