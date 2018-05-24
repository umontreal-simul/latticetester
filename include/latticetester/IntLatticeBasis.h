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
#include <string>

namespace LatticeTester {

/**
 * \class IntLatticeBasis
 *
 * \brief This class represents a Lattice and its basis
 *
 * This class offers tools to manipulate lattice bases. Each lattice is
 * at least represented by a basis \f$V\f$, a dimension and a Norm.
 * Users can beside precise the dual lattice \f$W\f$ and the modulo.
 * In that case, the flag \f$m_withDual\f$ is set to \f$true\f$.
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
     * Constructor. The primal basis is initialize with \f$primalbasis\f$,
     * the dual basis is initualize with \f$dualbasis\f$, the
     * dimension of the lattice with \f$dim\f$ and the Norm used for
     * reduction with norm.
     */
    IntLatticeBasis (
        const BMat primalbasis,
        const BMat dualbasis,
        const MScal modulo,
        const int dim,
        NormType norm = L2NORM);


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
     * Cleans and releases all the memory allocated for this lattice.
     */
    void kill ();

    /**
     * Copy the lattice
     */
    void copyBasis (const IntLatticeBasis & lat);

    /**
     * Copy the n first elements of the lattice lat
     */
    void copyBasis (const IntLatticeBasis & lat, int n);

    /**
     * Set all the norm of the matrix to -1
     */
    void initVecNorm ();

    /**
     * Return the basis, a BMat-type matrix
     */
    BMat & getBasis () { return m_basis; }

    /**
     * Return the dual basis, a BMat-type Matrix
     */
    BMat & getDualBasis () { return m_dualbasis; }

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
     * Return the norm of the i-th line of the basis
     */
    NScal getDualVecNorm (const int & i) { return m_dualvecNorm[i]; }

    /**
     * Return the norm of each line of the basis
     */
    NVect getDualVecNorm () const { return m_dualvecNorm; }

    /**
     * Return the modulo used for the Dual
     */
    MScal getModulo () const { return m_modulo; }

    /**
     * Set the dimension of the basis
     */
    void setDim (const int & d) { if(d>0) m_dim = d;}

    /**
     * Set the norm used by the lattice
     */
    void setNorm (const NormType & norm) { m_norm = norm; }

    /**
     * Set the norm i egal to the value of NScal
     */
    void setVecNorm ( const NScal & value, const int & i){ m_vecNorm[i] = value; }

    /**
     * Set the dual norm i egal to the value of NScal
     */
    void setDualVecNorm ( const NScal & value, const int & i){ m_dualvecNorm[i] = value; }

    /**
     * Return True if we use Dual.
     */
    bool withDual() { return m_withDual; }

    /**
     * Set the Dual Flag.
     */
    void setDualFlag(bool flag) { m_withDual = flag; }

    /**
     * Get and det m_xx but I don't now what means m_xx
     */
    bool getXX (int i) const { return m_xx[i]; } //???
    void setXX (bool val, int i) { m_xx[i] = val; } //??

     /**
     * Set the norm of all vectors to -1
     */
    void setNegativeNorm ();

    /**
     * Set the norm of the i-eme vector to -1
     */
    void setNegativeNorm (const int & i){ m_vecNorm[i] = -1; }

    /**
     * Set the norm of all vectors of the dual basis to -1
     */
    void setDualNegativeNorm ();

    /**
     * Set the norm of the i-eme vector in the dual Basis to -1
     */
    void setDualNegativeNorm (const int & i){ m_dualvecNorm[i] = -1; }

    /**
     * Recalculates the norm of each vector in the basis of
     * the lattice
     */
    void updateVecNorm ();

    /**
     * Recalculates the norm of each vector in the basis of
     * the lattice from d to dim
     */
    void updateVecNorm (const int & d);

    /**
     * Recalculates the norm of each vector in the basis of
     * the lattice
     */
    void updateDualVecNorm ();

    /**
     * Recalculates the norm of each vector in the dual basis of
     * the lattice from d to dim
     */
    void updateDualVecNorm (const int & d);

    /**
     * Updates the norm of vector at dimension `d` using the `L2NORM`.
     */
    void updateScalL2Norm (const int i);

    /**
     * Updates the norm of all basis vectors from dimensions `d1` to `d2`
     * (inclusive) using the `L2NORM`.
     */
    void updateScalL2Norm (const int k1, const int k2);

    /**
     * Updates the norm of vector in the Dual Basis at dimension `d`
     * using the `L2NORM`.
     */
    void updateDualScalL2Norm (const int i);

    /**
     * Updates the norm of all dual basis vectors from dimensions `d1`
     * to `d2` (inclusive) using the `L2NORM`.
     */
    void updateDualScalL2Norm (const int k1, const int k2);

    /**
     * Exchanges vectors \f$i\f$ and \f$j\f$ in the basis.
     */
    void permute (int i, int j);

    /**
     * Check duality
     */
    bool checkDuality();

    /**
     * Sorts the basis vectors with indices from \f$d\f$ to the dimension
     * of the basis by increasing length. The dual vectors are permuted
     * accordingly. Assumes that the lengths of the corresponding basis
     * vectors are up to date.
     */
    void sort (int d);

    /**
     * Return a string with the primal basis and its norms
     */
    std::string toStringBasis() const;

    /**
     * Return a string with the dual basis and its norms
     */
    std::string toStringDualBasis() const;

    /**
     * Writes the lattice and the parameters on standard output
     */
    void write () const;

protected:

    /**
    * Represent the basis of the lattice. BMat is defined in Types.h
    * This is a matrix where the type of the number can change.
    */
    BMat m_basis;

    /**
    * Represent the dual basis of the lattice. BMat is defined in Types.h
    * This is a matrix where the type of the number can change.
    */
    BMat m_dualbasis;


    /**
    * The dimension of the space.
    */
    int m_dim;

    /**
    * The norm used in the reduction and for this lattice
    */
    NormType m_norm;

    /**
     * The norm of each vector in the basis.
     */
    NVect m_vecNorm;

    /**
     * The norm of each vector in the dual basis.
     */
    NVect m_dualvecNorm;

    /**
     * The modulo linked to the m-dual.
     */
    MScal m_modulo;

    /**
     * If m_withDual is true, we can use Dual Basis.
     */
    bool m_withDual;

    /**
     * This table is used in the Minkowski reduction.
     */
    bool *m_xx;


};

} // namespace LatticeTester

#endif


