// Defining constants for the execution of the algorithms
BasisConstruction<BScal> constr; // The basis constructor we will use
BMat bas_mat, dua_mat; // basis matrix and dual basis matrix
// Creating a lattice basis
IntLatticeBasis<MScal, BScal, NScal, RScal> lattice(bas_mat, numlines);
// Constructing a basis matrix with GCDConstruction
constr.GCDConstruction(bas_mat);
// Constructing the matrix of a dual to bas_mat
MScal modulo(1);
constr.DualConstruction(bas_mat, dua_mat, modulo);
// Constructing a basis for lattice with LLLConstruction
constr.LLLConstruction(lattice.getBasis());
modulo = MScal(1);
constr.DualConstruction(lattice.getBasis(), lattice.getDualBasis(), modulo);
// The preceding line works to compute a dual for lattice but does not set all
// the properties of lattice to properly work with a dual.
// This next line sets the lattice to know it has a dual. Computing the norm of
// the vectors in the lattice would also be wise.
lattice.setDualFlag(true);
