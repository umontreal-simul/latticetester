/**
 * This example showcases the usage of the BasisConstruction module. This reads
 * matrices from files and builds a basis and a dual for IntLatticeBasis objects.
 *
 * This should also compare the speed of methods, namely, the speed of LLLConstruction + GCDConstruction vs the only usage of GCDConstruction.
 * */
/**
 * For now, this example benchmarks the different implementations we have
 * available for basis reduction.
 *
 * Add an error counter to see if both methods work.
 *
 * As of now the results look like this:
 * ma_basis: 21739192
 * cou_basis: 249924718
 * ma_basis/cou_basis: 0.086983
 * Dimension: 5
 * ma_basis: 300
 * cou_basis: 764
 * ma_basis/cou_basis: 0.39267
 * Dimension: 10
 * ma_basis: 1892
 * cou_basis: 3883
 * ma_basis/cou_basis: 0.487252
 * Dimension: 15
 * ma_basis: 7809
 * cou_basis: 18444
 * ma_basis/cou_basis: 0.42339
 * Dimension: 20
 * ma_basis: 48107
 * cou_basis: 186311
 * ma_basis/cou_basis: 0.258208
 * Dimension: 25
 * ma_basis: 70023
 * cou_basis: 329619
 * ma_basis/cou_basis: 0.212436
 * Dimension: 30
 * ma_basis: 449102
 * cou_basis: 6603505
 * ma_basis/cou_basis: 0.0680096
 * Dimension: 35
 * ma_basis: 21162057
 * cou_basis: 242782302
 * ma_basis/cou_basis: 0.0871647
 * ma_dual: 179732000
 * cou_dual: 5370
 * ma_dual/cou_dual: 33469.6
 *
 * This means that 1) Obviously the dual basis construction from solving a
 * linear system is crap. Doing the Euclid algorithm on the lines of the matrix
 * is about 10 times faster for the instances studied here but we should verify
 * if there is a difference in this ratio depending on the size of the instance.
 * Also, it is important to note that basis construction takes way less memory
 * in my method because it is done in place.
 *
 * Although it is not definitive, the triangularization seems faster in my
 * implementation, especially in larger dimensions.
 *
 * Aussi, clairement la construction de base par LLL est juste beaucoup trop
 * mieux que la triangularisation justement parce que les vecteurs n'escaladent
 * pas vers des tailles démoniaques.
 * */

// This should always use Types 2 or 3, because we get too big numbers with GCD
// elimination.
#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"

using namespace LatticeTester;

/*
 * Ce qu'on veut c'est de sonder les fichiers, puis
 * construire la base puis le dual pour chacun et le mettre dans un IntLatticeBasis.
 * On veut comparer la vitesse pour la construction de la base pour plusieurs méthodes différentes.
 * pour comparer les exécutions il suffit de mesurer le temps pris par chaque tour de boucle.
 * */
int main() {
  // The different clocks we will use for benchmarking
  int max_dim = 7; // Actual max dim is 5*max_dim
  clock_t tmp;
  clock_t gcd_time[max_dim], lll_time[max_dim];
  for (int i = 0; i < max_dim; i++) {
    gcd_time[i] = 0;
    lll_time[i] = 0;
  }

  std::string prime = "1021";
  ParamReader<MScal, BScal, RScal> reader;

  // The constructor we will use
  BasisConstruction<BScal> constr;
  BMat bas_mat, dua_mat;
  int numlines;
  unsigned int ln;
  std::string name;
  // This loop builds the basis with GCDConstruction straight up.
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      std::cout << name << std::endl;
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      bas_mat.kill();
      bas_mat.SetDims(numlines, numlines);
      dua_mat.kill();
      dua_mat.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(bas_mat, ln, 0, numlines);
      // We want to avoid singular matrix because we can't compute the dual, and
      // IntLatticeBasis only really supports square matrices.
      if (NTL::determinant(bas_mat) == 0) {
        std::cout << name << " is singular\n";
        continue;
      }

      // Timing ma first
      constr.GCDConstruction(bas_mat);
      // If you don't need the dual basis, the following line is sufficient
      // basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, numlines);
      MScal modulo(1);

      constr.DualConstruction(bas_mat, dua_mat, modulo);
      IntLatticeBasis <MScal, BScal, NScal, RScal> basis =
        IntLatticeBasis<MScal, BScal, NScal, RScal>(bas_mat, dua_mat, modulo, numlines);
      gcd_time[j] += clock() - tmp;
    }
  }
  // This loop builds the basis with LLLConstruction first.
  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      tmp = clock();
      // Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      reader = ParamReader<MScal, BScal, RScal>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      bas_mat.kill();
      bas_mat.SetDims(numlines, numlines);
      dua_mat.kill();
      dua_mat.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(bas_mat, ln, 0, numlines);
      // We want to avoid singular matrix because we can't compute the dual, and
      // IntLatticeBasis only really supports square matrices.
      if (NTL::determinant(bas_mat) == 0) {
        std::cout << name << " is singular\n";
        continue;
      }

      // Timing ma first
      constr.GCDConstruction(bas_mat);
      // If you don't need the dual basis, the following line is sufficient
      // basis = IntLatticeBasis<MScal, BScal, NScal, RScal>(matrix, numlines);
      MScal modulo(1);

      constr.DualConstruction(bas_mat, dua_mat, modulo);
      IntLatticeBasis <MScal, BScal, NScal, RScal> basis =
        IntLatticeBasis<MScal, BScal, NScal, RScal>(bas_mat, dua_mat, modulo, numlines);
      lll_time[j] += clock() - tmp;
    }
  }
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += gcd_time[i];
  }
  int width1 = log10(tmp)+2;
  std::cout << "           GCD";
  for (int i = 0; i<width1-3; i++) {
    std::cout << " ";
  }
  std::cout << "LLL\n";
  std::cout << "Total time " << tmp << " ";
  tmp = 0;
  for (int i = 0; i < max_dim; i++) {
    tmp += lll_time[i];
  }
  int width2 = log10(tmp)+2;
  std::cout << tmp << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 << std::setw(width1)
      << gcd_time[i] << std::setw(width2) << lll_time[i] << std::endl;
  }

  return 0;
}
