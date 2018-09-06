# This file contains a few notes on what work has to be done in this repo

## Urgent work
- Have examples (usage of most functions on a simple basis)
  - One on basis/dual basis construction
  - One on Weights???
- Update the doc to be publishable
  - Finish documentation on Weight classes, Lacunary and Rank1Lattice/IntLattice
  - Re-read what has been done to fit with the pdf guide.
- Implement upper bounds as in Cohn and Elkies and look for higher dimensions
- Build projections basis and dual basis
- Benchmark our LLL implementation
- Check if LLL can be done with all norms.

## New features
- Choose a license
- Implement the P_alpha calculation
- Test current implementation for supported types combinations
- Make sure everything that is advertised in the manual is implemented
- Implement a test suite to look for bugs

## Design changes

## Minor changes

# Questions
J'ai pensé à ça, et c'est possible de retirer les switch. Il suffit en fait de
réimplémenter les enum en des classes et créer des sous-classes pour chaque cas
qui elles-même implémentent les fonctions nécessaires, ce qui nous permet de ne
pas avoir de switch et qui fait en sorte que tout ce que les switch font se
décide à la compilation!

# Notes personnelles
- Est-ce que LatticeTester devrait avoir des classes pour supporter les réseaux
qui ne sont pas générés par une bases de n vecteur en n dimensions (présentement
le logiciel supporte seulement les matrices de base carrées).
  - R: Probablement
- Est-ce que Lacunary est pertinent dans LatticeTester? Je ne pense pas parce
que c'est principalement un truc qui est utile pour étudier les MRGs.
  - R: Oui. On peut faire la projection d'un réseau sur certaines de ses
    coordonnées de manière générale en prenant simplement un sous-ensemble des
    colonnes pour constituer l'ensemble générateur.
- Le calcul des bornes de Rogers ne vient pas de l'article de 1959 parce que le
terme avec les log donne quelque chose de négatif pour les valeurs que l'on
considère. C'est vraiment de l'article de 58, mais il faut arriver à calculer le truc...
