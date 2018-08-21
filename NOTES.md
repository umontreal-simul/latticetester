# This file contains a few notes on what work has to be done in this repo

## Urgent work
- Have examples (usage of most functions on a simple basis)
  - One on lacunary indices
  - One on basis/dual basis construction
- Update the doc to be publishable
  - Finish documentation on Weight classes, Lacunary and Rank1Lattice
  - Re-read what has been done to fit with the pdf guide.
- Choose a license
- Implement upper bounds as in cohn and elkies

## New features
- Implement the P_alpha calculation
- Test current implementation for supported types combinations
- Make sure everything that is advertised in the manual is implemented
- Implement a test suite to look for bugs

## Design changes

## Minor changes

# Questions
- D'où viennent les algorithmes pour calculer le dual? Ou bien les bases avec
des indices manquants?
- Il faut que je m'inscrive
- Le calcul des bornes de Rogers ne vient pas de l'article de 1959 parce que le
terme avec les log donne quelque chose de négatif pour les valeurs que l'on
considère. C'est vraiment de l'article de 58, mais il faut arriver à calculer le truc...

# Notes personnelles
- Est-ce que LatticeTester devrait avoir des classes pour supporter les réseaux
qui ne sont pas générés par une bases de n vecteur en n dimensions (présentement
le logiciel supporte seulement les matrices de base carrées).
  - R: Probablement
