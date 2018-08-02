# This file contains a few notes on what work has to be done in this repo

## Urgent work
- Have examples (usage of most functions on a simple basis)
- Update the doc to be publishable
- Choose a license

## New features
- Test current implementation for supported types combinations
- Make sure everything that is advertised in the manual is implemented
- Implement a test suite to look for bugs

## Design changes

## Minor changes
- Add the `.bib` files as a submodule to keep them up to date.

# Questions
- Quelle est la pertinence d'avoir les normalisations dans LatticeTester sachant
que
  - LatNetBuilder n'utilise pas les normalisations de LatticeTester
  - Aucun test n'est implémenté directement dans LatticeTester?
- Cela soulève également la question: est-ce que le test spectral devrait plutôt
être implémenté dans LatticeTester plutôt que dans LatMRG?
  - Si oui, alors LatticeTester devient un vrai petit programme en soit
  - Si non, LatticeTester est vraiment simplement une librairie avec des
  algorithmes de réduction et de manipulation de lattices.
- D'où viennent les algorithmes pour calculer le dual? Ou bien les bases avec
des indices manquants?
