# Using the radii

## Which set to use

Use `data/radii/density_pbe0_anorcc.csv` (PBE0, DKH2, ANO-RCC-TZP) unless you have a reason not
to. The other three combinations exist to show how little the radii depend on that choice, not
because any of them is preferable. For ions use `density_cations_pbe0_anorcc.csv` and
`density_anions_pbe0_anorcc.csv`; anion radii are reported with ANO-RCC only, because the
x2c-TZVPPall basis has no diffuse functions and shrinks them by up to 0.09 A.

## File format

One row per element and kind:

| column | meaning |
|---|---|
| `set_id` | identifier of the set |
| `symbol` | element symbol |
| `variant` | spin state or other qualifier, empty when there is only one |
| `kind` | `v` van der Waals-like, `r` valence, `c` core (orbital-based, c_orb), `c_count` sphere enclosing the core electrons (the `c` of version 1.0), `c_elf` outermost ELF minimum inside r_outer (where defined) |
| `value_A` | radius in angstrom |
| `esd_A` | estimated standard deviation, where the source provides one |
| `source` | publication or, for our values, the level of theory |
| `table_ref` | table, figure or criterion the value comes from |
| `note` | caveats: flagged configuration, weakly bound anion, misprint in the source |

Read the `note` column. It is where the exceptions live.

## Domain of validity

**Heavy alkali metals.** K, Rb, Cs and Fr come out 0.67-0.81 A smaller than crystallographic and
isolated-atom radii. The outer density of a group 1 atom comes from a single, spin-isolated s
electron whose Pauli repulsion is stronger than the same density elsewhere, so a single density
threshold underestimates their effective size. Hydrogen behaves the same way. If you need van der
Waals radii for these elements, use an experimental set.

**Weakly bound anions.** 54 of the 81 anions bind the extra electron at the PBE0 level. Whether
one does is predicted almost entirely by its experimental electron affinity: 49 of the 53 anions
above 0.35 eV are bound, only 5 of the 28 below it. Anions whose computed electron affinity
differs from experiment by more than 0.3 eV are flagged in the `note` column; the heavy 6p and 5d
anions (Pb, Po, At, W, Pt, Au) carry errors from the neglect of spin-orbit coupling.

**Free-atom configurations.** A valence radius describes the atom's ground state, not the atom
prepared for bonding. Palladium is the clear case: its free atom is 4d10 5s0, so its outermost
orbital is a compact 4d orbital, while the atom bonds through the empty 5s. Its valence radius
(0.58 A) is therefore not a covalent radius. Recomputed in 4d9 5s1 it becomes 1.38 A, close to
the covalent radius of Cordero et al., and Pr and Nd in 4f(n-1) 5d1 6s2 contract by 0.13 A
(`results/revision/bonding_configurations.csv`). Even so, sums of valence radii do not predict bond
lengths (mean absolute error 0.24 A against 0.07 A for Cordero et al., `scripts/bond_length_test.py`):
the valence radius follows the trends down each group, but it is not a bonding radius.

**Symmetry-broken solutions.** La, Ce and Th (neutral) and Sc+, Ti+, Hf+ have valence orbitals
that are s-d hybrids in a single-determinant description. Their radii are usable - none deviates
from published density radii by more than 0.033 A - and a state-averaged multiconfigurational
calculation shifts them by about +0.04 A in `v` and +0.07 A in `r`, a change of the same size as
the one caused by the functional. For nine ordinary open-shell atoms (B, C, O, F, Al, Cl, Sc, Ti,
Fe) the Hartree-Fock single determinants differ from state-averaged CASSCF by at most 0.011 A in
`v`, 0.055 A in `r` and 0.002 A in `c` (`results/revision/casscf_comparison_open_shell.csv`).
Ce+ is not reported.

**Reference configuration.** Every species is computed in the configuration of the lowest
experimental (NIST ASD) level. For Ni, W and Tb the lowest J-averaged term belongs to another
configuration; the radii in that configuration are in `results/revision/SI_javg_configurations.csv`.

**Which radius for what.** `v` can replace crystallographic van der Waals radii in the
penetration index of covalent bonds (`scripts/penetration_index_test.py`). `r` and `c` describe
the shell structure of the free atom; they should not replace empirical covalent radii.

**Density threshold.** The threshold is a calibration parameter, not a constant of nature. Ours
is 0.001 e bohr-3, following Boyd and Rahm et al.; a calibration that reproduces Bondi's helium
radius gives 0.00072 e bohr-3 instead, which moves the radius of hydrogen by 0.13 A. If you mix
our radii with another density-based set, check that both use the same threshold.

## Precision

Values are tabulated to 0.001 A because that is what the code produces, not because they are
accurate to that. A realistic figure for the spread between reasonable computational choices is
0.01-0.03 A for atoms and cations, and up to 0.05 A for anions; the threshold convention moves
everything by far more. Quoting these radii to two decimals is usually honest enough.
