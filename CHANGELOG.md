# Changelog

## 1.1.0 (2026-11)

Revision of the article (J. Phys. Chem. A, transferred from J. Chem. Theory Comput.).

### Changed
- **Core radius.** `c` in `data/radii/density_*.csv` is now the orbital-based core radius c_orb,
  the most probable radius (outermost maximum of r^2|phi|^2) of the outermost core orbital. The
  sphere enclosing the chemical-core electrons, which was `c` in 1.0, is kept as `c_count`, and the
  outermost ELF minimum inside r_outer is added as `c_elf` where defined.
- Cations with a bare noble-gas configuration (Rb+, Cs+, Fr+) take the core of the isoelectronic
  noble gas.
- The SLURM template loads `orca/6.1.1` and gives every array task a private TMPDIR, which avoids
  ORCA/OpenMPI jobs killing each other when they share a node.
- `make_casscf_tests.py` includes nine open-shell atoms (B, C, O, F, Al, Cl, Sc, Ti, Fe).

### Added
- Reference-configuration analysis: `scripts/fetch_nist_levels.py`, `scripts/javg_nist.py`,
  `workflow/make_javg_alt_inputs.py`, `workflow/fix_javg.py`, `scripts/audit_javg.py`,
  `data/atomic_states_javg.csv` and the NIST ASD levels (`data/reference_data/nist_asd_levels.zip`).
- Core radii: `scripts/core_radii_all.py`, `scripts/core_radii_elf.py`,
  `scripts/update_core_radii_tables.py`, `results/core_radii_*.csv`.
- Bonding configurations (Pd, Pr, Nd): `data/atomic_states_bond.csv`, `scripts/audit_bond.py`,
  `workflow/examples/bond_fixes`.
- State-averaged CASSCF comparison: `scripts/analyse_casscf.py`.
- Tests of the radii: `scripts/compare_desclaux.py` (Dirac-Fock orbital radii),
  `scripts/bond_length_test.py` (403 single bonds of Pyykko and Atsumi),
  `scripts/penetration_index_test.py` (14 sets of Pardos et al.), `scripts/plot_crust.py`.
- `results/revision/` with the outputs of all the above, and `data/reference_data/` with their
  transcribed inputs.

### Unchanged
- The van der Waals-like radius `v` and the valence radius `r` of every species.

## 1.0.0 (2026-09)

First release, with the submitted manuscript.
