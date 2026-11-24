Restart inputs used for the bonding configurations of Ce, Pr and Nd (version 1.1).
The first runs (make_javg_alt_inputs.py with data/atomic_states_bond.csv) collapsed to other
configurations; these restarts read their orbitals (%moinp) and move one electron:
  A/  DeltaSCF (maximum overlap method) with ALPHACONF/BETACONF
  B/  90-degree rotation of the two orbitals (%scf Rotate)
Orbital indices count from 0 within each spin, in the order of the gbw/Molden file.
Outcome: B converged to the target configuration for Pr and Nd; Ce (5d2 6s2) fell back to
4f1 5d2 6s1 with both A and B. Adjust the %moinp paths before reuse.
