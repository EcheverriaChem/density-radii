Ground-state configurations used for the calculations.

atomic_states.csv      96 neutral atoms (H-Cm). Configuration and spin multiplicity from the
                       NIST Atomic Spectra Database.
atomic_states_ions.csv 95 cations M+ and 81 anions X-. Cation configurations from the NIST
                       database; anion configurations and experimental electron affinities from
                       the compilation of Ning and Lu, J. Phys. Chem. Ref. Data 2022, 51, 021502
                       (values as tabulated on the Wikipedia electron-affinity data page, which
                       reproduces that review). Anions are listed for every element with a
                       positive experimental electron affinity; whether the anion is bound at the
                       PBE0 level is decided by the calculation, not by this file.

Columns
  Z, symbol            atomic number and element symbol
  charge               0, +1 or -1 (ions file only)
  configuration        e.g. [Ar]3d6 4s2; the bracketed noble gas is expanded by the code
  ground_term          spectroscopic term, where recorded
  multiplicity         2S+1 imposed in the calculation
  status               only rows marked APPROVED are used by the input generators
  note                 provenance and any caveat
  EA_eV, EA_quality    anions only: experimental electron affinity and whether it is measured,
                       calculated or estimated

The generators read only rows whose status begins with APPROVED, so a configuration can be put
on hold by editing that column, without deleting the row.
