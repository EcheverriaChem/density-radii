# Data files

## data/radii/*.csv

Radius tables, ours and the transcribed published sets, in the format described in
`docs/radii_usage.md`. Load them with `density_radii.load_set`, which resolves variants (spin
states, coordination numbers) and raises rather than guessing when a choice is ambiguous.

Transcribed sets and the check applied to each:

| file | source | transcription check |
|---|---|---|
| `bondi1964_vdw.csv` | Bondi 1964, Table I | values compared with the Ar correction published by Vogt and Alvarez |
| `mantina2009_vdw.csv` | Mantina et al. 2009 | - |
| `batsanov2001_vdw_crystallographic.csv`, `batsanov2001_vdw_equilibrium.csv` | Batsanov 2001, Table 9 | all 65 crystallographic values match the column quoted by Alvarez 2013 |
| `alvarez2013_vdw.csv` | Alvarez 2013, Table 1 | crust widths reproduce those quoted in Echeverria and Alvarez 2024 |
| `vogt2014_noble_gas_vdw.csv` | Vogt and Alvarez 2014 + corrigendum | Ar value from the corrigendum |
| `cordero2008_covalent.csv` | Cordero et al. 2008 | - |
| `pyykko2009_covalent.csv` | Pyykko and Atsumi 2009, Figure 2 | three bond-length sums quoted in their text |
| `rahm2017_density.csv` | Rahm et al. 2016 with the 2017 corrigendum | published statistics against Alvarez recovered |
| `shannon1976_ionic.csv` | Shannon 1976, Table 1, CN VI | crystal radius minus ionic radius = 0.14 A for every entry |
| `charry2024_vdw.csv` | Charry and Tkatchenko 2024, Figure 7 | every value reproduces their own radius-polarizability relation |

## results/*.csv

| file | content |
|---|---|
| `atomic_radii_raw_all.csv` | every neutral calculation, before corrections |
| `atomic_radii_selected.csv` | the calculation chosen for each element, level and quantity, with its provenance |
| `ion_radii_selected.csv` | the same for the ions, with ionization energies and electron affinities |
| `reference_set_comparison.csv` | regression statistics against every published set |
| `ion_sensitivity.csv` | effect of method, Hamiltonian and basis on the ion radii |
| `casscf_comparison.csv` | state-averaged CASSCF test for La, Ce and Th |

Columns in the per-calculation tables: the job name encodes `Z_symbol[charge]_method_scheme`;
`config_ok` and `max_pop_dev` report the configuration audit; `origin` names the calculation the
value came from; `usable` and `flag` carry the decisions described in `docs/radii_usage.md`.

## results/revision/*.csv (version 1.1)

| file | content | produced by |
|---|---|---|
| `javg_summary.csv`, `javg_terms.csv` | J-averaged NIST terms of 191 atoms and cations, lowest level vs lowest J-averaged term | `scripts/javg_nist.py` |
| `SI_javg_configurations.csv` | radii of Ni, W and Tb in their lowest J-averaged configuration | `workflow/make_javg_alt_inputs.py`, `scripts/audit_javg.py` |
| `sensitivity_three_radii.csv` | Table 3: mean and maximum change of v, r and c with the computational choices | `scripts/core_radii_all.py` |
| `desclaux_comparison.csv` | r_outer against Dirac-Fock orbital radii | `scripts/compare_desclaux.py` |
| `casscf_comparison_open_shell.csv` | state-averaged CASSCF vs single determinants | `scripts/analyse_casscf.py` |
| `bonding_configurations.csv` | radii of Pd, Pr, Nd in configurations for bonding | `scripts/audit_bond.py`, `scripts/analyse_casscf.py` (radii) |
| `crust_comparison.csv`, `Fig_crust.*` | v - r against r_vdW(Alvarez) - r_cov(Cordero), Figure 5 | `scripts/plot_crust.py` |
| `bond_length_test.csv` | 403 single bonds of Pyykko and Atsumi, classified, with predicted lengths | `scripts/bond_length_test.py` |
| `penetration_index_sets_p.csv`, `penetration_index_R2.csv` | Table 6: penetration indices and r2 with four combinations of radii | `scripts/penetration_index_test.py` |
| `FigS_desclaux.*` | Figure S1 | |

`results/core_radii_pbe0_anorcc.csv` and `results/core_radii_other_levels.csv` hold c_orb, c_count,
the ELF minima and the electrons they enclose for every species, written by `scripts/core_radii_all.py`;
`scripts/update_core_radii_tables.py` copies them into `data/radii/`.

`data/reference_data/` holds transcribed inputs of these tests: Desclaux's orbital radii, Table 4 of
Pyykko and Atsumi, the 14 sets of the penetration-index test (Pardos et al.), and the NIST ASD
energy levels used for the J-averaging (`nist_asd_levels.zip`, downloaded with `scripts/fetch_nist_levels.py`).
