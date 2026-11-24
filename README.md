# density-radii

Van der Waals-like, valence and core radii of the free atoms and singly charged ions of
elements 1–96, derived from spherically averaged all-electron relativistic electron densities,
together with the code that produces and validates them.

Companion repository to:

> A. Urriolabeitia and J. Echeverría, *Van der Waals-Like, Valence and Core Radii of the Atoms and
> Singly Charged Ions of Elements 1–96 from the Same Relativistic Calculations*, submitted to
> J. Phys. Chem. A (2026).

Version 1.1.0 accompanies the revised manuscript; see `CHANGELOG.md` for what changed.

## What is here

Three radii are derived from one density by one recipe, so that they are mutually consistent and
available for the same elements:

| symbol | quantity | definition |
|---|---|---|
| `v` | van der Waals-like radius | outermost radius at which the spherically averaged density equals 0.001 e bohr⁻³ (also given at 0.0005 and 0.002) |
| `r` | valence radius | most probable radius of the outermost occupied orbital (outermost maximum of *r*²\|φ\|²) |
| `c` | core radius | most probable radius of the outermost core orbital (c_orb); the sphere enclosing the core electrons (`c_count`) and the outermost ELF minimum (`c_elf`) are reported alongside it |

Coverage: **96 neutral atoms, 94 cations M⁺ and 54 bound anions X⁻** at the primary level of
theory (PBE0, DKH2, ANO-RCC-TZP), plus three further method/basis combinations for the neutral
atoms and a sensitivity subset for the ions.

## Using the radii

The tables are plain CSV and need no code:

```
data/radii/density_pbe0_anorcc.csv        primary set: neutral atoms, kinds v, r, c
data/radii/density_pbe0_x2c.csv           X2C / x2c-TZVPPall
data/radii/density_hf_anorcc.csv          Hartree-Fock variants
data/radii/density_hf_x2c.csv
data/radii/density_cations_pbe0_anorcc.csv
data/radii/density_anions_pbe0_anorcc.csv
```

Each row carries the value, its source, the criterion it came from, and a note flagging anything
the user should know (a symmetry-broken solution, an anion whose computed electron affinity
differs from experiment by more than 0.3 eV). The same folder holds the published reference sets
transcribed for this work: Bondi, Mantina, Batsanov (crystallographic and isolated-atom),
Alvarez, Vogt and Alvarez (noble gases), Cordero, Pyykkö and Atsumi, Rahm (with the 2017
corrigendum), Shannon, and Charry and Tkatchenko.

From Python:

```python
from density_radii import load_set

load_set("density_pbe0_anorcc.csv", "v").value("Fe")   # 2.378 A, van der Waals-like
load_set("density_pbe0_anorcc.csv", "r").value("Fe")   # 1.333 A, valence
load_set("density_pbe0_anorcc.csv", "c").value("Fe")   # 0.379 A, core (c_orb)
```

`load_set` refuses to guess: asking for an element with several variants (a spin state, a
coordination number in the reference sets) without naming one raises instead of picking.

**Read `docs/radii_usage.md` before using the radii in anger.** It states the domain of validity:
heavy alkali metals are undersized by any single density threshold, weakly bound anions are not
bound at the density-functional level, and elements whose free-atom ground state differs from
their bonding configuration (Pd is the clear case) give valence radii that do not describe the
bonded atom.

**What each radius is for.** The van der Waals-like radius `v` can replace crystallographic van der
Waals radii in the penetration index of covalent bonds without loss of chemical information
(`scripts/penetration_index_test.py`). The valence and core radii describe the shell structure of
the free atom and follow the variation of size down each group, but they are not bonding radii:
sums of valence radii miss single-bond lengths by 0.24 Å on average, against 0.07 Å for the
covalent radii of Cordero et al. (`scripts/bond_length_test.py`).

## Reproducing the calculations

The quantum-chemistry outputs are large and are not distributed; everything needed to regenerate
them is. The workflow assumes ORCA 6.1.1 and a SLURM cluster.

```bash
cp workflow/config.example.yaml workflow/config.yaml     # edit paths, account, module name
python workflow/make_free_atom_inputs.py --outdir $SCRATCH/free_atoms
python workflow/make_free_atom_inputs.py --outdir $SCRATCH/ions --charge 1
```

f-block species need first-order convergence with a level shift:

```bash
python workflow/make_free_atom_inputs.py --outdir $SCRATCH/fblock --elements 57-71,89-96 \
    --scf-keywords "NormalSCF SlowConv NoSOSCF NoTRAH" \
    --scf-extra "Shift Shift 0.40 ErrOff 0.0 end" --maxiter 2000
```

Then check, audit, correct and analyse:

```bash
python workflow/check_orca.py $SCRATCH/free_atoms            # convergence, electron counts, files
python scripts/audit_configurations.py $SCRATCH/free_atoms   # populations vs intended configuration
python workflow/make_config_fixes.py --results ... --outdir ...   # DeltaSCF / orbital rotation
python scripts/extract_atomic_radii.py $SCRATCH/free_atoms results/raw.csv
python scripts/build_density_radii.py results/raw.csv results/fixes.csv
```

`docs/workflow.md` documents each generator and the order in which they are used, including the
analyses added in version 1.1 (J-averaged configurations, core radii, CASSCF tests, bonding
configurations, bond-length and penetration-index tests). These run from the CSV files in the
repository and need no quantum-chemistry output:

```bash
python scripts/compare_desclaux.py        # Table S8
python scripts/bond_length_test.py        # Section 5.1, Table 5
python scripts/penetration_index_test.py  # Section 5.2, Table 6
python scripts/plot_crust.py              # Figure 5
```

## Why the configuration audit matters

Of the 617 self-consistent solutions behind this work, 93 converged to a configuration other than
the experimental ground state, which changes the van der Waals-like radius by up to 0.48 Å and the
valence radius by up to 1.02 Å — far more than the functional, the Hamiltonian or the basis set.
The failure is silent: those calculations converge normally and report plausible energies. The
audit compares Löwdin populations with the intended configuration, classifies the deviation as an
occupation error or a symmetry-broken hybrid, and generates the corresponding correction. It is
cheap, and worth running for any published set of calculated atomic properties.

## Layout

```
density_radii/     package: radius tables, density analysis, ORCA output parsing
scripts/           analysis: extraction, audit, set building, comparisons, figures
workflow/          input generators and the SLURM templates
data/radii/        our radii and the transcribed published sets (CSV)
data/reference_data/  transcribed inputs of the tests (Desclaux, Pyykkö-Atsumi bonds, penetration-index sets, NIST levels)
data/              ground-state, J-averaged and bonding configurations of atoms and ions
results/           per-calculation tables, comparison statistics, figures
results/revision/  outputs of the analyses added in version 1.1
tests/             transcription and density-code tests
docs/              usage, workflow and data-format documentation
```

## Requirements

Python 3.10 or newer with NumPy, pandas, Matplotlib and PyYAML (`pip install -r requirements.txt`, or
`conda env create -f environment.yml`). ORCA 6.1.1 is needed only to regenerate the calculations.

```bash
pytest tests            # or: python tests/test_radius_tables.py
```

The tests check each transcribed table against a quantity stated in its own source (for example
the bond-length sums quoted by Pyykkö and Atsumi, or Shannon's constant 0.14 Å offset between
crystal and ionic radii), and check the density code against an analytic Gaussian density.

## License

Code: MIT (see `LICENSE`). The radii computed here are released with the article; the transcribed
reference sets remain the work of their original authors and are included with full citation so
that the comparisons are reproducible.
