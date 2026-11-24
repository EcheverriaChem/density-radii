# Workflow

The calculations were run on a SLURM cluster with ORCA 6.1.1. Every generator writes a folder of
job directories, a `joblist.txt` and a job-array submission script; nothing needs editing by hand.

## 0. Configuration

    cp workflow/config.example.yaml workflow/config.yaml

Set the scratch path, the SLURM partition and account, the ORCA module name, the cores per job
and the memory. `config.yaml` is not tracked by git, so updates never overwrite local settings.

## 1. Inputs

`make_free_atom_inputs.py` builds the free-atom and free-ion inputs from the approved
configurations in `data/`.

| option | purpose |
|---|---|
| `--charge {-1,0,1}` | neutral atoms, cations or anions |
| `--elements 57-71,89-96` | restrict to atomic-number ranges |
| `--only Fe Co` | restrict to elements |
| `--methods PBE0 HF`, `--schemes anorcc x2c` | level of theory |
| `--scf-keywords`, `--scf-extra`, `--maxiter` | convergence settings |
| `--basis`, `--extra-keywords`, `--tag` | basis overrides, e.g. `DecontractBas` |
| `--nprocs`, `--walltime` | resources; the job script matches the input's `%pal` |

Job lists accumulate, so several calls can fill one folder. Keep species that need different SCF
settings in different folders: the submission script is written for the inputs beside it.

## 2. Checking

`workflow/check_orca.py` reports, per job: termination status, SCF convergence, spin
contamination, integrated electron count against Z, and whether the Molden file was written.
`scripts/audit_configurations.py` compares Loewdin populations with the intended configuration
and writes `configuration_audit.csv`.

## 3. Correcting

| situation | tool |
|---|---|
| SCF did not converge | `make_scf_retries.py` - ladder of three strategies (level shift and first-order steps; tight restart from those orbitals; model-potential guess), run in sequence until one converges |
| wrong occupation | `make_config_fixes.py` - DeltaSCF with the maximum overlap method (`--mode deltascf`) or a 90-degree orbital rotation (`--mode rotate`); handles up to two electrons per spin channel; `--jobs` restricts to named jobs |
| symmetry-broken s-d hybrid | `make_hybrid_retries.py` - restart PBE0 from the purer Hartree-Fock orbitals with symmetry enabled (works for the d block; fails for the f block) |
| no single-determinant solution | `make_casscf_tests.py` - state-averaged CASSCF over the components of the ground term |

A practical warning: with a level shift, ORCA writes the occupied orbitals first and the virtuals
after, which is **not** energy order. Orbital indices for `ALPHACONF`, `BETACONF` and `Rotate`
must follow the file order. The code does this; if you write your own, check it.

## 4. Analysis

    python scripts/extract_atomic_radii.py <folder> results/raw.csv
    python scripts/build_density_radii.py results/raw.csv results/fixes.csv ...
    python scripts/build_ion_radii.py <base folders> --fixes <folders> --neutral-dirs <folders>

`build_*` merges the base calculations with the corrections, preferring a solution that restores
the configuration or reduces the deviation, and records in every row which calculation supplied
the value. Then:

    python scripts/compare_with_reference_sets.py    # regressions against the published sets
    python scripts/plot_density_radii.py             # Figure 1
    python scripts/plot_router_vs_cordero.py         # Figure 2
    python scripts/plot_ion_radii.py                 # Figure 3
    python scripts/plot_shannon_comparison.py        # Figure 4

## 5. Revision analyses (version 1.1)

    python scripts/fetch_nist_levels.py --outdir nist_levels          # NIST ASD levels, atoms and cations
    python scripts/javg_nist.py nist_levels --out results/revision    # lowest level vs lowest J-averaged term
    python workflow/make_javg_alt_inputs.py --outdir $SCRATCH/javg_alt   # Ni, W, Tb in the J-averaged configuration
    python scripts/audit_javg.py $SCRATCH/javg_alt
    python scripts/core_radii_all.py --root $SCRATCH --out results/core_radii_pbe0_anorcc.csv
    python scripts/core_radii_all.py --root $SCRATCH --levels HF/anorcc HF/x2c PBE0/x2c --out results/core_radii_other_levels.csv
    python scripts/update_core_radii_tables.py                       # c = c_orb in data/radii
    python workflow/make_javg_alt_inputs.py --states data/atomic_states_bond.csv --tag _bond \
        --purpose "bonding configuration" --submit-name submit_bond.sh --outdir $SCRATCH/bond_configs
    python scripts/audit_bond.py $SCRATCH/bond_configs
    python workflow/make_casscf_tests.py --outdir $SCRATCH/casscf_open_shell   # B C O F Al Cl Sc Ti Fe
    python scripts/analyse_casscf.py --casscf $SCRATCH/casscf_open_shell --sd-root $SCRATCH
    python scripts/compare_desclaux.py
    python scripts/bond_length_test.py
    python scripts/penetration_index_test.py
    python scripts/plot_crust.py

The bonding configurations of the f elements usually collapse to the ground configuration; restart
them from the converged orbitals with DeltaSCF or a 90-degree rotation (the inputs used for Ce, Pr and
Nd are in `workflow/examples/bond_fixes`, with orbital indices taken in file order) and audit again. The tetravalent Ce 5d2 6s2
could not be stabilized.

**Several ORCA jobs on one node.** ORCA's OpenMPI keeps its session files in /tmp. When two array
tasks of the same user share a node, the first to finish removes them and the others die without
an error message. The SLURM template therefore gives every task its own `TMPDIR`
(`OMPI_MCA_orte_tmpdir_base`). Keep this if you write your own scripts, and redirect the standard
input of ORCA from /dev/null when you run several jobs in a `while read` loop.
