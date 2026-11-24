"""State-averaged CASSCF tests for the atoms that stay symmetry-broken in a single determinant.

La ([Xe]5d1 6s2), Ce ([Xe]4f1 5d1 6s2) and Th ([Rn]6d2 7s2) have no spherically symmetric
single-determinant solution of their ground configuration: the SCF mixes ns with (n-1)d_z2,
which belong to the same irreducible representation of every Abelian point group. Averaging the
density over the degenerate components of the atomic term removes the artefact. ORCA stores the
state-averaged natural orbitals with their fractional occupations in the gbw file, which
orca_2mkl converts to a Molden file that the analysis code reads unchanged.

The integral transformation is done exactly by default: "TrafoStep RI" needs an auxiliary basis,
which ANO-RCC and x2c-TZVPPall do not provide, and would require AutoAux (use --trafostep RI to
get both).

The active spaces hold the valence shells of each configuration, and the number of roots equals
the number of components of the ground term (D: 5, F: 7, G: 9). Check the root energies afterwards:
the average is only spherical if the averaged roots are degenerate. For La the five roots agree to
1e-9 Eh; for Th and Ce the default counts do not close a term (splittings of 400 and 1400 cm-1), so
use --nroots to include the remaining components (for example 10 for Th: 3F + 3P).

Usage:
    python workflow/make_casscf_tests.py --outdir $SCRATCH/density-radii/casscf_step11   # B C O F Al Cl Sc Ti Fe
    python workflow/make_casscf_tests.py --outdir ... --only La Th Ce                     # the original tests
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow"))
from make_free_atom_inputs import load_config  # noqa: E402

# symbol: (Z, multiplicity, active electrons, active orbitals, roots, term, active space)
TARGETS = {
    "La": ("057", 2, 3, 6, 5, "2D", "6s + 5d"),
    "Th": ("090", 3, 4, 6, 7, "3F", "7s + 6d"),
    "Ce": ("058", 3, 4, 13, 9, "1G/3H manifold", "6s + 5d + 4f"),
    # Step 11 (R1.4): open-shell p and d atoms whose single determinant is not spherical.
    # The active space holds the whole valence shell, so the averaged roots span the full term.
    "B": ("005", 2, 3, 4, 3, "2P", "2s + 2p"),
    "C": ("006", 3, 4, 4, 3, "3P", "2s + 2p"),
    "O": ("008", 3, 6, 4, 3, "3P", "2s + 2p"),
    "F": ("009", 2, 7, 4, 3, "2P", "2s + 2p"),
    "Al": ("013", 2, 3, 4, 3, "2P", "3s + 3p"),
    "Cl": ("017", 2, 7, 4, 3, "2P", "3s + 3p"),
    "Sc": ("021", 2, 3, 6, 5, "2D", "4s + 3d"),
    "Ti": ("022", 3, 4, 6, 7, "3F", "4s + 3d"),
    "Fe": ("026", 5, 8, 6, 5, "5D", "4s + 3d"),
}
STEP11 = ["B", "C", "O", "F", "Al", "Cl", "Sc", "Ti", "Fe"]
SCHEMES = {"anorcc": "DKH ANO-RCC-TZP", "x2c": "X2C x2c-TZVPPall"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--only", nargs="*", default=STEP11)
    ap.add_argument("--schemes", nargs="*", default=["anorcc"])
    ap.add_argument("--nroots", type=int, help="override the number of averaged roots")
    ap.add_argument("--trafostep", default="Exact", choices=["Exact", "RI"],
                    help="Exact (default; cheap for a single atom) or RI, which also adds AutoAux")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--walltime", default="24:00:00")
    ap.add_argument("--submit-name", default="submit_free_atoms.sh")
    args = ap.parse_args()
    cfg = load_config(args.config)
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    out = Path(args.outdir)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    jobs = []
    for el in args.only:
        z, mult, nel, norb, nroots, term, space = TARGETS[el]
        nroots = args.nroots or nroots
        for sch in args.schemes:
            name = f"{z}_{el}_CASSCF_{sch}"
            (out / name).mkdir(exist_ok=True)
            (out / name / "job.inp").write_text(
                f"# state-averaged CASSCF({nel},{norb}) for {el}, {nroots} roots of the {term} term ({space})\n"
                f"! CASSCF {SCHEMES[sch]} NoRI TightSCF"
                f"{' AutoAux' if args.trafostep == 'RI' else ''}\n"
                f"%pal nprocs {args.nprocs} end\n%maxcore {fa['maxcore_mb']}\n"
                f"%casscf\n  nel {nel}\n  norb {norb}\n  mult {mult}\n  nroots {nroots}\n"
                f"  maxiter 300\n  TrafoStep {args.trafostep}\nend\n"
                f"* xyz 0 {mult}\n{el}  0.0  0.0  0.0\n*\n")
            jobs.append(name)
    (out / "joblist.txt").write_text("\n".join(jobs) + "\n")
    slurm = (ROOT / "workflow" / "templates" / "slurm_orca_array_agustina.sh").read_text()
    extra = ag.get("extra_modules", "NONE")
    (out / args.submit_name).write_text(slurm.format(
        jobname="pen_casscf", partition=ag["partition"], account=ag["account"], nprocs=args.nprocs,
        walltime=args.walltime, njobs=len(jobs), mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3),
        orca_module=ag["orca_module"],
        extra_module_line="" if extra == "NONE" else f"module load {extra}"))
    print(f"{len(jobs)} CASSCF inputs in {out}\nsubmit with:  cd {out} && sbatch {args.submit_name}")


if __name__ == "__main__":
    main()
