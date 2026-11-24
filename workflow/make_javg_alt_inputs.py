"""Generate ORCA inputs for the J-averaged alternative configurations (R1.7, Step 1).

The species and configurations are read from data/atomic_states_javg.csv
(Ni 3d9 4s1, W 5d5 6s1, Tb 4f8 5d1 6s2). The inputs use exactly the same template,
SLURM script and settings as make_free_atom_inputs.py. f-block species get the
f-block recipe of the paper (first-order steps, level shift 0.40 Eh, NormalSCF,
up to 2000 cycles).

Usage (from the package root, with workflow/config.yaml filled in):
    python workflow/make_javg_alt_inputs.py --outdir $SCRATCH/density-radii/javg_alt
    python workflow/make_javg_alt_inputs.py --outdir ... --methods PBE0 HF --schemes anorcc x2c
Job folders are named like the primary set plus the tag _javg (e.g. 074_W_PBE0_anorcc_javg),
so the existing parsing and radius scripts read them unchanged.
"""
import argparse
import csv
from pathlib import Path

from make_free_atom_inputs import load_config, in_range

ROOT = Path(__file__).resolve().parent.parent
F_BLOCK = set(range(58, 72)) | set(range(90, 104))
FBLOCK_SCF = "NormalSCF SlowConv NoSOSCF NoTRAH"
FBLOCK_EXTRA = "Shift Shift 0.40 ErrOff 0.0 end"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--states", default=ROOT / "data" / "atomic_states_javg.csv")
    ap.add_argument("--methods", nargs="*", default=["PBE0"])
    ap.add_argument("--schemes", nargs="*", default=["anorcc"])
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--tag", default="_javg")
    ap.add_argument("--walltime", default="24:00:00")
    ap.add_argument("--purpose", default="J-averaged alternative configuration (R1.7)")
    ap.add_argument("--submit-name", default="submit_javg.sh")
    ap.add_argument("--jobname", default="javg_alt")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    tmpl = (ROOT / "workflow" / "templates" / "orca_free_atom.inp").read_text()
    slurm = (ROOT / "workflow" / "templates" / "slurm_orca_array_agustina.sh").read_text()
    out = Path(args.outdir)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    jobs = []
    for row in csv.DictReader(open(args.states, newline="")):
        sym, z, charge = row["symbol"], int(row["Z"]), int(row["charge"])
        if args.only and sym not in args.only:
            continue
        fblock = z in F_BLOCK
        for sch in args.schemes:
            sc = cfg[f"scheme_{sch}"]
            if not in_range(sym, sc["elements"]):
                print(f"skipped {sym}/{sch}: outside scheme range")
                continue
            for method in args.methods:
                sign = {1: "+", -1: "-", 0: ""}[charge]
                name = f"{z:03d}_{sym}{sign}_{method}_{sch}{args.tag}"
                d = out / name
                d.mkdir(exist_ok=True)
                (d / "job.inp").write_text(tmpl.format(
                    label=f"{sym}{sign} {row['ground_term']} ({row['configuration']}) 2S+1={row['multiplicity']} "
                          f"{method}/{sch} - {args.purpose}",
                    scf_type="UHF" if method == "HF" else "UKS",
                    method="" if method == "HF" else method,
                    hamiltonian=sc["hamiltonian"], basis=sc["basis"],
                    charge=charge, nprocs=fa["nprocs"], maxcore=fa["maxcore_mb"],
                    scf_keywords=FBLOCK_SCF if fblock else "TightSCF SlowConv",
                    maxiter=2000 if fblock else 500,
                    scf_extra=("  " + FBLOCK_EXTRA + "\n") if fblock else "",
                    mult=int(row["multiplicity"]), symbol=sym))
                jobs.append(name)
    (out / "joblist.txt").write_text("\n".join(jobs) + "\n")
    extra = ag.get("extra_modules", "NONE")
    (out / args.submit_name).write_text(slurm.format(
        jobname=args.jobname, partition=ag["partition"], account=ag["account"],
        nprocs=fa["nprocs"], walltime=args.walltime, njobs=len(jobs),
        mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3), orca_module=ag["orca_module"],
        extra_module_line="" if extra == "NONE" else f"module load {extra}"))
    print(f"{len(jobs)} inputs written to {out}:", " ".join(jobs))
    print(f"submit with:  cd {out} && sbatch {args.submit_name}")


if __name__ == "__main__":
    main()
