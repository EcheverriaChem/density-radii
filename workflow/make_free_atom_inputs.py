"""Generate ORCA 6.1.1 free-atom inputs and a SLURM array script for Agustina.

Usage:
    python workflow/make_free_atom_inputs.py --outdir $SCRATCH/density-radii/free_atoms
    python workflow/make_free_atom_inputs.py --outdir ... --only H C O Zn Te I Fr --schemes x2c anorcc   # pilot
One job per (element, method, scheme); elements outside a scheme's range are skipped.
Only rows of data/atomic_states.csv whose status starts with "APPROVED" are used.
"""
import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
           "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
           "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
           "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm").split()


def load_config(path):
    """Minimal YAML reader (nested mappings, inline lists) - no PyYAML needed."""
    cfg, stack = {}, []
    for raw in open(path):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        key, _, val = line.strip().partition(":")
        val = val.strip()
        while stack and stack[-1][0] >= indent:
            stack.pop()
        parent = stack[-1][1] if stack else cfg
        if val == "":
            parent[key] = {}
            stack.append((indent, parent[key]))
        elif val.startswith("["):
            parent[key] = [v.strip() for v in val.strip("[]").split(",")]
        else:
            parent[key] = val.strip('"')
    return cfg


def in_range(symbol, rng):
    lo, hi = rng.split("-")
    z = SYMBOLS.index(symbol)
    return SYMBOLS.index(lo) <= z <= SYMBOLS.index(hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--only", nargs="*", help="restrict to these symbols (pilot)")
    ap.add_argument("--schemes", nargs="*", help="override schemes from config")
    ap.add_argument("--methods", nargs="*", help="override methods from config")
    ap.add_argument("--elements", help="restrict to Z ranges, e.g. 58-71,89-96 or 64")
    ap.add_argument("--charge", type=int, default=0, choices=[-1, 0, 1],
                    help="0: neutral atoms (atomic_states.csv); +1/-1: ions (atomic_states_ions.csv)")
    ap.add_argument("--basis", help="override the scheme basis (e.g. ANO-RCC-QZP)")
    ap.add_argument("--extra-keywords", default="", help="appended to the ! line (e.g. DecontractBas)")
    ap.add_argument("--tag", default="", help="suffix for job folder names (e.g. _decon)")
    ap.add_argument("--scf-keywords", default="TightSCF SlowConv",
                    help='SCF keywords; f-block recipe: "NormalSCF SlowConv NoSOSCF NoTRAH"')
    ap.add_argument("--scf-extra", default="", help='extra %%scf lines, e.g. "Shift Shift 0.40 ErrOff 0.0 end"')
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--nprocs", type=int, help="override cores per job")
    ap.add_argument("--walltime", help="override walltime, e.g. 12:00:00")
    args = ap.parse_args()

    if not Path(args.config).exists():
        raise SystemExit(f"{args.config} not found: copy workflow/config.example.yaml to "
                         "workflow/config.yaml and fill in account and module name")
    cfg = load_config(args.config)
    if "CHANGE_ME" in cfg["agustina"]["account"]:
        raise SystemExit("set agustina.account in workflow/config.yaml before generating jobs")
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    schemes = args.schemes or fa["schemes"]
    methods = args.methods or fa["methods"]
    tmpl = (ROOT / "workflow" / "templates" / "orca_free_atom.inp").read_text()
    slurm = (ROOT / "workflow" / "templates" / "slurm_orca_array_agustina.sh").read_text()

    out = Path(args.outdir)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    jobs, skipped = [], []
    states_file = "atomic_states.csv" if args.charge == 0 else "atomic_states_ions.csv"
    with open(ROOT / "data" / states_file, newline="") as fh:
        states = [r for r in csv.DictReader(fh) if r["status"].startswith("APPROVED")
                  and int(r.get("charge", 0) or 0) == args.charge]
    for row in states:
        sym = row["symbol"]
        if args.only and sym not in args.only:
            continue
        if args.elements:
            ranges = [tuple(map(int, (part.split("-") * 2)[:2])) for part in args.elements.split(",")]
            if not any(lo <= int(row["Z"]) <= hi for lo, hi in ranges):
                continue
        for sch in schemes:
            sc = cfg[f"scheme_{sch}"]
            if not in_range(sym, sc["elements"]):
                skipped.append(f"{sym}/{sch}")
                continue
            sign = {1: "+", -1: "-", 0: ""}[args.charge]
            for method in methods:
                scf = "UHF" if method == "HF" else "UKS"
                meth = "" if method == "HF" else method
                name = f"{int(row['Z']):03d}_{sym}{sign}_{method}_{sch}{args.tag}"
                d = out / name
                d.mkdir(exist_ok=True)
                (d / "job.inp").write_text(tmpl.format(
                    label=f"{sym}{sign} {row.get('ground_term', '')} ({row['configuration']}) 2S+1={row['multiplicity']} {method}/{sch}".replace("  ", " "),
                    scf_type=scf, method=meth, hamiltonian=sc["hamiltonian"],
                    basis=(args.basis or sc["basis"]) + (" " + args.extra_keywords if args.extra_keywords else ""),
                    charge=args.charge, nprocs=args.nprocs or fa["nprocs"], maxcore=fa["maxcore_mb"],
                    scf_keywords=args.scf_keywords, maxiter=args.maxiter,
                    scf_extra=("  " + args.scf_extra + "\n") if args.scf_extra else "",
                    mult=int(row["multiplicity"]), symbol=sym))
                jobs.append(name)
    # append to an existing job list, so several generator calls can fill one folder
    jl = out / "joblist.txt"
    old = [l.strip() for l in jl.read_text().splitlines() if l.strip()] if jl.exists() else []
    jobs = old + [j for j in jobs if j not in old]
    jl.write_text("\n".join(jobs) + "\n")
    extra = ag.get("extra_modules", "NONE")
    (out / "submit_free_atoms.sh").write_text(slurm.format(
        jobname="pen_atoms", partition=ag["partition"], account=ag["account"],
        nprocs=args.nprocs or fa["nprocs"], walltime=args.walltime or fa["walltime"], njobs=len(jobs),
        mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3),
        orca_module=ag["orca_module"],
        extra_module_line="" if extra == "NONE" else f"module load {extra}"))
    print(f"{len(jobs)} inputs written to {out}")
    if skipped:
        print("outside scheme range (skipped):", " ".join(skipped))
    print(f"submit with:  cd {out} && sbatch submit_free_atoms.sh")


if __name__ == "__main__":
    main()
