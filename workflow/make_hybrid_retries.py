"""PBE0 retries for s/d- or f/d-hybridised free atoms, starting from the purer HF orbitals.

Variant A: plain DIIS (no SOSCF/TRAH), UseSym, MORead from the HF solution with the
           correct configuration. A symmetric (unhybridised) start can stay symmetric.
Variant B: as A plus a level shift, to keep the SCF close to the starting orbitals.

Usage (on Agustina):
    python workflow/make_hybrid_retries.py \
        --hf-dirs $SCRATCH/density-radii/config_fixes $SCRATCH/density-radii/free_atoms \
        --outdir $SCRATCH/density-radii/hybrid_retries
Creates <outdir>/A and <outdir>/B, each with its own joblist.txt and submit script.
The first --hf-dirs entry that contains <Z>_<El>_HF_<scheme>/job.gbw is used.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow"))
from make_free_atom_inputs import load_config  # noqa: E402

DEFAULT_TARGETS = ["Zr", "Tc", "W"]
SCHEMES = {"x2c": "X2C x2c-TZVPPall", "anorcc": "DKH ANO-RCC-TZP"}
def variants(maxiter, extra):
    base = f"  MaxIter {maxiter}\n" + (f"  {extra}\n" if extra else "")
    return {"A": base, "B": base + "  Shift Shift 0.30 ErrOff 0.0 end\n" if "Shift" not in extra else base}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hf-dirs", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--only", nargs="*", default=DEFAULT_TARGETS, help="element symbols to retry")
    ap.add_argument("--schemes", nargs="*", default=["x2c", "anorcc"])
    ap.add_argument("--scf-keywords", default="TightSCF", help='f block: "NormalSCF SlowConv"')
    ap.add_argument("--scf-extra", default="", help='extra %%scf line, e.g. "Shift Shift 0.40 ErrOff 0.0 end"')
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--nprocs", type=int)
    ap.add_argument("--walltime")
    args = ap.parse_args()
    import csv as _csv
    states = {r["symbol"]: r for r in _csv.DictReader(open(ROOT / "data" / "atomic_states.csv"))
              if r["status"].startswith("APPROVED")}
    targets = [(f"{int(states[s]['Z']):03d}", s, int(states[s]["multiplicity"])) for s in args.only]
    cfg = load_config(args.config)
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    slurm = (ROOT / "workflow" / "templates" / "slurm_orca_array_agustina.sh").read_text()
    extra = ag.get("extra_modules", "NONE")
    for var, scf in variants(args.maxiter, args.scf_extra).items():
        out = Path(args.outdir) / var
        (out / "logs").mkdir(parents=True, exist_ok=True)
        jobs = []
        for z, el, mult in targets:
            for sch, kw in {k: v for k, v in SCHEMES.items() if k in args.schemes}.items():
                src = f"{z}_{el}_HF_{sch}"
                gbw = next((f"{d.rstrip('/')}/{src}/job.gbw" for d in args.hf_dirs
                            if Path(f"{d.rstrip('/')}/{src}/job.gbw").exists()), None)
                if gbw is None:
                    gbw = f"{args.hf_dirs[0].rstrip('/')}/{src}/job.gbw"
                    print(f"warning: {src}/job.gbw not found in --hf-dirs (expected when generating off-cluster)")
                name = f"{z}_{el}_PBE0_{sch}"
                (out / name).mkdir(exist_ok=True)
                (out / name / "job.inp").write_text(
                    f"# hybrid retry {var} for {name}, start from {src}\n"
                    f"! UKS PBE0 {kw} NoRI {args.scf_keywords} UseSym MORead NoTRAH NoSOSCF\n"
                    f'%moinp "{gbw}"\n'
                    f"%pal nprocs {args.nprocs or fa['nprocs']} end\n%maxcore {fa['maxcore_mb']}\n"
                    f"%scf\n{scf}end\n* xyz 0 {mult}\n{el}  0.0  0.0  0.0\n*\n")
                jobs.append(name)
        (out / "joblist.txt").write_text("\n".join(jobs) + "\n")
        (out / "submit_free_atoms.sh").write_text(slurm.format(
            jobname=f"pen_hyb{var}", partition=ag["partition"], account=ag["account"],
            nprocs=args.nprocs or fa["nprocs"], walltime=args.walltime or fa["walltime"], njobs=len(jobs),
            mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3), orca_module=ag["orca_module"],
            extra_module_line="" if extra == "NONE" else f"module load {extra}"))
        print(f"variant {var}: {len(jobs)} inputs in {out}")


if __name__ == "__main__":
    main()
