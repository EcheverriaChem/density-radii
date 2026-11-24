"""Retry ladder for free-atom jobs whose SCF did not converge (open f shells).

Each element gets one directory containing three inputs that the job script runs in
sequence, stopping at the first one that converges:

  a: NormalSCF, first-order only (NoSOSCF, NoTRAH), large level shift, MaxIter 2000
  b: TightSCF restart from a.gbw, smaller shift, still first-order
  c: NormalSCF with a PModel guess and an even larger shift (independent attempt)

The converged result is copied to job.out / job.gbw / job.molden.input, so the usual
check_orca.py, audit_configurations.py and extract_atomic_radii.py all work unchanged.

Usage:
    python workflow/make_scf_retries.py --jobs 058_Ce_PBE0_anorcc 058_Ce_PBE0_x2c 092_U_PBE0_anorcc \
        --src $SCRATCH/density-radii/fblock_pilot \
        --outdir $SCRATCH/density-radii/fblock_retry
`--src` supplies the original job.inp of each job (its method, basis and multiplicity are reused).
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow"))
from make_free_atom_inputs import load_config  # noqa: E402

STRATEGIES = [
    ("a", "NormalSCF SlowConv NoSOSCF NoTRAH", "  MaxIter 2000\n  Shift Shift 0.40 ErrOff 0.0 end\n", None),
    ("b", "TightSCF SlowConv NoSOSCF NoTRAH MORead", "  MaxIter 1000\n  Shift Shift 0.15 ErrOff 0.0 end\n", "a"),
    ("c", "NormalSCF SlowConv NoSOSCF NoTRAH", "  MaxIter 3000\n  Shift Shift 0.70 ErrOff 0.0 end\n  Guess PModel\n", None),
]

RUNNER = """#!/bin/bash -l
#SBATCH -J {jobname}
#SBATCH -o logs/{jobname}.o%A_%a
#SBATCH -e logs/{jobname}.e%A_%a
#SBATCH -p {partition}
#SBATCH --account={account}
#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={mem_per_cpu}M
#SBATCH --time={walltime}
#SBATCH --array=1-{njobs}%50

echo "host: $(hostname)  task: $SLURM_ARRAY_TASK_ID"
cd "$SLURM_SUBMIT_DIR" || exit 4
JOBDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" joblist.txt)
cd "$JOBDIR" || exit 5
module load {orca_module} || {{ echo "ERROR: module load failed" > status.txt; exit 2; }}
{extra_module_line}
ORCA_EXE=$(command -v orca); ORCADIR=$(dirname "$ORCA_EXE")
[ -z "$ORCA_EXE" ] && {{ echo "ERROR: orca not found" > status.txt; exit 3; }}

: > status.txt
for s in a b c; do
  [ -f "$s.inp" ] || continue
  echo "--- strategy $s" | tee -a status.txt
  "$ORCA_EXE" "$s.inp" > "$s.out" 2> "$s.err"
  if grep -q "SCF CONVERGED" "$s.out"; then
     echo "converged with strategy $s" | tee -a status.txt
     cp "$s.inp" job.inp; cp "$s.out" job.out; cp "$s.gbw" job.gbw
     "$ORCADIR/orca_2mkl" "$s" -molden > /dev/null 2>&1 && cp "$s.molden.input" job.molden.input
     rm -f ?.tmp* *.tmp
     exit 0
  fi
  echo "strategy $s did not converge" | tee -a status.txt
done
echo "ALL STRATEGIES FAILED" | tee -a status.txt
exit 1
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", nargs="+", required=True, help="job directory names, e.g. 058_Ce_PBE0_anorcc")
    ap.add_argument("--src", required=True, help="folder holding the original job directories")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--walltime", default="48:00:00")
    args = ap.parse_args()
    cfg = load_config(args.config)
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    out = Path(args.outdir); (out / "logs").mkdir(parents=True, exist_ok=True)
    names = []
    for job in args.jobs:
        src = Path(args.src) / job / "job.inp"
        lines = src.read_text().splitlines()
        bang = next(l for l in lines if l.startswith("!"))
        head = bang.split("TightSCF")[0].replace("!", "").strip()          # method, Hamiltonian, basis, NoRI
        geom = lines[lines.index(next(l for l in lines if l.startswith("* xyz"))):]
        d = out / job; d.mkdir(exist_ok=True)
        for tag, keywords, scf, restart in STRATEGIES:
            txt = [f"# {job}: SCF retry strategy {tag}", f"! {head} {keywords}"]
            if restart:
                txt.append(f'%moinp "{restart}.gbw"')
            txt += [f"%pal nprocs {args.nprocs} end", f"%maxcore {fa['maxcore_mb']}", "%scf", scf.rstrip("\n"), "end"] + geom
            (d / f"{tag}.inp").write_text("\n".join(txt) + "\n")
        names.append(job)
    (out / "joblist.txt").write_text("\n".join(names) + "\n")
    extra = ag.get("extra_modules", "NONE")
    (out / "submit_retries.sh").write_text(RUNNER.format(
        jobname="pen_scfretry", partition=ag["partition"], account=ag["account"], nprocs=args.nprocs,
        walltime=args.walltime, njobs=len(names), mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3),
        orca_module=ag["orca_module"], extra_module_line="" if extra == "NONE" else f"module load {extra}"))
    print(f"{len(names)} retry directories in {out}\nsubmit with:  cd {out} && sbatch submit_retries.sh")


if __name__ == "__main__":
    main()
