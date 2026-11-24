"""Generate DeltaSCF (MOM) restart inputs for free-atom jobs that converged to the wrong
orbital occupation (s/d, f/d or f/s swaps, one electron per spin channel).

Usage (run locally or on Agustina, where the original job.gbw files live):
    python workflow/make_config_fixes.py --results /path/to/free_atoms \
        --cluster-dir $SCRATCH/density-radii/free_atoms \
        --outdir $SCRATCH/density-radii/config_fixes

--results     folder with job.out / job.molden.input (can be a local copy)
--cluster-dir folder on Agustina holding the original job.gbw files (MORead source)
Jobs with s-d *hybridised* frontier orbitals (symmetry breaking, e.g. Zr, La, W) are
reported but not changed: MOM cannot remove the hybridisation.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "workflow"))
from density_radii.density import AtomDensity  # noqa: E402
from density_radii.orca_parse import (read, reduced_populations, l_totals, valence_shells,  # noqa: E402
                                 expected_spin_populations, observed_spin_populations,
                                 load_states, parse_job_name)
from make_free_atom_inputs import load_config  # noqa: E402

TOL = 0.15


def conf_string(i_empty, i_fill, i_homo):
    """ORCA ALPHACONF/BETACONF: occupations listed from the lowest to the highest
    orbital involved. Works whether the orbital to fill lies above or below the one
    to empty (the latter happens when a virtual sits below the occupied one)."""
    lo, hi = min(i_empty, i_fill), max(i_empty, i_fill)
    occ = []
    for i in range(lo, hi + 1):
        if i == i_empty:
            occ.append(0)
        elif i == i_fill:
            occ.append(1)
        else:
            occ.append(1 if i <= i_homo else 0)
    return ",".join(map(str, occ))


def conf_string_multi(empties, fills, i_homo):
    """Occupation string for several simultaneous moves in one spin channel."""
    lo, hi = min(empties + fills), max(empties + fills)
    occ = []
    for i in range(lo, hi + 1):
        if i in empties:
            occ.append(0)
        elif i in fills:
            occ.append(1)
        else:
            occ.append(1 if i <= i_homo else 0)
    return ",".join(map(str, occ))


def plan_fix(a, config, pop, tol=0.4):
    """Per spin, work out which orbital must be emptied and which filled.

    Deviations of the per-spin l populations from the configuration (Hund's rule within
    each open shell) give the excess l (electron to remove) and the deficient l (electron
    to add). Orbitals are then identified by their dominant l character: the emptied one
    is the highest occupied orbital of the excess l, the filled one the lowest virtual of
    the deficient l (its valence orbital).

    Returns ({spin: (i_empty, i_fill, i_homo, info)}, notes).
    """
    exp = expected_spin_populations(config)
    obs = observed_spin_populations(pop)
    plan, notes = {}, []
    for spin in ("alpha", "beta"):
        dev = {l: obs[spin][l] - exp[spin][l] for l in "spdfg"}
        # electrons to remove per l (excess) and to add per l (deficit), rounded
        remove = {l: int(round(v)) for l, v in dev.items() if v > tol}
        add = {l: int(round(-v)) for l, v in dev.items() if v < -tol}
        remove = {l: n for l, n in remove.items() if n > 0}
        add = {l: n for l, n in add.items() if n > 0}
        if not remove and not add:
            continue
        n_rm, n_add = sum(remove.values()), sum(add.values())
        if n_rm != n_add or n_rm > 2:
            notes.append(f"{spin}: deviations do not balance or exceed two electrons "
                         f"{({l: round(v, 2) for l, v in dev.items() if abs(v) > 0.15})}")
            continue
        ch, nocc = a.mo_characters(spin, n_virtual=16)
        empties, fills, labels = [], [], []
        for l, n in remove.items():
            cand = sorted([c for c in ch if c["occ"] > 0 and c["l"] == l and c["energy"] > -3.0 and c["w"] > 0.5],
                          key=lambda c: -c["energy"])[:n]
            if len(cand) < n:
                break
            empties += cand
        for l, n in add.items():
            cand = sorted([c for c in ch if c["occ"] == 0 and c["l"] == l and c["w"] > 0.5],
                          key=lambda c: c["energy"])[:n]
            if len(cand) < n:
                break
            fills += cand
        if len(empties) != n_rm or len(fills) != n_add:
            notes.append(f"{spin}: cannot identify the orbitals for {remove} -> {add}")
            continue
        hyb = [c for c in empties + fills if c["w"] < 0.9]
        if hyb:
            notes.append(f"{spin}: orbitals involved are hybridised (weight {min(c['w'] for c in hyb):.2f})")
        info = f"{spin}: empty " + ", ".join(f"MO {c['i']+1} ({c['l']})" for c in empties) + \
               "; fill " + ", ".join(f"MO {c['i']+1} ({c['l']})" for c in fills)
        plan[spin] = ([c["i"] for c in empties], [c["i"] for c in fills], nocc - 1, info)
    return plan, notes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--cluster-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=ROOT / "workflow" / "config.yaml")
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--walltime", help="override walltime, e.g. 24:00:00")
    ap.add_argument("--jobs", nargs="*", help="only these job folders (e.g. 062_Sm+_PBE0_anorcc)")
    ap.add_argument("--nprocs", type=int,
                    help="cores per job; default: the %%pal value of the source input (script and input must agree)")
    ap.add_argument("--scf-extra", default="", help='extra %%scf line, e.g. "Shift Shift 0.40 ErrOff 0.0 end"')
    ap.add_argument("--mode", default="deltascf",
                    choices=["deltascf", "deltascf-soscf", "deltascf-fresh", "rotate"],
                    help="deltascf: DeltaSCF+MORead (default); deltascf-soscf: same but keep SOSCF/TRAH; "
                         "deltascf-fresh: DeltaSCF without MORead; rotate: MORead + %%scf Rotate (no DeltaSCF)")
    ap.add_argument("--diagnose", action="store_true",
                    help="write all four modes as subfolders A-D of --outdir, to find what an ORCA build accepts")
    args = ap.parse_args()
    res, out = Path(args.results), Path(args.outdir)
    cfg_early = load_config(args.config)
    fa_nprocs = cfg_early["free_atoms"]["nprocs"]
    job_nprocs = []
    states = load_states()
    jobs, report = [], []
    for d in sorted(p for p in res.iterdir() if (p / "job.molden.input").exists()):
        if args.jobs and d.name not in args.jobs:
            continue
        z, sym, charge, method, scheme = parse_job_name(d.name)
        conf = states[(sym, charge)]["configuration"]
        try:
            pop = reduced_populations(read(d / "job.out"))
        except (ValueError, FileNotFoundError):
            report.append(f"{d.name:22s} SKIPPED: no population block (SCF did not converge or job unfinished)")
            continue
        exp = l_totals(conf)
        dev = max((pop["charge"][l] - exp[l] for l in "spdf"), key=abs)
        if abs(dev) <= TOL:
            continue
        plan, notes = plan_fix(AtomDensity(d / "job.molden.input"), conf, pop)
        if not plan:
            report.append(f"{d.name:22s} dev={dev:+.2f}  NOT FIXED  " + "; ".join(notes))
            continue
        inp_lines = (d / "job.inp").read_text().splitlines()
        bang = next(l for l in inp_lines if l.startswith("!"))
        src_pal = next((int(m.group(1)) for l in inp_lines
                        for m in [re.search(r"%pal\s+nprocs\s+(\d+)", l)] if m), None)
        nprocs = args.nprocs or src_pal or int(fa_nprocs)
        job_nprocs.append(nprocs)
        geom = inp_lines[inp_lines.index(next(l for l in inp_lines if l.startswith("* xyz"))):]
        def build(mode):
            scf = ["%scf", f"  MaxIter {args.maxiter}"]
            if args.scf_extra:
                scf.append("  " + args.scf_extra)
            head = bang
            if mode.startswith("deltascf"):
                head = head + " DeltaSCF"
                if mode == "deltascf-soscf":      # let ORCA use its default second-order steps
                    head = head.replace(" NoSOSCF", "").replace(" NoTRAH", "")
                for spin, (ie, jf, ih, _info) in plan.items():
                    scf.append(f"  {'ALPHACONF' if spin == 'alpha' else 'BETACONF'} {conf_string_multi(ie, jf, ih)}")
            else:   # rotate: swap the two orbitals by 90 degrees, ORCA counts MOs from 0
                rot = " ".join(f"{{{e},{f},90,{0 if spin == 'alpha' else 1},{0 if spin == 'alpha' else 1}}}"
                               for spin, (ie, jf, _ih, _i) in plan.items() for e, f in zip(ie, jf))
                scf.append(f"  Rotate {rot} end")
            scf.append("end")
            lines = [f"# config fix of {d.name} [{mode}]: " + " | ".join(p[3] for p in plan.values())]
            if mode == "deltascf-fresh":
                lines.append(head)
            else:
                lines += [head + " MORead", f'%moinp "{gbw}"']
            lines += [f"%pal nprocs {nprocs} end",
                      next(l for l in inp_lines if l.startswith("%maxcore"))] + scf + geom
            return "\n".join(lines) + "\n"

        gbw = f"{args.cluster_dir.rstrip('/')}/{d.name}/job.gbw"
        modes = ["deltascf", "deltascf-soscf", "deltascf-fresh", "rotate"] if args.diagnose else [args.mode]
        for mi, mode in enumerate(modes):
            base = out / "ABCD"[mi] if args.diagnose else out
            (base / d.name).mkdir(parents=True, exist_ok=True)
            (base / d.name / "job.inp").write_text(build(mode))
        jobs.append(d.name)
        report.append(f"{d.name:22s} dev={dev:+.2f}  FIX  " + " | ".join(p[3] for p in plan.values())
                      + (("  [" + "; ".join(notes) + "]") if notes else ""))
    print("\n".join(report))
    if not jobs:
        return
    out.mkdir(parents=True, exist_ok=True)
    folders = [out / c for c in "ABCD"] if args.diagnose else [out]
    for fo in folders:
        (fo / "logs").mkdir(parents=True, exist_ok=True)
        (fo / "joblist.txt").write_text("\n".join(jobs) + "\n")
    cfg = load_config(args.config)
    ag, fa = cfg["agustina"], cfg["free_atoms"]
    extra = ag.get("extra_modules", "NONE")
    slurm = (ROOT / "workflow" / "templates" / "slurm_orca_array_agustina.sh").read_text()
    for fo in folders:
      (fo / "submit_free_atoms.sh").write_text(slurm.format(
        jobname="pen_fix", partition=ag["partition"], account=ag["account"], nprocs=max(job_nprocs),
        walltime=fa["walltime"], njobs=len(jobs), mem_per_cpu=int(int(fa["maxcore_mb"]) * 1.3),
        orca_module=ag["orca_module"],
        extra_module_line="" if extra == "NONE" else f"module load {extra}"))
    print(f"\n{len(jobs)} fix inputs written to {', '.join(str(f) for f in folders)}")


if __name__ == "__main__":
    main()
