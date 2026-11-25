"""FHI-aims (all-electron, PBE, scalar ZORA) volume scans for the compression test.

Each element is computed in its ambient structure and in the bcc/fcc/hcp alternatives
(1 atom per cell, hcp 2), at 24 volumes from 1.05 V0 down to the compression given below.
Energies and the analytical stress (pressure) are written by aims; collect_aims_eos.py reads both.

Usage (on Agustina):
    module load gnu12 openmpi4/4.1.4 hwloc/2.7.0 fhi-aims/260331
    export AIMS_SPECIES=<.../species_defaults/defaults_2020/tight>   # see README
    python3 make_aims_eos.py --out eos_aims
    cd eos_aims && sbatch submit_eos.sh
"""
import argparse, math, os, sys
from pathlib import Path

ELEMENTS = {   # Z, structures, V0 per atom (A^3, experimental, sets the grid only), min V/V0
    "Li": (3, ["bcc", "fcc"], 21.3, 0.22),
    "Na": (11, ["bcc", "fcc"], 39.5, 0.17),
    "K": (19, ["bcc", "fcc"], 75.0, 0.20),
    "Rb": (37, ["bcc", "fcc"], 92.0, 0.22),
    "Cs": (55, ["bcc", "fcc"], 115.0, 0.25),
    "Mg": (12, ["hcp", "fcc", "bcc"], 23.2, 0.15),
    "Ca": (20, ["fcc", "bcc"], 43.6, 0.25),
    "Sr": (38, ["fcc", "bcc"], 56.3, 0.25),
    "Ba": (56, ["bcc"], 63.4, 0.25),
    "Al": (13, ["fcc", "hcp", "bcc"], 16.6, 0.15),
}
NPTS = 24
KMAX = 28             # cap per direction; 0.01 Ha smearing makes denser grids unnecessary
KSPACING = 0.12        # 1/A, VASP-like: n_i = |b_i| / KSPACING


def cell(struct, v):
    if struct == "fcc":
        a = (4 * v) ** (1 / 3)
        return [[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]], [[0, 0, 0]]
    if struct == "bcc":
        a = (2 * v) ** (1 / 3)
        return [[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]], [[0, 0, 0]]
    if struct == "hcp":
        ca = math.sqrt(8 / 3)
        a = (2 * v / (math.sqrt(3) / 2 * ca)) ** (1 / 3)
        return [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, ca * a]], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    raise ValueError(struct)


def kgrid(lat):
    a, b, c = lat
    cross = lambda u, w: [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    dot = lambda u, w: sum(x * y for x, y in zip(u, w))
    vol = dot(a, cross(b, c))
    rec = [cross(b, c), cross(c, a), cross(a, b)]
    return [min(KMAX, max(4, math.ceil(2 * math.pi * math.sqrt(dot(r, r)) / vol / KSPACING))) for r in rec]   # n_i = |b_i| / KSPACING


CONTROL = """# {tag}
xc                  pbe
relativistic        atomic_zora scalar
spin                none
occupation_type     gaussian 0.01
k_grid              {k[0]} {k[1]} {k[2]}
sc_accuracy_rho     1E-6
sc_accuracy_eev     1E-4
sc_accuracy_etot    1E-6
sc_iter_limit       200
compute_analytical_stress .true.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="eos_aims")
    ap.add_argument("--species", default=os.environ.get("AIMS_SPECIES"))
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--ntasks", type=int, default=16)
    args = ap.parse_args()
    if not args.species:
        sys.exit("set AIMS_SPECIES or --species")
    sdir, out = Path(args.species), Path(args.out)
    out.mkdir(exist_ok=True)
    jobs = []
    for el, (z, structs, v0, xmin) in ELEMENTS.items():
        if args.only and el not in args.only:
            continue
        sp = sorted(sdir.glob(f"{z:02d}_{el}_*"))
        if not sp:
            print(f"{el}: no species file {z:02d}_{el}_* in {sdir}, skipped"); continue
        species = sp[0].read_text()
        for st in structs:
            for i in range(NPTS):
                x = 1.05 * (xmin / 1.05) ** (i / (NPTS - 1))
                v = v0 * x
                tag = f"{el}_{st}_{i:02d}"
                d = out / tag
                d.mkdir(exist_ok=True)
                lat, pos = cell(st, v)
                with open(d / "geometry.in", "w") as f:
                    f.write(f"# {tag}  V/atom = {v:.5f} A^3\n")
                    for r in lat:
                        f.write("lattice_vector " + " ".join(f"{c:16.10f}" for c in r) + "\n")
                    for p in pos:
                        f.write("atom_frac " + " ".join(f"{c:14.10f}" for c in p) + f" {el}\n")
                (d / "control.in").write_text(CONTROL.format(tag=tag, k=kgrid(lat)) + "\n" + species)
                jobs.append(tag)
        print(f"{el}: species {sp[0].name}, {len(structs)} structures x {NPTS} volumes")
    groups = sorted({j.rsplit('_', 1)[0] for j in jobs})
    (out / "groups.txt").write_text("\n".join(groups) + "\n")
    (out / "submit_eos.sh").write_text(SUBMIT.format(n=len(groups), nt=args.ntasks))
    print(f"{len(jobs)} calculations in {len(groups)} groups -> {out}/submit_eos.sh")


SUBMIT = """#!/bin/bash -l
#SBATCH -J aims_eos
#SBATCH -o logs/eos.o%A_%a
#SBATCH -e logs/eos.e%A_%a
#SBATCH -p agustina_thin
#SBATCH --account=noncovnano
#SBATCH --nodes=1
#SBATCH --ntasks={nt}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=24:00:00
#SBATCH --array=1-{n}
# One element/structure per array task; its volumes run one after another.
cd "$SLURM_SUBMIT_DIR" || exit 4
mkdir -p logs
module load gnu12 openmpi4/4.1.4 hwloc/2.7.0 fhi-aims/260331 || exit 2
ulimit -s unlimited
export OMP_NUM_THREADS=1
export TMPDIR=/tmp/aims_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}; export OMPI_MCA_orte_tmpdir_base=$TMPDIR; mkdir -p $TMPDIR
G=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" groups.txt)
for d in ${{G}}_*; do
  [ -f "$d/aims.out" ] && grep -q "Have a nice day" "$d/aims.out" && continue
  ( cd "$d" && mpirun -np $SLURM_NTASKS $FHIAIMS_BINARY > aims.out 2> aims.err < /dev/null )
done
rm -rf $TMPDIR
"""

if __name__ == "__main__":
    main()
