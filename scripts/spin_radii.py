"""Spin-resolved isodensity radii of the free atoms and ions.

v_tot   outermost r where the total density rho = rho_c (0.001 e/bohr^3), the v_rho of the paper
v_maj   outermost r where the majority-spin density rho_sigma = rho_c/2
v_min   the same for the minority spin (empty for H and one-electron species)
zeta    spin polarization (rho_a - rho_b)/rho at v_tot

For a closed shell rho_a = rho_b = rho/2, so v_maj = v_min = v_tot. In H and the alkali metals the
tail is made of one spin only, and v_maj lies outside v_tot.

Usage (on the cluster, from the package root):
    python scripts/spin_radii.py --root /fs/agustina/jorge/density_radii --out spin_radii.csv
It picks for each species the folder whose v_rho matches the selected value of the paper (as
scripts/core_radii_all.py does), so run it over the same root.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from density_radii.density import AtomDensity, BOHR_TO_ANG  # noqa: E402
from core_radii_all import index_folders, targets  # noqa: E402


def iso(r, rho, c):
    lr = np.log(np.maximum(rho, 1e-300))
    above = np.where(lr >= math.log(c))[0]
    if len(above) == 0 or above[-1] == len(r) - 1:
        return float("nan")
    j = above[-1]
    x = r[j] + (math.log(c) - lr[j]) * (r[j + 1] - r[j]) / (lr[j + 1] - lr[j])
    return x * BOHR_TO_ANG


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="spin_radii.csv")
    ap.add_argument("--rho", type=float, default=0.001)
    ap.add_argument("--tol", type=float, default=0.003)
    args = ap.parse_args()
    idx = index_folders(args.root)
    print(f"{sum(len(v) for v in idx.values())} job folders with Molden files under {args.root}", flush=True)
    r = np.linspace(0.5, 25.0, 5000)
    rows = []
    for t in targets([("PBE0", "anorcc")]):
        job = t["job"]
        try:
            vsel = float(t["v_rho_0.001"])
        except ValueError:
            continue
        best = None
        for c in idx.get(job, []):
            try:
                a = AtomDensity(c / "job.molden.input")
                ra, rb = a.rho(r, "alpha"), a.rho(r, "beta")
            except Exception as err:
                print(f"  {job}: cannot read {c} ({err})")
                continue
            v = iso(r, ra + rb, args.rho)
            dv = abs(v - vsel)
            if best is None or dv < best[0]:
                best = (dv, c, ra, rb, v)
        if best is None:
            print(f"{job:24s} no Molden file found")
            continue
        dv, folder, ra, rb, v = best
        maj, mnr = (ra, rb) if ra[-500:].sum() >= rb[-500:].sum() else (rb, ra)
        v_maj = iso(r, maj, args.rho / 2)
        v_min = iso(r, mnr, args.rho / 2) if mnr.max() > args.rho / 2 else float("nan")
        k = np.argmin(abs(r * BOHR_TO_ANG - v))
        zeta = (ra[k] - rb[k]) / (ra[k] + rb[k])
        row = {"job": job, "Z": t["Z"], "symbol": t["symbol"], "charge": t["charge"],
               "configuration": t["configuration"], "status": "ok" if dv <= args.tol else f"mismatch {dv:.4f}",
               "v_tot": round(v, 4), "v_maj": round(v_maj, 4),
               "v_min": "" if math.isnan(v_min) else round(v_min, 4), "zeta": round(float(zeta), 3),
               "folder": str(folder)}
        rows.append(row)
        print(f"{job:24s} {row['status']:8s} v_tot={row['v_tot']:.3f} v_maj={row['v_maj']:.3f} "
              f"v_min={row['v_min']} zeta={row['zeta']:+.2f}", flush=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)} species -> {args.out}")


if __name__ == "__main__":
    main()
