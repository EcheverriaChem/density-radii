"""Orbital radii per angular momentum for the free atoms and singly charged cations.

For every species (PBE0/DKH2/ANO-RCC-TZP, the folder that gives the published v_rho) it writes
  r_s, r_p, r_d, r_f      outermost maximum of r^2|phi|^2 among the OCCUPIED spin orbitals of each l
  e_s, e_p, e_d, e_f      their orbital energies (Eh)
  rv_s, rv_p, rv_d        the same for the LOWEST VIRTUAL orbital of each l (alpha spin), with its
  ev_s, ev_p, ev_d        energy; Kohn-Sham virtuals depend on the basis, so treat these as indicative
together with v_rho (0.001 e/bohr^3), r_outer and c_orb as a check.

Usage (on the cluster, from the package root; no pandas needed):
    python scripts/orbital_radii_by_l.py --root /fs/agustina/jorge/density_radii --out orbital_radii_by_l.csv
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from density_radii.density import AtomDensity  # noqa: E402
from core_radii_all import index_folders, targets  # noqa: E402

LS = "spdf"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="orbital_radii_by_l.csv")
    ap.add_argument("--tol", type=float, default=0.003)
    ap.add_argument("--n-virtual", type=int, default=40)
    args = ap.parse_args()
    idx = index_folders(args.root)
    print(f"{sum(len(v) for v in idx.values())} job folders with Molden files under {args.root}", flush=True)
    rows = []
    for t in targets([("PBE0", "anorcc")]):
        if int(t["charge"]) not in (0, 1):
            continue
        job = t["job"]
        try:
            vsel = float(t["v_rho_0.001"])
        except ValueError:
            continue
        best = None
        for c in idx.get(job, []):
            try:
                a = AtomDensity(c / "job.molden.input")
                v = a.isodensity_radius(0.001)
            except Exception as err:
                print(f"  {job}: cannot read {c} ({err})")
                continue
            if best is None or abs(v - vsel) < best[0]:
                best = (abs(v - vsel), c, a, v)
        if best is None:
            print(f"{job:24s} no Molden file found")
            continue
        dv, folder, a, v = best
        row = {"job": job, "Z": t["Z"], "symbol": t["symbol"], "charge": t["charge"],
               "configuration": t["configuration"], "status": "ok" if dv <= args.tol else f"mismatch {dv:.4f}",
               "v_rho": round(v, 4)}
        occ = a.orbital_radii()                       # (rmax, l, energy, spin) of occupied spin orbitals
        for l in LS:
            mine = [o for o in occ if o[1] == l]
            if mine:
                o = max(mine, key=lambda x: x[0])
                row[f"r_{l}"], row[f"e_{l}"] = round(o[0], 4), round(o[2], 5)
            else:
                row[f"r_{l}"], row[f"e_{l}"] = "", ""
        row["r_outer"] = round(max(o[0] for o in occ), 4)
        try:
            mos, nocc = a.mo_characters("alpha", n_virtual=args.n_virtual)
            virt = mos[nocc:]
            for l in "spd":
                cand = [m for m in virt if m["l"] == l and m["w"] > 0.8]
                if cand:
                    m = min(cand, key=lambda x: x["energy"])
                    row[f"rv_{l}"], row[f"ev_{l}"] = round(m["rmax"], 4), round(m["energy"], 5)
                else:
                    row[f"rv_{l}"], row[f"ev_{l}"] = "", ""
        except ValueError as err:                     # level-shifted file: virtual order unusable
            for l in "spd":
                row[f"rv_{l}"], row[f"ev_{l}"] = "", ""
            row["status"] += f"; virtuals skipped ({err})"
        row["folder"] = str(folder)
        rows.append(row)
        print(f"{job:24s} {row['status'][:8]:8s} r_s={row['r_s']} r_p={row['r_p']} r_d={row['r_d']} "
              f"rv_p={row['rv_p']} rv_d={row['rv_d']}", flush=True)
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)} species -> {args.out}")


if __name__ == "__main__":
    main()
