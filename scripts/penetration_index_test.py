"""Penetration index with density-derived radii (Section 5.2, Table 6).

p_AB = 100 (v_A + v_B - d_AB) / (v_A + v_B - r_A - r_B)

The 14 sets of systems (distance d, covalency indicator y and the penetration index reported in
the original work, p_paper) are in data/reference_data/penetration_index_sets.csv. For each set
the script reports the coefficient of determination of the linear correlation of y with d and
with p computed from four combinations of radii:
  Alvarez + Cordero     crystallographic van der Waals radii (Vogt and Alvarez for the noble
                        gases) and covalent radii (sp2 C for C2, sp C for CN and BC, sp3 C otherwise,
                        low-spin Mn/Fe/Co)
  v_rho + Cordero       van der Waals-like radius of this work with the covalent radii
  Alvarez + r_outer     crystallographic van der Waals radii with the valence radius of this work
  v_rho + r_outer       both radii from this work
  Rahm + Cordero        density radii of Rahm et al. (2016/2017) with the covalent radii
  Charry + Cordero      polarizability-derived radii of Charry and Tkatchenko (2024) with the covalent radii
The first combination reproduces p_paper (largest deviation printed as a check).

Writes results/revision/penetration_index_R2.csv (the r2 of Table 6, p_paper being the published
index) and, with --out, the per-system penetration indices.

Usage (from the package root):  python scripts/penetration_index_test.py [--out file.csv]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAD = ROOT / "data" / "radii"


def radius_map(fn, kind, variants=()):
    d = pd.read_csv(RAD / fn)
    d = d[(d.kind == kind) & (d.variant.isna() | d.variant.isin(variants))]
    return d.drop_duplicates("symbol").set_index("symbol").value_A


def r2(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return np.corrcoef(x[ok], y[ok])[0, 1] ** 2 if ok.sum() > 2 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="write the per-system penetration indices here")
    args = ap.parse_args()
    D = pd.read_csv(ROOT / "data" / "reference_data" / "penetration_index_sets.csv")
    alv = radius_map("alvarez2013_vdw.csv", "v")
    alv.update(radius_map("vogt2014_noble_gas_vdw.csv", "v"))
    cor3 = radius_map("cordero2008_covalent.csv", "r", ("sp3", "ls"))
    cor2 = radius_map("cordero2008_covalent.csv", "r", ("sp2", "ls"))
    cor1 = radius_map("cordero2008_covalent.csv", "r", ("sp", "ls"))
    vd = radius_map("density_pbe0_anorcc.csv", "v")
    rd = radius_map("density_pbe0_anorcc.csv", "r")
    rahm = radius_map("rahm2017_density.csv", "v")
    cha = radius_map("charry2024_vdw.csv", "v")

    def cordero(row, sym):
        if row.A == row.B == "C":
            return cor2[sym]                  # C2: sp2 carbon, as in the original work
        if row.set == "diatomic ions":
            return cor1[sym]                  # CN+, CN-, BC-: sp carbon, as in the original work
        return cor3[sym]

    def p(row, v, r):
        vA, vB = v.get(row.A, np.nan), v.get(row.B, np.nan)
        rA, rB = r(row, row.A), r(row, row.B)
        return 100 * (vA + vB - row.d) / (vA + vB - rA - rB)

    combos = {"Alvarez+Cordero": (alv, cordero),
              "v_rho+Cordero": (vd, cordero),
              "Alvarez+r_outer": (alv, lambda row, s: rd.get(s, np.nan)),
              "v_rho+r_outer": (vd, lambda row, s: rd.get(s, np.nan)),
              "Rahm+Cordero": (rahm, cordero),
              "Charry+Cordero": (cha, cordero)}
    for k, (v, r) in combos.items():
        D[f"p_{k}"] = D.apply(lambda row: p(row, v, r), axis=1)
    dev = (D["p_Alvarez+Cordero"] - D.p_paper).abs().max()
    print(f"check: largest |p(Alvarez+Cordero) - p_paper| = {dev:.3f}\n")
    rows = []
    for (s, lab), g in D.groupby(["set", "ylabel"], sort=False):
        y = g.y.values.astype(float)
        row = {"set": s, "n": len(g), "descriptor": lab, "distance": r2(g.d.values.astype(float), y),
               "p_paper": r2(g.p_paper.values.astype(float), y)}
        row.update({k: r2(g[f"p_{k}"].values.astype(float), y) for k in combos})
        rows.append(row)
        if s == "X2 diatomics":               # second row of Table 6: without the alkali metals
            g = g[~g.A.isin(["Li", "Na", "K", "Rb", "Cs", "Fr"])]
            y = g.y.values.astype(float)
            row = {"set": s + " without alkali metals", "n": len(g), "descriptor": lab,
                   "distance": r2(g.d.values.astype(float), y), "p_paper": r2(g.p_paper.values.astype(float), y)}
            row.update({k: r2(g[f"p_{k}"].values.astype(float), y) for k in combos})
            rows.append(row)
    T = pd.DataFrame(rows)
    print(T.round(3).to_string(index=False))
    T.round(3).to_csv(ROOT / "results" / "revision" / "penetration_index_R2.csv", index=False)
    if args.out:
        D.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
