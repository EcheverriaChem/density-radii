"""Bond-length test of the valence radius (Section 5.1, Table 5).

Bond lengths: the 403 single bonds of Table 4 of Pyykko and Atsumi (Chem. Eur. J. 2009, 15, 186),
classified in results/revision/bond_length_test.csv (experimental or calculated, ion, solid,
elements 1-96). Predictions are sums of radii, R(A-B) = r(A) + r(B), recomputed here from
data/radii/: r_outer (this work, PBE0/anorcc), Cordero et al. (sp3 C, low-spin Mn/Fe/Co) and
Pyykko and Atsumi.

Usage (from the package root):  python scripts/bond_length_test.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAD = ROOT / "data" / "radii"
GROUPS = {1: ["Li", "Na", "K", "Rb", "Cs", "Fr"], 2: ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
          13: ["B", "Al", "Ga", "In", "Tl"], 14: ["C", "Si", "Ge", "Sn", "Pb"],
          15: ["N", "P", "As", "Sb", "Bi"], 16: ["O", "S", "Se", "Te", "Po"],
          17: ["F", "Cl", "Br", "I", "At"]}


def radius_map(fn, kind="r", variants=("sp3", "ls")):
    d = pd.read_csv(RAD / fn)
    d = d[(d.kind == kind) & (d.variant.isna() | d.variant.isin(variants))]
    return d.drop_duplicates("symbol").set_index("symbol").value_A


def main():
    d = pd.read_csv(ROOT / "results" / "revision" / "bond_length_test.csv")
    sets = {"r_outer": radius_map("density_pbe0_anorcc.csv"),
            "Cordero": radius_map("cordero2008_covalent.csv"),
            "Pyykko": radius_map("pyykko2009_covalent.csv")}
    for k, r in sets.items():
        d[f"pred_{k}"] = d.sym1.map(r) + d.sym2.map(r)
    inr = d[d.in_range_Z_le_96]
    subsets = {"all bonds, Z <= 96": inr,
               "experimental, neutral, Z <= 96": inr[(inr.source_type == "experimental") & (~inr.ion)],
               "experimental, neutral, no hydrogen": inr[(inr.source_type == "experimental") & (~inr.ion)
                                                         & (inr.sym1 != "H") & (inr.sym2 != "H")]}
    print("Mean absolute error (A)")
    for name, s in subsets.items():
        print(f"  {name:34s} n={len(s):3d}  " + "  ".join(
            f"{k} {(s[f'pred_{k}'] - s.R_A).abs().mean():.3f}" for k in sets))
    # linear rescaling R = a + b (rA + rB) for r_outer, without hydrogen
    s = inr[(inr.sym1 != "H") & (inr.sym2 != "H")].dropna(subset=["pred_r_outer"])
    b, a = np.polyfit(s.pred_r_outer, s.R_A, 1)
    print(f"  r_outer rescaled (all bonds without H, Z <= 96): R = {a:.3f} + {b:.3f} x sum, MAE "
          f"{(a + b * s.pred_r_outer - s.R_A).abs().mean():.3f}")
    # Table 5: mean A-H bond length of each element against its radius, down each group
    ah = inr[(inr["class"] == "A–H") & (~inr.ion)].copy()
    ah["A"] = np.where(ah.sym1 == "H", ah.sym2, ah.sym1)
    mean = ah.groupby("A").R_A.mean()
    print("\nTable 5: r2 of the mean A-H bond length against the radius of A, down each group")
    for g, els in GROUPS.items():
        e = [x for x in els if x in mean.index]
        line = f"  group {g:2d} ({len(e)} el.)"
        for k in ("r_outer", "Cordero"):
            x = sets[k].reindex(e).values
            line += f"  {k} {np.corrcoef(x, mean[e].values)[0, 1] ** 2:.3f}"
        print(line)


if __name__ == "__main__":
    main()
