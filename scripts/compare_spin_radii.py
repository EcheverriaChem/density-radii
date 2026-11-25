"""Test the spin-resolved isodensity radius against polarizability radii and van der Waals dimers.

Input: the CSV written by scripts/spin_radii.py.
Reference data (data/reference_data/):
  polarizability_gould2016_atoms.csv / _cations.csv   Gould and Bucko, J. Chem. Theory Comput. 2016
  vdw_dimers.csv                                       equilibrium distances of van der Waals dimers
R_alpha = 2.54 alpha^(1/7) bohr (Fedorov, Sadhukhan, Stoehr, Tkatchenko, PRL 2018).

Usage:  python scripts/compare_spin_radii.py spin_radii.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data" / "reference_data"
BOHR = 0.529177210903
ALK = {"H", "Li", "Na", "K", "Rb", "Cs", "Fr"}


def stats(x, t):
    q = x / t
    return f"ratio {q.mean():.3f} +- {q.std():.3f}  MAD {abs(x - t).mean():.3f}  r2 {np.corrcoef(x, t)[0, 1] ** 2:.3f}"


def main():
    S = pd.read_csv(sys.argv[1])
    S = S[S.status == "ok"]
    for charge, fn in ((0, "polarizability_gould2016_atoms.csv"), (1, "polarizability_gould2016_cations.csv")):
        a = pd.read_csv(REF / fn).set_index("symbol").alpha
        s = S[S.charge == charge].set_index("symbol")
        el = [x for x in s.index if x in a.index]
        t = 2.54 * a[el] ** (1 / 7) * BOHR
        print(f"\n{'atoms' if charge == 0 else 'cations'}: radius vs 2.54 alpha^(1/7)")
        for sub, e in (("all", el), ("H and alkali metals", [x for x in el if x in ALK]),
                       ("others", [x for x in el if x not in ALK])):
            if len(e) < 3:
                continue
            print(f"  {sub:20s} n={len(e):3d}")
            for k in ("v_tot", "v_maj"):
                print(f"     {k:6s} {stats(s.loc[e, k].astype(float), t[e])}")
        if charge == 0:
            out = s.loc[el, ["v_tot", "v_maj", "zeta"]].assign(R_alpha=t.round(4))
            out.to_csv(Path(sys.argv[1]).with_name("spin_radii_vs_alpha.csv"))
            worst = (out.v_maj - out.R_alpha).abs().sort_values(ascending=False).head(10)
            print("  largest |v_maj - R_alpha|:", worst.round(2).to_dict())
    D = pd.read_csv(REF / "vdw_dimers.csv")
    s = S[S.charge == 0].set_index("symbol")
    print("\nvan der Waals dimers, R_e / (rA + rB)")
    for cls in [None] + list(D["class"].unique()):
        d = D if cls is None else D[D["class"] == cls]
        line = f"  {str(cls or 'all'):22s} n={len(d):2d}"
        for k in ("v_tot", "v_maj"):
            sm = d.A.map(s[k]) + d.B.map(s[k])
            q = d.Re_A / sm
            line += f"   {k}: {q.mean():.3f} +- {q.std():.3f} (MAE {abs(sm - d.Re_A).mean():.2f})"
        print(line)


if __name__ == "__main__":
    main()
