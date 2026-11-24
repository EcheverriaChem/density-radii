"""Valence radius against the Dirac-Fock orbital radii of Desclaux (Table S8, Figure S1).

data/reference_data/desclaux1973_rmax.csv holds the radius of maximum radial density (bohr) of
the outermost orbital of each atom, transcribed from J. P. Desclaux, At. Data Nucl. Data Tables
1973, 12, 311. In Desclaux's notation an orbital without an asterisk is j = l + 1/2 and one with
an asterisk is j = l - 1/2; the two are averaged with weights 2l + 2 and 2l, i.e. their
occupations in a closed subshell. For C-Ne the outermost maximum of r^2|phi|^2 belongs to the 2s
orbital, so the 2s radii of Desclaux are used for these elements.

Usage (from the package root):  python scripts/compare_desclaux.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
BOHR = 0.529177210903
L = {"s": 0, "p": 1, "d": 2, "f": 3}


def r_outer(fn):
    d = pd.read_csv(ROOT / "data" / "radii" / fn)
    return d[d.kind == "r"].set_index("symbol").value_A


def main():
    D = pd.read_csv(ROOT / "data" / "reference_data" / "desclaux1973_rmax.csv")
    l = D.orbital.str[-1].map(L)
    hi, lo = D.rmax_j_l_plus_half_bohr, D.rmax_j_l_minus_half_bohr
    w_hi, w_lo = 2 * l + 2, (2 * l).where(lo.notna(), 0)
    D["r_DF_avg_A"] = ((w_hi * hi + w_lo * lo.fillna(0)) / (w_hi + w_lo) * BOHR).round(4)
    for tag, fn in (("PBE0", "density_pbe0_anorcc.csv"), ("HF", "density_hf_anorcc.csv")):
        D[f"r_outer_{tag}_A"] = D.symbol.map(r_outer(fn))
        D[f"dev_{tag}"] = (D[f"r_outer_{tag}_A"] - D.r_DF_avg_A).round(4)
        print(f"{tag:4s} r_outer - r_DF:  MAD {D[f'dev_{tag}'].abs().mean():.3f} A, "
              f"max {D[f'dev_{tag}'].abs().max():.3f} A ({D.loc[D[f'dev_{tag}'].abs().idxmax(), 'symbol']}), "
              f"n = {D[f'dev_{tag}'].notna().sum()}")
    D.to_csv(ROOT / "results" / "revision" / "desclaux_comparison.csv", index=False)


if __name__ == "__main__":
    main()
