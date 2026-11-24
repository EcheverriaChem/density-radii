"""Crust thickness, van der Waals-type minus covalent-type radius (Section 4, Figure 5).

Compares v_rho - r_outer (this work, PBE0/anorcc) with r_vdW(Alvarez; Vogt and Alvarez for the noble
gases) - r_cov(Cordero et al., sp3 C, low-spin Mn/Fe/Co)
and writes results/revision/crust_comparison.csv and Fig_crust.png/.eps.
H and the alkali metals are left out of the figure.

Usage (from the package root):  python scripts/plot_crust.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RAD = ROOT / "data" / "radii"
OUT = ROOT / "results" / "revision"
ALKALI = {"H", "Li", "Na", "K", "Rb", "Cs", "Fr"}


def radius_map(fn, kind, variants=("sp3", "ls")):
    d = pd.read_csv(RAD / fn)
    d = d[(d.kind == kind) & (d.variant.isna() | d.variant.isin(variants))]
    return d.drop_duplicates("symbol").set_index("symbol")


def main():
    dens = radius_map("density_pbe0_anorcc.csv", "v")
    d = pd.DataFrame({"symbol": dens.index, "v_rho": dens.value_A.values})
    zmap = pd.read_csv(ROOT / "results" / "atomic_radii_selected.csv").drop_duplicates("symbol").set_index("symbol").Z
    d["Z"] = d.symbol.map(zmap)
    d["r_outer"] = d.symbol.map(radius_map("density_pbe0_anorcc.csv", "r").value_A)
    alv = radius_map("alvarez2013_vdw.csv", "v").value_A.copy()
    alv.update(radius_map("vogt2014_noble_gas_vdw.csv", "v").value_A)     # noble gases: Vogt and Alvarez
    d["v_alv"] = d.symbol.map(alv)
    d["r_cor"] = d.symbol.map(radius_map("cordero2008_covalent.csv", "r").value_A)
    d["crust_calc"] = (d.v_rho - d.r_outer).round(3)
    d["crust_emp"] = (d.v_alv - d.r_cor).round(3)
    d["diff"] = (d.crust_calc - d.crust_emp).round(3)
    d = d.sort_values("Z")
    d[["symbol", "Z", "v_rho", "r_outer", "v_alv", "r_cor", "crust_calc", "crust_emp", "diff"]].to_csv(
        OUT / "crust_comparison.csv", index=False)
    p = d[~d.symbol.isin(ALKALI)]
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(p.Z, p.crust_calc, "-", color="0.6", lw=0.8, zorder=1)
    ax.scatter(p.Z, p.crust_calc, s=18, color="tab:blue", zorder=2, label=r"this work, $v_\rho - r_{\rm outer}$")
    e = p.dropna(subset=["crust_emp"])
    ax.scatter(e.Z, e.crust_emp, s=18, marker="s", facecolors="none", edgecolors="tab:red", lw=1, zorder=3,
               label="Alvarez 2013 $-$ Cordero 2008")
    ax.set_xlabel("Z")
    ax.set_ylabel("crust thickness (Å)")
    ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
    ax.set_xlim(-3, 100)
    ax.set_ylim(0.55, 1.75)
    fig.tight_layout()
    fig.savefig(OUT / "Fig_crust.png", dpi=200)
    fig.savefig(OUT / "Fig_crust.eps")
    print(f"{len(d)} elements -> {OUT / 'crust_comparison.csv'}; r2(calc, emp) without H and alkali metals = "
          f"{p.dropna(subset=['crust_emp'])[['crust_calc', 'crust_emp']].corr().iloc[0, 1] ** 2:.3f}")


if __name__ == "__main__":
    main()
