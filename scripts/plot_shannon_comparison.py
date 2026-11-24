"""Figure 4: free-ion radii (this work) against Shannon's CN VI crystal radii for singly charged ions."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sh = {}
for r in csv.DictReader(open(ROOT / "data" / "radii" / "shannon1976_ionic.csv")):
    if r["variant"].endswith("_CR"):
        sh[(r["symbol"], 1 if r["variant"].startswith("+") else -1)] = float(r["value_A"])
ions = {(r["symbol"], int(r["charge"])): r for r in csv.DictReader(open(ROOT / "results" / "ion_radii_selected.csv"))
        if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"}
groups = {"alkali cations": ["Li", "Na", "K", "Rb", "Cs", "Fr"], "d$^{10}$ and 6s cations": ["Cu", "Ag", "Au", "Hg", "Tl"],
          "halide anions": ["F", "Cl", "Br", "I"]}
style = {"alkali cations": ("tab:red", "o"), "d$^{10}$ and 6s cations": ("tab:orange", "s"), "halide anions": ("tab:blue", "^")}
fig, ax = plt.subplots(figsize=(5.2, 4.6))
for g, syms in groups.items():
    q = -1 if "anion" in g else 1
    pts = [(sh[(s, q)], float(ions[(s, q)]["v_rho_0.001"]), s) for s in syms if (s, q) in sh and (s, q) in ions]
    c, m = style[g]
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], color=c, marker=m, s=30, label=g)
    for x, y, s in pts:
        ax.annotate(s + ("$^-$" if q == -1 else "$^+$"), (x, y), textcoords="offset points", xytext=(4, -3), fontsize=7)
ax.plot([0.8, 2.7], [0.8, 2.7], "-", color="0.6", lw=0.8, label="y = x")
ax.set_xlabel("Shannon crystal radius, CN VI (Å)")
ax.set_ylabel("$v_\\rho$ of the free ion, this work (Å)")
ax.set_xlim(0.8, 2.2); ax.set_ylim(0.8, 2.8)
ax.legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(ROOT / "results" / "shannon_comparison.png", dpi=200)
print("saved")
