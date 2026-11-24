"""Figure 2: r_outer (most probable radius of the outermost occupied orbital) vs Cordero covalent radii."""
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii import load_set  # noqa: E402

cor = load_set("cordero2008_covalent.csv")
def cv(s):
    try:
        return cor.value(s)
    except ValueError:
        return cor.value(s, "ls")

sel = [r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"
       and r["symbol"] in cor.symbols()]
sel.sort(key=lambda r: int(r["Z"]))
x = np.array([cv(r["symbol"]) for r in sel])
y = np.array([float(r["r_outer"]) for r in sel])
Z = np.array([int(r["Z"]) for r in sel])
block = []
for r in sel:
    z = int(r["Z"])
    if r["outer_l"] == "p":
        block.append("p block")
    elif 57 <= z <= 71 or 89 <= z <= 96:
        block.append("f block")
    elif (21 <= z <= 30) or (39 <= z <= 48) or (72 <= z <= 80):
        block.append("d block")
    else:
        block.append("s block")
block = np.array(block)
slope, intercept = np.polyfit(x, y, 1)
r2 = np.corrcoef(x, y)[0, 1] ** 2

fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.3), width_ratios=[1.15, 1])
colors = {"s block": "tab:red", "p block": "tab:blue", "d block": "tab:green", "f block": "tab:purple"}
for b, c in colors.items():
    m = block == b
    ax[0].scatter(x[m], y[m], s=18, color=c, label=b, alpha=0.85)
lim = [0.2, 2.8]
ax[0].plot(lim, lim, "-", color="0.6", lw=0.8, label="y = x")
xs = np.linspace(*lim, 10)
ax[0].plot(xs, slope * xs + intercept, "--", color="k", lw=0.9,
           label=f"fit: {slope:.2f}x − {abs(intercept):.2f}  (r² = {r2:.3f})")
for r, xi, yi in zip(sel, x, y):
    if abs(yi - xi) > 0.35:
        ax[0].annotate(r["symbol"], (xi, yi), textcoords="offset points", xytext=(4, -2), fontsize=7)
ax[0].set_xlabel("covalent radius, Cordero 2008 (Å)")
ax[0].set_ylabel("$r_{outer}$, this work (Å)")
ax[0].set_xlim(lim); ax[0].set_ylim(lim)
ax[0].legend(fontsize=7, frameon=False, loc="upper left")

ax[1].axhline(0, color="0.5", lw=0.6)
for b, c in colors.items():
    m = block == b
    ax[1].scatter(Z[m], y[m] - x[m], s=16, color=c)
for r, zi, d in zip(sel, Z, y - x):
    if abs(d) > 0.35:
        ax[1].annotate(r["symbol"], (zi, d), textcoords="offset points", xytext=(3, -2), fontsize=7)
ax[1].set_xlabel("Z"); ax[1].set_ylabel("$r_{outer}$ − $r_{Cordero}$ (Å)")
fig.tight_layout()
fig.savefig(ROOT / "results" / "router_vs_cordero.png", dpi=200)
print(f"n={len(sel)} slope={slope:.3f} intercept={intercept:.3f} r2={r2:.3f} "
      f"mean diff={np.mean(y-x):+.3f} sd={np.std(y-x, ddof=1):.3f}")
for b in colors:
    m = block == b
    print(f"  {b}: n={m.sum():2d} mean {np.mean(y[m]-x[m]):+.3f} sd {np.std(y[m]-x[m], ddof=1):.3f} r2 {np.corrcoef(x[m],y[m])[0,1]**2:.3f}")
print("largest deviations:", [(r["symbol"], round(d, 2)) for r, d in sorted(zip(sel, y - x), key=lambda t: -abs(t[1]))[:8]])
