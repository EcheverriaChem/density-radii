"""Table of contents graphic: periodic table of the van der Waals-like radius v_rho.

Sized for the ACS TOC box (3.25 x 1.75 inch at 300 dpi and above).
"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

ROOT = Path(__file__).resolve().parent.parent
sel = {r["symbol"]: r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"}

# (symbol, column, row) on the standard 18-column layout; f block on two separate rows
MAIN = [("H",1,1),("He",18,1)]
ROW2 = ["Li","Be","B","C","N","O","F","Ne"]
ROW3 = ["Na","Mg","Al","Si","P","S","Cl","Ar"]
for i, s in enumerate(ROW2):
    MAIN.append((s, 1 + i if i < 2 else 13 + (i - 2), 2))
for i, s in enumerate(ROW3):
    MAIN.append((s, 1 + i if i < 2 else 13 + (i - 2), 3))
ROW4 = ["K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"]
ROW5 = ["Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"]
ROW6 = ["Cs","Ba","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"]
ROW7 = ["Fr","Ra"]
for row, series in [(4, ROW4), (5, ROW5), (6, ROW6), (7, ROW7)]:
    for i, s in enumerate(series):
        MAIN.append((s, i + 1, row))
LAN = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb"]
ACT = ["Ac","Th","Pa","U","Np","Pu","Am","Cm"]
for i, s in enumerate(LAN):
    MAIN.append((s, 3 + i, 8.5))
for i, s in enumerate(ACT):
    MAIN.append((s, 3 + i, 9.5))

vals = {s: float(r["v_rho_0.001"]) for s, r in sel.items()}
vmin, vmax = min(vals.values()), max(vals.values())
norm = colors.Normalize(vmin, vmax)
cmap = plt.get_cmap("viridis")

fig, ax = plt.subplots(figsize=(6.5, 3.5))
for sym, col, row in MAIN:
    if sym not in vals:
        continue
    v = vals[sym]
    x, y = col, -row
    ax.add_patch(plt.Rectangle((x, y), 0.94, 0.94, facecolor=cmap(norm(v)), edgecolor="white", lw=0.6))
    txt = "white" if norm(v) < 0.55 else "black"
    ax.text(x + 0.47, y + 0.60, sym, ha="center", va="center", fontsize=6.4, color=txt, fontweight="bold")
    ax.text(x + 0.47, y + 0.26, f"{v:.2f}", ha="center", va="center", fontsize=4.9, color=txt)

ax.set_xlim(0.6, 19.6); ax.set_ylim(-10.7, -0.2)
ax.set_aspect("equal"); ax.axis("off")
ax.text(3.0, -1.15, "Van der Waals-like, valence and core radii", fontsize=9.5, fontweight="bold", va="center")
ax.text(3.0, -1.95, "from one electron density, for elements 1–96", fontsize=8.2, va="center", color="0.25")
ax.text(3.0, -2.75, "atoms, cations and anions", fontsize=8.2, va="center", color="0.25")
ax.text(3.0, -9.9, "colour: $v_\\rho$ at 0.001 e bohr$^{-3}$, PBE0/DKH2/ANO-RCC-TZP", fontsize=6.6, va="center", color="0.35")

cax = fig.add_axes([0.62, 0.17, 0.30, 0.035])
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label("$v_\\rho$ (Å)", fontsize=7.5, labelpad=1)
cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.02)
fig.savefig(ROOT / "results" / "toc_periodic_table.png", dpi=600)
fig.savefig(ROOT / "results" / "toc_periodic_table.pdf")
print(f"{len(vals)} elements, v_rho from {vmin:.2f} to {vmax:.2f} A")
