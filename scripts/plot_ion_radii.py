"""Figure 3: v_rho of cations, neutral atoms and anions against Z, and the change on ionization/attachment."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
neu = {r["symbol"]: r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"}
ions = [r for r in csv.DictReader(open(ROOT / "results" / "ion_radii_selected.csv"))
        if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"]
cat = sorted([r for r in ions if r["charge"] == "1"], key=lambda r: int(r["Z"]))
an = sorted([r for r in ions if r["charge"] == "-1"], key=lambda r: int(r["Z"]))
nz = sorted(neu.values(), key=lambda r: int(r["Z"]))

fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.8), height_ratios=[2, 1.1], sharex=True)
ax[0].plot([int(r["Z"]) for r in nz], [float(r["v_rho_0.001"]) for r in nz], "-o", ms=3, lw=0.8, color="0.35", label="neutral atoms")
ax[0].plot([int(r["Z"]) for r in cat], [float(r["v_rho_0.001"]) for r in cat], "-s", ms=3, lw=0.8, color="tab:red", label="cations M$^+$")
good = [r for r in an if r["EA_within_0.3eV"] == "True"]
off = [r for r in an if r["EA_within_0.3eV"] != "True"]
ax[0].scatter([int(r["Z"]) for r in good], [float(r["v_rho_0.001"]) for r in good], s=16, marker="^", color="tab:blue", label="anions X$^-$ (EA within 0.3 eV)", zorder=3)
ax[0].scatter([int(r["Z"]) for r in off], [float(r["v_rho_0.001"]) for r in off], s=18, marker="^", facecolors="none", edgecolors="tab:blue", label="anions X$^-$ (EA error > 0.3 eV)", zorder=3)
ax[0].set_ylabel("$v_\\rho$ (0.001 e bohr$^{-3}$), Å")
ax[0].legend(fontsize=7, frameon=False, ncol=2, loc="upper left")

ax[1].axhline(0, color="0.5", lw=0.6)
dc = [(int(r["Z"]), float(r["v_rho_0.001"]) - float(neu[r["symbol"]]["v_rho_0.001"])) for r in cat if r["symbol"] in neu]
da = [(int(r["Z"]), float(r["v_rho_0.001"]) - float(neu[r["symbol"]]["v_rho_0.001"]), r["EA_within_0.3eV"] == "True") for r in an if r["symbol"] in neu]
ax[1].bar([z for z, _ in dc], [d for _, d in dc], width=0.8, color="tab:red", alpha=0.75, label="M$^+$ − M")
ax[1].bar([z for z, _, _ in da], [d for _, d, _ in da], width=0.8, color="tab:blue", alpha=0.75, label="X$^-$ − X")
# label the alkali cations just below the tip of their bar, clear of every bar
for sym, z in [("Li", 3), ("Na", 11), ("K", 19), ("Rb", 37), ("Cs", 55), ("Fr", 87)]:
    v = [d for zz, d in dc if zz == z]
    if v:
        ax[1].annotate(sym, (z, v[0]), xytext=(0, -9), textcoords="offset points",
                       ha="center", va="top", fontsize=7, color="tab:red")
ax[1].set_ylim(min(d for _, d in dc) - 0.22, max(d for _, d, _ in da) * 1.15)
ax[1].set_ylabel("change in $v_\\rho$ (Å)"); ax[1].set_xlabel("Z")
ax[1].legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(ROOT / "results" / "ion_radii.png", dpi=200)
print(f"cations {len(cat)}, anions {len(an)} ({len(good)} within 0.3 eV)")
dcv = [d for _, d in dc]; dav = [d for _, d, _ in da]
print(f"contraction on ionization: mean {sum(dcv)/len(dcv):+.3f}, range {min(dcv):+.3f} to {max(dcv):+.3f}")
print(f"expansion on attachment: mean {sum(dav)/len(dav):+.3f}, range {min(dav):+.3f} to {max(dav):+.3f}")
