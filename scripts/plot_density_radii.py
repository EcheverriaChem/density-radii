"""Figure: v_rho (this work) vs Rahm 2017 and Alvarez 2013 across the periodic table."""
import csv
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii import load_set  # noqa: E402

sel = [r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc"]
sel.sort(key=lambda r: int(r["Z"]))
rahm = load_set("rahm2017_density.csv"); alv = load_set("alvarez2013_vdw.csv")
fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), height_ratios=[2, 1], sharex=True)
Z = [int(r["Z"]) for r in sel]
ok = [r["config_ok"] == "True" for r in sel]
v = [float(r["v_rho_0.001"]) for r in sel]
ax[0].plot(Z, v, "-", color="0.6", lw=0.8)
ax[0].scatter([z for z, o in zip(Z, ok) if o], [x for x, o in zip(v, ok) if o], s=14, color="tab:blue", label="this work, PBE0/ANO-RCC")
ax[0].scatter([z for z, o in zip(Z, ok) if not o], [x for x, o in zip(v, ok) if not o], s=30, marker="x", color="tab:red", label="wrong/hybrid configuration")
ax[0].scatter(Z, [rahm.value(r["symbol"]) for r in sel], s=10, facecolors="none", edgecolors="k", label="Rahm 2016/2017")
za = [int(r["Z"]) for r in sel if r["symbol"] in alv.symbols()]
ax[0].scatter(za, [alv.value(r["symbol"]) for r in sel if r["symbol"] in alv.symbols()], s=10, marker="^", color="tab:green", label="Alvarez 2013")
ax[0].set_ylabel("radius (Å)"); ax[0].legend(fontsize=7, frameon=False, ncol=2)
ax[1].axhline(0, color="0.5", lw=0.6)
ax[1].scatter(Z, [x - rahm.value(r["symbol"]) for x, r in zip(v, sel)], s=12,
              c=["tab:blue" if o else "tab:red" for o in ok])
ax[1].set_ylabel("this work − Rahm (Å)"); ax[1].set_xlabel("Z")
for s_, z_ in [("Ni", 28), ("Tc", 43), ("Ir", 77)]:
    y_ = [x - rahm.value(r["symbol"]) for x, r in zip(v, sel) if r["symbol"] == s_]
    if y_:
        ax[1].annotate(s_, (z_, y_[0]), textcoords="offset points", xytext=(4, -3), fontsize=7)
fig.tight_layout()
fig.savefig(ROOT / "results" / "density_radii_vs_rahm_alvarez.png", dpi=200)
print("saved")
