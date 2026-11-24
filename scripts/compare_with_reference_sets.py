"""Compare the density-derived radii (PBE0/ANO-RCC) with the empirical and literature sets.

v_rho is compared with van der Waals-type sets (Alvarez 2013, Batsanov 2001 crystallographic and
equilibrium, Rahm 2017); r_outer with covalent sets (Cordero 2008, Pyykko 2009). Writes
results/reference_set_comparison.csv and prints a summary.
"""
import csv
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii import load_set  # noqa: E402

NOBLE_PD = ("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Pd")
sel = [r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc" and r["usable"] == "True"]
ours = {"v_rho": {r["symbol"]: float(r["v_rho_0.001"]) for r in sel},
        "r_outer": {r["symbol"]: float(r["r_outer"]) for r in sel}}


def val(S, s):
    try:
        return S.value(s)
    except ValueError:
        return S.value(s, "ls")


cases = [("v_rho", "alvarez2013_vdw.csv", ()), ("v_rho", "rahm2017_density.csv", ()),
         ("v_rho", "charry2024_vdw.csv", ()), ("v_rho", "charry2024_vdw.csv", ("Li", "Na", "K", "Rb", "Cs", "Fr")),
         ("v_rho", "vogt2014_noble_gas_vdw.csv", ()),
         ("v_rho", "batsanov2001_vdw_crystallographic.csv", ()),
         ("v_rho", "batsanov2001_vdw_equilibrium.csv", ("Bi",)),   # Bi: misprint in the source
         ("r_outer", "cordero2008_covalent.csv", ()), ("r_outer", "pyykko2009_covalent.csv", ()),
         ("r_outer", "cordero2008_covalent.csv", NOBLE_PD), ("r_outer", "pyykko2009_covalent.csv", NOBLE_PD)]
rows = []
for quantity, fn, excl in cases:
    ref = load_set(fn)
    k = [s for s in ours[quantity] if s in ref.symbols() and s not in excl]
    x = np.array([val(ref, s) for s in k]); y = np.array([ours[quantity][s] for s in k])
    slope, icpt = np.polyfit(x, y, 1)
    rows.append({"quantity": quantity, "reference": ref.set_id, "excluded": " ".join(excl), "n": len(k),
                 "r2": round(np.corrcoef(x, y)[0, 1] ** 2, 3), "slope": round(slope, 3), "intercept": round(icpt, 3),
                 "mean_diff": round(float(np.mean(y - x)), 3), "sd_diff": round(float(np.std(y - x, ddof=1)), 3),
                 "mad": round(float(np.mean(abs(y - x))), 3)})
out = ROOT / "results" / "reference_set_comparison.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['quantity']:8s} vs {r['reference']:22s} excl[{r['excluded']:>20s}] n={r['n']:2d} r2={r['r2']:.3f} "
          f"slope={r['slope']:.2f} mean={r['mean_diff']:+.3f} MAD={r['mad']:.3f}")
