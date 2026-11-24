"""Crust widths w = v - r (Cordero covalent radii) for every van der Waals set.

Writes results/crust_widths_by_set.csv and prints the elements whose crust
width is non-positive or unusually small (< 0.4 A), for which p_AB is
ill-defined or extremely sensitive.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from density_radii import load_set  # noqa: E402

cov = load_set("cordero2008_covalent.csv")
sets = {name: load_set(f) for name, f in [("alvarez2013", "alvarez2013_vdw.csv"),
                                          ("bondi1964", "bondi1964_vdw.csv"),
                                          ("mantina2009", "mantina2009_vdw.csv"),
                                          ("rahm2017", "rahm2017_density.csv")]}
out = Path(__file__).resolve().parent.parent / "results" / "crust_widths_by_set.csv"
rows, flags = [], []
for name, vdw in sets.items():
    for sym in vdw.symbols():
        if sym not in cov.symbols():
            flags.append(f"{name:12s} {sym:3s} no Cordero covalent radius")
            continue
        v = vdw.get(sym, "xray" if sym == "Pt" and name == "bondi1964" else None).value
        for var in cov.variants(sym):
            r = cov.get(sym, var).value
            w = round(v - r, 3)
            rows.append({"set": name, "symbol": sym, "cov_variant": var, "r": r, "v": v, "w": w})
            if w < 0.4:
                flags.append(f"{name:12s} {sym:3s}{('/' + var) if var else '':5s} r={r:.2f} v={v:.2f} w={w:+.2f}")
with open(out, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
    wr.writeheader()
    wr.writerows(rows)
print(f"{len(rows)} rows -> {out}\nProblematic crusts (w < 0.4 A):")
print("\n".join(flags))
