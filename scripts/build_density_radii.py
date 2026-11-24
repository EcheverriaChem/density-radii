"""Merge original and configuration-fixed free-atom results into the density radius sets.

Usage: python scripts/build_density_radii.py results/atomic_radii_raw.csv [fixes1.csv fixes2.csv ...]

Outputs
  results/atomic_radii_selected.csv   one row per (element, method, scheme): the job used
  data/radii/density_<method>_<scheme>.csv   radius sets in density_radii format (kinds v, r, c)
  results/density_radii_summary.csv   per element: PBE0/ANO-RCC values + spread over the 4 variants
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
raw = list(csv.DictReader(open(sys.argv[1])))
# later files take precedence; only fixes with the correct configuration are kept
fixmap = {}
for fname in sys.argv[2:]:
    for f in csv.DictReader(open(fname)):
        key = (f["symbol"], f["method"], f["scheme"])
        if f["config_ok"] == "True" or key not in fixmap:
            fixmap[key] = dict(f, source_file=Path(fname).name)

selected = []
raw.sort(key=lambda r: (int(r["Z"]), r["method"], r["scheme"]))
for x in raw:
    key = (x["symbol"], x["method"], x["scheme"])
    y, origin = x, "original"
    if x["config_ok"] != "True" and key in fixmap:
        f = fixmap[key]
        # take the fix if it restored the configuration, or at least reduced the deviation
        if f["config_ok"] == "True" or abs(float(f["max_pop_dev"])) < abs(float(x["max_pop_dev"])):
            y, origin = f, "fix:" + f["source_file"]
    y = {k: v for k, v in y.items() if k != "source_file"}
    y["origin"] = origin
    y["usable"] = y["config_ok"] == "True" or abs(float(y["max_pop_dev"])) < 0.5
    y["flag"] = "" if y["config_ok"] == "True" else (
        "s/d hybridised (symmetry broken)" if abs(float(y["max_pop_dev"])) < 0.5 else "wrong configuration")
    selected.append(y)

out = ROOT / "results" / "atomic_radii_selected.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(selected[0]))
    w.writeheader(); w.writerows(selected)

by_set = defaultdict(list)
for y in selected:
    by_set[(y["method"], y["scheme"])].append(y)
for (m, sch), rows in by_set.items():
    path = ROOT / "data" / "radii" / f"density_{m.lower()}_{sch}.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set_id", "symbol", "variant", "kind", "value_A", "esd_A", "source", "table_ref", "note"])
        sid = f"density_{m.lower()}_{sch}"
        for y in sorted(rows, key=lambda r: int(r["Z"])):
            if not y["usable"]:
                continue
            note = y["flag"] + ("; " if y["flag"] else "") + y["origin"]
            src = f"this work, {m}/{sch}, free atom, spherical average"
            w.writerow([sid, y["symbol"], "", "v", y["v_rho_0.001"], "", src, "rho=0.001 e/bohr3", note])
            w.writerow([sid, y["symbol"], "", "r", y["r_outer"], "", src, "most probable radius, outermost occupied orbital", note])
            if y["c_count"] not in ("", "nan"):
                w.writerow([sid, y["symbol"], "", "c", y["c_count"], "", src, "sphere with chemical-core electrons", note])

summ = []
per_el = defaultdict(dict)
for y in selected:
    per_el[(int(y["Z"]), y["symbol"])][(y["method"], y["scheme"])] = y
for (z, sym), d in sorted(per_el.items()):
    ref = d.get(("PBE0", "anorcc"))
    usable = [float(v["v_rho_0.001"]) for v in d.values() if v["usable"]]
    summ.append({"Z": z, "symbol": sym,
                 "v_rho": ref["v_rho_0.001"] if ref and ref["usable"] else "",
                 "v_rho_spread": round(max(usable) - min(usable), 3) if len(usable) > 1 else "",
                 "n_variants": len(usable),
                 "r_outer": ref["r_outer"] if ref else "", "c_count": ref["c_count"] if ref else "",
                 "c_elf": ref["c_elf"] if ref else "",
                 "flags": "; ".join(sorted({f"{k[0]}/{k[1]}: {v['flag']}" for k, v in d.items() if v["flag"]}))})
with open(ROOT / "results" / "density_radii_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
print(f"{len(selected)} selected jobs; {sum(1 for y in selected if y['origin'].startswith('fix'))} from fixes; "
      f"{sum(1 for y in selected if not y['usable'])} unusable (need fix)")
print("radius sets:", ", ".join(f"density_{m.lower()}_{s}" for m, s in by_set))
