"""Write the orbital-based core radius into the radius tables (version 1.1).

In version 1.0 the kind "c" of data/radii/density_*.csv was the radius of the sphere enclosing
the chemical-core electrons. From version 1.1 "c" is the orbital-based core radius c_orb (most
probable radius of the outermost core orbital), and the former value is kept as kind "c_count".
The ELF core radius (outermost ELF minimum inside r_outer) is added as kind "c_elf" where defined.

Input: results/core_radii_pbe0_anorcc.csv (atoms and ions, primary level) and
results/core_radii_other_levels.csv (atoms, HF/anorcc, HF/x2c, PBE0/x2c), both written by
scripts/core_radii_all.py.

Usage (from the package root):  python scripts/update_core_radii_tables.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RADII = ROOT / "data" / "radii"
TABLES = {  # file: (core-radii csv, method, scheme, charge)
    "density_pbe0_anorcc.csv": ("core_radii_pbe0_anorcc.csv", None, None, 0),
    "density_cations_pbe0_anorcc.csv": ("core_radii_pbe0_anorcc.csv", None, None, 1),
    "density_anions_pbe0_anorcc.csv": ("core_radii_pbe0_anorcc.csv", None, None, -1),
    "density_hf_anorcc.csv": ("core_radii_other_levels.csv", "HF", "anorcc", 0),
    "density_hf_x2c.csv": ("core_radii_other_levels.csv", "HF", "x2c", 0),
    "density_pbe0_x2c.csv": ("core_radii_other_levels.csv", "PBE0", "x2c", 0),
}
REF_ORB = "most probable radius, outermost core orbital (c_orb)"
REF_COUNT = "sphere with chemical-core electrons (c_count, kind c in version 1.0)"
REF_ELF = "outermost ELF minimum inside r_outer (c_elf)"


def fmt(x):
    return "" if x in ("", None) else f"{float(x):.4f}"


def main():
    for fn, (src, method, scheme, charge) in TABLES.items():
        core = {}
        for r in csv.DictReader(open(ROOT / "results" / src)):
            if int(r["charge"]) != charge:
                continue
            if method and (r["method"], r["scheme"]) != (method, scheme):
                continue
            core[r["symbol"]] = r
        rows = list(csv.DictReader(open(RADII / fn)))
        fields = list(rows[0].keys())
        out, missing = [], []
        for row in rows:
            if row["kind"] in ("c_count", "c_elf"):
                continue                      # rerun-safe
            if row["kind"] != "c":
                out.append(row)
                continue
            c = core.get(row["symbol"])
            if c is None:
                missing.append(row["symbol"])
                out.append(row)
                continue
            n_core = int(c["n_core"]) if c["n_core"] else 0
            orb = dict(row, value_A=fmt(c["c_orb"]) if n_core else "0.0", table_ref=REF_ORB)
            cnt = dict(row, kind="c_count", value_A=fmt(c["c_count"]) if n_core else "0.0", table_ref=REF_COUNT)
            out += [orb, cnt]
            if c.get("c_elf"):
                out.append(dict(row, kind="c_elf", value_A=fmt(c["c_elf"]), table_ref=REF_ELF))
        with open(RADII / fn, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(out)
        print(f"{fn:36s} {sum(r['kind'] == 'c' for r in out)} c_orb, "
              f"{sum(r['kind'] == 'c_elf' for r in out)} c_elf" + (f", no core data for {missing}" if missing else ""))


if __name__ == "__main__":
    main()
