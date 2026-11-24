"""Merge all ion calculations, check anion binding, and write the ion radius sets.

Usage: python scripts/build_ion_radii.py base1 [base2 ...] --fixes dirA [dirB ...] --neutral-dirs d1 [d2 ...]
Later folders take precedence for the same job when they restore the configuration or reduce the
deviation. Electron affinities and ionization energies use the neutral PBE0/ANO-RCC energies found in
--neutral-dirs (the folder whose calculation is selected in results/atomic_radii_selected.csv).
Outputs: results/ion_radii_selected.csv, data/radii/density_{cations,anions}_pbe0_anorcc.csv
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii.orca_parse import read, final_energy, load_states  # noqa: E402

HARTREE_EV = 27.211386
ap = argparse.ArgumentParser()
ap.add_argument("bases", nargs="+")
ap.add_argument("--fixes", nargs="*", default=[])
ap.add_argument("--neutral-dirs", nargs="+", required=True)
args = ap.parse_args()


def extract(folder):
    out = ROOT / "results" / f"_tmp_{Path(folder).name}.csv"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "extract_atomic_radii.py"), folder, str(out)],
                   check=True, capture_output=True)
    rows = list(csv.DictReader(open(out)))
    out.unlink()
    for r in rows:
        r["folder"] = str(folder)
    return rows


best = {}
for folder in args.bases + args.fixes:
    for r in extract(folder):
        key = (r["symbol"], int(r["charge"]), r["method"], r["scheme"])
        old = best.get(key)
        better = old is None or (r["config_ok"] == "True" and old["config_ok"] != "True") or \
                 (old["config_ok"] != "True" and abs(float(r["max_pop_dev"])) < abs(float(old["max_pop_dev"])))
        if better:
            best[key] = r

# neutral energies: the calculation actually selected for each element
sel = {r["symbol"]: r for r in csv.DictReader(open(ROOT / "results" / "atomic_radii_selected.csv"))
       if r["method"] == "PBE0" and r["scheme"] == "anorcc"}
def neutral_energy(sym):
    job = sel[sym]["job"]
    for d in args.neutral_dirs:
        f = Path(d) / job / "job.out"
        if f.exists():
            e = final_energy(read(f))
            if e is not None:
                return e
    return None

states = load_states()
rows = []
for (sym, q, m, s), r in sorted(best.items(), key=lambda kv: (kv[0][1], int(kv[1]["Z"]))):
    e_ion = final_energy(read(Path(r["folder"]) / r["job"] / "job.out"))
    e_neu = neutral_energy(sym) if sym in sel else None
    r = dict(r)
    r["E_ion_Eh"] = e_ion
    if e_neu is not None and e_ion is not None:
        de = (e_ion - e_neu) * HARTREE_EV if q == 1 else (e_neu - e_ion) * HARTREE_EV
        r["IE_calc_eV" if q == 1 else "EA_calc_eV"] = round(de, 3)
    if q == -1:
        st = states[(sym, -1)]
        ea_exp = float(st["EA_eV"])
        r["EA_exp_eV"] = ea_exp
        r["EA_quality"] = st["EA_quality"]
        ea = r.get("EA_calc_eV")
        r["bound"] = ea is not None and ea > 0
        r["EA_within_0.3eV"] = ea is not None and abs(ea - ea_exp) <= 0.3
    r["usable"] = (r["config_ok"] == "True" or abs(float(r["max_pop_dev"])) < 0.7) and \
                  (q == 1 or r.get("bound", False))
    r["flag"] = "" if r["config_ok"] == "True" else "configuration deviates by %s e" % r["max_pop_dev"]
    rows.append(r)

fields = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
with open(ROOT / "results" / "ion_radii_selected.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)

for q, name in [(1, "cations"), (-1, "anions")]:
    with open(ROOT / "data" / "radii" / f"density_{name}_pbe0_anorcc.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set_id", "symbol", "variant", "kind", "value_A", "esd_A", "source", "table_ref", "note"])
        sid = f"density_{name}_pbe0_anorcc"
        for r in rows:
            if int(r["charge"]) != q or r["method"] != "PBE0" or r["scheme"] != "anorcc" or not r["usable"]:
                continue
            note = r["flag"]
            if q == -1 and not r.get("EA_within_0.3eV"):
                note = (note + "; " if note else "") + "computed EA differs from experiment by more than 0.3 eV"
            src = "this work, PBE0/ANO-RCC, free ion, spherical average"
            w.writerow([sid, r["symbol"], "", "v", r["v_rho_0.001"], "", src, "rho=0.001 e/bohr3", note])
            w.writerow([sid, r["symbol"], "", "r", r["r_outer"], "", src, "outermost occupied orbital", note])
            w.writerow([sid, r["symbol"], "", "c", r["c_count"], "", src, "sphere with chemical-core electrons", note])
print(f"{len(rows)} ions merged -> results/ion_radii_selected.csv")
