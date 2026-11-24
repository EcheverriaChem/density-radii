"""Compare the state-averaged CASSCF densities with the single-determinant (SD) ones (Step 11, R1.4).

Usage (from the package root, on Agustina):
    python scripts/analyse_casscf.py --casscf /fs/agustina/jorge/density_radii/casscf_step11 \
        --sd-root /fs/agustina/jorge/density_radii --out casscf_comparison_step11.csv

For every folder named ZZZ_Sym_CASSCF_scheme under --casscf (searched recursively, so the old
La/Th/Ce folders are picked up too if they sit below it) the script reads the Molden file of the
state-averaged natural orbitals and computes
  v_rho      radius of the 0.001 e/bohr^3 isodensity sphere
  r_outer    outermost maximum of r^2|phi|^2 among the natural orbitals with occupation > --occ-min
  c_orb      the same for the outermost core orbital (core = the n_core/2 doubly occupied
             orbitals lowest in energy, n_core as in the paper)
  c_count    sphere enclosing n_core electrons
It then looks for the SD folder of the same atom (ZZZ_Sym_PBE0_anorcc, and HF if present) under
--sd-root, picking the one that reproduces the selected v_rho of the paper, and prints the
differences. Root energies and natural-orbital occupations of the active space are read from
job.out so that the degeneracy of the averaged roots can be checked (a spherical average needs a
spread close to zero).
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from density_radii.density import AtomDensity  # noqa: E402
from density_radii.orca_parse import load_states, core_electrons_ion  # noqa: E402
from core_radii_all import orbital_table, index_folders, targets  # noqa: E402

HARTREE_TO_CM = 219474.63


def roots(text):
    """Energies (Eh) of the averaged roots in the final CAS-SCF STATES block."""
    i = text.rfind("CAS-SCF STATES FOR BLOCK")
    if i < 0:
        return []
    m = re.search(r"NROOTS=\s*(\d+)", text[i:i + 200])
    n = int(m.group(1)) if m else 1
    return [float(x) for x in re.findall(r"ROOT\s+\d+:\s+E=\s+(-?\d+\.\d+)", text[i:])[:n]]


def active_occupations(text):
    m = re.findall(r"N\[\s*\d+\]\s*=\s*(\d\.\d+)", text)
    return m


def radii(a, z, conf, occ_min):
    orbs = [o for o in orbital_table(a) if o["occ"] > occ_min]
    n_core = core_electrons_ion(z, conf)
    res = {"v_rho": round(a.isodensity_radius(0.001), 4),
           "r_outer": round(max(o["rmax"] for o in orbs), 4),
           "r_outer_l": max(orbs, key=lambda o: o["rmax"])["l"], "n_core": n_core}
    if n_core:
        core = []
        for spin in ("alpha", "beta"):
            mine = sorted((o for o in orbs if o["spin"] == spin), key=lambda o: o["energy"])
            if not mine:
                continue
            # n_core/2 orbitals per spin (unrestricted SD) or n_core/2 doubly occupied ones (restricted CASSCF)
            core += mine[: n_core // 2]
        top = max(core, key=lambda o: o["rmax"])
        res.update(c_orb=round(top["rmax"], 4), c_orb_l=top["l"], c_count=round(a.enclosing_radius(n_core), 4))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--casscf", required=True)
    ap.add_argument("--sd-root", required=True)
    ap.add_argument("--out", default="casscf_comparison_step11.csv")
    ap.add_argument("--occ-min", type=float, default=0.05)
    args = ap.parse_args()

    states = load_states()
    sel = {(t["symbol"], t["method"], t["scheme"]): t for t in targets([("PBE0", "anorcc"), ("HF", "anorcc")])
           if int(t["charge"]) == 0}
    sd_idx = index_folders(args.sd_root)
    rows = []
    for dirpath, _d, files in sorted(os.walk(args.casscf)):
        d = Path(dirpath)
        if "_CASSCF_" not in d.name or "job.molden.input" not in files:
            continue
        z, sym = int(d.name.split("_")[0]), d.name.split("_")[1]
        conf = states[(sym, 0)]["configuration"]
        text = (d / "job.out").read_text(errors="replace") if (d / "job.out").exists() else ""
        e = roots(text)
        spread = (max(e) - min(e)) * HARTREE_TO_CM if len(e) > 1 else 0.0
        cas = radii(AtomDensity(d / "job.molden.input"), z, conf, args.occ_min)
        row = {"symbol": sym, "Z": z, "configuration": conf, "nroots": len(e),
               "root_spread_cm-1": round(spread, 1), "active_NO_occ": " ".join(active_occupations(text)[-12:]),
               **{f"{k}_CAS": v for k, v in cas.items()}}
        line = (f"{sym:3s} roots={len(e)} spread={spread:8.1f} cm-1  CAS: v={cas['v_rho']:.4f} "
                f"r_out={cas['r_outer']:.4f}({cas['r_outer_l']}) c_orb={cas.get('c_orb', '')} c_cnt={cas.get('c_count', '')}")
        for method in ("PBE0", "HF"):
            t = sel.get((sym, method, "anorcc"))
            if not t:
                continue
            best = None
            for c in sd_idx.get(f"{z:03d}_{sym}_{method}_anorcc", []):
                try:
                    a = AtomDensity(c / "job.molden.input")
                    dv = abs(a.isodensity_radius(0.001) - float(t["v_rho_0.001"]))
                except Exception:
                    continue
                if best is None or dv < best[0]:
                    best = (dv, a)
            if best is None:
                continue
            sd = radii(best[1], z, conf, args.occ_min)
            for k in ("v_rho", "r_outer", "c_orb", "c_count"):
                if k in sd and k in cas:
                    row[f"{k}_{method}"] = sd[k]
                    row[f"d_{k}_{method}"] = round(sd[k] - cas[k], 4)
            line += (f"\n      {method:4s} SD - CAS: dv={sd['v_rho'] - cas['v_rho']:+.4f} "
                     f"dr_out={sd['r_outer'] - cas['r_outer']:+.4f} "
                     + (f"dc_orb={sd['c_orb'] - cas['c_orb']:+.4f} dc_cnt={sd['c_count'] - cas['c_count']:+.4f}"
                        if 'c_orb' in sd and 'c_orb' in cas else ""))
        print(line, flush=True)
        rows.append(row)
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)} CASSCF folders -> {args.out}")


if __name__ == "__main__":
    main()
