"""Orbital-based core radius c_orb and ELF shell populations for the whole primary set
(atoms and ions, PBE0/DKH2/ANO-RCC-TZP). Steps 2 and 3 of the revision (R1.8, R2.1, R2.2).

Usage (from the package root, on Agustina):
    python scripts/core_radii_all.py --root /fs/agustina/jorge/density_radii --out core_radii_all.csv

The species and their selected values come from results/atomic_radii_selected.csv and
results/ion_radii_selected.csv. For each species the script looks for every folder under
--root named like the job (e.g. 028_Ni_PBE0_anorcc) that holds a job.molden.input, and
picks the one whose v_rho(0.001) reproduces the selected value, so the solution used in the
paper is analysed whatever folder (original run, ΔSCF fix, retry) it came from.

Columns:
  c_orb      most probable radius (outermost maximum of r^2|phi|^2) of the outermost core
             orbital. Core orbitals: per spin, the n_core/2 lowest-energy occupied orbitals,
             with n_core the chemical core of the paper (orca_parse.core_electrons_ion).
  core_l_ok  whether the selected core orbitals have the expected l composition.
  c_count    sphere enclosing n_core electrons (previous definition).
  c_elf, N_inside_c_elf   outermost ELF minimum inside r_outer (when the number of minima
             is the expected one) and the number of electrons it encloses.
"""
import argparse
import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii.density import AtomDensity, BOHR_TO_ANG, LCHAR, _trapz  # noqa: E402
from density_radii.orca_parse import (CORES, expand, load_states, core_electrons_ion,  # noqa: E402
                                      max_principal_n)

NOBLE = [(0, ""), (2, "[He]"), (10, "[Ne]"), (18, "[Ar]"), (36, "[Kr]"), (54, "[Xe]"), (86, "[Rn]")]


def expected_core_l(n_core):
    """Per-spin orbital counts per l of a chemical core with n_core electrons."""
    base, tag = max((n, t) for n, t in NOBLE if n <= n_core)
    cnt = {l: 0 for l in "spdf"}
    if tag:
        for _n, l, _o in expand(CORES[tag]):
            cnt[l] += 2 * "spdf".index(l) + 1
    extra = {0: {}, 10: {"d": 5}, 14: {"f": 7}, 24: {"d": 5, "f": 7}}.get(n_core - base)
    if extra is None:
        return None
    for l, n in extra.items():
        cnt[l] += n
    return cnt


def orbital_table(a, rmin=0.01, rmax=20.0, npts=6000):
    r = np.linspace(rmin, rmax, npts)
    R = {k: a._radial(k, r) for k in range(len(a.atom.shells))}
    out = []
    for mo in a.atom.mos:
        if mo.occ <= 0:
            continue
        dens = np.zeros_like(r); lw = {}
        for (l, _m), members in a._groups.items():
            phi = sum(mo.coef[idx] * R[k] for k, idx in members)
            lw[l] = lw.get(l, 0.0) + float(_trapz(phi * phi * r * r, r))
            dens += phi * phi
        f = r * r * dens
        ldom = max(lw, key=lw.get)
        out.append({"spin": mo.spin, "energy": mo.energy, "occ": mo.occ, "l": LCHAR[ldom],
                    "w": lw[ldom] / (sum(lw.values()) or 1), "rmax": float(r[int(np.argmax(f))] * BOHR_TO_ANG)})
    return out


def electrons_inside(a, r_ang, npts=40000):
    r = np.linspace(0, r_ang / BOHR_TO_ANG, npts)
    return float(_trapz(4 * math.pi * r ** 2 * a.rho(r), r))


def index_folders(root, maxdepth=4):
    idx = {}
    root = Path(root)
    for dirpath, dirnames, filenames in os.walk(root):
        depth = len(Path(dirpath).relative_to(root).parts)
        if depth >= maxdepth:
            dirnames[:] = []
        if "job.molden.input" in filenames:
            idx.setdefault(Path(dirpath).name, []).append(Path(dirpath))
    return idx


def targets(levels):
    out = []
    for fn in ("atomic_radii_selected.csv", "ion_radii_selected.csv"):
        for r in csv.DictReader(open(ROOT / "results" / fn)):
            if (r["method"], r["scheme"]) in levels:
                out.append(r)
    return out


# Cations with a bare noble-gas configuration take the core of the isoelectronic noble-gas
# atom, so that e.g. Fr+ and Rn have the same core ([Xe]4f14 5d10 = 78 electrons).
ISO_CORE = {"[He]": 0, "[Ne]": 2, "[Ar]": 10, "[Kr]": 28, "[Xe]": 46, "[Rn]": 78}


def analyse(a, z, charge, conf):
    n_el = z - charge
    n_core = ISO_CORE[conf.strip()] if conf.strip() in ISO_CORE else core_electrons_ion(n_el, conf)
    orbs = orbital_table(a)
    r_outer = max(o["rmax"] for o in orbs)
    res = {"n_core": n_core, "r_outer": round(r_outer, 4)}
    if n_core:
        core = []
        for spin in ("alpha", "beta"):
            mine = sorted((o for o in orbs if o["spin"] == spin), key=lambda o: o["energy"])
            core += mine[: n_core // 2]
        top = max(core, key=lambda o: o["rmax"])
        got = {l: sum(1 for o in core if o["l"] == l) // 2 for l in "spdf"}
        exp = expected_core_l(n_core)
        res.update(c_orb=round(top["rmax"], 4), c_orb_l=top["l"], c_orb_l_weight=round(top["w"], 3),
                   core_l_ok=(exp == got) if exp else "unknown",
                   c_count=round(a.enclosing_radius(n_core), 4))
    else:
        res.update(c_orb="", c_orb_l="", c_orb_l_weight="", core_l_ok="", c_count="")
    mins = [m for m in a.elf_shell_minima() if m < r_outer]
    nexp = max_principal_n(conf) - 1
    ok = len(mins) == nexp
    res.update(elf_minima=" ".join(f"{m:.4f}" for m in mins), elf_minima_ok=ok,
               N_inside_minima=" ".join(f"{electrons_inside(a, m):.2f}" for m in mins),
               c_elf=round(mins[-1], 4) if (mins and ok) else "",
               N_inside_c_elf=round(electrons_inside(a, mins[-1]), 2) if (mins and ok) else "")
    if res["c_orb"] != "" and res["c_elf"] != "":
        res["order_corb_celf_router"] = res["c_orb"] < res["c_elf"] < res["r_outer"]
    else:
        res["order_corb_celf_router"] = ""
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="folder holding all the job folders (searched recursively)")
    ap.add_argument("--out", default="core_radii_all.csv")
    ap.add_argument("--only", nargs="*", help="restrict to these job names")
    ap.add_argument("--levels", nargs="+", default=["PBE0/anorcc"],
                    help="method/scheme pairs, e.g. PBE0/anorcc HF/anorcc PBE0/x2c HF/x2c")
    ap.add_argument("--tol", type=float, default=0.003, help="max |v_rho - selected| (angstrom) to accept a folder")
    args = ap.parse_args()

    states = load_states()
    idx = index_folders(args.root)
    print(f"{sum(len(v) for v in idx.values())} job folders with Molden files found under {args.root}")
    rows = []
    levels = [tuple(x.split("/")) for x in args.levels]
    for t in targets(levels):
        job = t["job"]
        if args.only and job not in args.only:
            continue
        z, sym, charge = int(t["Z"]), t["symbol"], int(t["charge"])
        base = {"job": job, "method": t["method"], "scheme": t["scheme"], "Z": z, "symbol": sym, "charge": charge,
                "configuration": t["configuration"], "v_rho_selected": t["v_rho_0.001"]}
        cands = idx.get(job, [])
        try:
            vsel = float(t["v_rho_0.001"])
        except ValueError:
            vsel = None
        best = None
        for c in cands:
            try:
                a = AtomDensity(c / "job.molden.input")
                v = a.isodensity_radius(0.001)
            except Exception as err:  # unreadable or truncated Molden file
                print(f"  {job}: cannot read {c} ({err})")
                continue
            dv = abs(v - vsel) if vsel is not None else 0.0
            if best is None or dv < best[0]:
                best = (dv, c, a, v)
        if best is None:
            rows.append({**base, "status": "no Molden file found"})
            print(f"{job:24s} no Molden file found")
            continue
        dv, folder, a, v = best
        status = "ok" if dv <= args.tol else f"v_rho mismatch {dv:.4f}"
        conf = states[(sym, charge)]["configuration"]
        res = analyse(a, z, charge, conf)
        rows.append({**base, "status": status, "folder": str(folder), "v_rho_recomputed": round(v, 4), **res})
        print(f"{job:24s} {status:8s} c_orb={res['c_orb']} c_count={res['c_count']} c_elf={res['c_elf']} "
              f"N_in={res['N_inside_c_elf']} r_outer={res['r_outer']}  [{folder}]", flush=True)

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    bad = [r["job"] for r in rows if r.get("status") != "ok"]
    print(f"\n{len(rows)} species -> {args.out};  {len(bad)} not ok: {' '.join(bad)}")


if __name__ == "__main__":
    main()
