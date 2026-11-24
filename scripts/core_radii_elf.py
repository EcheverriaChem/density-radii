"""Orbital-based core radius c_orb and ELF shell populations (Steps 2 and 3, R1.8, R2.1, R2.2).

Usage: python scripts/core_radii_elf.py /path/to/jobs [out.csv]
       (job folders with job.molden.input and job.out, named like 026_Fe_PBE0_anorcc)

Per job:
  c_orb      outermost maximum of r^2|phi|^2 among the core spin orbitals, i.e. the
             most probable radius of the outermost core orbital (same construction as r_outer).
             Core orbitals: per spin, the lowest-energy occupied orbitals of each l, in the
             number required by the chemical core (noble-gas core, plus f14 in the d block,
             plus d10 and f14 in the p block, as in orca_parse.core_electrons).
  c_count    radius of the sphere enclosing the core electrons (previous definition).
  elf_minima ELF shell minima inside r_outer and N_inside, the number of electrons enclosed
             by each minimum (4 pi int_0^r rho r'^2 dr').
"""
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii.density import AtomDensity, BOHR_TO_ANG, LCHAR, _trapz  # noqa: E402
from density_radii.orca_parse import (CORES, expand, valence_shells, core_electrons,  # noqa: E402
                                      load_states, parse_job_name, D_BLOCK)


def core_orbital_counts(Z, config):
    """Number of core spatial orbitals per l (per spin), e.g. Fe -> {'s': 3, 'p': 6}."""
    core_cfg = ""
    for core in CORES:
        if config.startswith(core):
            core_cfg = core
    counts = {l: 0 for l in "spdfg"}
    for _n, l, _occ in (expand(core_cfg + " ") if core_cfg else []):
        counts[l] += 2 * "spdfg".index(l) + 1
    for tok in valence_shells(config):
        m = re.fullmatch(r"(\d)([spdf])(\d+)", tok)
        l, n = m.group(2), int(m.group(3))
        if (l == "f" and n == 14) or (l == "d" and n == 10 and Z not in D_BLOCK):
            counts[l] += 2 * "spdfg".index(l) + 1
    return counts


def orbital_table(a, rmin=0.01, rmax=20.0, npts=8000):
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
        out.append({"spin": mo.spin, "energy": mo.energy, "l": LCHAR[ldom],
                    "w": lw[ldom] / (sum(lw.values()) or 1), "rmax": float(r[int(np.argmax(f))] * BOHR_TO_ANG)})
    return out


def electrons_inside(a, r_ang, npts=40000):
    r = np.linspace(0, r_ang / BOHR_TO_ANG, npts)
    f = 4 * math.pi * r ** 2 * a.rho(r)
    return float(_trapz(f, r))


def main(src, dst):
    states = load_states()
    rows = []
    for d in sorted(p for p in Path(src).iterdir() if (p / "job.molden.input").exists()):
        z, sym, charge, method, scheme = parse_job_name(d.name)
        conf = states[(sym, charge)]["configuration"]
        a = AtomDensity(d / "job.molden.input")
        orbs = orbital_table(a)
        r_outer = max(o["rmax"] for o in orbs)
        need = core_orbital_counts(z, conf)
        core = []
        for spin in ("alpha", "beta"):
            mine = sorted((o for o in orbs if o["spin"] == spin), key=lambda o: o["energy"])
            for l, n in need.items():
                core += [o for o in mine if o["l"] == l][:n]
        c_orb_o = max(core, key=lambda o: o["rmax"])
        n_core = core_electrons(z, conf)
        mins = [m for m in a.elf_shell_minima() if m < r_outer]
        n_in = [electrons_inside(a, m) for m in mins]
        rows.append({
            "job": d.name, "Z": z, "symbol": sym, "configuration": conf,
            "n_core": n_core, "core_orbitals_per_spin": " ".join(f"{l}{n}" for l, n in need.items() if n),
            "c_orb": round(c_orb_o["rmax"], 4), "c_orb_orbital_l": c_orb_o["l"],
            "c_orb_l_weight": round(c_orb_o["w"], 3),
            "c_count": round(a.enclosing_radius(n_core), 4),
            "c_elf": round(mins[-1], 4) if mins else "",
            "N_inside_c_elf": round(n_in[-1], 2) if mins else "",
            "elf_minima": " ".join(f"{m:.4f}" for m in mins),
            "N_inside_minima": " ".join(f"{n:.2f}" for n in n_in),
            "r_outer": round(r_outer, 4),
            "order_corb_celf_router": (c_orb_o["rmax"] < mins[-1] < r_outer) if mins else "",
        })
        r = rows[-1]
        print(f"{d.name:24s} c_orb={r['c_orb']:.3f}({r['c_orb_orbital_l']}) c_count={r['c_count']:.3f} "
              f"c_elf={r['c_elf']} N_in={r['N_inside_c_elf']} n_core={n_core} r_outer={r['r_outer']:.3f}")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} jobs -> {dst}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else Path(sys.argv[1]) / "core_radii_elf.csv")
