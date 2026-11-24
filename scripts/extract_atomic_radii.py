"""Density-derived atomic radii from ORCA free-atom Molden files (milestone M2).

Usage: python scripts/extract_atomic_radii.py /path/to/free_atoms [results.csv]

Per job: electron-count check, v_rho at 0.0005/0.001/0.002 e/bohr^3, most probable
radius of the HOMO (r_homo) and of the outermost occupied orbital (r_outer), ELF shell minima and core radius c_elf (outermost
shell boundary inside r_outer; empty when ELF misses a shell), core radius c_count
(sphere containing the chemical core electrons, see orca_parse.core_electrons), plus the configuration audit.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii.density import AtomDensity  # noqa: E402
from density_radii.orca_parse import (read, reduced_populations, l_totals, s2, core_electrons_ion,  # noqa: E402
                                 load_states, parse_job_name, max_principal_n)

PERIOD = [0] + [1] * 2 + [2] * 8 + [3] * 8 + [4] * 18 + [5] * 18 + [6] * 32 + [7] * 32
TOL = 0.15


def main(src, dst):
    src = Path(src)
    states = load_states()
    rows = []
    for d in sorted(p for p in src.iterdir() if (p / "job.molden.input").exists()):
        z, sym, charge, method, scheme = parse_job_name(d.name)
        st = states[(sym, charge)]
        txt = read(d / "job.out")
        try:
            pop = reduced_populations(txt)
        except ValueError:
            print(f"{d.name:24s} SKIPPED: no population block (SCF failed?)")
            continue
        a = AtomDensity(d / "job.molden.input")
        exp = l_totals(st["configuration"])
        dev = max((pop["charge"][l] - exp[l] for l in "spdf"), key=abs)
        r_homo, l_homo, w_homo, _, e_homo = a.homo_radius()
        orbs = a.orbital_radii()
        r_outer, l_outer, e_outer, _ = max(orbs, key=lambda o: o[0])
        mins = [m for m in a.elf_shell_minima() if m < r_outer]
        conf = st["configuration"]
        n_el = z - charge
        n_core = core_electrons_ion(n_el, conf)      # isoelectronic rule for cations with no valence shell
        c_count = a.enclosing_radius(n_core) if n_core else 0.0
        nexp = max_principal_n(st["configuration"]) - 1   # occupied shells, not the period
        rows.append({
            "job": d.name, "Z": z, "symbol": sym, "charge": charge, "method": method, "scheme": scheme,
            "configuration": st["configuration"],
            "config_ok": abs(dev) <= TOL, "max_pop_dev": round(dev, 3),
            "S2": s2(txt), "norm_convention": a.convention,
            "n_electrons": round(a.electron_count(), 5),
            "v_rho_0.0005": round(a.isodensity_radius(0.0005), 4),
            "v_rho_0.001": round(a.isodensity_radius(0.001), 4),
            "v_rho_0.002": round(a.isodensity_radius(0.002), 4),
            "r_homo": round(r_homo, 4), "homo_l": l_homo, "homo_l_weight": round(w_homo, 3),
            "homo_energy_Eh": round(e_homo, 5),
            "r_outer": round(r_outer, 4), "outer_l": l_outer, "outer_energy_Eh": round(e_outer, 5),
            "n_core": n_core, "c_count": round(c_count, 4),
            "elf_minima": " ".join(f"{m:.4f}" for m in mins),
            "elf_minima_ok": len(mins) == nexp,
            "c_elf": round(mins[-1], 4) if mins and len(mins) == nexp else "",
        })
        print(f"{d.name:24s} N={rows[-1]['n_electrons']:9.4f} v={rows[-1]['v_rho_0.001']:.3f} "
              f"rO={r_outer:.3f}{l_outer} cE={rows[-1]['c_elf']} cN={c_count:.3f} cfg={'ok' if rows[-1]['config_ok'] else 'WRONG'}")
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} jobs -> {dst}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else Path(sys.argv[1]) / "atomic_radii_raw.csv")
