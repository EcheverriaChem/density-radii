"""Compare Loewdin s/p/d/f populations of each free-atom SCF with the approved configuration.

Usage: python scripts/audit_configurations.py /path/to/free_atoms
Writes <dir>/configuration_audit.csv; prints jobs whose populations deviate by > tol.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from density_radii.orca_parse import read, reduced_populations, l_totals, s2, final_energy, load_states, parse_job_name  # noqa: E402

TOL = 0.15
root = Path(sys.argv[1])
states = load_states()
rows = []
for d in sorted(p for p in root.iterdir() if (p / "job.out").exists()):
    z, sym, charge, method, scheme = parse_job_name(d.name)
    st = states[(sym, charge)]
    exp = l_totals(st["configuration"])
    txt = read(d / "job.out")
    try:
        pop = reduced_populations(txt)
    except ValueError:
        print(f"  {d.name:22s} SKIPPED: no population block (SCF failed?)")
        continue
    dev = {l: round(pop["charge"][l] - exp[l], 3) for l in "spdf"}
    mult = int(st["multiplicity"]); S = (mult - 1) / 2
    s2v = s2(txt)
    worst = max(dev.values(), key=abs)
    rows.append({"job": d.name, "symbol": sym, "method": method, "scheme": scheme,
                 "config": st["configuration"], "E_Eh": final_energy(txt),
                 "S2": s2v, "S2_exp": round(S * (S + 1), 4),
                 "dev_s": dev["s"], "dev_p": dev["p"], "dev_d": dev["d"], "dev_f": dev["f"],
                 "spin_s": pop["spin"]["s"], "spin_p": pop["spin"]["p"],
                 "spin_d": pop["spin"]["d"], "spin_f": pop["spin"]["f"],
                 "max_abs_dev": abs(worst),
                 "ok": abs(worst) <= TOL})
with open(root / "configuration_audit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
bad = [r for r in rows if not r["ok"]]
print(f"{len(rows)} jobs audited, {len(bad)} with |population deviation| > {TOL} e")
for r in bad:
    print(f"  {r['job']:22s} {r['config']:22s} ds={r['dev_s']:+.2f} dp={r['dev_p']:+.2f} dd={r['dev_d']:+.2f} df={r['dev_f']:+.2f}  S2={r['S2']}/{r['S2_exp']}")
