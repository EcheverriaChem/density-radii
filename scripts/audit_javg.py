"""Audit the J-averaged alternative-configuration jobs against data/atomic_states_javg.csv.

Usage: python scripts/audit_javg.py /path/to/javg_alt
Same criterion as audit_configurations.py (Loewdin s/p/d/f populations within 0.15 e).
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from density_radii.orca_parse import read, reduced_populations, l_totals, s2, final_energy, parse_job_name  # noqa: E402

TOL = 0.15
states = {(r["symbol"], int(r["charge"])): r
          for r in csv.DictReader(open(ROOT / "data" / "atomic_states_javg.csv"))}
root = Path(sys.argv[1])
for d in sorted(p for p in root.iterdir() if (p / "job.out").exists()):
    z, sym, charge, method, scheme = parse_job_name(d.name)
    st = states[(sym, charge)]
    txt = read(d / "job.out")
    try:
        pop = reduced_populations(txt)
    except ValueError:
        print(f"{d.name:28s} no population block (SCF not converged or job unfinished)")
        continue
    exp = l_totals(st["configuration"])
    dev = {l: pop["charge"][l] - exp[l] for l in "spdf"}
    worst = max(dev.values(), key=abs)
    S = (int(st["multiplicity"]) - 1) / 2
    flag = "OK" if abs(worst) <= TOL else "WRONG CONFIGURATION"
    print(f"{d.name:28s} {st['configuration']:20s} E={final_energy(txt)}  <S2>={s2(txt)} (exp {S*(S+1):.2f})  "
          + " ".join(f"d{l}={v:+.2f}" for l, v in dev.items()) + f"  {flag}")
