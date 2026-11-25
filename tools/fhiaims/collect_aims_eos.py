"""Collect volume, energy and pressure from the FHI-aims scans (no pandas needed).
Usage:  python3 collect_aims_eos.py eos_aims > eos_pv.csv
"""
import re, sys
from pathlib import Path
NUM = r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)"
print("element,structure,i,V_atom_A3,E_eV_atom,P_GPa,converged")
for d in sorted(Path(sys.argv[1]).iterdir()):
    o = d / "aims.out"
    if not o.is_file():
        continue
    t = o.read_text(errors="replace")
    el, st, i = d.name.split("_")
    nat = 2 if st == "hcp" else 1
    v = re.findall(r"Unit cell volume\s*:\s*" + NUM, t)
    e = re.findall(r"Total energy corrected\s*:\s*(?:" + NUM + r"\s*Ha\s*)?" + NUM + r"\s*eV", t)
    p = re.findall(r"Pressure\s*:\s*" + NUM + r"\s*\[eV/A", t)
    ok = "Have a nice day" in t
    if not (v and e):
        print(f"# {d.name}: volume {'found' if v else 'MISSING'}, energy {'found' if e else 'MISSING'}, finished={ok}", file=sys.stderr)
        continue
    pg = f"{float(p[-1]) * 160.21766208:.3f}" if p else ""
    print(f"{el},{st},{i},{float(v[-1])/nat:.5f},{float(e[-1][1])/nat:.6f},{pg},{ok}")
