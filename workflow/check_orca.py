"""Quality checks for ORCA single-point / optimization outputs.

Usage: python workflow/check_orca.py /path/to/free_atoms  -> writes checks.csv and flags.csv
Only jobs listed in joblist.txt are checked (if the file exists).
"""
import csv
import re
import sys
from pathlib import Path

MULT_RE = re.compile(r"^\*\s*xyz\s+-?\d+\s+(\d+)", re.I | re.M)


def check(jobdir: Path) -> dict:
    res = {"job": jobdir.name, "status": "ok", "energy_Eh": "", "S2": "",
           "S2_expected": "", "scf_iterations": "", "issues": ""}
    out, inp = jobdir / "job.out", jobdir / "job.inp"
    issues = []
    if not out.exists():
        st = jobdir / "status.txt"
        if st.exists():
            res.update(status="failed", issues="failed before ORCA: " + st.read_text().strip().replace("\n", " | "))
        else:
            res.update(status="not started", issues="no job.out and no status.txt (job not run yet, or old job script)")
        return res
    text = out.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        issues.append("abnormal termination")
    if re.search(r"SCF NOT CONVERGED", text):
        issues.append("SCF not converged")
    m = re.findall(r"SCF CONVERGED AFTER\s+(\d+)", text)
    if m:
        res["scf_iterations"] = m[-1]
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if m:
        res["energy_Eh"] = m[-1]
    else:
        issues.append("no final energy")
    for pat in ("There are no basis functions", "Unknown basis", "INPUT ERROR"):
        if pat.lower() in text.lower():
            issues.append(f"input/basis problem: {pat}")
    m = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)", text)
    mm = MULT_RE.search(inp.read_text()) if inp.exists() else None
    if mm:
        s = (int(mm.group(1)) - 1) / 2
        res["S2_expected"] = f"{s * (s + 1):.4f}"
        if m:
            res["S2"] = m[-1]
            if abs(float(m[-1]) - s * (s + 1)) > 0.05:
                issues.append("spin contamination > 0.05")
    if not (jobdir / "job.molden.input").exists():
        issues.append("no molden file")
    if issues:
        res["status"] = "flag"
        res["issues"] = "; ".join(issues)
    return res


def main(root):
    root = Path(root)
    jl = root / "joblist.txt"
    if jl.exists():
        names = [n.strip() for n in jl.read_text().splitlines() if n.strip()]
        dirs = [root / n for n in names]
        extra = sorted(d.name for d in root.iterdir() if (d / "job.inp").exists() and d.name not in names)
        if extra:
            print(f"note: {len(extra)} job folders not in joblist.txt were ignored (e.g. {extra[0]})")
    else:
        dirs = sorted(d for d in root.iterdir() if (d / "job.inp").exists())
    rows = [check(d) for d in dirs]
    fields = list(rows[0]) if rows else []
    for name, subset in (("checks.csv", rows), ("flags.csv", [r for r in rows if r["status"] != "ok"])):
        with open(root / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(subset)
    n_flag = sum(r["status"] != "ok" for r in rows)
    print(f"{len(rows)} jobs checked, {n_flag} flagged -> {root / 'flags.csv'}")


if __name__ == "__main__":
    main(sys.argv[1])
