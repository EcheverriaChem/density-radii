#!/usr/bin/env python3
"""Standalone audit of the bonding-configuration jobs (Step 8b).

Needs only Python 3 (no density_radii package). Run it from anywhere:
    python3 audit_bond.py /fs/agustina/jorge/density_radii/bond_configs

For every job folder with a job.out it compares the Loewdin reduced s/p/d/f
populations with the target configuration, using the same criterion as
audit_configurations.py (all deviations within 0.15 e), and also prints the
per-spin deviations, which show which shell holds an extra electron.
"""
import re
import sys
from pathlib import Path

TARGETS = {   # symbol: (configuration, multiplicity)
    "Pd": ("[Kr]4d9 5s1", 3),
    "Ce": ("[Xe]5d2 6s2", 3),
    "Pr": ("[Xe]4f2 5d1 6s2", 4),
    "Nd": ("[Xe]4f3 5d1 6s2", 5),
}
TOL = 0.15
L = "spdfg"
CORES = {
    "[He]": "1s2",
    "[Ne]": "1s2 2s2 2p6",
    "[Ar]": "1s2 2s2 2p6 3s2 3p6",
    "[Kr]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
    "[Xe]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6",
    "[Rn]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6",
}


def expand(config):
    for core, full in CORES.items():
        if config.startswith(core):
            config = full + " " + config[len(core):]
    out = []
    for tok in config.split():
        m = re.fullmatch(r"(\d)([spdf])(\d+)", tok)
        out.append((int(m.group(1)), m.group(2), int(m.group(3))))
    return out


def l_totals(config):
    tot = {l: 0 for l in L}
    for _, l, n in expand(config):
        tot[l] += n
    return tot


def expected_spin(config):
    out = {"alpha": {l: 0.0 for l in L}, "beta": {l: 0.0 for l in L}}
    for _n, l, occ in expand(config):
        a = min(occ, 2 * L.index(l) + 1)
        out["alpha"][l] += a
        out["beta"][l] += occ - a
    return out


def reduced_populations(text):
    head = "LOEWDIN REDUCED ORBITAL CHARGES AND SPIN POPULATIONS"
    idx = text.rfind(head)
    if idx < 0:
        raise ValueError
    res, cur = {"charge": {}, "spin": {}}, None
    for line in text[idx:].split("\n")[2:200]:
        s = line.strip()
        if s in ("CHARGE", "SPIN"):
            cur = s.lower()
            continue
        if cur is None:
            continue
        for m in re.finditer(r"\b([spdfg])\s*:\s*(-?\d+\.\d+)", line):
            res[cur][m.group(1)] = float(m.group(2))
        if cur == "spin" and (s.startswith("---") or s.startswith("*")):
            break
    for part in res.values():
        for l in L:
            part.setdefault(l, 0.0)
    return res


def last(pattern, text):
    m = re.findall(pattern, text)
    return float(m[-1]) if m else None


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    root = Path(sys.argv[1])
    found = False
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        sym = d.name.split("_")[1].rstrip("+-") if "_" in d.name else ""
        if sym not in TARGETS:
            continue
        found = True
        conf, mult = TARGETS[sym]
        out = d / "job.out"
        if not out.exists():
            print(f"{d.name:28s} no job.out yet (queued or running?)")
            continue
        text = out.read_text(errors="replace")
        conv = "SCF CONVERGED" in text or "SUCCESS" in text
        try:
            pop = reduced_populations(text)
        except ValueError:
            print(f"{d.name:28s} no population block (SCF not converged or job unfinished)")
            continue
        exp = l_totals(conf)
        dev = {l: pop["charge"][l] - exp[l] for l in "spdf"}
        worst = max(dev.values(), key=abs)
        es = expected_spin(conf)
        obs = {"alpha": {l: (pop["charge"][l] + pop["spin"][l]) / 2 for l in L},
               "beta": {l: (pop["charge"][l] - pop["spin"][l]) / 2 for l in L}}
        parts = []
        for sp in ("alpha", "beta"):
            d_ = [f"{l}{obs[sp][l] - es[sp][l]:+.2f}" for l in "spdf" if abs(obs[sp][l] - es[sp][l]) > 0.05]
            parts.append(f"{sp}:" + (",".join(d_) if d_ else "ok"))
        spin_dev = " ".join(parts)
        S = (mult - 1) / 2
        e = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
        s2 = last(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)", text)
        flag = "OK" if abs(worst) <= TOL else "WRONG CONFIGURATION"
        print(f"{d.name:28s} target {conf:18s} E={e}  <S2>={s2} (exp {S*(S+1):.2f})  "
              + " ".join(f"d{l}={v:+.2f}" for l, v in dev.items())
              + f"  [{spin_dev}]  {flag}" + ("" if conv else "  (convergence message not found)"))
    if not found:
        print(f"no Pd/Ce/Pr/Nd job folders found in {root}")


if __name__ == "__main__":
    main()
