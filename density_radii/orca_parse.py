"""Parsers for ORCA 6 outputs of free-atom calculations."""
from __future__ import annotations

import re
from pathlib import Path

CORES = {
    "[He]": "1s2",
    "[Ne]": "1s2 2s2 2p6",
    "[Ar]": "1s2 2s2 2p6 3s2 3p6",
    "[Kr]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
    "[Xe]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6",
    "[Rn]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6",
}
L = "spdfg"


def expand(config: str) -> list[tuple[int, str, int]]:
    """'[Ar]3d5 4s1' -> [(1,'s',2), ..., (3,'d',5), (4,'s',1)]"""
    for core, full in CORES.items():
        if config.startswith(core):
            config = full + " " + config[len(core):]
    out = []
    for tok in config.split():
        m = re.fullmatch(r"(\d)([spdf])(\d+)", tok)
        if not m:
            raise ValueError(f"cannot parse '{tok}' in configuration")
        out.append((int(m.group(1)), m.group(2), int(m.group(3))))
    return out


def l_totals(config: str) -> dict[str, int]:
    tot = {l: 0 for l in L}
    for _, l, n in expand(config):
        tot[l] += n
    return tot


def valence_shells(config: str) -> list[str]:
    """Shells written explicitly after the noble-gas core (e.g. ['3d5', '4s1'])."""
    for core in CORES:
        if config.startswith(core):
            config = config[len(core):]
    return config.split()


def reduced_populations(out_text: str, kind: str = "LOEWDIN") -> dict:
    """Last '<kind> REDUCED ORBITAL CHARGES AND SPIN POPULATIONS' block of a one-atom output.

    Returns {'charge': {l: pop}, 'spin': {l: pop}}.
    """
    head = f"{kind} REDUCED ORBITAL CHARGES AND SPIN POPULATIONS"
    idx = out_text.rfind(head)
    if idx < 0:
        raise ValueError(f"no {head} block")
    block = out_text[idx:].split("\n")
    res, cur = {"charge": {}, "spin": {}}, None
    for line in block[2:200]:
        s = line.strip()
        if s == "CHARGE":
            cur = "charge"
            continue
        if s == "SPIN":
            cur = "spin"
            continue
        if cur is None:
            continue
        for m in re.finditer(r"\b([spdfg])\s*:\s*(-?\d+\.\d+)", line):
            res[cur][m.group(1)] = float(m.group(2))
        if cur == "spin" and s.startswith("---") or (cur == "spin" and s.startswith("*")):
            break
    for part in res.values():
        for l in L:
            part.setdefault(l, 0.0)
    return res


def final_energy(out_text: str) -> float | None:
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out_text)
    return float(m[-1]) if m else None


def s2(out_text: str) -> float | None:
    m = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)", out_text)
    return float(m[-1]) if m else None


def read(path: Path) -> str:
    return Path(path).read_text(errors="replace")


D_BLOCK = set(range(21, 31)) | set(range(39, 49)) | {57} | set(range(72, 81)) | {89}


def core_electrons(Z: int, config: str) -> int:
    """Chemical core: noble-gas core plus filled f14 shells (d block) or filled
    d10 and f14 shells (p block) that are written explicitly in the configuration."""
    val = 0
    for tok in valence_shells(config):
        m = re.fullmatch(r"(\d)([spdf])(\d+)", tok)
        l, n = m.group(2), int(m.group(3))
        if l == "f" and n == 14:
            continue
        if l == "d" and n == 10 and Z not in D_BLOCK:
            continue
        val += n
    return Z - val


def expected_spin_populations(config: str) -> dict:
    """Per-spin l populations implied by a configuration with maximum multiplicity
    (Hund's rule within each open shell). Returns {'alpha': {l: n}, 'beta': {l: n}}."""
    out = {"alpha": {l: 0.0 for l in L}, "beta": {l: 0.0 for l in L}}
    for _n, l, occ in expand(config):
        cap = 2 * L.index(l) + 1
        a = min(occ, cap)
        out["alpha"][l] += a
        out["beta"][l] += occ - a
    return out


def observed_spin_populations(pop: dict) -> dict:
    """Split Loewdin reduced charges/spins into alpha and beta l populations."""
    return {"alpha": {l: (pop["charge"][l] + pop["spin"][l]) / 2 for l in L},
            "beta": {l: (pop["charge"][l] - pop["spin"][l]) / 2 for l in L}}


import csv as _csv
from pathlib import Path as _Path
_DATA = _Path(__file__).resolve().parent.parent / "data"


def load_states(include_ions: bool = True) -> dict:
    """{(symbol, charge): row} from atomic_states.csv and, if present, atomic_states_ions.csv."""
    out = {}
    for fn in ["atomic_states.csv"] + (["atomic_states_ions.csv"] if include_ions else []):
        path = _DATA / fn
        if not path.exists():
            continue
        for r in _csv.DictReader(open(path)):
            out[(r["symbol"], int(r.get("charge", 0) or 0))] = r
    return out


def parse_job_name(name: str):
    """'011_Na+_PBE0_anorcc' -> (11, 'Na', +1, 'PBE0', 'anorcc'); neutral names have no sign."""
    z, sym, method, scheme = name.split("_")[:4]
    charge = 1 if sym.endswith("+") else -1 if sym.endswith("-") else 0
    return int(z), sym.rstrip("+-"), charge, method, scheme


def max_principal_n(config: str) -> int:
    """Highest principal quantum number among occupied shells (number of electron shells)."""
    return max((n for n, _l, occ in expand(config) if occ > 0), default=0)


_NOBLE_CORES = [("[He]", 2), ("[Ne]", 10), ("[Ar]", 18), ("[Kr]", 36), ("[Xe]", 54), ("[Rn]", 86)]


def core_electrons_ion(n_electrons: int, config: str) -> int:
    """Chemical core for any species. When nothing is left outside the noble-gas core
    (alkali-metal cations such as Na+ = [Ne]), the species is treated like its isoelectronic
    noble gas: the outermost closed shell becomes the valence shell and the next inner
    noble-gas core is the core (Li+ -> 0, Na+ -> 2, K+ -> 10, ...)."""
    n_core = core_electrons(n_electrons, config)
    if n_core < n_electrons:
        return n_core
    inner = [n for tag, n in _NOBLE_CORES if n < n_electrons]
    return inner[-1] if inner else 0
