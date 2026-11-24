"""Loading and querying atomic radius sets.

Each radius set is a CSV in data/radii/ with at least the columns
    set_id, symbol, variant, kind, value_A, esd_A
where kind is "r" (covalent or valence), "v" (van der Waals), "c" (core; c_orb for the density sets
from version 1.1), and, for the density sets, "c_count" and "c_elf" (alternative core radii).

Rules enforced here:
  * values are never invented: a missing element raises KeyError;
  * elements with several variants and no declared default (Mn, Fe, Co
    spin states in Cordero 2008) raise an error unless a variant is given;
  * C defaults to sp3 in Cordero 2008 (declared in the CSV note).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "radii"


@dataclass(frozen=True)
class Radius:
    symbol: str
    kind: str
    value: float
    esd: float | None
    variant: str
    set_id: str


class RadiusSet:
    def __init__(self, set_id: str, kind: str, rows: list[dict]):
        self.set_id = set_id
        self.kind = kind
        self._by_symbol: dict[str, dict[str, Radius]] = {}
        self._default: dict[str, str] = {}
        for row in rows:
            sym = row["symbol"].strip()
            var = (row.get("variant") or "").strip()
            esd = row.get("esd_A", "").strip()
            rad = Radius(sym, kind, float(row["value_A"]),
                         float(esd) if esd else None, var, set_id)
            self._by_symbol.setdefault(sym, {})[var] = rad
            if (row.get("note") or "").strip().lower().startswith("default variant"):
                self._default[sym] = var

    def symbols(self) -> list[str]:
        return list(self._by_symbol)

    def variants(self, symbol: str) -> list[str]:
        return list(self._by_symbol[symbol])

    def get(self, symbol: str, variant: str | None = None) -> Radius:
        if symbol not in self._by_symbol:
            raise KeyError(f"{symbol} not in radius set {self.set_id}")
        options = self._by_symbol[symbol]
        if variant is not None:
            if variant not in options:
                raise KeyError(f"{symbol} has no variant '{variant}' in "
                               f"{self.set_id}; available: {list(options)}")
            return options[variant]
        if len(options) == 1:
            return next(iter(options.values()))
        if symbol in self._default:
            return options[self._default[symbol]]
        raise ValueError(f"{symbol} has several variants in {self.set_id} "
                         f"({list(options)}) and no default: specify one")

    def value(self, symbol: str, variant: str | None = None) -> float:
        return self.get(symbol, variant).value


def load_set(filename: str, kind: str | None = None,
             data_dir: Path = DATA_DIR) -> RadiusSet:
    path = Path(data_dir) / filename
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"empty radius file {path}")
    kinds = {r["kind"] for r in rows}
    if kind is None:
        if len(kinds) != 1:
            raise ValueError(f"{path} mixes kinds {kinds}; pass kind=")
        kind = kinds.pop()
    rows = [r for r in rows if r["kind"] == kind]
    return RadiusSet(rows[0]["set_id"], kind, rows)


def sanity_check(cov: RadiusSet, vdw: RadiusSet,
                 w_min: float = 0.4, w_max: float = 1.5) -> list[str]:
    """Return a list of problems for elements present in both sets."""
    problems = []
    for sym in vdw.symbols():
        if sym not in cov.symbols():
            continue
        for var in cov.variants(sym):
            r = cov.get(sym, var).value
            v = vdw.get(sym).value if len(vdw.variants(sym)) == 1 else math.nan
            if math.isnan(v):
                continue
            w = v - r
            if w <= 0:
                problems.append(f"{sym}{'/' + var if var else ''}: v <= r")
            elif not (w_min <= w <= w_max):
                problems.append(f"{sym}{'/' + var if var else ''}: "
                                f"crust width {w:.2f} A outside [{w_min}, {w_max}]")
    return problems
