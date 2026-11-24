"""Tests of the transcribed reference radius tables.

Each test checks a transcription against a quantity stated in the source itself, so that
a typo in the tables is caught rather than propagated. Run with pytest, or directly."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from density_radii import load_set  # noqa: E402



def test_cordero_table():
    cov = load_set("cordero2008_covalent.csv")
    assert len(cov.symbols()) == 96
    assert cov.value("C") == 0.76                      # sp3 default
    assert cov.value("C", "sp2") == 0.73
    assert cov.value("Fe", "hs") == 1.52
    for sym in ("Mn", "Fe", "Co"):
        try:
            cov.value(sym)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{sym} should require a spin-state variant")
    # spot checks against Table 2
    for sym, val, esd in [("H", 0.31, 0.05), ("K", 2.03, 0.12), ("Te", 1.38, 0.04),
                          ("Zn", 1.22, 0.04), ("Ta", 1.70, 0.08), ("O", 0.66, 0.02),
                          ("Ra", 2.21, 0.02), ("Cm", 1.69, 0.03)]:
        rad = cov.get(sym)
        assert rad.value == val and rad.esd == esd, sym
    assert cov.get("Pm").esd is None


def test_vdw_tables():
    alv = load_set("alvarez2013_vdw.csv")
    bon = load_set("bondi1964_vdw.csv")
    man = load_set("mantina2009_vdw.csv")
    assert len(alv.symbols()) == 93 and len(bon.symbols()) == 38 and len(man.symbols()) == 44
    for sym, val in [("H", 1.20), ("O", 1.50), ("Zn", 2.39), ("Hg", 2.45), ("Cm", 3.05)]:
        assert alv.value(sym) == val, sym
    assert bon.value("Pt") == 1.72 and bon.value("Pt", "uv") == 1.75
    assert man.value("H") == 1.10 and man.value("Ra") == 2.83


def test_rahm_table():
    rahm = load_set("rahm2017_density.csv")
    assert len(rahm.symbols()) == 96
    # corrigendum values replace the 2016 ones
    for sym, val in [("Fe", 2.37), ("Co", 2.33), ("Ni", 2.29), ("Zr", 2.69),
                     ("Tc", 2.52), ("Th", 2.88), ("Ir", 2.40), ("Np", 2.81)]:
        assert rahm.value(sym) == val, sym
    # values quoted in the text of the 2016 paper
    for sym, val in [("H", 1.54), ("Be", 2.19), ("Ne", 1.56), ("C", 1.90), ("Pd", 2.15),
                     ("Rh", 2.33), ("Ag", 2.25), ("Pt", 2.30), ("Cr", 2.33), ("V", 2.52),
                     ("Mn", 2.42), ("Ba", 2.93), ("Ra", 2.92), ("Cd", 2.38), ("Hg", 2.29)]:
        assert rahm.value(sym) == val, sym


def test_pyykko_batsanov_tables():
    py = load_set("pyykko2009_covalent.csv")
    assert len(py.symbols()) == 118
    # bond-length sums quoted in the text of Pyykko & Atsumi (pm)
    assert round(100 * (py.value("Xe") + py.value("H"))) == 163
    assert round(100 * (py.value("H") + py.value("Ar"))) == 128
    assert round(100 * (py.value("Ar") + py.value("F"))) == 160
    bc = load_set("batsanov2001_vdw_crystallographic.csv")
    be = load_set("batsanov2001_vdw_equilibrium.csv")
    assert len(bc.symbols()) == 65 and len(be.symbols()) == 65
    assert bc.value("Cs") == 3.0 and be.value("Cs") == 3.30 and be.value("Bi") == 3.52   # Bi as printed


def _read_reply_widths():
    import csv
    path = Path(__file__).resolve().parent.parent / "data" / "radii" / "reply2024_crust_widths.csv"
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests passed")
