"""Tests for the Molden density code and configuration helpers."""
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "workflow"))
from density_radii.density import AtomDensity, BOHR_TO_ANG  # noqa: E402
from density_radii.orca_parse import core_electrons, l_totals  # noqa: E402

MOLDEN = """[Molden Format]
[Atoms] AU
He 1 2 0.0 0.0 0.0
[GTO]
  1 0
s   1 1.0
  {a}  {c}
p   1 1.0
  0.5  1.0

[5D]
[MO]
 Sym= 1a
 Ene= -0.9
 Spin= Alpha
 Occup= 1.0
  1  1.0
  2  0.0
  3  0.0
  4  0.0
 Sym= 1a
 Ene= -0.9
 Spin= Beta
 Occup= 1.0
  1  1.0
  2  0.0
  3  0.0
  4  0.0
"""


def _density(a, c):
    tmp = Path(tempfile.mkdtemp()) / "t.molden"
    tmp.write_text(MOLDEN.format(a=a, c=c))
    return AtomDensity(tmp)


def test_single_gaussian_density():
    a = 0.8
    for c in (1.0, 37.5):                         # arbitrary contraction scale must not matter
        d = _density(a, c)
        assert d.mo_norm_error < 1e-10
        assert abs(d.electron_count() - 2.0) < 1e-6
        # analytic: rho = 2 * (2a/pi)^(3/2) exp(-2 a r^2)
        rc = 0.001
        r_exact = math.sqrt(math.log(2 * (2 * a / math.pi) ** 1.5 / rc) / (2 * a)) * BOHR_TO_ANG
        assert abs(d.isodensity_radius(rc) - r_exact) < 1e-3
        # enclosing radius of 1 electron: half the density
        r = d.enclosing_radius(1.0) / BOHR_TO_ANG
        from math import erf, exp, pi, sqrt
        x = sqrt(2 * a) * r
        frac = erf(x) - 2 * x * exp(-x * x) / sqrt(pi)
        assert abs(frac - 0.5) < 1e-4


def test_core_electrons():
    assert core_electrons(31, "[Ar]3d10 4s2 4p1") == 28
    assert core_electrons(46, "[Kr]4d10") == 36
    assert core_electrons(74, "[Xe]4f14 5d4 6s2") == 68
    assert core_electrons(81, "[Xe]4f14 5d10 6s2 6p1") == 78
    assert core_electrons(8, "[He]2s2 2p4") == 2
    assert core_electrons(1, "1s1") == 0
    assert l_totals("[Xe]5d1 6s2") == {"s": 12, "p": 24, "d": 21, "f": 0, "g": 0}


def test_ion_helpers():
    from density_radii.orca_parse import parse_job_name, max_principal_n, core_electrons
    assert parse_job_name("009_F-_PBE0_anorcc_decon") == (9, "F", -1, "PBE0", "anorcc")
    assert parse_job_name("026_Fe+_HF_x2c") == (26, "Fe", 1, "HF", "x2c")
    assert parse_job_name("026_Fe_HF_x2c") == (26, "Fe", 0, "HF", "x2c")
    assert max_principal_n("[Kr]4d10") == 4          # Pd: four shells, not five
    assert core_electrons(10, "[Ne]") == 10          # Na+: no valence left ...
    from density_radii.orca_parse import core_electrons_ion
    assert core_electrons_ion(10, "[Ne]") == 2       # ... so it is treated like Ne: [He] core
    assert core_electrons_ion(2, "[He]") == 0        # Li+ like He
    assert core_electrons_ion(54, "[Xe]") == 36      # Cs+ like Xe: [Kr] core
    assert core_electrons_ion(28, "[Ar]3d10") == 18  # Cu+ keeps its 3d10 valence shell
    assert core_electrons(18, "[Ne]3s2 3p6") == 10   # Cl-


def test_conf_string():
    from make_config_fixes import conf_string
    assert conf_string(12, 13, 12) == "0,1"        # HOMO -> LUMO
    assert conf_string(11, 13, 12) == "0,1,1"      # HOMO-1 -> LUMO (ORCA manual)
    assert conf_string(12, 14, 12) == "0,0,1"      # HOMO -> LUMO+1 (ORCA manual)
    assert conf_string(12, 10, 12) == "1,1,0"      # fill an orbital below the emptied one
    assert conf_string(14, 12, 13) == "1,1,0"      # same, with the hole above the HOMO
    from make_config_fixes import conf_string_multi
    assert conf_string_multi([12], [13], 12) == conf_string(12, 13, 12)
    assert conf_string_multi([11, 12], [13, 14], 12) == "0,0,1,1"   # two electrons moved in one spin


if __name__ == "__main__":
    for k, v in dict(globals()).items():
        if k.startswith("test_"):
            v(); print("PASS", k)
