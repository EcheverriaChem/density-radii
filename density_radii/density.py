"""Spherically averaged electron density of a single atom from a Molden file.

For one atom at the origin, basis functions are chi_{k,m}(r) = R_k(r) Y_{l,m}(Omega)
with real, orthonormal spherical harmonics (Molden [5D]/[7F]/[9G] = pure functions).
The angular average kills all cross terms between different (l, m), so

    rho_bar(r) = 1/(4 pi) * sum_i n_i * sum_{l,m} ( sum_{k in l} C_{i,k,m} R_k(r) )^2

and only the radial parts are needed: no angular conventions, no m ordering beyond
"slot j of every shell with the same l is the same m", which Molden guarantees.

The radial primitives' normalization convention in the file is detected, not assumed:
both conventions are tried and the one for which every occupied MO is normalized
(C^T S C = 1) is kept. The electron count integral is reported as a check.

Units: bohr internally; radii are returned in angstrom.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz   # NumPy < 2.0 compatibility

BOHR_TO_ANG = 0.529177210903
LCHAR = "spdfghi"
PURE_TAGS = {"[5D]", "[7F]", "[9G]", "[5D7F]", "[5D10F]"}


@dataclass
class Shell:
    l: int
    exps: np.ndarray
    coefs: np.ndarray


@dataclass
class MO:
    spin: str
    occ: float
    energy: float
    coef: np.ndarray


@dataclass
class MoldenAtom:
    symbol: str
    Z: int
    shells: list[Shell]
    mos: list[MO]
    pure: bool
    notes: list[str] = field(default_factory=list)

    @property
    def nbf(self) -> int:
        return sum((2 * s.l + 1) if self.pure else (s.l + 1) * (s.l + 2) // 2 for s in self.shells)


def read_molden(path) -> MoldenAtom:
    lines = Path(path).read_text(errors="replace").splitlines()
    sec, shells, mos, atoms = None, [], [], []
    tags = set()
    i = 0
    cur_mo = None
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("["):
            tag = line.split("]")[0] + "]"
            tags.add(tag.upper())
            sec = tag.upper()
            i += 1
            continue
        if sec == "[ATOMS]" and line:
            p = line.split()
            atoms.append((p[0], int(p[2]), [float(x) for x in p[3:6]]))
        elif sec == "[GTO]" and line:
            p = line.split()
            if len(p) >= 2 and p[0].lower() in LCHAR and p[1].isdigit():
                l = LCHAR.index(p[0].lower())
                n = int(p[1])
                ex, co = [], []
                for j in range(1, n + 1):
                    a, c = lines[i + j].replace("D", "E").split()[:2]
                    ex.append(float(a))
                    co.append(float(c))
                shells.append(Shell(l, np.array(ex), np.array(co)))
                i += n
        elif sec == "[MO]" and line:
            if line.startswith("Sym="):
                cur_mo = {"coef": {}}
            elif line.startswith("Ene="):
                cur_mo["energy"] = float(line.split("=")[1])
            elif line.startswith("Spin="):
                cur_mo["spin"] = line.split("=")[1].strip().lower()
            elif line.startswith("Occup="):
                cur_mo["occ"] = float(line.split("=")[1])
                mos.append(cur_mo)
            else:
                k, v = line.split()[:2]
                cur_mo["coef"][int(k)] = float(v)
        i += 1
    if len(atoms) != 1:
        raise ValueError(f"{path}: expected one atom, found {len(atoms)}")
    if any(np.linalg.norm(atoms[0][2]) > 1e-8 for _ in [0]):
        raise ValueError(f"{path}: atom not at the origin")
    pure = bool(tags & PURE_TAGS) or any(s.l >= 2 for s in shells) and "[5D]" in tags
    if any(s.l >= 2 for s in shells) and not (tags & PURE_TAGS):
        raise ValueError(f"{path}: Cartesian d/f functions are not supported")
    atom = MoldenAtom(atoms[0][0], atoms[0][1], shells, [], True)
    nbf = atom.nbf
    for m in mos:
        c = np.zeros(nbf)
        for k, v in m["coef"].items():
            c[k - 1] = v
        atom.mos.append(MO(m.get("spin", "alpha"), m["occ"], m.get("energy", float("nan")), c))
    return atom


def _prim_norm(a: np.ndarray, l: int) -> np.ndarray:
    return np.sqrt(2.0 * (2.0 * a) ** (l + 1.5) / math.gamma(l + 1.5))


def _radial_overlap(e1, c1, e2, c2, l) -> float:
    s = e1[:, None] + e2[None, :]
    return float(c1 @ (math.gamma(l + 1.5) / (2.0 * s ** (l + 1.5))) @ c2)


class AtomDensity:
    """Spherically averaged density of an atom read from Molden."""

    def __init__(self, path, convention: str = "auto"):
        self.atom = read_molden(path)
        self.path = str(path)
        convs = ["normalized", "unnormalized"] if convention == "auto" else [convention]
        best = None
        for conv in convs:
            eff = self._effective_coefs(conv)
            err, nel = self._mo_norm_error(eff)
            if best is None or err < best[1]:
                best = (conv, err, nel, eff)
        self.convention, self.mo_norm_error, self.n_electrons, self._eff = best
        # group basis-function indices by l and m-slot
        self._groups: dict[tuple[int, int], list[tuple[int, int]]] = {}
        idx = 0
        for k, sh in enumerate(self.atom.shells):
            for m in range(2 * sh.l + 1):
                self._groups.setdefault((sh.l, m), []).append((k, idx))
                idx += 1

    # -- normalization ----------------------------------------------------------
    def _effective_coefs(self, conv):
        out = []
        for sh in self.atom.shells:
            c = sh.coefs * (_prim_norm(sh.exps, sh.l) if conv == "normalized" else 1.0)
            n = math.sqrt(_radial_overlap(sh.exps, c, sh.exps, c, sh.l))
            out.append(c / n)
        return out

    def _mo_norm_error(self, eff):
        shells = self.atom.shells
        # overlap matrix block-diagonal in (l, m)
        idx, blocks = 0, {}
        for k, sh in enumerate(shells):
            for m in range(2 * sh.l + 1):
                blocks.setdefault((sh.l, m), []).append((k, idx))
                idx += 1
        nbf = idx
        S = np.zeros((nbf, nbf))
        for members in blocks.values():
            for k1, i1 in members:
                for k2, i2 in members:
                    l = shells[k1].l
                    S[i1, i2] = _radial_overlap(shells[k1].exps, eff[k1], shells[k2].exps, eff[k2], l)
        worst, nel = 0.0, 0.0
        for mo in self.atom.mos:
            if mo.occ <= 0:
                continue
            nrm = float(mo.coef @ S @ mo.coef)
            worst = max(worst, abs(nrm - 1.0))
            nel += mo.occ * nrm
        return worst, nel

    # -- density ------------------------------------------------------------------
    def _radial(self, k, r):
        sh = self.atom.shells[k]
        return (r[:, None] ** sh.l * np.exp(-np.outer(r * r, sh.exps))) @ self._eff[k]

    def rho(self, r_bohr, spin: str | None = None) -> np.ndarray:
        """Spherically averaged density (e/bohr^3). spin: None (total), 'alpha' or 'beta'."""
        r = np.atleast_1d(np.asarray(r_bohr, float))
        radials = {k: self._radial(k, r) for k in range(len(self.atom.shells))}
        occ = [mo for mo in self.atom.mos if mo.occ > 0 and (spin is None or mo.spin == spin)]
        C = np.array([mo.coef for mo in occ])
        n = np.array([mo.occ for mo in occ])
        tot = np.zeros_like(r)
        for members in self._groups.values():
            # phi[i, r] = sum_k C[i, idx_k] R_k(r)
            phi = sum(np.outer(C[:, idx], radials[k]) for k, idx in members)
            tot += n @ (phi ** 2)
        return tot / (4.0 * math.pi)

    def electron_count(self, rmax=30.0, npts=6000) -> float:
        x = np.linspace(0, 1, npts)
        r = rmax * x ** 3            # dense near the nucleus
        f = 4 * math.pi * r ** 2 * self.rho(r)
        return float(_trapz(f, r))

    # -- radii ----------------------------------------------------------------------
    def isodensity_radius(self, rho_c=0.001, rmin=0.5, rmax=25.0, npts=4000) -> float:
        """Outermost r (angstrom) where rho_bar(r) = rho_c (log-interpolated)."""
        r = np.linspace(rmin, rmax, npts)
        lr = np.log(np.maximum(self.rho(r), 1e-300))
        above = np.where(lr >= math.log(rho_c))[0]
        if len(above) == 0 or above[-1] == len(r) - 1:
            return float("nan")
        j = above[-1]
        x0, x1, y0, y1 = r[j], r[j + 1], lr[j], lr[j + 1]
        rc = x0 + (math.log(rho_c) - y0) * (x1 - x0) / (y1 - y0)
        return rc * BOHR_TO_ANG

    def rdf_outer_max(self, rmin=0.05, rmax=15.0, npts=6000) -> float:
        """Outermost maximum of r^2 rho_bar(r), in angstrom."""
        r = np.linspace(rmin, rmax, npts)
        f = r * r * self.rho(r)
        d = np.diff(f)
        idx = np.where((d[:-1] > 0) & (d[1:] <= 0))[0]
        if len(idx) == 0:
            return float("nan")
        j = idx[-1] + 1
        # parabolic refinement
        y0, y1, y2 = f[j - 1], f[j], f[j + 1]
        den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return (r[j] + off * (r[1] - r[0])) * BOHR_TO_ANG


# ---------------------------------------------------------------------------------
# Orbital-resolved and ELF quantities (spherical averages)
# ---------------------------------------------------------------------------------
def _radial_deriv(self, k, r):
    sh = self.atom.shells[k]
    e = np.exp(-np.outer(r * r, sh.exps))
    lpart = (sh.l * r[:, None] ** (sh.l - 1)) if sh.l > 0 else np.zeros((len(r), 1))
    d = (lpart - 2.0 * np.outer(r ** (sh.l + 1), sh.exps)) * e
    return d @ self._eff[k]


AtomDensity._radial_deriv = _radial_deriv


def _spin_quantities(self, r, spin):
    """rho_s, |grad rho_s|, tau_s (tau = sum_i n_i |grad phi_i|^2, no 1/2), angularly averaged."""
    R = {k: self._radial(k, r) for k in range(len(self.atom.shells))}
    dR = {k: self._radial_deriv(k, r) for k in range(len(self.atom.shells))}
    occ = [mo for mo in self.atom.mos if mo.occ > 0 and mo.spin == spin]
    rho = np.zeros_like(r); drho = np.zeros_like(r); tau = np.zeros_like(r)
    if not occ:
        return rho, drho, tau
    C = np.array([mo.coef for mo in occ]); n = np.array([mo.occ for mo in occ])
    for (l, _m), members in self._groups.items():
        phi = sum(np.outer(C[:, idx], R[k]) for k, idx in members)
        dphi = sum(np.outer(C[:, idx], dR[k]) for k, idx in members)
        rho += n @ (phi ** 2)
        drho += n @ (2 * phi * dphi)
        tau += n @ (dphi ** 2 + l * (l + 1) * phi ** 2 / r ** 2)
    f = 1.0 / (4 * math.pi)
    return rho * f, np.abs(drho) * f, tau * f


AtomDensity._spin_quantities = _spin_quantities


def elf(self, r_bohr):
    """Becke-Edgecombe ELF from spherically averaged spin densities (approximation for open shells)."""
    r = np.asarray(r_bohr, float)
    cf = 0.6 * (6 * math.pi ** 2) ** (2.0 / 3.0)
    D = np.zeros_like(r); D0 = np.zeros_like(r)
    for spin in ("alpha", "beta"):
        rho, g, tau = self._spin_quantities(r, spin)
        with np.errstate(divide="ignore", invalid="ignore"):
            D += np.where(rho > 1e-30, tau - 0.25 * g * g / rho, 0.0)
        D0 += cf * np.maximum(rho, 0) ** (5.0 / 3.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = np.where(D0 > 1e-300, D / D0, 1e30)
    return 1.0 / (1.0 + chi ** 2)


AtomDensity.elf = elf


def elf_shell_minima(self, rmin=0.005, rmax=8.0, npts=8000):
    """Radii (angstrom) of ELF minima between shells, innermost first."""
    r = np.geomspace(rmin, rmax, npts)
    e = self.elf(r)
    idx = np.where((e[1:-1] < e[:-2]) & (e[1:-1] <= e[2:]))[0] + 1
    idx = [i for i in idx if e[i] < 0.95]           # ignore numerical ripples near 1
    return [float(r[i] * BOHR_TO_ANG) for i in idx]


AtomDensity.elf_shell_minima = elf_shell_minima


def homo_radius(self, rmin=0.05, rmax=20.0, npts=6000):
    """Most probable radius (angstrom) of the highest occupied spin orbital:
    outermost maximum of r^2 |phi|^2 (angularly averaged)."""
    occ = [mo for mo in self.atom.mos if mo.occ > 0]
    homo = max(occ, key=lambda m: m.energy)
    r = np.linspace(rmin, rmax, npts)
    R = {k: self._radial(k, r) for k in range(len(self.atom.shells))}
    dens = np.zeros_like(r); lw = {}
    for (l, _m), members in self._groups.items():
        phi = sum(homo.coef[idx] * R[k] for k, idx in members)
        w = float(_trapz(phi * phi * r * r, r))
        lw[l] = lw.get(l, 0.0) + w
        dens += phi * phi
    f = r * r * dens
    j = int(np.argmax(np.where(r > 0.3 * r[np.argmax(f)], f, 0)))
    # outermost local maximum
    d = np.diff(f)
    loc = np.where((d[:-1] > 0) & (d[1:] <= 0))[0] + 1
    loc = [i for i in loc if f[i] > 0.05 * f.max()]
    j = loc[-1] if loc else j
    tot = sum(lw.values()) or 1.0
    dominant_l = max(lw, key=lw.get)
    return float(r[j] * BOHR_TO_ANG), LCHAR[dominant_l], lw[dominant_l] / tot, homo.spin, homo.energy


AtomDensity.homo_radius = homo_radius


def orbital_radii(self, rmin=0.02, rmax=20.0, npts=5000):
    """Most probable radius (angstrom), dominant l and energy of every occupied spin orbital."""
    r = np.linspace(rmin, rmax, npts)
    R = {k: self._radial(k, r) for k in range(len(self.atom.shells))}
    out = []
    for mo in self.atom.mos:
        if mo.occ <= 0:
            continue
        dens = np.zeros_like(r); lw = {}
        for (l, _m), members in self._groups.items():
            phi = sum(mo.coef[idx] * R[k] for k, idx in members)
            lw[l] = lw.get(l, 0.0) + float(_trapz(phi * phi * r * r, r))
            dens += phi * phi
        f = r * r * dens
        out.append((float(r[int(np.argmax(f))] * BOHR_TO_ANG), LCHAR[max(lw, key=lw.get)], mo.energy, mo.spin))
    return out


AtomDensity.orbital_radii = orbital_radii


def enclosing_radius(self, n_el, rmax=15.0, npts=20000):
    """Radius (angstrom) of the sphere containing n_el electrons."""
    x = np.linspace(0, 1, npts)
    r = rmax * x ** 3
    f = 4 * math.pi * r ** 2 * self.rho(r)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])
    j = int(np.searchsorted(cum, n_el))
    if j == 0 or j >= len(r):
        return float("nan")
    t = (n_el - cum[j - 1]) / (cum[j] - cum[j - 1])
    return float((r[j - 1] + t * (r[j] - r[j - 1])) * BOHR_TO_ANG)


AtomDensity.enclosing_radius = enclosing_radius


def mo_characters(self, spin, n_virtual=15, rmin=0.02, rmax=20.0, npts=4000):
    """Per spin: list of dicts (index within spin, energy, occ, dominant l, l-weight, r_max in angstrom)
    for all occupied MOs and the first n_virtual virtuals, in ORCA's file order."""
    r = np.linspace(rmin, rmax, npts)
    R = {k: self._radial(k, r) for k in range(len(self.atom.shells))}
    # Keep ORCA's own order: it is the numbering used by Rotate and ALPHACONF/BETACONF.
    # With a level shift ORCA writes occupied orbitals first, which is NOT energy order.
    mos = [m for m in self.atom.mos if m.spin == spin]
    nocc = sum(1 for m in mos if m.occ > 0)
    if [i for i, m in enumerate(mos) if m.occ > 0] != list(range(nocc)):
        raise ValueError("occupied orbitals are not the first ones in the file: cannot index for Rotate/CONF")
    out = []
    for i, mo in enumerate(mos[: nocc + n_virtual]):
        dens = np.zeros_like(r); lw = {}
        for (l, _m), members in self._groups.items():
            phi = sum(mo.coef[idx] * R[k] for k, idx in members)
            lw[l] = lw.get(l, 0.0) + float(_trapz(phi * phi * r * r, r))
            dens += phi * phi
        tot = sum(lw.values()) or 1.0
        ldom = max(lw, key=lw.get)
        out.append({"i": i, "energy": mo.energy, "occ": mo.occ, "l": LCHAR[ldom],
                    "w": lw[ldom] / tot, "rmax": float(r[int(np.argmax(r * r * dens))] * BOHR_TO_ANG)})
    return out, nocc


AtomDensity.mo_characters = mo_characters
