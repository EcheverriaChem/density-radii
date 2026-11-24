#!/usr/bin/env python3
"""
J-averaged analysis of NIST ASD levels (Step 1 of the revision plan, R1.7).

For every spectrum CSV produced by fetch_nist_levels.py this script
  1. reads the levels (configuration, term, J, energy),
  2. groups them into LS terms (contiguous blocks with the same configuration
     and term label, as given by the term-ordered NIST output),
  3. computes the J-averaged energy of every complete LS term,
         E(term) = sum_J (2J+1) E_J / sum_J (2J+1),
  4. compares the configuration of the ground level with the configuration
     of the lowest complete J-averaged term.

Classification of each spectrum
  Terms with missing J levels (at least two known) get a centre of gravity
  estimated with the Lande interval rule, flagged as *_estimated.

  same        the lowest J-averaged term belongs to the ground-level configuration
  reordered   a complete term of another configuration lies below the
              J-averaged ground term (the scalar-relativistic reference differs)
  possible    no complete term proves a reordering, but some level of another
              configuration lies below the J-averaged ground term
  same (config. average) / reordered (config. average)
              decided from complete configuration averages when no LS term
              of the ground configuration can be averaged (jj-coupled spectra)
  incomplete  no complete term or configuration (J-average undefined)
  ground level only  the ASD lists a single level for this spectrum

Usage
  python3 javg_nist.py nist_levels/            -> javg_summary.csv, javg_terms.csv
  python3 javg_nist.py nist_levels/ --emax 3.0 (eV window for the term table)
"""
import argparse
import csv
import glob
import os
import re
from fractions import Fraction

CM_TO_EV = 1.0 / 8065.543937
LS_RE = re.compile(r'^(?:[a-z]\'*\s+)?(\d+)([SPDFGHIKLMNOQRTUVWXYZ])(\*?)$')
L_LETTERS = 'SPDFGHIKLMNOQRTUVWXYZ'   # L = 0, 1, 2, ... (J is skipped in spectroscopy)
ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
            "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
            "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
            "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm").split()


def clean(x):
    x = x.strip()
    if x.startswith('='):
        x = x[1:]
    return x.strip('"').strip()


def norm_config(conf):
    """Remove parent terms and couplings: '3d8.(3F).4s2' -> '3d8.4s2'."""
    c = conf
    prev = None
    while prev != c:
        prev = c
        c = re.sub(r'\([^()]*\)', '', c)
    c = re.sub(r'<[^<>]*>', '', c).replace(' ', '')
    parts = [p for p in c.split('.') if p and not re.fullmatch(r'\dp6', p)]
    return '.'.join(parts)


def parse_j(s):
    s = s.strip().rstrip('?')
    if not re.fullmatch(r'\d+(/2)?', s):
        return None
    return float(Fraction(s))


def expected_js(mult, lett):
    S = Fraction(mult - 1, 2)
    L = L_LETTERS.index(lett)
    lo, hi = abs(L - S), L + S
    js, j = [], lo
    while j <= hi:
        js.append(float(j))
        j += 1
    return js


def read_levels(path):
    with open(path, encoding='utf-8') as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    levels = []
    for r in rows[1:]:
        d = {k: clean(v) for k, v in zip(head, r)}
        term = d.get('Term', '')
        if term == 'Limit' or 'Limit' in d.get('Configuration', ''):
            continue
        e = d.get('Level (cm-1)', '')
        suffix = d.get('Suffix', '')
        if not re.fullmatch(r'-?\d+(\.\d*)?', e) or 'x' in suffix:
            continue   # unknown or relative (+x) energies cannot be used
        levels.append({
            'conf': d.get('Configuration', ''), 'term': term,
            'J': parse_j(d.get('J', '')), 'E': float(e),
            'flag': (d.get('Prefix', '') + suffix).strip(),
        })
    return levels


def group_terms(levels):
    """Contiguous blocks of (configuration, term); a repeated J starts a new block."""
    groups, cur, key = [], None, None
    for lv in levels:
        k = (lv['conf'], lv['term'])
        if k != key or (lv['J'] is not None and lv['J'] in [x['J'] for x in cur['levels']]):
            cur = {'conf': lv['conf'], 'term': lv['term'], 'levels': []}
            groups.append(cur)
            key = k
        cur['levels'].append(lv)
    for g in groups:
        g['nconf'] = norm_config(g['conf'])
        g['m'] = LS_RE.match(g['term'])

    # attach unlabelled levels (term blank or '*') of the same configuration
    # to an incomplete LS term that lacks exactly that J
    unl = [g for g in groups if g['term'] in ('', '*') and len(g['levels']) == 1]
    used = set()
    for g in groups:
        if not g['m']:
            continue
        exp = expected_js(int(g['m'].group(1)), g['m'].group(2))
        have = {lv['J'] for lv in g['levels']}
        for j in sorted(set(exp) - have):
            cands = [u for u in unl if id(u) not in used and u['nconf'] == g['nconf']
                     and u['levels'][0]['J'] == j]
            if len(cands) == 1:   # attach only when unambiguous
                g['levels'].append(dict(cands[0]['levels'][0], attached=True))
                used.add(id(cands[0]))
    groups = [g for g in groups if id(g) not in used]

    out = []
    for g in groups:
        m = g['m']
        js = [lv['J'] for lv in g['levels']]
        g.update(emin=min(lv['E'] for lv in g['levels']), ls=bool(m), complete=False,
                 estimated=False, eavg=None, missing=None)
        if m and all(j is not None for j in js):
            mult, lett = int(m.group(1)), m.group(2)
            exp = expected_js(mult, lett)
            missing = sorted(set(exp) - set(js))
            g['missing'] = missing
            if not missing and sorted(js) == sorted(exp):
                g['complete'] = True
                w = [2 * lv['J'] + 1 for lv in g['levels']]
                g['eavg'] = sum(wi * lv['E'] for wi, lv in zip(w, g['levels'])) / sum(w)
            elif len(js) >= 2 and set(js) <= set(exp):
                # Lande interval rule: E_J = Ec + (A/2) [J(J+1) - L(L+1) - S(S+1)].
                # The (2J+1)-weighted mean of the bracket is zero, so the
                # centre of gravity is Ec, obtained from a linear fit.
                L = L_LETTERS.index(lett)
                S = (mult - 1) / 2
                x = [j * (j + 1) - L * (L + 1) - S * (S + 1) for j in js]
                y = [lv['E'] for lv in g['levels']]
                n = len(x)
                mx, my = sum(x) / n, sum(y) / n
                sxx = sum((xi - mx) ** 2 for xi in x)
                if sxx > 0:
                    a = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / sxx
                    g['eavg'] = my - a * mx
                    g['estimated'] = True
        out.append(g)
    return out


from math import comb


def config_states(nconf):
    """Number of states of a configuration, e.g. '6s2.6p2' -> 1 * 15."""
    total = 1
    for part in nconf.split('.'):
        m = re.fullmatch(r'\d+([spdfg])(\d*)', part)
        if not m:
            return None
        l = 'spdfg'.index(m.group(1))
        occ = int(m.group(2) or 1)
        total *= comb(2 * (2 * l + 1), occ)
    return total


def config_average(levels, nconf):
    """g-weighted mean over all levels of a configuration, if all its states are known."""
    lv = [x for x in levels if norm_config(x['conf']) == nconf and x['J'] is not None]
    n = config_states(nconf)
    if not lv or n is None or sum(2 * x['J'] + 1 for x in lv) != n:
        return None
    return sum((2 * x['J'] + 1) * x['E'] for x in lv) / n


def analyse(path):
    levels = read_levels(path)
    if not levels:
        return None, []
    terms = group_terms(levels)
    ground_lv = min(levels, key=lambda lv: lv['E'])
    gterm = next(t for t in terms if ground_lv in t['levels'])
    gconf = gterm['nconf']
    res = {'ground_conf': gconf, 'ground_term': gterm['term'],
           'ground_term_complete': gterm['complete'],
           'E_ground_term_avg': gterm['eavg'],
           'ground_term_estimated': gterm['estimated']}
    complete = [t for t in terms if t['eavg'] is not None]
    # lowest complete term of the ground configuration and of any other configuration
    same = [t for t in complete if t['nconf'] == gconf]
    other = [t for t in complete if t['nconf'] != gconf]
    best_same = min(same, key=lambda t: t['eavg']) if same else None
    best_other = min(other, key=lambda t: t['eavg']) if other else None
    ref = best_same['eavg'] if best_same else None
    res.update(
        best_same_term=best_same['term'] if best_same else '',
        E_best_same=best_same['eavg'] if best_same else None,
        alt_conf=best_other['nconf'] if best_other else '',
        alt_estimated=best_other['estimated'] if best_other else '',
        same_estimated=best_same['estimated'] if best_same else '',
        alt_term=best_other['term'] if best_other else '',
        E_alt=best_other['eavg'] if best_other else None,
    )
    lowest_other_level = min((lv for lv in levels if norm_config(lv['conf']) != gconf),
                             key=lambda lv: lv['E'], default=None)
    res['lowest_other_level_conf'] = norm_config(lowest_other_level['conf']) if lowest_other_level else ''
    res['E_lowest_other_level'] = lowest_other_level['E'] if lowest_other_level else None

    if len(levels) == 1:
        status = 'ground level only'
    elif ref is None:
        # fallback: complete configuration averages (useful for jj-labelled spectra)
        cav = config_average(levels, gconf)
        others = sorted({norm_config(lv['conf']) for lv in levels} - {gconf})
        oav = [(config_average(levels, c), c) for c in others]
        oav = [t for t in oav if t[0] is not None]
        res['E_ground_config_avg'] = cav
        if cav is not None and lowest_other_level and lowest_other_level['E'] > cav:
            status = 'same (config. average)'
        elif cav is not None and oav and min(oav)[0] < cav:
            status = 'reordered (config. average)'
        else:
            status = 'incomplete'
    elif best_other and best_other['eavg'] < ref:
        status = 'reordered'
    elif lowest_other_level and lowest_other_level['E'] < ref:
        status = 'possible'
    else:
        status = 'same'
    res['status'] = status
    res['gap_cm'] = (best_other['eavg'] - ref) if (best_other and ref is not None) else None
    return res, terms


def fmt(x, nd=1):
    return '' if x is None else f'{x:.{nd}f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder')
    ap.add_argument('--emax', type=float, default=3.0, help='eV window for javg_terms.csv')
    ap.add_argument('--out', default='.')
    a = ap.parse_args()

    summ, termrows = [], []
    for z, sym in enumerate(ELEMENTS, start=1):
        for ion in ('I', 'II'):
            path = os.path.join(a.folder, f'{sym}_{ion}.csv')
            if not os.path.exists(path):
                continue
            res, terms = analyse(path)
            if res is None:
                summ.append({'Z': z, 'species': f'{sym} {ion}', 'status': 'no data'})
                continue
            row = {'Z': z, 'species': f'{sym} {ion}', **res}
            summ.append(row)
            for t in terms:
                if t['emin'] * CM_TO_EV <= a.emax:
                    termrows.append({'Z': z, 'species': f'{sym} {ion}', 'config': t['nconf'],
                                     'config_full': t['conf'], 'term': t['term'],
                                     'complete': t['complete'], 'lande_estimate': t['estimated'],
                                     'missing_J': '' if not t['missing'] else ' '.join(str(j) for j in t['missing']),
                                     'n_levels': len(t['levels']),
                                     'E_min_cm': fmt(t['emin'], 3), 'E_avg_cm': fmt(t['eavg'], 3)})

    cols = ['Z', 'species', 'status', 'ground_conf', 'ground_term', 'ground_term_complete',
            'ground_term_estimated', 'E_ground_term_avg', 'same_estimated', 'alt_estimated', 'best_same_term', 'E_best_same', 'alt_conf', 'alt_term', 'E_alt',
            'gap_cm', 'gap_eV', 'E_ground_config_avg', 'lowest_other_level_conf', 'E_lowest_other_level']
    with open(os.path.join(a.out, 'javg_summary.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in summ:
            r = dict(r)
            if r.get('gap_cm') is not None:
                r['gap_eV'] = fmt(r['gap_cm'] * CM_TO_EV, 4)
            for k in ('E_ground_term_avg', 'E_best_same', 'E_alt', 'gap_cm', 'E_ground_config_avg', 'E_lowest_other_level'):
                if isinstance(r.get(k), float):
                    r[k] = fmt(r[k], 1)
            w.writerow(r)
    with open(os.path.join(a.out, 'javg_terms.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(termrows[0].keys()))
        w.writeheader()
        w.writerows(termrows)

    from collections import Counter
    print(Counter(r['status'] for r in summ))
    for r in summ:
        if r['status'] != 'same':
            gap = r.get('gap_cm')
            print(f"{r['species']:7s} {r['status']:10s} ground {r.get('ground_conf',''):22s} {r.get('ground_term',''):8s}"
                  f" alt {r.get('alt_conf',''):22s} {r.get('alt_term',''):8s} gap {fmt(gap,0):>7s} cm-1"
                  f"  lowest other lvl {r.get('lowest_other_level_conf','')} {fmt(r.get('E_lowest_other_level'),0)}")


if __name__ == '__main__':
    main()
