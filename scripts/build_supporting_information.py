"""Build the Supporting Information (docx) from the density_radii package data."""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

PKG = Path(__file__).resolve().parent.parent
R, REV, DATA = PKG / 'results', PKG / 'results' / 'revision', PKG / 'data'
import sys
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else PKG / 'results' / 'Supporting_Information.docx'
FONT = 'Times New Roman'

doc = Document()
st = doc.styles['Normal']; st.font.name = FONT; st.font.size = Pt(11)
st.element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
for s in ('Heading 1', 'Heading 2'):
    h = doc.styles[s]; h.font.name = FONT; h.font.color.rgb = None
    h.element.rPr.rFonts.set(qn('w:eastAsia'), FONT)
doc.styles['Heading 1'].font.size = Pt(13); doc.styles['Heading 2'].font.size = Pt(12)
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for side in ('left_margin', 'right_margin'): setattr(sec, side, Cm(2.0))
sec.top_margin = sec.bottom_margin = Cm(2.0)

SUB = re.compile(r'([A-Za-zΔρ]+)_(\([^)]*\)|[A-Za-zρσπ0-9]+)')
SUP = re.compile(r'\^([0-9+\-−]+)')


def runs(p, text, size=None, bold=False, italic=False):
    """Add text with X_sub and ^sup markup."""
    text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", "", str(text))
    pos = 0
    tokens = []
    for m in re.finditer(r'([A-Za-zΔρ]+)_(\([^)]*\)|[A-Za-zρσπ0-9]+)|\^([0-9+\-−]+)', text):
        tokens.append((m.start(), m.end(), m))
    for a, b, m in tokens:
        if a > pos:
            r = p.add_run(text[pos:a]); r.bold, r.italic = bold, italic
            if size: r.font.size = Pt(size)
        if m.group(3) is not None:
            r = p.add_run(m.group(3)); r.font.superscript = True
        else:
            r = p.add_run(m.group(1)); r.bold, r.italic = bold, italic
            if size: r.font.size = Pt(size)
            r = p.add_run(m.group(2).strip('()')); r.font.subscript = True
        r.bold = bold
        if size: r.font.size = Pt(size)
        pos = b
    if pos < len(text):
        r = p.add_run(text[pos:]); r.bold, r.italic = bold, italic
        if size: r.font.size = Pt(size)


def para(text, bold_lead=None, size=None, align=WD_ALIGN_PARAGRAPH.JUSTIFY, space=6):
    p = doc.add_paragraph(); p.alignment = align
    p.paragraph_format.space_after = Pt(space)
    if bold_lead:
        r = p.add_run(bold_lead); r.bold = True
        if size: r.font.size = Pt(size)
    runs(p, text, size)
    return p


def set_landscape(on):
    s = doc.add_section()
    if on:
        s.orientation = WD_ORIENT.LANDSCAPE; s.page_width, s.page_height = Cm(29.7), Cm(21.0)
    else:
        s.orientation = WD_ORIENT.PORTRAIT; s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    for side in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'): setattr(s, side, Cm(2.0))


def border(cell, top=False, bottom=False):
    tcPr = cell._tc.get_or_add_tcPr(); b = OxmlElement('w:tcBorders')
    for side, on in (('top', top), ('bottom', bottom)):
        e = OxmlElement(f'w:{side}'); e.set(qn('w:val'), 'single' if on else 'nil')
        e.set(qn('w:sz'), '6'); e.set(qn('w:color'), '000000'); b.append(e)
    tcPr.append(b)


def table(caption, header, rows, size=8, widths=None, note=None):
    p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(10); p.paragraph_format.keep_with_next = True
    lead, rest = caption.split('.', 1)
    r = p.add_run(lead + '.'); r.bold = True; r.font.size = Pt(10)
    runs(p, rest, 10)
    t = doc.add_table(rows=1 + len(rows), cols=len(header)); t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = widths is None
    for j, h in enumerate(header):
        c = t.rows[0].cells[j]; c.text = ''
        runs(c.paragraphs[0], str(h), size, bold=True); border(c, top=True, bottom=True)
    for i, row in enumerate(rows, 1):
        for j, v in enumerate(row):
            c = t.rows[i].cells[j]; c.text = ''
            runs(c.paragraphs[0], '' if v is None else str(v), size)
            border(c, bottom=(i == len(rows)))
    for row in t.rows:
        trPr = row._tr.get_or_add_trPr(); h = OxmlElement('w:cantSplit'); trPr.append(h)
        for c in row.cells:
            for pp in c.paragraphs:
                pp.paragraph_format.space_after = Pt(0); pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    hdr = t.rows[0]._tr.get_or_add_trPr(); e = OxmlElement('w:tblHeader'); hdr.append(e)
    if widths:
        for row in t.rows:
            for c, w in zip(row.cells, widths): c.width = Cm(w)
    if note:
        q = para(note, size=8.5); q.paragraph_format.space_before = Pt(3)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def f2(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '–'
    t = f'{x:.{d}f}'
    return t[1:] if t.startswith('-') and float(t) == 0 else t


def figure(path, caption, width=15):
    doc.add_picture(str(path), width=Cm(width)); doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = doc.add_paragraph(); lead, rest = caption.split('.', 1)
    r = p.add_run(lead + '.'); r.bold = True; r.font.size = Pt(10); runs(p, rest, 10)


# ------------------------------------------------------------------ title
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run('Supporting Information'); r.bold = True; r.font.size = Pt(14)
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run('Van der Waals-Like, Valence and Core Radii of the Atoms and Singly Charged Ions of Elements 1–96 from the Same Relativistic Calculations'); r.bold = True; r.font.size = Pt(13)
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.add_run('Asier Urriolabeitia and Jorge Echeverría*')
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run('Instituto de Síntesis Química y Catálisis Homogénea (ISQCH), Facultad de Ciencias, Universidad de Zaragoza, Pedro Cerbuna 12, 50009 Zaragoza, Spain. E-mail: jorge.echeverria@unizar.es'); r.font.size = Pt(10)

doc.add_heading('Contents', level=1)
CONTENTS = [
    'S1. Configuration-control protocol',
    'Table S1. Radii of the neutral atoms at the four levels of theory',
    'Table S2. Radii of the cations and anions, ionization energies and electron affinities',
    'Table S3. Solutions corrected by the configuration audit',
    'Table S4. Regression statistics against published sets of radii',
    'Table S5. Sensitivity of the ion radii to the computational choices',
    'Table S6. J-averaged configurations and radii in the alternative configurations',
    'Table S7. Core radii, ELF minima and ELF shell populations',
    'Table S8 and Figure S1. Comparison with the Dirac–Fock orbital radii of Desclaux',
    'Table S9. State-averaged CASSCF tests',
    'Table S10. Bond-length data set',
    'Table S11. Penetration indices of the systems of Table 7',
    'Table S12. Valence radii in configurations for bonding',
    'S2. Structure maps and Table S13',
    'S3. FHI-aims calculations, Table S14 and Table S15',
]
for c in CONTENTS:
    q = doc.add_paragraph(c); q.paragraph_format.space_after = Pt(0)

# ------------------------------------------------------------------ S1 protocol
doc.add_heading('S1. Configuration-control protocol', level=1)
for t in [
    'All calculations were run with ORCA 6.1.1 on a SLURM cluster. The inputs were generated from the approved configurations of each species (lowest experimental level from the NIST Atomic Spectra Database for neutral atoms and cations, and the compilation of Ning and Lu for anions) and include the occupation pattern, the spin multiplicity and the SCF settings. Every calculation was checked for normal termination, SCF convergence, spin contamination (⟨S²⟩), the integrated electron count of the spherically averaged density, and the presence of the Molden file used for the density analysis.',
    'The configuration of each solution was audited with Löwdin reduced orbital populations. The total s, p, d and f populations were compared with those of the intended configuration, and a deviation larger than ±0.15 e in any of them was flagged. Two kinds of deviations were found. Integer deviations correspond to an electron that has moved between two shells, for instance 4f^n6s^2 converging to 4f^(n−1)5d^16s^2, or to a wrong occupation of the d shell in the transition metals. Fractional deviations reveal symmetry-broken s–d or s–d–f hybrids, in which the unrestricted solution mixes orbitals of different angular momentum to lower its energy.',
    'Integer deviations were corrected with two complementary approaches. In the first one, a ΔSCF calculation with the maximum overlap method was started from the orbitals of the faulty solution, with the intended occupation imposed through the ALPHACONF and BETACONF keywords. In the second one, the two orbitals involved were exchanged by a 90° rotation (Rotate block) before the SCF was restarted. Fractional deviations were treated by restarting PBE0 from the Hartree–Fock orbitals of the same species, which are purer in angular momentum, with point-group symmetry enabled. This procedure removes the hybrids of the d-block atoms, but it fails for some f-block atoms, which are flagged in Table S3 and marked with an asterisk in Table 4 of the main text. When the SCF did not converge, a ladder of three increasingly conservative strategies was applied in sequence (level shift and first-order steps, a tight restart from the resulting orbitals, and a guess from the model potential), and open 4f and 5f shells were always converged with first-order steps only, a level shift of 0.40 E_h and up to 2000 cycles.',
    'A practical aspect of the orbital rotations deserves to be mentioned. When a level shift is applied, ORCA writes the occupied orbitals first and the virtual orbitals afterwards, which does not correspond to the energy order. Therefore, the orbital indices used in ALPHACONF, BETACONF and Rotate must follow the order in the orbital file, and not the energy order. The scripts provided in the density_radii package handle this automatically.',
    'After the corrections, the values of each radius were taken from the solution that restores the intended configuration or, when none does, from the one with the smallest deviation, and every value in the tables records which calculation supplied it. Of the 617 self-consistent solutions, 93 had converged to a configuration different from the intended one. The corrected solutions of the neutral atoms are listed in Table S3.',
]:
    para(t)

# ------------------------------------------------------------------ Table S1
sel = pd.read_csv(R / 'atomic_radii_selected.csv')
core_p = pd.read_csv(R / 'core_radii_pbe0_anorcc.csv'); core_p = core_p[core_p.charge == 0]
core_o = pd.read_csv(R / 'core_radii_other_levels.csv'); core_o = core_o[core_o.charge == 0]
core_p['method'], core_p['scheme'] = 'PBE0', 'anorcc'
core = pd.concat([core_p, core_o])
LV = [('PBE0', 'anorcc'), ('PBE0', 'x2c'), ('HF', 'anorcc'), ('HF', 'x2c')]
set_landscape(True)
rows = []
for z in range(1, 97):
    a = sel[(sel.Z == z) & (sel.method == 'PBE0') & (sel.scheme == 'anorcc')].iloc[0]
    row = [z, a.symbol]
    for m, sc in LV:
        s = sel[(sel.Z == z) & (sel.method == m) & (sel.scheme == sc)]
        c = core[(core.Z == z) & (core.method == m) & (core.scheme == sc)]
        if len(s):
            s = s.iloc[0]; cc = c.c_orb.iloc[0] if len(c) else np.nan
            row += [f2(s['v_rho_0.001']), f2(s.r_outer), f2(cc)]
        else:
            row += ['–'] * 3
    row += [f2(a['v_rho_0.0005']), f2(a['v_rho_0.002'])]
    rows.append(row)
hdr = ['Z', 'El']
for lab in ('PBE0/DKH2', 'PBE0/X2C', 'HF/DKH2', 'HF/X2C'):
    hdr += [f'v_ρ {lab}', f'r_outer {lab}', f'c_orb {lab}']
hdr += ['v_ρ(0.0005)', 'v_ρ(0.002)']
table('Table S1. Van der Waals-like radius v_ρ (0.001 e bohr^−3), valence radius r_outer and core radius c_orb (Å) of the neutral atoms at the four levels of theory, namely PBE0 and HF with the DKH2 Hamiltonian and the ANO-RCC-TZP basis (H–Cm), and with the X2C Hamiltonian and the x2c-TZVPPall basis (H–Rn). The last two columns give v_ρ at 0.0005 and 0.002 e bohr^−3 for the primary set (PBE0/DKH2).',
      hdr, rows, size=6.5)

# ------------------------------------------------------------------ Table S2
ions = pd.read_csv(R / 'ion_radii_selected.csv')
corei = pd.read_csv(R / 'core_radii_pbe0_anorcc.csv')
rows = []
for z in range(1, 97):
    cat = ions[(ions.Z == z) & (ions.charge == 1)]
    an = ions[(ions.Z == z) & (ions.charge == -1)]
    sym = sel[sel.Z == z].symbol.iloc[0]
    if not len(cat) and not len(an):
        continue
    row = [z, sym]
    if len(cat):
        c = cat.iloc[0]; cc = corei[(corei.Z == z) & (corei.charge == 1)]
        row += [c.configuration, f2(c['v_rho_0.001']), f2(c.r_outer) if c.n_electrons > 0 else '–',
                f2(cc.c_orb.iloc[0]) if len(cc) else '–', f2(c.IE_calc_eV)]
    else:
        row += ['–'] * 5
    if len(an):
        a = an.iloc[0]
        status = 'reliable' if bool(a['EA_within_0.3eV']) else ('bound' if bool(a.bound) else 'not bound')
        row += [a.configuration, f2(a['v_rho_0.001']), f2(a.r_outer), f2(a.EA_calc_eV), f2(a.EA_exp_eV), status]
    else:
        row += ['–'] * 6
    rows.append(row)
table('Table S2. Radii of the singly charged cations and anions (Å) from PBE0/DKH2/ANO-RCC-TZP calculations, together with the computed adiabatic ionization energies (IE) and electron affinities (EA, eV) and the experimental electron affinities. Anions are labelled reliable when the computed EA lies within 0.3 eV of experiment, bound when it is positive but outside this range, and not bound otherwise.',
      ['Z', 'El', 'M^+ configuration', 'v_ρ(M^+)', 'r_outer(M^+)', 'c_orb(M^+)', 'IE calc', 'X^− configuration', 'v_ρ(X^−)', 'r_outer(X^−)', 'EA calc', 'EA exp', 'X^− status'],
      rows, size=6.5)

# ------------------------------------------------------------------ Table S3
set_landscape(False)
ROUTE = {'fix:fblock_fixes2.csv': 'f block, second round', 'fix:fblock_fixes.csv': 'f block, first round',
         'fix:atomic_radii_fixes.csv': 'ΔSCF/rotation or HF restart', 'fix:atomic_radii_retryA.csv': 'SCF retry ladder'}
fx = sel[(sel.origin != 'original') | sel.flag.notna()].sort_values(['Z', 'method', 'scheme'])
rows = []
for _, a in fx.iterrows():
    lev = f"{a.method}/{'DKH2' if a.scheme == 'anorcc' else 'X2C'}"
    rows.append([a.Z, a.symbol, lev, a.configuration, ROUTE.get(a.origin, 'none'),
                 f2(a.max_pop_dev), 'yes' if a.config_ok else 'no', a.flag if isinstance(a.flag, str) else ''])
table('Table S3. Neutral-atom solutions that were replaced by a corrected calculation or remain flagged after the configuration audit. The maximum deviation (e) of the Löwdin s, p, d or f populations from the intended configuration and its acceptance (±0.15 e) refer to the final solution.',
      ['Z', 'El', 'Level', 'Intended configuration', 'Correction', 'Max. dev.', 'Accepted', 'Remaining flag'], rows, size=7.5)

# ------------------------------------------------------------------ Table S4, S5
ref = pd.read_csv(R / 'reference_set_comparison.csv')
rows = [[r.quantity, r.reference, '' if pd.isna(r.excluded) else r.excluded, r.n, f2(r.r2, 3), f2(r.slope, 3), f2(r.intercept, 3), f2(r.mean_diff, 3), f2(r.sd_diff, 3), f2(r.mad, 3)] for r in ref.itertuples()]
table('Table S4. Linear regression of the radii of this work (PBE0/DKH2) against published sets of radii (y = slope·x + intercept, with x the published radius), with the mean, standard deviation and mean absolute value of the differences (Å).',
      ['Radius', 'Reference', 'Excluded', 'n', 'r^2', 'Slope', 'Intercept', 'Mean diff.', 'SD diff.', 'MAD'], rows, size=7.5)
ion_s = pd.read_csv(R / 'ion_sensitivity.csv')
rows = [[r.ion, r.variant, r.n, f2(r.mean, 3), f2(r.mad, 3), f2(r.max, 3)] for r in ion_s.itertuples()]
table('Table S5. Changes of v_ρ of the ions (Å) with respect to the primary set (PBE0/DKH2/ANO-RCC-TZP) when the method, Hamiltonian or basis set is changed.',
      ['Ions', 'Variant', 'n', 'Mean', 'MAD', 'Max'], rows, size=8)
s3 = pd.read_csv(REV / 'sensitivity_three_radii.csv')
rows = [list(r) for r in s3.itertuples(index=False)]
table('Table S5b. Changes of the three radii of the neutral atoms (Å) with respect to the primary set, given as mean absolute deviation (maximum).',
      list(s3.columns), rows, size=8)

# ------------------------------------------------------------------ Table S6
jv = pd.read_csv(REV / 'SI_javg_configurations.csv')
rows = [list(r) for r in jv.itertuples(index=False)]
set_landscape(True)
table('Table S6. Species whose lowest J-averaged term belongs to a configuration different from that of the lowest experimental level, with the experimental and PBE0 energies of the J-averaged configuration relative to the primary one (eV) and the radii (Å) in both configurations.',
      list(jv.columns), rows, size=7)
js = pd.read_csv(REV / 'javg_summary.csv')
js = js[js.status != 'same']
rows = [[r.Z, r.species, r.status, r.ground_conf if isinstance(r.ground_conf, str) else '', r.ground_term if isinstance(r.ground_term, str) else '',
         r.alt_conf if isinstance(r.alt_conf, str) else '', r.alt_term if isinstance(r.alt_term, str) else '', f2(r.gap_eV, 3)] for r in js.itertuples()]
table('Table S6b. Spectra of the NIST database (neutral atoms, I, and singly charged cations, II) for which the J-averaged analysis does not confirm the configuration of the ground level. Reordered, a complete term of another configuration lies below the J-averaged ground term. Possible, no complete term proves a reordering, but a level of another configuration lies below the J-averaged ground term. Same (config. average), decided from complete configuration averages in jj-coupled spectra. Incomplete, no complete term or configuration. Ground level only, the database lists a single level. The gap is the energy of the lowest alternative term or level above the J-averaged ground term (eV). The other 160 spectra keep the configuration of the ground level.',
      ['Z', 'Spectrum', 'Status', 'Ground configuration', 'Term', 'Alternative configuration', 'Term', 'Gap'], rows, size=7)

# ------------------------------------------------------------------ Table S7
set_landscape(False)
cp = pd.read_csv(R / 'core_radii_pbe0_anorcc.csv'); cp = cp[(cp.charge == 0) & (cp.n_core > 0)]
rows = []
for r in cp.itertuples():
    rows.append([r.Z, r.symbol, int(r.n_core), f2(r.c_orb), r.c_orb_l if isinstance(r.c_orb_l, str) else '', f2(r.c_count),
                 f2(r.c_elf), f2(r.N_inside_c_elf, 1), f2(r.r_outer)])
table('Table S7. Core radii of the neutral atoms (Å, PBE0/DKH2), namely the orbital radius c_orb with the angular momentum of the outermost core orbital, the electron-count radius c_count, defined by 4π∫₀^c ρ̄(r)r^2dr = n_core, and the outermost ELF minimum inside r_outer, c_ELF, together with the number of electrons enclosed by c_ELF and the valence radius r_outer.',
      ['Z', 'El', 'n_core', 'c_orb', 'l', 'c_count', 'c_ELF', 'N(c_ELF)', 'r_outer'], rows, size=7.5)

# ------------------------------------------------------------------ Table S8, Figure S1
de = pd.read_csv(REV / 'desclaux_comparison.csv')
rows = [[r.Z, r.symbol, r.orbital, f2(r.r_DF_avg_A, 3), f2(r.r_outer_HF_A, 3), f2(r.dev_HF, 3), f2(r.r_outer_PBE0_A, 3), f2(r.dev_PBE0, 3)] for r in de.itertuples()]
table('Table S8. Most probable radius of the outermost orbital (Å) from the Dirac–Fock calculations of Desclaux, as the (2j + 1)-weighted average of the two j components for p and d orbitals, compared with r_outer from HF/DKH2 and PBE0/DKH2 calculations, and their deviations from the Dirac–Fock value.',
      ['Z', 'El', 'Orbital', 'r_DF', 'r_outer HF', 'Dev. HF', 'r_outer PBE0', 'Dev. PBE0'], rows, size=7.5)
figure(REV / 'FigS_desclaux.png', 'Figure S1. Valence radius r_outer from HF/DKH2 (and PBE0/DKH2) calculations against the Dirac–Fock radius of maximum radial density of the outermost orbital reported by Desclaux.', 13)

# ------------------------------------------------------------------ Table S9
cas = pd.read_csv(REV / 'casscf_comparison_open_shell.csv')
rows = [[r.symbol, r.configuration, int(r.nroots), f2(r['root_spread_cm-1'], 1), f2(r.v_rho_CAS, 3), f2(r.d_v_rho_HF, 3), f2(r.d_v_rho_PBE0, 3),
         f2(r.r_outer_CAS, 3), f2(r.d_r_outer_HF, 3), f2(r.d_r_outer_PBE0, 3), f2(r.c_orb_CAS, 3), f2(r.d_c_orb_HF, 3), f2(r.d_c_orb_PBE0, 3)]
        for _, r in cas.iterrows()]
casL = pd.read_csv(R / 'casscf_comparison.csv')
for r in casL.itertuples():
    rows.append([r.symbol, '', r.roots, '0', f2(r.v_rho_CAS, 3), '–', f2(r.v_rho_SD - r.v_rho_CAS, 3), f2(r.r_outer_CAS, 3), '–', f2(r.r_outer_SD - r.r_outer_CAS, 3), '–', '–', '–'])
table('Table S9. Radii (Å) from state-averaged CASSCF densities averaged over the components of the ground term, and differences (single determinant − CASSCF) for Hartree–Fock (HF) and PBE0 single determinants. The spread of the averaged roots is given in cm^−1. For La, Ce and Th the single-determinant values are from PBE0.',
      ['El', 'Configuration', 'Roots', 'Spread', 'v_ρ CAS', 'Δv_ρ HF', 'Δv_ρ PBE0', 'r_outer CAS', 'Δr_outer HF', 'Δr_outer PBE0', 'c_orb CAS', 'Δc_orb HF', 'Δc_orb PBE0'], rows, size=7)

# ------------------------------------------------------------------ Table S10
bl = pd.read_csv(REV / 'bond_length_test.csv')
rows = [[r.bond, r.species_and_ref if isinstance(r.species_and_ref, str) else r.species, r.source_type, f2(r.R_A, 3), f2(r.pred_router, 3), f2(r.pred_cordero, 3), 'yes' if r.in_range_Z_le_96 and not r.ion else 'no']
        for r in bl.itertuples()]
table('Table S10. Single-bond data set of Pyykkö and Atsumi used in Section 5.2.1, with the reference bond length R (Å), the bond lengths predicted from the sums of r_outer and of the covalent radii of Cordero et al. (Å), and whether the bond enters the test (neutral species of elements 1–96).',
      ['Bond', 'Species', 'Source', 'R', 'r_outer(A) + r_outer(B)', 'Cordero', 'Used'], rows, size=7)

# ------------------------------------------------------------------ Table S11
pi = pd.read_csv(REV / 'penetration_index_sets_p.csv')
set_landscape(True)
rows = [[r.set, r.system, r.A, r.B, f2(r.d, 3), f2(r.y, 3), r.ylabel, f2(r.p_paper, 1), f2(r['p_v_rho+Cordero'], 1), f2(r['p_Rahm+Cordero'], 1), f2(r['p_Charry+Cordero'], 1), f2(r['p_v_rho+r_outer'], 1)] for _, r in pi.iterrows()]
table('Table S11. Penetration indices (%) of the systems of Table 7 computed with the published radii (Alvarez + Cordero, as reported by Pardos et al.), with v_ρ + Cordero, with the radii of Rahm et al. + Cordero, with the radii of Charry and Tkatchenko + Cordero, and with v_ρ + r_outer, together with the interatomic distance d (Å) and the covalency indicator y.',
      ['Set', 'System', 'A', 'B', 'd', 'y', 'Indicator', 'Alvarez + Cordero', 'v_ρ + Cordero', 'Rahm + Cordero', 'Charry + Cordero', 'v_ρ + r_outer'], rows, size=6.5)

# ------------------------------------------------------------------ Table S12
set_landscape(False)
bc = pd.read_csv(REV / 'bonding_configurations.csv')
rows = [[r.symbol, r.ground_config, r.bonding_config, f2(r.dE_eV_PBE0), f2(r.v_rho_ground), f2(r.v_rho_bond), f2(r.r_outer_ground), f2(r.r_outer_bond), f2(r.c_orb_ground), f2(r.c_orb_bond), f2(r.r_cov_Cordero), f2(r.r_cov_Pyykko)] for r in bc.itertuples()]
table('Table S12. Radii (Å) of Pd, Pr and Nd in the ground configuration and in the configuration used for bonding, with the PBE0 excitation energy (eV) and the covalent radii of Cordero et al. and of Pyykkö and Atsumi. For Ce, the tetravalent 5d^26s^2 configuration could not be stabilized in a single determinant.',
      ['El', 'Ground', 'Bonding', 'ΔE', 'v_ρ ground', 'v_ρ bond', 'r_outer ground', 'r_outer bond', 'c_orb ground', 'c_orb bond', 'Cordero', 'Pyykkö'], rows, size=7)

# ------------------------------------------------------------------ S2 structure maps
doc.add_heading('S2. Structure maps', level=1)
para('The valence s and p orbital radii were extracted from the PBE0/DKH2/ANO-RCC-TZP Molden files of the neutral atoms as the outermost maxima of r^2|φ|^2 of the outermost occupied s and p orbitals of majority spin. For the elements of groups 1, 2, 11 and 12, whose valence p shell is empty, r_p is the most probable radius of the lowest virtual p orbital of majority spin. The radii of Zunger were transcribed from Table I of Phys. Rev. B 1980, 22, 5839 and converted from bohr to Å. The structure-map coordinates are r_σ = |r_p(A) + r_s(A) − r_p(B) − r_s(B)| and r_π = |r_p(A) − r_s(A)| + |r_p(B) − r_s(B)|. Both sets of coordinates were standardized before the classification and the kernel ridge regression. The linear boundary of Figure 6 is that of a logistic regression with weak regularization (C = 10^4), and the leave-one-out validation used a linear support-vector classifier (C = 100). The kernel ridge regression used a Gaussian kernel, with the regularization parameter (10^−6 to 10^−1) and the kernel width (γ = 0.01 to 100) chosen to minimize the error of 150 random partitions into 90% training and 10% test compounds. The two-dimensional descriptor of Ghiringhelli et al. was evaluated with our radii and with our PBE0 ionization energies of the cations and electron affinities of the anions (EA = −0.07 eV for N, which is not bound).')
rc = pd.read_csv(REV / 'structure_maps' / 'radii_comparison_octet.csv')
rows = [[r.symbol, f2(r.rs_ours, 3), f2(r.rp_ours, 3), f2(r.rs_Zunger, 3), f2(r.rp_Zunger, 3)] for r in rc.itertuples()]
table('Table S13a. Valence s and p orbital radii (Å) of the 34 elements of the octet binary compounds, from this work and from Zunger.',
      ['El', 'r_s this work', 'r_p this work', 'r_s Zunger', 'r_p Zunger'], rows, size=8)
oc = pd.read_csv(REV / 'structure_maps' / 'octet_ghiringhelli2015.csv')
ZS = {3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 19: 'K', 20: 'Ca', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 37: 'Rb', 38: 'Sr', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 55: 'Cs', 56: 'Ba'}
R2 = rc.set_index('symbol')
rows = []
for r in oc.itertuples():
    A, B = ZS[r.ZA], ZS[r.ZB]
    def co(sfx):
        a, b = R2.loc[A], R2.loc[B]
        rs = abs(a['rp_' + sfx] + a['rs_' + sfx] - b['rp_' + sfx] - b['rs_' + sfx])
        rp = abs(a['rp_' + sfx] - a['rs_' + sfx]) + abs(b['rp_' + sfx] - b['rs_' + sfx])
        return rs, rp
    o, z = co('ours'), co('Zunger')
    rows.append([A + B, f2(r.dE_RS_ZB_eV, 3), 'ZB' if r.dE_RS_ZB_eV > 0 else 'RS', f2(o[0], 3), f2(o[1], 3), f2(z[0], 3), f2(z[1], 3)])
table('Table S13b. Energy difference ΔE = E(RS) − E(ZB) (eV, LDA, from Ghiringhelli et al.), more stable structure and structure-map coordinates (Å) of the 82 octet binary compounds computed with the radii of this work and with those of Zunger.',
      ['AB', 'ΔE', 'Stable', 'r_σ this work', 'r_π this work', 'r_σ Zunger', 'r_π Zunger'], rows, size=7.5)

# ------------------------------------------------------------------ S3 FHI-aims
doc.add_heading('S3. FHI-aims calculations', level=1)
para('The equations of state were computed with FHI-aims (version 260331) with the PBE functional, the atomic ZORA scalar-relativistic treatment, tight species defaults, no spin polarization, Gaussian smearing of 0.01 eV and the analytical stress tensor. The SCF convergence criteria were 10^−6 for the density, 10^−4 eV for the sum of eigenvalues and 10^−6 eV for the total energy. The cells were the conventional primitive cells of the bcc, fcc and ideal hcp (c/a = √(8/3)) structures, with one atom (two for hcp), and 24 volumes were spaced logarithmically between 1.05 V₀ and x_min·V₀, with x_min = 0.15 for Mg and Al and 0.20–0.25 for the other elements. The k-point grids were set as n_i = ⌈|b_i|/0.12 Å^−1⌉, with at least 4 and at most 28 points along each reciprocal vector b_i (including the factor 2π). At the highest compressions and for a few volumes of Na and Cs, the overlap matrix of the basis became ill-conditioned, and the corresponding calculations were run with basis_threshold = 10^−5 and override_illconditioning. The pressures from the stress tensor agree with −dE/dV obtained from a cubic spline of E(ln V) within 0.2–4 GPa for all structures (0.2–1.3% above 100 GPa).')
eos = pd.read_csv(REV / 'compression' / 'eos_pv_fhiaims.csv')
from scipy.interpolate import CubicSpline
rows = []
for (el, st), g in eos.groupby(['element', 'structure']):
    g = g.sort_values('V_atom_A3')
    lnV, P = np.log(g.V_atom_A3.values), g.P_GPa.values
    cs = CubicSpline(lnV, P)
    from scipy.optimize import brentq
    x0 = brentq(cs, lnV[-4], lnV[-1]) if cs(lnV[-1]) < 0 else np.nan
    V0 = np.exp(x0); K0 = -cs(x0, 1)
    rows.append([el, st, f2(V0), f2(K0, 1), f2(g.V_atom_A3.min()), f2(g.P_GPa.max(), 0)])
table('Table S14a. Equilibrium volume V₀ (Å^3 per atom) and bulk modulus K₀ = −V(dP/dV) (GPa) obtained from a cubic spline of the PBE pressures against ln V, and the smallest volume and highest pressure computed for each element and structure.',
      ['El', 'Structure', 'V₀', 'K₀', 'V_min', 'P_max'], rows, size=8)
rows = [[r.element, r.structure, r.i, f2(r.V_atom_A3, 4), f2(r.E_eV_atom, 6), f2(r.P_GPa, 3)] for r in eos.itertuples()]
table('Table S14b. FHI-aims PBE volumes (Å^3 per atom), total energies (eV per atom) and pressures from the stress tensor (GPa) of the elemental metals.',
      ['El', 'Structure', 'Point', 'V', 'E', 'P'], rows, size=7)
cpnt = pd.read_csv(REV / 'compression' / 'compression_points.csv')
set_landscape(True)
rows = [[r.el, r.state, f2(r.P_GPa, 1), f2(r.V_A3, 2), r.source, f2(r.r_WS, 3), f2(r.r_outer, 3), f2(r.v_rho, 3), f2(r['r_WS/r_outer'], 3), f2(r['r_WS/v_rho'], 3)] for _, r in cpnt.iterrows()]
table('Table S15. Volumes per atom (Å^3) and Wigner–Seitz radii (Å) of the elemental metals at ambient pressure and at their phase transitions, with their ratios to the valence radius r_outer and to the van der Waals-like radius v_ρ. The cI16 phase of Na is listed for completeness but does not enter the averages of Section 5.1.2.',
      ['El', 'State', 'P (GPa)', 'V', 'Source', 'r_WS', 'r_outer', 'v_ρ', 'r_WS/r_outer', 'r_WS/v_ρ'], rows, size=7)

doc.save(OUT)
print('written', OUT)
