#!/usr/bin/env python3
"""
Download NIST ASD energy levels (CSV) for the neutral atoms (I) and singly
charged cations (II) of elements 1-96.

Usage examples
  python3 fetch_nist_levels.py                 # all spectra, into ./nist_levels
  python3 fetch_nist_levels.py --only Fe Ni Pd # a subset
  python3 fetch_nist_levels.py --dry-run       # print the URLs only
  python3 fetch_nist_levels.py --force         # re-download existing files

Each spectrum is saved as <outdir>/<Symbol>_<I|II>.csv, and a summary of every
request is written to <outdir>/download_log.csv. Files that already exist are
skipped, so the script can be re-run after an interruption.
Only the Python standard library is used.
"""
import argparse
import csv
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm",
]

BASE_URL = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"

# Parameters taken from the Fe I query. units=0 -> cm-1, format=2 -> CSV,
# output=0 -> whole table (no paging), multiplet_ordered=0 -> term ordered.
PARAMS = {
    "de": "0",
    "units": "0",
    "format": "2",
    "output": "0",
    "page_size": "15",
    "multiplet_ordered": "0",
    "conf_out": "on",
    "term_out": "on",
    "level_out": "on",
    "unc_out": "1",
    "j_out": "on",
    "g_out": "on",
    "lande_out": "on",
    "perc_out": "on",
    "biblio": "on",
    "temp": "",
    "submit": "Retrieve Data",
}

HEADERS = {"User-Agent": "Mozilla/5.0 (density-radii project; research use)"}


def spectra(only=None):
    """Yield (Z, symbol, charge label). H II has no electrons and is skipped."""
    for z, sym in enumerate(ELEMENTS, start=1):
        if only and sym not in only:
            continue
        yield z, sym, "I"
        if sym != "H":
            yield z, sym, "II"


def build_url(sym, ion):
    params = dict(PARAMS, spectrum=f"{sym} {ion}")
    return BASE_URL + "?" + urllib.parse.urlencode(params)


def looks_like_csv(text):
    """NIST returns an HTML page (not CSV) when there is no data or an error."""
    head = text.lstrip()[:500].lower()
    if head.startswith("<") or "<html" in head or "<!doctype" in head:
        return False
    return "j" in head and "level" in head


def fetch(url, retries=3, timeout=60):
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError) as err:
            last_err = err
            time.sleep(5 * attempt)
    raise RuntimeError(f"failed after {retries} attempts: {last_err}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="nist_levels")
    ap.add_argument("--only", nargs="+", metavar="SYMBOL",
                    help="restrict to these element symbols")
    ap.add_argument("--delay", type=float, default=2.0,
                    help="seconds between requests (be gentle with the NIST server)")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    ap.add_argument("--dry-run", action="store_true", help="print URLs and exit")
    args = ap.parse_args()

    if args.only:
        unknown = set(args.only) - set(ELEMENTS)
        if unknown:
            sys.exit(f"Unknown symbols: {', '.join(sorted(unknown))}")

    todo = list(spectra(args.only))
    if args.dry_run:
        for z, sym, ion in todo:
            print(f"{sym} {ion}\t{build_url(sym, ion)}")
        print(f"{len(todo)} spectra")
        return

    os.makedirs(args.outdir, exist_ok=True)
    log_path = os.path.join(args.outdir, "download_log.csv")
    log_rows = []
    n_ok = n_skip = n_fail = 0

    for i, (z, sym, ion) in enumerate(todo, start=1):
        fname = os.path.join(args.outdir, f"{sym}_{ion}.csv")
        tag = f"[{i:3d}/{len(todo)}] {sym} {ion}"
        if os.path.exists(fname) and not args.force:
            print(f"{tag}: exists, skipped")
            n_skip += 1
            log_rows.append([z, sym, ion, "skipped", fname, ""])
            continue
        url = build_url(sym, ion)
        try:
            text = fetch(url)
        except RuntimeError as err:
            print(f"{tag}: ERROR {err}")
            n_fail += 1
            log_rows.append([z, sym, ion, "error", "", str(err)])
            time.sleep(args.delay)
            continue

        if looks_like_csv(text):
            with open(fname, "w", encoding="utf-8", newline="") as fh:
                fh.write(text)
            n_lines = max(text.count("\n") - 1, 0)
            print(f"{tag}: OK ({n_lines} levels)")
            n_ok += 1
            log_rows.append([z, sym, ion, "ok", fname, n_lines])
        else:
            # keep the HTML answer to see what NIST said (no data, error, ...)
            bad = os.path.join(args.outdir, f"{sym}_{ion}_NOT_CSV.html")
            with open(bad, "w", encoding="utf-8") as fh:
                fh.write(text)
            print(f"{tag}: no CSV returned, saved {bad}")
            n_fail += 1
            log_rows.append([z, sym, ion, "not_csv", bad, ""])
        time.sleep(args.delay)

    with open(log_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Z", "symbol", "spectrum", "status", "file", "n_levels_or_error"])
        w.writerows(log_rows)

    print(f"\nDone: {n_ok} downloaded, {n_skip} skipped, {n_fail} failed. Log: {log_path}")
    if n_fail:
        print("Re-run the same command to retry the failed spectra.")


if __name__ == "__main__":
    main()
