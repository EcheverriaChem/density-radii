"""Run make_config_fixes.py for the _javg jobs, using the alternative configurations.

Usage (same arguments as make_config_fixes.py):
    python workflow/fix_javg.py --results $SCRATCH/density-radii/javg_alt \
        --cluster-dir $SCRATCH/density-radii/javg_alt --outdir $SCRATCH/density-radii/javg_fix
For the f block add the f-block recipe, as for the primary set:
        --scf-extra "Shift Shift 0.40 ErrOff 0.0 end"
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow"))
import make_config_fixes as mcf  # noqa: E402

_JAVG = {(r["symbol"], int(r["charge"])): r
         for r in csv.DictReader(open(ROOT / "data" / "atomic_states_javg.csv"))}
mcf.load_states = lambda *a, **k: _JAVG

if __name__ == "__main__":
    mcf.main()
