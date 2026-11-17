"""Run a copied step-4 script (ag_proposals_table_s2b.py / ag_test_summary_s2b.py: frozen code, output root
changed only) with the verified frozen step-4 src on sys.path.   python runf.py SCRIPT ARGS..."""
import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen4               # noqa: E402
frozen4.use()
script = str(Path(__file__).resolve().parent / sys.argv[1])
sys.argv = [script] + sys.argv[2:]
runpy.run_path(script, run_name='__main__')
