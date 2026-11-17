"""Shared helpers for stage thin-data-20261114 (read-only imports of the mosaic-typer-20261113 code).

Set PYTHONDONTWRITEBYTECODE=1 so that importing the other stage never writes into it.
"""
from __future__ import annotations
import hashlib, json, os, re, subprocess, sys
from pathlib import Path

sys.dont_write_bytecode = True
MT_SRC = '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src'
if MT_SRC not in sys.path:
    sys.path.insert(0, MT_SRC)
import common as C            # noqa: E402  (frozen shared library, read-only)
import run_donor as RD        # noqa: E402

TS = C.E / 'thin-data-20261114'
SEED = 20261114
LEVELS = {2: 8, 4: 4, 8: 2, 16: 1}           # denominator -> keep if hash mod 16 < value
ARMS_ALL = ('A', 'B', 'N', 'M')
ARM_DIR = {'A': 'genomic', 'B': 'genomic_core_hprc_asian', 'N': 'genomic_core_nonasian',
           'M': 'genomic_core_arab_saudi'}
N_LEVELS = (4, 16)                            # depth denominators at which arm N is run
SIF = RD.SIF
RUN_T1K = RD.RUN_T1K
THREADS = 4


def readsets():
    """{'D-real':{donor:...}, 'T-real':..., 'T-sim':...} -> list of donors, plus run folders."""
    runs = C.readset_runs()
    dsim = set(runs['D-sim'])
    return {'D-real': sorted(d for d in runs['real'] if d in dsim),
            'T-real': sorted(d for d in runs['real'] if d not in dsim),
            'T-sim': sorted(runs['T-sim'])}


def arm_run(rs, donor, arm):
    """Existing full-depth T1K run folder of an arm for a donor (or None)."""
    d = ARM_DIR[arm]
    if arm == 'A':
        base = {'T-sim': C.REGIONAL / 't1k/runs', 'T-real': C.ARMS_REAL / 'runs',
                'D-real': C.ARMS_REAL / 'runs'}[rs]
    elif rs == 'T-sim':
        base = C.AUGMENTED / 'runs'
    else:
        base = C.ARMS_REAL / 'runs'
        if arm == 'M':
            return None
    p = base / donor / d
    return p if (p / 'COMPLETE.json').exists() else None


def arm_info(rs, donor, arm):
    """Reference path and reads of the existing run of an arm (read from its receipt)."""
    run = arm_run(rs, donor, arm)
    if run is None:
        return None
    rec = C.load_json(run / 'COMPLETE.json')
    cmd = RD.parse_cmd(rec['command'])
    return dict(run=run, ref=RD.to_host(cmd['-f']), r1=RD.to_host(cmd['-1']), r2=RD.to_host(cmd['-2']),
                threads=int(cmd['-t']), ref_sha=rec.get('reference_sha256'),
                table_sha=rec.get('native_table_sha256'), status=rec['status'],
                seconds=rec['finished'] - rec['started'])


def levels_for_arm(arm):
    return {'A': (2, 4, 8, 16), 'B': (2, 4, 8, 16), 'M': (2, 4, 8, 16), 'N': N_LEVELS}[arm]


def read_key(name: str) -> int:
    """Stable hash of a read name (mate suffix and comment removed), seed 20261114."""
    n = name.split()[0]
    if n.endswith('/1') or n.endswith('/2'):
        n = n[:-2]
    h = hashlib.sha256(f'{SEED}:{n}'.encode()).digest()
    return int.from_bytes(h[:8], 'big') % 16


def t1k(r1, r2, ref, outdir: Path, threads: int, log: Path):
    """Same command as run_donor.t1k, with this stage bound as well (outputs live here)."""
    binds = ','.join(f'{v}:{k}' for k, v in RD.BIND.items()) + f',{TS}:{TS}'
    cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--bind', binds, str(SIF), RUN_T1K,
           '-1', RD.to_container(Path(r1)), '-2', RD.to_container(Path(r2)), '--preset', 'hla-wgs',
           '-f', RD.to_container(Path(ref)), '-t', str(threads),
           '--alleleDigitUnits', '4', '--alleleDelimiter', ':', '-o', RD.to_container(outdir / 'calls')]
    outdir.mkdir(parents=True, exist_ok=True)
    with open(log, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    t = outdir / 'calls_genotype.tsv'
    assert t.exists() and t.stat().st_size > 0, 'empty genotype table'
    return cmd


def log_stats(path: Path):
    """Truth-free quantities from a T1K run.log."""
    txt = Path(path).read_text()
    out = {}
    m = re.findall(r'Found (\d+) read fragments', txt)
    if m:
        out['candidate_fragments'] = int(m[0])
        out['analyzer_fragments'] = int(m[-1])
    m = re.findall(r'(\d+) read fragments can be assigned', txt)
    if m:
        out['assigned_fragments_analyzer'] = int(m[-1])
    return out
