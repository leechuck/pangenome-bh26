"""Library for stage dose-response-20261114 (step 3 of the 2026-11-14 route).

Fixed T1K, changing reference content only: IPD genomic alleles (common.GENOMIC_FA, bytes unchanged)
plus a panel of cohort gene sequences drawn from whole family components of a source population.
Everything here runs on IBEX under Slurm. Sequences never leave IBEX.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
MT = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113')
sys.path.insert(0, str(MT / 'src'))
import common as C  # noqa: E402  (read-only import of the frozen shared library)

E = Path('/encrypted/e3008/hohndor-dogohla-20261020')
STAGE = E / 'dose-response-20261114'
SEED = 20261114
GENES = C.GENES
NONASIAN_DIR = C.ZEN / 'nonasian-panel-control-20261018'
SIF = C.BENCH / 'containers/ubuntu24.sif'
BIND = {'/benchmark': C.ZEN, '/baseline': C.REGIONAL, '/control': C.HPRC_CONTROL,
        '/augmented': C.AUGMENTED}
RUN_T1K = '/benchmark/baselines/assets/t1k/run-t1k'

GROUPS = ('ME', 'EAS', 'SAS')
NONASIAN = 'NONASIAN'
MATCHED_N = (5, 10, 20, 40)
CROSS_N = (10, 40)


def group_of_population(pop: str) -> str:
    if pop.startswith('Arab') or pop.startswith('Saudi'):
        return 'ME'
    if pop in ('CDX', 'CHB', 'CHS', 'KHV') or 'JPT' in pop:
        return 'EAS'
    if pop in ('BEB', 'GIH', 'ITU', 'PJL', 'STU'):
        return 'SAS'
    raise ValueError('unknown population ' + pop)


def cell_names_for(group: str, with_nonasian: bool = True):
    """(cell name, source group, n requested) for a test donor of `group`."""
    out = [(f'match_n{n:02d}', group, n) for n in MATCHED_N]
    for src in GROUPS:
        if src != group:
            out += [(f'x_{src}_n{n:02d}', src, n) for n in CROSS_N]
    if with_nonasian:
        out += [(f'x_{NONASIAN}_n{n:02d}', NONASIAN, n) for n in CROSS_N]
    return out


def seed_of(donor: str, cell: str) -> int:
    return int(hashlib.sha256(f'{donor}|{cell}|{SEED}'.encode()).hexdigest(), 16)


def draw_panel(donor: str, cell: str, own_comp, comps: dict, n: int):
    """Whole components drawn at random without replacement, never the donor's own component,
    until at least n donors are drawn (or the source is exhausted). comps: {comp_key: [donors]}."""
    rng = random.Random(seed_of(donor, cell))
    keys = sorted(k for k in comps if k != own_comp)
    rng.shuffle(keys)
    drawn, used = [], []
    for k in keys:
        if len(drawn) >= n:
            break
        drawn.extend(sorted(comps[k]))
        used.append(k)
    return drawn, used


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def gene_body(seq: str, exons) -> str:
    return seq[exons[0][0]:exons[-1][1]]


def load_cohort():
    """Population map, group map, components, truth records by donor (this needs truth: IBEX only)."""
    folds = C.load_json(MT / 'inputs/FOLDS.json')
    pop = {d: v['population'] for d, v in folds['donors'].items()}
    grp = {d: group_of_population(p) for d, p in pop.items()}
    # connected components of FOLDS groups (union-find), via the frozen exclusion map
    excl = C.exclusion_map()
    comp_of = {d: '+'.join(sorted(excl[d])) for d in pop}
    return folds, pop, grp, comp_of


def components_by_group(grp, comp_of):
    """{source group: {comp_key: [member donors that belong to that group]}}"""
    out = {g: {} for g in GROUPS}
    for d, g in grp.items():
        out[g].setdefault(comp_of[d], []).append(d)
    return out


def load_nonasian():
    """Per non-Asian donor: list of dict(gene, haplotype, core_sequence, core_exons), rebuilt from the
    `nonasian` panel sidecars of the five folds of nonasian-panel-control-20261018 (the panels behind the
    genomic_core_nonasian arm). Core = sequence[2000:-2000], exons shifted by -2000 (build_core.py)."""
    FL = 2000
    donors = {}
    seen = set()
    for fold in range(5):
        base = NONASIAN_DIR / f'fold{fold}/references/nonasian'
        for g in GENES:
            seqs = {}
            name, parts = None, []
            for line in open(base / f'HLA-{g}.fa'):
                line = line.strip()
                if line.startswith('>'):
                    if name:
                        seqs[name] = ''.join(parts).upper()
                    name, parts = line[1:], []
                else:
                    parts.append(line)
            if name:
                seqs[name] = ''.join(parts).upper()
            for r in json.load(open(base / f'HLA-{g}.json')):
                full = seqs[r['path']]
                core = full[FL:-FL]
                ex = [[x - FL, y - FL] for x, y in r['exons']]
                for src in r['sources']:
                    d, hap, _ = src.split('#')
                    key = (d, hap, g)
                    if key in seen:
                        continue
                    seen.add(key)
                    donors.setdefault(d, []).append(dict(gene=g, haplotype=hap, core_sequence=core,
                                                         core_exons=ex, source_id=src))
    return donors


def ipd_sequence_set():
    s = set()
    for name, ints, seq in C.read_fasta_headers_and_seqs(C.GENOMIC_FA):
        s.add(seq)
    return s


def build_additions(recs, ipd_seqs):
    """recs: truth-like records (gene, core_sequence, core_exons). Returns (additions list, stats).
    Same rules as t1k_reference.build: skip missing exons, sequences with N, sequences already in IPD;
    dedupe by sequence; alias HLA-<gene>*PG<sha[:12]>."""
    adds, seen, stats = [], {}, dict(skipped_ipd=0, skipped_N=0, skipped_noexon=0, dup=0)
    for r in sorted(recs, key=lambda r: (r['gene'], r['core_sequence'])):
        seq = r['core_sequence'].upper()
        if not r['core_exons']:
            stats['skipped_noexon'] += 1
            continue
        if 'N' in seq:
            stats['skipped_N'] += 1
            continue
        if seq in ipd_seqs:
            stats['skipped_ipd'] += 1
            continue
        key = (r['gene'], seq)
        h = sha(seq)
        if key in seen:
            stats['dup'] += 1
            if seen[key]['exons'] != [list(e) for e in r['core_exons']]:
                raise ValueError('identical core with conflicting exon annotation')
            continue
        name = f"HLA-{r['gene']}*PG{h[:12]}"
        ex = [tuple(e) for e in r['core_exons']]
        seen[key] = dict(name=name, gene=r['gene'], sha=h, exons=[list(e) for e in ex],
                         body_sha=sha(gene_body(seq, ex)), length=len(seq),
                         fa=C.t1k_header(name, ex, len(seq)) + seq + '\n')
        adds.append(seen[key])
    names = [a['name'] for a in adds]
    if len(set(names)) != len(names):
        raise ValueError('PG alias collision')
    stats['added'] = len(adds)
    return adds, stats


# ---------------- T1K (copied from mosaic-typer run_donor.py: same container, preset, options)

def to_host(p: str) -> Path:
    for k, v in BIND.items():
        if p == k or p.startswith(k + '/'):
            return C.guard_path(Path(str(v) + p[len(k):]))
    return C.guard_path(Path(p))


def to_container(p: Path) -> str:
    s = str(p)
    for k, v in BIND.items():
        if s.startswith(str(v) + '/'):
            return k + s[len(str(v)):]
    return s


def parse_cmd(cmd):
    out = {}
    for i, t in enumerate(cmd):
        if t in ('-1', '-2', '-f', '-t', '-o'):
            out[t] = cmd[i + 1]
    return out


def verify_runtime():
    rt = C.load_json(C.AUGMENTED / 'T1K_RUNTIME.json')
    for a in rt['artifacts']:
        got = C.sha256_file(to_host(a['path']))
        assert got == a['sha256'], 'T1K runtime artifact changed: ' + a['path']
    return rt['version']


def t1k(r1: str, r2: str, ref: Path, out: Path, threads: int, log: Path):
    binds = ','.join(f'{v}:{k}' for k, v in BIND.items()) + f',{STAGE}:{STAGE}'
    cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--bind', binds, str(SIF), RUN_T1K,
           '-1', r1, '-2', r2, '--preset', 'hla-wgs', '-f', to_container(ref), '-t', str(threads),
           '--alleleDigitUnits', '4', '--alleleDelimiter', ':', '-o', to_container(out / 'calls')]
    with open(log, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    table = out / 'calls_genotype.tsv'
    assert table.exists() and table.stat().st_size > 0, 'empty genotype table'
    return cmd


def donor_runs(dsim_donors=None):
    """{donor: dict(readset, armA_run, r1, r2, threads)} for the 147 benchmark donors, read type by donor:
    ME -> T-sim (simulated), EAS/SAS -> real (T-real / D-real; never the simulated HPRC copies)."""
    folds, pop, grp, comp_of = load_cohort()
    runs = C.readset_runs()
    dsim = set(runs['D-sim'])
    out = {}
    for d in sorted(pop):
        g = grp[d]
        if g == 'ME':
            rs, run = 'T-sim', runs['T-sim'].get(d)
        else:
            rs = 'D-real' if d in dsim else 'T-real'
            run = runs['real'].get(d)
        if run is None:
            out[d] = dict(readset=rs, armA_run=None)
            continue
        rec = C.load_json(run / 'COMPLETE.json')
        cmd = parse_cmd(rec['command'])
        out[d] = dict(readset=rs, armA_run=str(run), r1=cmd['-1'], r2=cmd['-2'],
                      threads=int(cmd['-t']), armA_ref=cmd['-f'],
                      armA_table_sha256=rec['native_table_sha256'])
    return out


def armB_spec(readset, donor):
    base = {'T-sim': C.AUGMENTED / 'runs', 'T-real': C.ARMS_REAL / 'runs',
            'D-real': C.ARMS_REAL / 'runs'}[readset]
    run = base / donor / 'genomic_core_hprc_asian'
    rec = C.load_json(run / 'COMPLETE.json')
    cmd = parse_cmd(rec['command'])
    ref = to_host(cmd['-f'])
    assert C.sha256_file(ref) == rec['reference_sha256'], 'arm-B reference changed'
    return dict(run=str(run), ref=str(ref), table_sha256=rec['native_table_sha256'],
                r1=cmd['-1'], r2=cmd['-2'], threads=int(cmd['-t']))
