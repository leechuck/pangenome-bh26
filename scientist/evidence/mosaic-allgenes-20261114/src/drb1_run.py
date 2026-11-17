"""DRB1 arm C3: T1K pass 1 (existing arm-A run) + local-version patch builder (drb1_patch2) +
pass-2 T1K + frozen acceptance rule. Panel = IPD genomic alleles only. No truth is read.

  python drb1_run.py READSET DONOR [--no-pass2]

Outputs under $S2/work/<readset>/<donor>/C3/ :
  proposals.json   decisions, patches, notes (no sequence)
  additions.fa     novel candidates (sequence; stays on IBEX)
  OVERRIDES.json   slots re-named to a known IPD allele whose gene body equals the patched allele
  pass2/, ACCEPT.json, COMPLETE.json | FAILURE.json
"""
from __future__ import annotations

import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drb1_lib as L          # noqa: E402
import drb1_patch2 as PT      # noqa: E402
import drb1_diag2 as D2       # noqa: E402
C, EV, PR = L.C, L.EV, L.PR
import run_donor as RD        # noqa: E402

S2 = D2.S2
G = 'DRB1'          # default gene; run() takes the gene as a parameter (stage mosaic-allgenes)
ARM = 'C3'
GUARD_KINDS = ('confounded_event', 'near_blind_region')
AMBIG_KINDS = ('two_supported_versions', 'pile_up')
# truth-free abstain rules (stage mosaic-allgenes, item 3): a haplotype is 'unresolved' under a rule
# when it still carries evidence of that kind; then its patched allele does not become the final call
# rule -> {status count: minimum that makes the haplotype unresolved}
ABSTAIN_RULES = {
    'R1': dict(guard=1, ambig=1),
    'R2': dict(guard=1, ambig=1, still=1, alt=1),
    'R3': dict(guard=1, ambig=1, still=1, alt=1, blind=1, untestable=1),
    # R4 was added after looking at the development outcome of R1-R3 (post hoc): a guard rejection,
    # or two or more family differences longer than a read can test
    'R4': dict(guard=1, untestable=2),
}


def unresolved_status(events, notes, keep, blind, max_len=120):
    """Counts of unresolved evidence for one haplotype, from statuses the builder already produces."""
    kinds = [n['kind'] for n in notes]
    return dict(guard=sum(k in GUARD_KINDS for k in kinds), ambig=sum(k in AMBIG_KINDS for k in kinds),
                still=kinds.count('own_version_still_supported'),
                alt=sum(p['kind'] != 'own_version_unsupported' for p in keep),
                blind=len(blind),
                untestable=sum(1 for e in events if not e['resolvable'] and e['cur'] is not None
                               and len(e['alleles']) >= 2 and e['longest'] > max_len),
                substring=sum(1 for e in events if not e['resolvable'] and e['cur'] is not None
                              and len(e['alleles']) >= 2 and e['longest'] <= max_len))
DEV = ('D-sim', 'D-real')
PARAMS = dict(max_edit=2, sup_share=0.25, min_sup=3, cur_ratio=0.25, het_share=0.3, pile=2.5,
              max_indel_cols=150, hom_ratio=1.5, blind_min_run=100, blind_near=150, blind_max_nm=2, n_max=4,
              family_theta=0.01)


def t1k(r1, r2, ref: Path, out: Path, threads: int, log: Path):
    """Same command as the frozen stage's run_donor.t1k, with this stage's folder bound too."""
    import subprocess
    binds = ','.join(f'{v}:{k}' for k, v in RD.BIND.items()) + f',{C.S}:{C.S},{S2}:{S2}'
    cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--bind', binds, str(RD.SIF), RD.RUN_T1K,
           '-1', r1, '-2', r2, '--preset', 'hla-wgs', '-f', RD.to_container(ref), '-t', str(threads),
           '--alleleDigitUnits', '4', '--alleleDelimiter', ':', '-o', RD.to_container(out / 'calls')]
    with open(log, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    table = out / 'calls_genotype.tsv'
    assert table.exists() and table.stat().st_size > 0, 'empty genotype table'
    return cmd


def run(readset, donor, out: Path, do_pass2=True, arm='C', gene=G, shared=None, confounded_path=None):
    G = gene
    tm = {}
    t0 = time.time()
    PT.WORK.clear()
    X = D2.setup(readset, donor, out, arm, gene=G, shared=shared)
    tm['setup'] = round(time.time() - t0, 2)
    panel, pu, msa, rows = X['panel'], X['pu'], X['msa'], X['rows']
    sa, sb = X['sa'], X['sb']
    frz = PR.Params(**C.load_json(RD.FROZEN_PARAMS))
    proj = L.project_all(pu, X['rd_u'], 6)
    ev1 = L.make_evidence(panel.C, proj, 1, X['bad_q'], True)
    D = L.tables_quant(ev1)[3]['D']
    reads = [x for x in X['rd'] if x['qname'] not in X['bad_q']]
    p1 = (sa[0], sb[0])
    slot_rows = (set(sa), set(sb))
    fams = (X['fam_a'], X['fam_b'])
    # reference alleles by gene body (first-exon start to last-exon end) and by full sequence
    body_names, full_names = {}, {}
    for i, r in enumerate(rows):
        ex = msa.exons(r)
        if ex:
            body_names.setdefault(panel.seq[i][ex[0][0]:ex[-1][1]], []).extend(X['all_names'][i])
        full_names.setdefault(panel.seq[i], []).extend(X['all_names'][i])
    tm['evidence'] = round(time.time() - t0 - tm['setup'], 2)
    rec = dict(donor=donor, readset=readset, gene=G, arm=arm + '3', panel_arm=arm, params=PARAMS, expected_depth=D,
               pass1=[X['names'][p1[0]], X['names'][p1[1]]], same_family=X['same'],
               reads_assigned=X['n_assigned'], reads_used=len(reads), haplotypes=[])
    additions, overrides, seen = [], [], set()
    homozygous = p1[0] == p1[1]
    cf = Path(confounded_path) if confounded_path else S2 / f'inputs/confounded_{G}.json'
    confounded = [tuple(w) for w in C.load_json(cf)['windows']] if cf.exists() else []
    rec['confounded_windows_loaded'] = len(confounded)
    rec['confounded_path'] = str(cf) if cf.exists() else None
    built = {}
    for h in (0, 1):
        if h == 1 and homozygous:
            continue
        own = slot_rows[h] | (slot_rows[1 - h] if X['same'] else set())
        spans = [(x['cs'], x['ce']) for x in reads
                 if x['prefer_row'] in own and x['nm'] <= PARAMS['blind_max_nm'] and x['clip'] <= 5]
        blind, _ = L.blind_intervals(panel, p1[h], spans, PARAMS['blind_min_run'])
        kw = {k: PARAMS[k] for k in ('max_edit', 'sup_share', 'min_sup', 'cur_ratio', 'het_share', 'pile',
                                     'max_indel_cols')}
        Hh = PT.build_haplotype(panel, fams[h], p1[h], p1[1 - h], reads, D, **kw)
        keep, notes = [], list(Hh['notes'])
        for p in Hh['patches']:
            if any(p['cl'] <= b_ and p['cr'] >= a_ for a_, b_ in confounded):
                notes.append(dict({k: v for k, v in p.items() if k != 'new'}, kind='confounded_event'))
            elif L.near_blind(panel, p1[h], p['cl'], p['cr'], blind, PARAMS['blind_near']):
                notes.append(dict({k: v for k, v in p.items() if k != 'new'}, kind='near_blind_region'))
            else:
                keep.append(p)
        built[h] = dict(Hh=Hh, keep=keep, notes=notes, blind=blind)
    tm['build'] = round(time.time() - t0 - tm['setup'] - tm['evidence'], 2)
    abst = {r: [] for r in ABSTAIN_RULES}
    # same family, two different pass-1 alleles: a version that both haplotypes would take but that
    # has reads for only one copy belongs to one of them, and short reads cannot say which
    if X['same'] and not homozygous:
        key = lambda p: (p['cl'], p['cr'], p['new_sha'])
        both = {key(p) for p in built[0]['keep']} & {key(p) for p in built[1]['keep']}
        for h in (0, 1):
            for p in built[h]['keep']:
                if key(p) in both and p['kind'] == 'own_version_unsupported' and \
                        p['n_new'] < PARAMS['hom_ratio'] * max(p['exp_h'], 1.0):
                    p['kind'] = 'same_family_phase_unknown'
    for h in (0, 1):
        if h == 1 and homozygous:
            rec['haplotypes'].append(dict(skipped='same pass-1 allele as haplotype 0'))
            continue
        Hh, keep, notes, blind = built[h]['Hh'], built[h]['keep'], built[h]['notes'], built[h]['blind']
        fp, events = Hh['fp'], Hh['events']
        path = PT.patched_path(fp, Hh['ra'], keep, events)
        info = dict(pass1_allele=X['names'][p1[h]], family_rows=len(Hh['fam']), events=len(events),
                    resolvable=sum(e['resolvable'] for e in events), blind_regions=blind[:20],
                    n_blind=len(blind), n_patches=len(keep),
                    patches=[{k: v for k, v in p.items() if k != 'new'} for p in keep],
                    notes=[{k: v for k, v in n.items() if k != 'new'} for n in notes],
                    decision='no_patch', pass1_msa_row=int(X['rows_msa'][p1[h]]))
        info['unresolved'] = unresolved_status(events, notes, keep, blind)
        info['unresolved_rules'] = sorted(r for r, ks in ABSTAIN_RULES.items()
                                          if any(info['unresolved'][k] >= m for k, m in ks.items()))
        n_ov = len(overrides)
        if keep:
            exr = {i: msa.exons(rows[r]) for i, r in enumerate(Hh['fam'])}
            seq = PR.path_sequence(fp, path)
            exo = PR.path_exons(fp, path, exr)
            sha = C.sha256_str(seq)
            sp = PT.spell(fp, path.codes(fp.P), frz)
            info.update(sha256=sha, length=len(seq), spell_mismatch=sp[0], spell_switches=sp[1],
                        kinds=sorted({p['kind'] for p in keep}))
            only_own = all(p['kind'] == 'own_version_unsupported' for p in keep)
            if exo is None:
                info['decision'] = 'exon_annotation_inconsistent'
            elif seq in full_names or seq[exo[0][0]:exo[-1][1]] in body_names:
                names = full_names.get(seq) or body_names[seq[exo[0][0]:exo[-1][1]]]
                names = sorted(set(names))
                if X['names'][p1[h]] in names:
                    info['decision'] = 'equals_pass1_allele'
                elif only_own:
                    info['decision'] = 'rename_to_known_allele'
                    info['rename_to'] = names[:10]
                    overrides.append(dict(kind='rename', pass1_label=X['names'][p1[h]], new_labels=names[:10],
                                          also_other_slot=homozygous))
                else:
                    info['decision'] = 'known_allele_but_phase_unknown'
            elif sha in seen:
                info['decision'] = 'duplicate_of_other_haplotype'
                overrides.append(dict(kind='novel', pass1_label=X['names'][p1[h]],
                                      new_labels=[X['names'][p1[h]], f'HLA-{G}*MOS_{sha[:12]}'],
                                      also_other_slot=False))
            elif len(additions) >= PARAMS['n_max']:
                info['decision'] = 'over_candidate_cap'
            else:
                seen.add(sha)
                name = f'HLA-{G}*MOS_{sha[:12]}'
                info.update(decision='proposed', name=name, exons=exo)
                additions.append(C.t1k_header(name, exo, len(seq)) + seq + '\n')
                # direct call: the patched allele replaces the pass-1 allele of this slot; when the
                # haplotype carrying a version is unknown, the named allele stays the call
                overrides.append(dict(kind='novel', pass1_label=X['names'][p1[h]],
                                      new_labels=[name] if only_own else [X['names'][p1[h]], name],
                                      also_other_slot=False))
        # overrides this haplotype produced (index range), so that abstain rules can drop them
        for o in overrides[n_ov:]:
            o['haplotype'] = h
            o['unresolved_rules'] = info['unresolved_rules']
        rec['haplotypes'].append(info)
    tm['decide'] = round(time.time() - t0 - sum(tm.values()), 2)
    rec['seconds'] = tm
    rec['work'] = dict(PT.WORK)
    out.mkdir(parents=True, exist_ok=True)
    C.save_json(out / 'proposals.json', rec)
    C.save_json(out / 'OVERRIDES.json', dict(gene=G, overrides=[{k: v for k, v in o.items()
                                                                  if k not in ('haplotype', 'unresolved_rules')}
                                                                 for o in overrides]))
    # abstain variants: 'Rk'  drops novel AND rename overrides of unresolved haplotypes (pass-1 name kept);
    #                   'Rkn' drops only novel final calls (renames to a known allele are kept)
    for r in ABSTAIN_RULES:
        for suf, kinds in (('', ('novel', 'rename')), ('n', ('novel',))):
            kept = [{k: v for k, v in o.items() if k not in ('haplotype', 'unresolved_rules')} for o in overrides
                    if not (r in o['unresolved_rules'] and o['kind'] in kinds
                            and (o['kind'] == 'rename' or o['new_labels'] == [o['new_labels'][-1]]
                                 and '*MOS_' in o['new_labels'][-1]))]
            C.save_json(out / f'OVERRIDES_{r}{suf}.json', dict(gene=G, rule=r + suf, overrides=kept))
    (out / 'additions.fa').write_text(''.join(additions))
    spec = X['spec']
    table = C.verified_table(spec['run']) if readset not in C.SYN_SETS else spec['run'] / 'calls_genotype.tsv'
    res = dict(n_candidates=len(additions), n_overrides=len(overrides), pass1_run=str(spec['run']),
               pass1_reference=str(spec['ref']))
    if not additions or not do_pass2:
        res.update(pass2='not_run_no_candidates' if not additions else 'not_run_proposal_only',
                   calls_table=str(table), calls_table_sha256=C.sha256_file(table))
        return res
    ref2 = out / 'pass2' / 'reference.fa'
    ref2.parent.mkdir(parents=True, exist_ok=True)
    with open(ref2, 'wb') as f:
        with open(spec['ref'], 'rb') as g0:
            shutil.copyfileobj(g0, f)
        f.write(''.join(additions).encode())
    ref_sha = C.sha256_file(ref2)
    cmd = t1k(spec['r1'], spec['r2'], ref2, out / 'pass2', spec['threads'], out / 'pass2/run.log')
    ref2.unlink()
    table2 = out / 'pass2/calls_genotype.tsv'
    res.update(pass2='complete', pass2_reference_sha256=ref_sha, command=cmd,
               calls_table=str(table2), calls_table_sha256=C.sha256_file(table2))
    import accept
    acc = accept.run(spec['run'], spec['ref'], table2, out / 'additions.fa', out, frz)
    res.update(mos_final_calls=sum(len(v) for v in acc.values()),
               mos_rejected=sum(1 for v in acc.values() for x in v if not x['accepted']))
    return res


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    do_pass2 = '--no-pass2' not in sys.argv
    arm = 'D' if '--arm-D' in sys.argv else 'C'
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    assert readset in DEV or (S2 / 'FREEZE.json').exists(), 'test read sets only after FREEZE.json'
    assert not (arm == 'D' and readset == 'D-sim'), 'no cohort-panel run exists for D-sim'
    RD.verify_runtime()
    sub = arm + '3' if do_pass2 else arm + '3_P'
    out = S2 / 'work' / readset / donor / sub
    if (out / 'COMPLETE.json').exists():
        return
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rec = dict(status='running', readset=readset, donor=donor, arm=sub, job_id=os.environ['SLURM_JOB_ID'],
               code_sha256={p.name: C.sha256_file(p) for p in sorted(Path(__file__).resolve().parent.glob('*.py'))},
               started=time.time())
    try:
        rec.update(run(readset, donor, out, do_pass2, arm), status='complete', finished=time.time())
        C.save_json(out / 'COMPLETE.json', rec)
    except BaseException as e:
        rec.update(status='failed', error=repr(e), traceback=traceback.format_exc(), finished=time.time())
        C.save_json(out / 'FAILURE.json', rec)
        raise


if __name__ == '__main__':
    main()
