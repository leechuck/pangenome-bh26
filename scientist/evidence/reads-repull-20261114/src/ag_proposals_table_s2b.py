"""Truth-free table of builder outcomes per haplotype (decisions, unresolved-evidence counts, rules,
sequence hashes, timings, work counters), and task/gene status. Numbers, names and hashes only.
  python ag_proposals_table.py TAG"""
import json
import sys
from pathlib import Path

S2 = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/ag')   # s2b: only change (output root)
tag = sys.argv[1]
H = ['readset', 'donor', 'arm', 'gene', 'gene_status', 'failure_kind', 'gene_seconds', 'same_family', 'hap',
     'pass1_allele', 'family_rows', 'events', 'resolvable', 'n_patches', 'n_notes', 'n_blind', 'decision',
     'sha256_16', 'kinds', 'guard', 'ambig', 'still', 'alt', 'blind', 'untestable', 'substring', 'rules',
     'work_family_rows', 'work_events', 'work_align_calls', 's_setup', 's_evidence', 's_build', 's_decide']
rows, tasks = ['\t'.join(H)], ['readset\tdonor\tarm\tstatus\tshared_seconds\tverify_seconds\ttotal_seconds\tjob']
for t in sorted((S2 / 'work' / tag).glob('*/*/*')):
    rs, d, arm = t.parts[-3:]
    c = t / 'COMPLETE.json'
    f = t / 'FAILURE.json'
    st = 'complete' if c.exists() else 'failed' if f.exists() else 'missing'
    rec = json.loads((c if c.exists() else f).read_text()) if st != 'missing' else {}
    tasks.append('\t'.join(map(str, (rs, d, arm, st, rec.get('shared_seconds'), rec.get('verify_runtime_seconds'),
                                     round(rec.get('finished', 0) - rec.get('started', 0), 1) if rec else None,
                                     f"{rec.get('job_id')}_{rec.get('task_id')}"))))
    for gd in sorted(p for p in t.iterdir() if p.is_dir()):
        g = gd.name
        gc, gf = gd / 'COMPLETE.json', gd / 'FAILURE.json'
        gst = 'complete' if gc.exists() else 'failed' if gf.exists() else 'missing'
        grec = json.loads((gc if gc.exists() else gf).read_text()) if gst != 'missing' else {}
        pp = gd / 'proposals.json'
        if not pp.exists():
            rows.append('\t'.join(map(str, [rs, d, arm, g, gst, grec.get('failure_kind'), grec.get('seconds')] +
                                      [''] * (len(H) - 7))))
            continue
        j = json.loads(pp.read_text())
        w, s = j.get('work', {}), j.get('seconds', {})
        for h, x in enumerate(j['haplotypes']):
            if 'skipped' in x:
                continue
            u = x.get('unresolved', {})
            rows.append('\t'.join(map(str, (rs, d, arm, g, gst, grec.get('failure_kind'), grec.get('seconds'), j['same_family'],
                                            h, x['pass1_allele'], x['family_rows'], x['events'], x['resolvable'],
                                            x['n_patches'], len(x['notes']), x['n_blind'], x['decision'],
                                            (x.get('sha256') or '')[:16], ','.join(x.get('kinds', [])),
                                            u.get('guard'), u.get('ambig'), u.get('still'), u.get('alt'), u.get('blind'),
                                            u.get('untestable'), u.get('substring'), ','.join(x.get('unresolved_rules', [])),
                                            w.get('family_rows'), w.get('events'), w.get('align_calls'),
                                            s.get('setup'), s.get('evidence'), s.get('build'), s.get('decide')))))
(S2 / 'results' / f'haplotypes_{tag}.tsv').write_text('\n'.join(rows) + '\n')
(S2 / 'results' / f'tasks_{tag}.tsv').write_text('\n'.join(tasks) + '\n')
print(len(rows) - 1, 'haplotype rows', len(tasks) - 1, 'tasks')
