"""Explain graph refinement outcomes on the fixed exposed development donors."""
import collections
import csv
import json
from pathlib import Path
from score_graph_pair_refinement import HERE, CONDITIONS, run
from score_linear_controls import sha, write_tsv


def audit():
    # Revalidate call hashes, baseline identity and preserved coarse calls first.
    run()
    root = HERE / 'development/graph-pair-refinement'
    scores = list(csv.DictReader((root / 'gene_scores.tsv').open(), delimiter='\t'))
    provenance = json.loads((root / 'provenance.json').read_text())
    native = {(r['donor'], r['gene']): r for r in scores
              if r['method'] == 'ipd_genome' and r['fields'] == '4'}
    decisions = {}
    for key in provenance['runs']:
        donor, condition = key.split('/', 1)
        if provenance['runs'][key]['status'] == 'complete':
            decisions[donor, condition] = json.loads(
                (HERE / 'work/graph-pair-refinement-snapshot' / key / 'decisions.json').read_text())
    rows = []
    for r in scores:
        key = r['donor'], r['method']
        if r['fields'] != '4' or r['method'] not in CONDITIONS or key not in decisions:
            continue
        base = native[r['donor'], r['gene']]
        decision = decisions[key][r['gene']]
        rows.append(dict(donor=r['donor'], condition=r['method'], gene=r['gene'],
                         eligible=int(r['eligible']), native_correct=int(base['correct']),
                         refined_correct=int(r['correct']), changed=int(decision['changed']),
                         reason=decision['reason'],
                         score_gain_over_native=decision.get('score_gain_over_native', ''),
                         score_gap=decision.get('score_gap', ''),
                         informative_fragments=decision.get('informative_fragments', '')))
    summary = {}
    for condition in CONDITIONS:
        subset = [r for r in rows if r['condition'] == condition]
        errors = [r for r in subset if r['eligible'] and not r['native_correct']]
        summary[condition] = dict(
            completed_donors=len({r['donor'] for r in subset}),
            eligible=sum(r['eligible'] for r in subset),
            changed_genes=sum(r['changed'] for r in subset),
            rescues=sum(r['refined_correct'] for r in errors),
            losses=sum(r['eligible'] and r['native_correct'] and not r['refined_correct'] for r in subset),
            residual_error_reasons=dict(collections.Counter(r['reason'] for r in errors
                                                           if not r['refined_correct'])),
            all_decision_reasons=dict(collections.Counter(r['reason'] for r in subset)))
    if rows:
        write_tsv(root / 'decision-audit.tsv', rows)
    report = dict(scope='Fixed exposed development cohort only; partial results are not cohort accuracy estimates',
                  conditions=summary, score_sha256=sha(root / 'gene_scores.tsv'),
                  provenance_sha256=sha(root / 'provenance.json'), driver_sha256=sha(Path(__file__)))
    (root / 'decision-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    audit()
