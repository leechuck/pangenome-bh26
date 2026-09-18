"""Audit additive versus replacement panels; summarize paired, completed outcomes."""
import collections
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def audit_membership(members, metadata, excluded):
    panels = {k: set(v) for k, v in members.items()}
    for name, paths in members.items():
        if len(paths) != len(panels[name]):
            raise ValueError('Duplicate panel member: ' + name)
        if panels[name] & excluded:
            raise ValueError('Excluded haplotype in ' + name)
        if panels[name] - metadata.keys():
            raise ValueError('Unknown panel member in ' + name)
    hprc, full, asian = (panels[k] for k in ('hprc', 'full', 'asian_matched'))
    if not hprc <= full or not asian <= full:
        raise ValueError('Full panel does not contain the comparator panels')
    added = full - hprc
    if any(metadata[p]['cohort'].startswith('HPRC') or metadata[p]['cohort'] == 'REF' for p in added):
        raise ValueError('Unexpected source among added haplotypes')
    return dict(counts={k: len(v) for k, v in panels.items()},
                additive_retains_all_hprc=True, added_haplotypes=len(added),
                added_sources=dict(collections.Counter(metadata[p]['cohort'] for p in added)),
                replacement_removes_hprc=len(hprc - asian), excluded_overlap=0)


def paired(rows, candidate, baseline, level, stratum='ALL'):
    subsets = []
    for method in (candidate, baseline):
        rs = [r for r in rows if r['method'] == method and r['level'] == level
              and (stratum == 'ALL' or r['stratum'] == stratum)]
        if not rs or any(r['state'] != 'complete' for r in rs):
            return None
        index = {(r['donor'], r['gene']): r for r in rs}
        if len(index) != len(rs):
            raise ValueError('Duplicate genotype row')
        subsets.append(index)
    a, b = subsets
    if a.keys() != b.keys() or any(a[k]['eligible'] != b[k]['eligible'] for k in a):
        raise ValueError('Unmatched truth denominators')
    keys = [k for k in a if int(a[k]['eligible'])]
    return dict(eligible=len(keys), candidate_correct=sum(int(a[k]['correct']) for k in keys),
                baseline_correct=sum(int(b[k]['correct']) for k in keys),
                gains=sum(int(a[k]['correct']) > int(b[k]['correct']) for k in keys),
                losses=sum(int(a[k]['correct']) < int(b[k]['correct']) for k in keys))


def main():
    repo = HERE.parents[1]
    metadata = {r['hap_id']: r for r in csv.DictReader(
        (repo / 'hla-typer/source/haplotype_donors.tsv').open(), delimiter='\t')}
    excluded = {r['hap_id'] for r in csv.DictReader((HERE / 'excluded_paths.tsv').open(), delimiter='\t')}
    audit = audit_membership(json.loads((HERE / 'panel-membership.json').read_text()), metadata, excluded)
    rows = list(csv.DictReader((HERE / 'analysis/gene_scores.tsv').open(), delimiter='\t'))
    report = ['# Additive and replacement panel comparisons', '',
              f"Full = HPRC/reference ({audit['counts']['hprc']} haplotypes) + "
              f"{audit['added_haplotypes']} additional haplotypes ({audit['added_sources']}). "
              'All HPRC members are retained; excluded donor/family paths are absent.', '',
              f"The size-matched Asian arm replaces {audit['replacement_removes_hprc']} HPRC haplotypes. "
              'It is a replacement experiment, not an additive experiment.', '',
              'Primary panel contrast: full versus HPRC. These exploratory paired comparisons change '
              'both upstream read collection and the DRB1 structural graph. They do not isolate graph alignment.', '',
              '| Contrast | Resolution | Stratum | Correct: candidate / baseline | Gains | Losses | Eligible |',
              '|---|---|---|---:|---:|---:|---:|']
    comparisons = []
    for arm in ('full', 'asian_matched'):
        for variant in ('DogoHLA', 'DogoHLA-no-graph'):
            for level in ('two_field', 'four_field'):
                for stratum in ('ALL', 'EAS', 'SAS', 'EUR', 'AFR'):
                    candidate, baseline = f'{arm}:{variant}', f'hprc:{variant}'
                    result = paired(rows, candidate, baseline, level, stratum)
                    comparisons.append(dict(candidate=candidate, baseline=baseline, level=level,
                                            stratum=stratum, result=result))
                    if result:
                        report.append(f"| {candidate} vs {baseline} | {level} | {stratum} | "
                                      f"{result['candidate_correct']} / {result['baseline_correct']} | "
                                      f"{result['gains']} | {result['losses']} | {result['eligible']} |")
    (HERE / 'analysis/panel_comparison.json').write_text(json.dumps(
        dict(membership=audit, comparisons=comparisons), indent=2) + '\n')
    (HERE / 'analysis/PANEL_COMPARISON.md').write_text('\n'.join(report) + '\n')


if __name__ == '__main__':
    main()
