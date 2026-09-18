"""Reserve assembly validation families using metadata and prior exposure IDs only.

Does not inspect genotypes, correctness, or sequence-derived truth labels.
Reservation is immutable; this is not permission to unblind its outcomes.
"""
import collections
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write(path, rows, fields):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    out = HERE / 'validation'
    if out.exists():
        raise FileExistsError('Validation reservation already exists; never silently reselect')
    sources = dict(cohorts=REPO/'hla-bench/source/cohorts.tsv',
                   matched=REPO/'hla-spechla-pg/series20260918/cohort.tsv',
                   legacy=REPO/'hla-bench/results/legacy_t1k_scores.tsv',
                   metadata=REPO/'hla-typer/source/haplotype_donors.tsv',
                   groups=REPO/'hla-spechla-pg/series20260918/validation_groups.tsv')
    cohorts = read(sources['cohorts'])
    exposed = {r['donor'] for r in read(sources['matched'])}
    # Only the donor identifiers are consumed from old result tables.
    exposed.update(r['donor'] for r in read(sources['legacy']))
    exposed.update(r['donor'] for r in cohorts if r['dev'] == '1' or r['cohort'] != 'test_A_loo')
    families = {r['family'] for r in cohorts if r['donor'] in exposed}
    reserved = sorted((r for r in cohorts if r['cohort'] == 'test_A_loo'
                       and r['donor'] not in exposed and r['family'] not in families), key=lambda r:r['donor'])
    assert reserved
    ids = {r['donor'] for r in reserved}
    assert len(ids) == len(reserved)
    groups = read(sources['groups'])
    meta = read(sources['metadata'])
    # Exclude the entire reserved set and its families, as well as development
    # donors/families, from every new inference panel, including alias assemblies.
    inference_excluded_donors = ids | exposed
    inference_excluded_families = families | {r['family'] for r in reserved}
    hit_groups = {r['group'] for r in groups if r['donor'] in inference_excluded_donors
                  or r['group'] in inference_excluded_families}
    excluded = {r['hap_id'] for r in groups if r['group'] in hit_groups}
    excluded |= {r['hap_id'] for r in meta if r['donor_id'] in inference_excluded_donors}
    assert not {r['hap_id'] for r in meta if r['donor_id'] in ids} - excluded
    out.mkdir()
    write(out/'cohort.tsv', [{k:r[k] for k in ('donor','family','population','superpopulation','cram')}
                            for r in reserved], ['donor','family','population','superpopulation','cram'])
    write(out/'excluded_paths.tsv', [r for r in meta if r['hap_id'] in excluded], list(meta[0]))
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    record = dict(status='reserved_not_unblinded',
                  rule='All test_A_loo families with no development, matched-64, Gourraud or indexed legacy-T1K exposure; no outcome-based selection.',
                  donors=len(ids), families=len({r['family'] for r in reserved}),
                  strata=dict(collections.Counter(r['superpopulation'] for r in reserved)),
                  excluded_panel_haplotypes=len(excluded),
                  source_sha256={k:digest(p) for k,p in sources.items()},
                  cohort_sha256=digest(out/'cohort.tsv'), excluded_paths_sha256=digest(out/'excluded_paths.tsv'),
                  limitation='Exposure audit covers repository cohort metadata and indexed legacy outputs; further discovered prior evaluations invalidate the affected families before unblinding, with an explicit audit amendment.')
    (out/'RESERVATION.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
