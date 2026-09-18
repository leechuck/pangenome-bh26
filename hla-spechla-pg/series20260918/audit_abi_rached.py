#!/usr/bin/env python3
"""Audit the published computational comparator without promoting it to laboratory truth."""
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = HERE / 'abi_rached2018'
GENES = ('A', 'B', 'C', 'DQB1', 'DRB1')


def read(path, delimiter='\t'):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def alternatives(raw):
    """Expand slash shorthand, preserving null-expression suffixes and source flags."""
    if raw in ('', 'None'):
        return [], 'missing'
    if '*' in raw:
        return [], 'source_asterisk_unresolved'
    result = []
    for group in raw.split():
        parts = group.split('/')
        if not re.fullmatch(r'\d+:\d+[A-Z]?', parts[0]):
            return [], 'invalid_source_value'
        first = parts[0].split(':')[0]
        for part in parts:
            allele = part if ':' in part else first + ':' + part
            if not re.fullmatch(r'\d+:\d+[A-Z]?', allele):
                return [], 'invalid_source_value'
            result.append(allele)
    return sorted(set(result)), 'ambiguous' if len(set(result)) > 1 else 'resolved_two_field'


def write(name, rows):
    with (DATA / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    raw = read(DATA / '20181129_HLA_types_full_1000_Genomes_Project_panel.txt')
    assert len(raw) == len({r['Sample ID'] for r in raw}) == 2693
    panel = read(REPO / 'hla-typer/source/haplotype_donors.tsv')
    pedigree = {r['SampleID']: r for r in read(REPO / 'hla-pilot/source/1000G.ped', ' ')}
    family = lambda d: pedigree.get(d, {}).get('FamilyID', d)
    panel_donors = {r['donor_id'] for r in panel if r['cohort'] != 'REF'}
    panel_families = {family(d) for d in panel_donors}
    gourraud = {r['id'] for r in read(REPO / '1000g_ground_truth/1000G_2014/20140702_hla_diversity.txt', ' ')}
    cohort = {r['donor'] for r in read(HERE / 'cohort.tsv')}
    urls = {}
    for line in (REPO / 'hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv').read_text().splitlines():
        donor, path = line.split('\t')[:2]
        urls[donor] = 'https://s3.amazonaws.com/1000genomes/' + path
    records, donors = [], []
    for r in raw:
        d = r['Sample ID']
        reason = ('graph_donor' if d in panel_donors else 'graph_family' if family(d) in panel_families
                  else 'pedigree_unverified' if d not in pedigree else 'no_cram_manifest' if d not in urls else 'eligible')
        donors.append(dict(donor=d, ancestry=r['Region'], population=r['Population'], family=family(d),
                           gourraud_overlap=int(d in gourraud), matched64=int(d in cohort),
                           eligibility=reason, cram=urls.get(d, 'UNAVAILABLE')))
        for gene in GENES:
            for copy in (1, 2):
                value = r[f'HLA-{gene} {copy}']
                alleles, status = alternatives(value)
                records.append(dict(donor=d, gene=gene, copy=copy, raw=value or 'EMPTY',
                                    alternatives=';'.join(gene+'*'+a for a in alleles) or 'UNRESOLVED',
                                    status=status, evidence='PolyPheMe2018_computational'))
    write('donor_audit.tsv', donors)
    write('calls.tsv', records)
    eligible = [r for r in donors if r['eligibility'] == 'eligible']
    write('eligible_donors.tsv', eligible)
    report = dict(samples=len(raw), gourraud_overlap=sum(r['gourraud_overlap'] for r in donors),
                  matched64_overlap=sum(r['matched64'] for r in donors),
                  eligibility=dict(collections.Counter(r['eligibility'] for r in donors)),
                  eligible_ancestry=dict(collections.Counter(r['ancestry'] for r in eligible)),
                  allele_status=dict(collections.Counter(r['status'] for r in records)),
                  source_sha256=hashlib.sha256((DATA / '20181129_HLA_types_full_1000_Genomes_Project_panel.txt').read_bytes()).hexdigest(),
                  policy='Secondary two-field concordance only. No four-field truth. Preserve ambiguity; unresolved annotations and malformed values excluded explicitly. Graph donors and known families excluded from independent cohort; matched64 evaluated only with their held-out panels. Eligibility is not a job submission.')
    (DATA / 'audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
