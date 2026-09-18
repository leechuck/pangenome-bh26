"""Check actual observed-panel sources against reserved HGSVC donors and families."""
import csv
import json
from pathlib import Path
from build_graph import sha

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def run():
    cohort_path = HERE/'validation-next/cohort.tsv'
    cohort = list(csv.DictReader(cohort_path.open(), delimiter='\t'))
    donors = {r['donor'] for r in cohort}
    families = {r['family'] for r in cohort}
    ped_path = REPO/'hla-pilot/source/1000G.ped'
    pedigree = {r['SampleID']:r['FamilyID'] for r in csv.DictReader(ped_path.open(), delimiter=' ')}
    metadata_path = REPO/'hla-typer/source/haplotype_donors.tsv'
    metadata = list(csv.DictReader(metadata_path.open(), delimiter='\t'))
    by_hap = {r['hap_id']:r['donor_id'] for r in metadata}
    assert len(by_hap) == len(metadata)
    reference = HERE/'references/observed-v2'
    manifest = json.loads((reference/'COMPLETE.json').read_text())
    assert manifest['status'] == 'complete'
    output_hashes = manifest['output_sha256']
    checked = 0
    source_donors = set()
    hits = []
    for filename, digest in output_hashes.items():
        assert sha(reference/filename) == digest, filename
        if not filename.endswith('.json'):
            continue
        for row in json.loads((reference/filename).read_text()):
            assert row['sources'], (filename, row['path'])
            for source in row['sources']:
                hap = '#'.join(source.split('#')[:2])
                donor = by_hap[hap]  # Unmapped source is an error, never silently excluded.
                checked += 1
                source_donors.add(donor)
                if donor in donors or pedigree.get(donor) in families:
                    hits.append(dict(source=source, donor=donor, panel=filename))
    unknown = sorted(d for d in source_donors if d not in pedigree)
    record = dict(scope='Actual observed reference source paths, both panels and all eight loci',
        cohort_sha256=sha(cohort_path), pedigree_sha256=sha(ped_path), metadata_sha256=sha(metadata_path),
        reference_manifest_sha256=sha(reference/'COMPLETE.json'),
        source_occurrences_checked=checked, source_donors=len(source_donors),
        reserved_donors=len(donors), reserved_families=len(families), overlaps=hits,
        source_donors_without_1000g_pedigree=unknown,
        limitation='Family relationships unavailable in the pinned pedigree cannot be inferred. Existing source aliases use the pinned haplotype-to-donor map.')
    assert not hits, hits
    output = HERE/'validation-next/GRAPH_SOURCE_AUDIT.json'
    if output.exists():raise FileExistsError(output)
    output.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':run()
