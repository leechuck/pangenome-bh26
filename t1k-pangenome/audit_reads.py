"""Check reserved donors' public CRAM/CRAI availability without reading outcomes."""
import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve(donor, recorded, public):
    if donor not in public:
        raise ValueError('No public CRAM manifest entry for ' + donor)
    url = 'https://s3.amazonaws.com/1000genomes/' + public[donor]
    if recorded.startswith('https://') and recorded != url:
        raise ValueError('Conflicting CRAM sources for ' + donor)
    if not url.endswith('/' + donor + '.final.cram'):
        raise ValueError('Manifest donor/path mismatch for ' + donor)
    return url


def probe(url):
    try:
        request = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(request, timeout=30) as response:
            size = int(response.headers.get('Content-Length', '0'))
            return dict(url=url, status=response.status, bytes=size,
                        accessible=response.status == 200 and size > 0,
                        etag=response.headers.get('ETag'))
    except (OSError, ValueError, urllib.error.URLError) as error:
        return dict(url=url, accessible=False, error=str(error))


def audit(reservation, manifest, output):
    if output.exists():
        raise FileExistsError(output)
    cohort = reservation / 'cohort.tsv'
    lock = json.loads((reservation / 'RESERVATION.json').read_text())
    if sha(cohort) != lock['cohort_sha256']:
        raise ValueError('Reserved cohort changed')
    with cohort.open() as stream:
        donors = list(csv.DictReader(stream, delimiter='\t'))
    public = {}
    for line in manifest.read_text().splitlines():
        donor, path = line.split('\t')[:2]
        if donor in public and public[donor] != path:
            raise ValueError('Conflicting public manifest entries: ' + donor)
        public[donor] = path
    entries = []
    # Resolve all IDs before network requests; no truth files are consumed.
    for row in donors:
        url = resolve(row['donor'], row['cram'], public)
        entries.append(dict(donor=row['donor'], recorded_source=row['cram'],
                            cram_url=url, crai_url=url + '.crai'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        futures = [(entry, kind, pool.submit(probe, entry[kind + '_url']))
                   for entry in entries for kind in ('cram', 'crai')]
        for entry, kind, future in futures:
            entry[kind] = future.result()
    report = dict(observed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  cohort_sha256=sha(cohort), public_manifest_sha256=sha(manifest),
                  auditor_sha256=sha(Path(__file__)),
                  scope='HTTP availability only; not CRAM decoding or index validation. '
                        'No donor is excluded by this audit; failed requests need investigation.',
                  available=sum(all(e[k]['accessible'] for k in ('cram', 'crai')) for e in entries),
                  planned=len(entries), donors=entries)
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'donors'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reservation', type=Path, required=True)
    parser.add_argument('--public-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.reservation, args.public_manifest, args.output)
