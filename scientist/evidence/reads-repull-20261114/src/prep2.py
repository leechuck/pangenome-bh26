"""Part B step 6: build read set D-real2 for one development donor.

Rule (truth-free): every read pair of the donor's whole-genome sequencing data (all primary paired records
of the 1000 Genomes/HPRC CRAM, mapped anywhere or unmapped, as prepared in baselines/prepared/wgs.<donor>;
for NA18944/45/70 the original paired FASTQ) that T1K's own k-mer read extractor (fastq-extractor,
`-s 0.97`, reference = IPD genomic.fa of all HLA and related genes, as called by run-t1k --preset hla-wgs)
selects. These selections already exist as calls_candidate_{1,2}.fq of the completed full-WGS T1K
baseline run baselines/results/wgs.<donor>/t1k/t1k_genomic_ipd (same T1K runtime, same genomic.fa);
this script copies them (gzip) and records provenance. No truth is read.

  python prep2.py DONOR     -> work/reads2/DONOR/{r1.fq.gz,r2.fq.gz,PAIRS.json}
"""
import gzip
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402
import run_donor as RD      # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')


def main():
    donor = sys.argv[1]
    assert os.environ.get('SLURM_JOB_ID')
    base = C.ZEN / 'baselines/results' / f'wgs.{donor}' / 't1k'
    rec = C.load_json(base / 'COMPLETE.json')
    assert rec['status'] == 'complete' and rec['donor'] == donor
    run = base / 't1k_genomic_ipd'
    cmd = [c for c in rec['commands'] if c[c.index('-f') + 1].endswith('/genomic.fa')][0]
    assert C.sha256_file(run / 'calls_genotype.tsv') == rec['output_sha256']['t1k_genomic_ipd/calls_genotype.tsv']
    assert C.sha256_file(C.GENOMIC_FA) == C.GENOMIC_SHA
    log = (run / 'run.log').read_text()
    assert 'fastq-extractor' in log and '-s 0.97' in log and log.rstrip().endswith('Finish.')
    out = ST / 'work/reads2' / donor
    out.mkdir(parents=True, exist_ok=True)
    names, counts = [], []
    for m in (1, 2):
        src = run / f'calls_candidate_{m}.fq'
        dst = out / f'r{m}.fq.gz'
        with open(dst, 'wb') as f:
            subprocess.run(['gzip', '-c', '-4', str(src)], stdout=f, check=True)
        n, nm = 0, []
        with gzip.open(dst, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    n += 1
                    if n <= 1000:
                        nm.append(line.split()[0])
        assert (i + 1) % 4 == 0
        names.append(nm)
        counts.append(n)
    assert counts[0] == counts[1] and names[0] == names[1], 'mates not paired'
    regional = RD.pass1_spec('D-real', donor, 'A')
    old_pairs = None
    try:
        old = RD.to_host(regional['r1']).parent.parent / 'PAIRS.json'
        old_pairs = C.load_json(old).get('pairs')
    except Exception:
        pass
    C.save_json(out / 'PAIRS.json', dict(
        status='paired_fastq_verified', donor=donor, readset='D-real2', pairs=counts[0],
        old_regional_pairs=old_pairs, old_regional_r1=regional['r1'],
        rule='T1K fastq-extractor (-s 0.97, IPD genomic.fa) over all primary read pairs of the whole-genome data',
        source_run=str(run), source_command=cmd, source_wgs_pairs_sha256=rec.get('pairs_sha256'),
        r1=str(out / 'r1.fq.gz'), r2=str(out / 'r2.fq.gz'),
        r1_sha256=C.sha256_file(out / 'r1.fq.gz'), r2_sha256=C.sha256_file(out / 'r2.fq.gz'),
        truth_loaded=False, job=os.environ['SLURM_JOB_ID']))


if __name__ == '__main__':
    main()
