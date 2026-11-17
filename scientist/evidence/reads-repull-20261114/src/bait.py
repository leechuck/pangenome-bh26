"""Part B, rule 2 (truth-free): build the bait k-mer set.

  python bait.py prep    -> work/bait/ipd_kmers.fa (canonical 31-mers of every IPD genomic allele in
                            genomic.fa, all genes) and work/bait/ref_outside_mhc.fa (GRCh38 analysis set:
                            primary chromosomes, unplaced/unlocalised, decoys and EBV; chr6:28,500,000-33,500,000
                            masked with N; *_alt and HLA-* contigs left out)
  (jellyfish counts which IPD k-mers occur in ref_outside_mhc.fa; see bait.sbatch)
  python bait.py final   -> work/bait/bait_kmers.txt = IPD k-mers minus those found outside the MHC
No donor data and no truth are read.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113/src')
import common as C          # noqa: E402

ST = Path('/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114')
OUT = ST / 'work/bait'
REF = C.BENCH / 'hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa'
K = 31
COMP = str.maketrans('ACGT', 'TGCA')


def canon(k):
    r = k[::-1].translate(COMP)
    return min(k, r)


def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    assert C.sha256_file(C.GENOMIC_FA) == C.GENOMIC_SHA
    ks = set()
    n = 0
    for name, ints, seq in C.read_fasta_headers_and_seqs(C.GENOMIC_FA):
        n += 1
        for i in range(len(seq) - K + 1):
            k = seq[i:i + K]
            if 'N' not in k:
                ks.add(canon(k))
    with open(OUT / 'ipd_kmers.fa', 'w') as f:
        for i, k in enumerate(sorted(ks)):
            f.write(f'>{i}\n{k}\n')
    with open(REF) as fi, open(OUT / 'ref_outside_mhc.fa', 'w') as fo:
        keep, name, pos = False, None, 0
        for line in fi:
            if line.startswith('>'):
                name = line[1:].split()[0]
                keep = not (name.startswith('HLA-') or name.endswith('_alt'))
                pos = 0
                if keep:
                    fo.write(line)
                continue
            if not keep:
                continue
            s = line.rstrip('\n')
            if name == 'chr6' and pos + len(s) > 28500000 and pos < 33500000:
                s = ''.join('N' if 28500000 <= pos + j < 33500000 else c for j, c in enumerate(s))
            pos += len(s)
            fo.write(s + '\n')
    C.save_json(OUT / 'PREP.json', dict(ipd_alleles=n, ipd_kmers=len(ks), k=K, reference=str(REF),
                                        masked='chr6:28500000-33500000 (0-based half-open)',
                                        excluded='HLA-* and *_alt contigs'))


def final():
    found = set()
    with open(OUT / 'found_outside_mhc.txt') as f:
        for line in f:
            k, c = line.split()
            if int(c) > 0:
                found.add(canon(k.upper()))
    n_all = n_keep = 0
    with open(OUT / 'ipd_kmers.fa') as f, open(OUT / 'bait_kmers.txt', 'w') as fo:
        for line in f:
            if line.startswith('>'):
                continue
            k = line.strip()
            n_all += 1
            if k not in found:
                fo.write(k + '\n')
                n_keep += 1
    info = C.load_json(OUT / 'PREP.json')
    info.update(ipd_kmers_found_outside_mhc=len(found), bait_kmers=n_keep)
    C.save_json(OUT / 'BAIT.json', info)


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID')
    {'prep': prep, 'final': final}[sys.argv[1]]()
