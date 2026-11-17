"""Step 2b item 1a (development donors only, reported only): apply the rule-2 bait test directly to every
read pair of the whole-genome data (no T1K prefilter) and count rule-2 pairs that the prefilter misses.

  python pref_check.py patterns            -> work/pref/bait_both.txt (bait k-mers, both orientations, for grep -F)
  python pref_check.py evaluate DONOR      -> work/pref/DONOR/PREF.json
     input: work/pref/DONOR/hits_r{1,2}.tsv = "name<TAB>sequence" of every read with >= 1 bait 31-mer anywhere
     (grep -F superset), from the full-WGS FASTQ baselines/prepared/wgs.DONOR (all primary pairs of the CRAM).
Counts only.
"""
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import select3 as S          # noqa: E402  (same hits() rule, same bait file)
C, RD = S.C, S.RD
ST = S.ST
COMP = str.maketrans('ACGT', 'TGCA')


def patterns():
    out = ST / 'work/pref'
    out.mkdir(parents=True, exist_ok=True)
    with open(ST / 'work/bait/bait_kmers.txt') as f, open(out / 'bait_both.txt', 'w') as fo:
        for line in f:
            k = line.strip()
            fo.write(k + '\n' + k[::-1].translate(COMP) + '\n')


def names_of(path):
    out = set()
    with gzip.open(path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                out.add(S.qname(line))
    return out


def evaluate(donor):
    w = ST / 'work/pref' / donor
    bait = set(l.strip() for l in open(ST / 'work/bait/bait_kmers.txt'))
    passing, n_hit = set(), 0
    for m in (1, 2):
        with open(w / f'hits_r{m}.tsv') as f:
            for line in f:
                nm, seq = line.rstrip('\n').split('\t')
                n_hit += 1
                if S.hits(seq, bait):
                    passing.add(nm.split()[0].lstrip('@').removesuffix('/1').removesuffix('/2'))
    r3 = ST / 'work/reads3' / donor
    sel = names_of(r3 / 'r1.fq.gz')
    old = RD.pass1_spec('D-real', donor, 'A')
    oldn = names_of(RD.to_host(old['r1']))
    added = sel - oldn
    full_rule = passing | oldn          # rule 2 applied without the prefilter
    C.save_json(w / 'PREF.json', dict(
        donor=donor, reads_with_any_bait_kmer=n_hit, pairs_passing_bait_rule_full_data=len(passing),
        pairs_passing_and_not_in_old_pull=len(passing - oldn),
        rule2_added_pairs_with_prefilter=len(added),
        missed_by_prefilter=len((passing - oldn) - added),
        added_with_prefilter_but_not_passing_full_scan=len(added - passing),
        rule2_total_with_prefilter=len(sel), rule2_total_without_prefilter=len(full_rule),
        job=os.environ.get('SLURM_JOB_ID')))


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID')
    if sys.argv[1] == 'patterns':
        patterns()
    else:
        assert sys.argv[2] in ('HG00558', 'HG00438'), 'development check donors only'
        evaluate(sys.argv[2])
