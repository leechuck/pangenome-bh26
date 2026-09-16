#!/usr/bin/env python3
"""Training paths only: dense diagnostic probes plus the original marker sketch."""
import csv,gzip,json
from pathlib import Path
import numpy as np
from Bio import SeqIO
ROOT=Path(__file__).resolve().parent
def main():
    families={r['family'] for r in csv.DictReader(open(ROOT/'source/validation_donors.tsv'),delimiter='\t')}
    rows=[r for r in csv.DictReader(open(ROOT.parent/'hla-structural/source/catalogue.tsv'),delimiter='\t') if r['locus']=='RCCX' and r['family'] not in families]
    keys=[int(k) for k in (ROOT/'source/markers.txt').read_text().splitlines()];lookup={k:i for i,k in enumerate(keys)}
    seq={}
    with gzip.open(ROOT.parent/'hla-structural/source/locus_sequences.fa.gz','rt') as f:
        for r in SeqIO.parse(f,'fasta'):
            if r.id.startswith('RCCX'):seq[r.id]=str(r.seq)
    X=np.zeros((len(rows),len(keys)),dtype=np.uint16);mask=(1<<62)-1;bases={'A':0,'C':1,'G':2,'T':3}
    for i,r in enumerate(rows):
        f=rev=n=0
        for c in seq[r['sequence_id']]:
            if c not in bases:f=rev=n=0;continue
            b=bases[c];f=((f<<2)|b)&mask;rev=(rev>>2)|((3-b)<<60);n+=1
            if n>=31:
                j=lookup.get(min(f,rev))
                if j is not None:X[i,j]+=1
    np.save(ROOT/'source/training_profiles.npy',X)
    (ROOT/'source/training_rows.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('Training profiles',X.shape,'held-out families excluded',len(families))
if __name__=='__main__':main()
