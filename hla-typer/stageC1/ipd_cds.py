#!/usr/bin/env python3
"""CDS catalogues from IPD-IMGT/HLA 3.65.0 nuc FASTAs, deduplicated by sequence.

Two outputs: complete CDS (graft donors for a whole coding sequence) and partial CDS
(most DRB1/DRB3/4/5 alleles are typed over exon 2 only; these are grafted into the matching
sub-span of a backbone CDS).

Keeps a CDS only if it starts with ATG, ends with a stop codon, is a whole number of codons and
has no internal stop (null alleles with internal stops are therefore excluded as graft donors).
Output: source/ipd_cds/<gene>.tsv  (cds_sha256, length, alleles ';'-joined, sequence).
"""
import hashlib,collections,csv
from pathlib import Path
H=Path(__file__).resolve().parent.parent
IMGT=H.parent/'hla-analysis/source/imgt/fasta'
GENES=['A','B','C','DRB1','DRB3','DRB4','DRB5','DQA1','DQB1','DPA1','DPB1']
STOPS={'TAA','TAG','TGA'}

def records(path):
    name=None;seq=[]
    for line in open(path):
        if line[0]=='>':
            if name:yield name,''.join(seq)
            name=line.split()[1];seq=[]
        else:seq.append(line.strip().upper())
    if name:yield name,''.join(seq)

def complete(s):
    if len(s)%3 or not s.startswith('ATG') or s[-3:] not in STOPS or set(s)-set('ACGT'):return False
    return not any(s[i:i+3] in STOPS for i in range(0,len(s)-3,3))

def main():
    out=H/'source/ipd_cds';out.mkdir(parents=True,exist_ok=True)
    part=H/'source/ipd_partial';part.mkdir(parents=True,exist_ok=True)
    for g in GENES:
        stem='DRB345' if g in ('DRB3','DRB4','DRB5') else g
        by=collections.defaultdict(list);partial=collections.defaultdict(list);total=0
        for name,s in records(IMGT/f'{stem}_nuc.fasta'):
            if name.split('*')[0]!=g:continue
            total+=1
            if complete(s):by[s].append(name)
            elif len(s)>=200 and not set(s)-set('ACGT'):partial[s].append(name)
        for d,table in [(out,by),(part,partial)]:
            with open(d/f'{g}.tsv','w') as f:
                w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['cds_sha256','length','alleles','sequence'])
                for s,names in sorted(table.items(),key=lambda x:x[1][0]):
                    w.writerow([hashlib.sha256(s.encode()).hexdigest(),len(s),';'.join(names),s])
        print(g,'alleles',total,'complete',sum(map(len,by.values())),'distinct',len(by),'| partial',sum(map(len,partial.values())),'distinct',len(partial))
if __name__=='__main__':main()
