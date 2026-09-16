#!/usr/bin/env python3
"""Post-freeze audit for diagnostic probe dropout in held-out gene sequences."""
import csv,json
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    probes={int(s.split()[0]):int(s.split()[1]) for s in (ROOT/'source/marker_groups.tsv').read_text().splitlines() if int(s.split()[1])}
    audit={r['gene_id']:r for r in read(ROOT/'results/validation_C4_diagnostic_audit.tsv')}
    groups=json.loads((ROOT/'source/marker_design.json').read_text())['groups'];bases={'A':0,'C':1,'G':2,'T':3};rows=[]
    for gene in SeqIO.parse(ROOT/'source/validation_C4_genes.fa','fasta'):
        f=rev=n=0;hits={g:set() for g in groups}
        for c in str(gene.seq).upper():
            if c not in bases:f=rev=n=0;continue
            b=bases[c];f=((f<<2)|b)&((1<<62)-1);rev=(rev>>2)|((3-b)<<60);n+=1
            if n<31:continue
            key=min(f,rev);mask=probes.get(key,0)
            for g,bit in groups.items():
                if mask&bit:hits[g].add(key)
        a=audit[gene.id];rows.append(dict(gene_id=gene.id,donor=a['donor'],form=a['inferred_form'],**{g+'_unique_probes':len(h) for g,h in hits.items()}))
    with open(ROOT/'results/validation_probe_coverage.tsv','w') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    missing=[r for r in rows if r[r['form'][2]+'_unique_probes']==0]
    print('Audited',len(rows),'genes; zero A/B probe coverage:',len(missing))
    for r in rows:
        if r['donor']=='HG04204' or r in missing:print(r)
if __name__=='__main__':main()
